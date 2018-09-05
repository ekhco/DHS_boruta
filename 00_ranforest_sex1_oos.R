# ------------------------------------------------------------
# A random forest prediction using Zimbabwe Standard DHS 2015
# Zimbabwe 2015 selected because it has IR, MR, and AR
# ky71_flattened.dta
#
# Eric Chow, 2018-05-14
# ------------------------------------------------------------

rm(list=ls())
gc()
# remove.packages("haven")
# devtools::install_version("haven", version = "1.1.0", repos = "http://cran.us.r-project.org")

library(tidyverse)
library(haven)
library(foreign)
library(caret)
library(stringr)
library(randomForest)
library(ggplot2)
library(plyr)
library(devtools)   # install.packages("devtools")
library(reprtree)   # install_github('araastat/reprtree')
library(glmnet)     # Library for glmnet ie: LASSO
library(doParallel) # Library for parallel ops
library(Boruta)     # feature selection heuristic
# setwd("D:/EricChow/DHS_ranforest")
setwd("~/QSU/DHS_boruta")
# load missinginess and murho fucntions
source("0_fn.R")
# # *****************************************************
# # A function that identifies the degree of missinginess
# # in a variable.
# missingness <- function(var) {
#     if (is.character(var)) {
#         count_miss <- sum(as.numeric(var == ""))
#     }
#     else {
#         count_miss <- sum(as.numeric(is.na(var)))
#     }
#     # count_blnk <- sum(as.numeric(ky[,var] == ""))
#     count_ttl <- length(var)
#     return(count_miss/count_ttl)
# }


# -------------------------------------------
# read in the data for Zimbabwe 9143 vars
# ky.raw <- read_dta("/Volumes/echow/DATA_DHS/Zimbabwe_2015/ky71_flattened_v12.dta")
ky.raw <- read_dta("/Users/echow/DHS_live/Kenya/Standard\ DHS\ 2003/ke42_flattened_v12.dta")

ky.raw <- data.frame(ky.raw)
head(ky.raw[,1:10])
dim(ky.raw) # 43,706 obs,   9143 vars

# meta data
meta <- read_dta("/Users/echow/DHS_live/Kenya/Standard\ DHS\ 2003/meta_data.dta")
meta <- data.frame(meta)
head(meta)
dim(meta)    # 5,809 vars,  5 metavars *** HMMMM why only 64% of vars have metavar?
# this meta definitely corresponds to ky71_flattened... which vars don't have meta data?


# which variables do not have a metadat?
meta_in_ky <- names(ky.raw) %in% meta$var_name
vars_not_used <- names(ky.raw)[!meta_in_ky]
# oh, turns out that all these variables are 100% missing variables! perfect.

# -------------------------------------------
# keep observations with HIV test result data
ky <- ky.raw[!is.na(ky.raw$hiv03), ]
dim(ky)                # 32,192 obs
nrow(ky)/nrow(ky.raw)  # 73.7% of obs had HIV
# keep only if HIV status == 0 or 1
ky <- ky[ky$hiv03 %in% c(0,1), ]
nrow(ky)               # still 32,192 obs - so they were all 0,1

# what percent had HIV data?
if (nrow(ky.raw) != nrow(ky)) {message( nrow(ky.raw) - nrow(ky)," (", round(100*(nrow(ky.raw)-nrow(ky))/nrow(ky.raw)), "%) observations were dropped." )}

# what is the prevalence of HIV in the original dataset?
table(ky$hiv03)           # 91% HIV (-), 8.7% HIV (+)
table(ky$hiv03)/nrow(ky)

# how many men/women?  HV104 is sex
table(ky$hv104)/nrow(ky)
table(ky$hiv03, ky$hv104)  # 6.8% of men are hiv[+], 10.3% of females are hiv[+]
table(ky$hiv03, ky$hv104)[2,]/table(ky$hv104) # proportion of men/women HIV+
prop.test(table(ky$hv104, ky$hiv03))

# put ky into another dataset to preserve ky for later
ky_ <- ky

# -------------------------------------------
# subset out a validation set
set.seed(314)
train <- createDataPartition(y = ky_$hiv03, p=0.80, list=FALSE)

# Subset to an out-of-sample validation set, and dataset for training
kyt <- ky_[train, ]    # used for fitting the model, and Mu and Rho etc
dim(kyt)               # 25,690 obs (~80%)
kyv  <- ky_[-train, ]  # THE OOS VALIDATION SET
dim(kyv)               # 6,422 obs  (~20%)





# ******************************************************************************
# try random forest once, with all genders 7-11-2018
table(ky$hiv03, ky$hv104)  # this is the N=32,112 population of known HIV status
dim(ky)                     # with 9143 vars (some empty)


# pre-process with mu and rho
kyt$train <- 1
kyv$train <- 0
# recombine the train and test sets
kya <- rbind(kyt, kyv)
# apply murho to both train and test
da <- murho(data = kya, mu = 0.15, rho = 0.2)
# resplit into train and test
kyt_2_2 <- da[da$train==1, ]
kyv_2_2 <- da[da$train==0, ]

dim(kyt_2_2)
dim(kyv_2_2)

# what is the prevalence of HIV in this dataset?
table(kyt_2_2$hiv03)           # 91% HIV (-), 8.7% HIV (+)
table(kyt_2_2$hiv03)/nrow(kyt_2_2)
#prevalence of HIV in VALIDATION
table(kyv_2_2$hiv03)
table(kyv$hiv03)/nrow(kyv)

# how many men/women?  HV104 is sex
table(kyt_2_2$hv104)/nrow(kyt_2_2)
table(kyt_2_2$hiv03, kyt_2_2$hv104)  # 6.8% of men are hiv[+], 10.3% of females are hiv[+]
table(kyt_2_2$hiv03, kyt_2_2$hv104)[2,]/table(kyt_2_2$hv104) # proportion of men/women HIV+
prop.test(table(kyt_2_2$hv104, kyt_2_2$hiv03))


# ---------------------------------------------
#  try to use the # try a random forest
str(kyt_2_2, list.len=10) # show factors and # levels

# just drop it
# names(kyt_2_2) != "shdist"
kyt_2_2 <- kyt_2_2[,!(names(kyt_2_2) %in% c("shdist","sdist"))]
registerDoParallel(4)

    # -----------------------------------------------------
    # 1. tweak the max # of vars - sqrt, % of vars, any number, mtry
    # 2. minium leaf size - avoid overfitting ie: node size
    # 3. number of trees, ntree
     # will auto center and scale
    nvars <- ncol(kyt_2_2) # sample about half of the variables (200/435)
    sqrt_nvars <- round(sqrt(nvars),0)
    nobs <- round(nrow(kyt_2_2),0)

    nvars # number of vars
    nobs # number of obseravations
    nodesize = 12 # min number in each leaf

    # set ntree and mtry
    ntree = round(nobs * nvars * 0.000131,0)  # eyeball around 150-200 trees for 3563*430 obs*vars
    mtry  = round(nvars / 2, 0) #

    ntree
    mtry
# ---------------------------------------------------------------
# function to get OOB sensitivity = 1 - class.error of hiv03==1
get_se <- function(model) {
    Se <- 1-median(model$err.rate[,3]) # gives error rate for each tree
    return(Se)
}

# -------------------------------------------------------------------
# trys a ran forest
try_rf <- function(ntree, mtry, nodesize) {
    modelFit <- NULL
    modelFit <- randomForest(hiv03 ~ . , data = kyt_2_2, importance = TRUE,
        ntree=ntree, strata = hv104, nodesize=nodesize, mtry=mtry)
    modelFit
    se <- get_se(modelFit)
    result <- c(ntree, mtry, nodesize, se)
    return(result)
}


# try Boruta
ntree=200; mtry=0.46;
b_fit <- Boruta(hiv03 ~. , data = kyt_2_2,
         ntree=ntree, mtry=mtry) # , nodesize=nodesize, strata=hv104)
b_fit
plot(b_fit)

# ------------------------------------------------------------------------------
# try to find best ntree, nodesize, and number of vars
#
# NTREES: try up to 150% of #obs / #nodesize (ie: try to grow enough single-node
# trees to use all the observations about once each)
#
# MTRY: try anywhere from 3 up to all the vars
#
# NODE SIZE: set at 12, so minimum of 12 observations in each leaf
# fileConn <- file("data/tuning_final.csv")

for (ntree in seq(60,1.5*nobs/nodesize,3)) {
        for (mtry in seq(3,nvars,3)) {
            cat(try_rf(ntree = ntree, mtry = mtry, nodesize = 12),"\n")
            # writeLines(str_c(sprintf("%10.9f",try_rf(ntree = ntree, mtry = mtry, nodesize = 12))), fileConn)
        }
}
# close(fileConn)

# get parametesr from random Forest
names(modelFit)
modelFit$confusion
modelFit$mtry
modelFit$ntree
modelFit$err.rate # gives error rate for each tree
get_se(modelFit)

# plot a big heat map


# SELECT BEST TUNING PARAMETERS HERE


# get the importance
imp_rf <- importance(modelFit)[,4]
imp_rf <- data.frame(imp_rf[order(-imp_rf)])
# list first 15 most important variables
head(imp_rf,15)

# get missingness of each vars
avar <- row.names(imp_rf)[234]

va <- row.names(imp_rf)[row.names(imp_rf) %in% names(ky)]
result <- NA
for (i in 1:length(va)) {
    result[i] <- missingness(ky[,va[i]])
}
miss_ky <- data.frame(var = va, miss = result)

miss_ky_s <- join(data.frame(var=row.names(imp_rf)), miss_ky, by="var", type="left", match="first")
miss_ky_s$completion <- 1 - miss_ky_s$miss
# barplot of variable importance
y <- miss_ky_s$completion + rnorm(nrow(miss_ky_s),0,0.01)
plot(imp_rf[,1], y, pch=20, xlab="Importance by mean decrease in Gini", ylab="Variable Completeness",
    main="Random Forest Variable Importance")
text(imp_rf[,1], y, row.names(imp_rf), cex=0.6, srt=45, pos=3)

# -------------------------
# plot a representative tree (but this sucks)
# reptree <- ReprTree(modelFit, kyt_2_2, metric='d2')
# plot(reptree)

# predict validation set
predictions <- predict(modelFit, newdata=kyv_2_2)
length(predictions)
# compare to actual
cm <- confusionMatrix(predictions, kyv_2_2$hiv03) # use confusion matrix to EVALUATE
cm





# ---------------------------------------------
#  try to use the vars that Kajal had identified from Kenya 2003, 2008
# read in Kajal's file

kenya_vars <- read.delim("data/comparison03-08.tsv", sep="\t", header=TRUE)
kenya_vars
kvars_l <- as.character(kenya_vars$Exposure) # the var names
# strip the level that Kajal added
kvars <- unique(substr(kvars_l, 1, str_length(kvars_l)-2))
length(kvars) # 105!

# get the variables in the trainign set that are in Kajal's
kvars_in <- (names(kyt_2_2) %in% kvars)  # the variables in kyt that are in kvars
kvi <- names(kyt_2_2)[kvars_in]
kvi # what are the vars that I have that are also in Kajal's?
table(kvars_in) # most of the kenya predictive variables didn't make it through murho processing

# which of the important vars is in the Kenya XWAS list?
rf_10 <- row.names(head(imp_rf,10))
rf_10 %in% kvars

# were all of the kenya variables in my data before murho?
kvars_in_orig <- (names(ky) %in% kvars)
table(kvars_in_orig)

# do predictions
rf_kenya <- randomForest(hiv03 ~ . , data = kyt_2_2[, c(kvi, "hiv03")], importance = TRUE)  # will auto center and scale
lg_kenya <- glm(hiv03 ~ . , data = kyt_2_2[, c(kvi, "hiv03")], family = "binomial")  # will auto center and scale

# predict validation set
prf_kenya <- predict(rf_kenya, newdata=kyv_2_2)
plg_kenya <- as.factor(round(predict(lg_kenya, newdata=kyv_2_2, type="response")))
table(plg_kenya)
# compare to actual
cmrf_kenya <- confusionMatrix(prf_kenya, kyv_2_2$hiv03) # use confusion matrix to EVALUATE
cmrf_kenya
cmlg_kenya <- confusionMatrix(plg_kenya, kyv_2_2$hiv03) # use confusion matrix to EVALUATE
cmlg_kenya








# ----------------------------------------------------------
# Try to do multiple imputation instead of mu?



# -----------
# STOP HERE!!
# -----------


if (1==0) {
    # ------------------------------------------------------------------------------
    # DO LASSO PREDICTION - code inspired by Don
    # Use sparse matrix for computational+space efficiency
      train_x <- sparse.model.matrix(hiv03~., kyt)[,-1]
      train_y <- as.factor(as.numeric(kyt$hiv03))  # Don did -1, why?
      test_x <- sparse.model.matrix(hiv03~., kyv)[,-1]
      test_y <- as.factor(as.numeric(kyv$hiv03))  # Don did -1, why?

      # Find the best lambda from our list via cross-validation (here, we didn't use a validation set). We do this for a fixed mu and rho
      registerDoParallel(8)
      cv.out <- cv.glmnet(train_x, train_y, family="binomial", alpha = 1, nfolds=10, type.measure="class", parallel=TRUE)
      best_lambda <- cv.out$lambda.min

      # Min. mean cross-validation error associated with best lambda
      min_CVE <- min(cv.out$cvm)

      # Construct lasso model of interest
      lasso.mod <- glmnet(train_x, train_y, family="binomial", alpha = 1) # alpha=1 => lasso
}








# *****************************************************
# cycle through a bunch of Mu and Rhos, and for each
# apply the restrictions and see how many are lost

# an array to store results
sumstats <- data.frame(mu=NA, rho=NA, vars=NA, cat_vars=NA, obs=NA, tn=NA, fn=NA, tp=NA, fp=NA, se=NA, sp=NA, ppv=NA, npv=NA)

# index var for storing meta results
i = 1
m = 5; r = 2
set.seed(314)
# takes about ... um ... some time to run! 5:26 -
# loop through Mu and Rho's
for (m in 1:9) {
    for (r in 1:4) {
        mu = m/20; rho = r/10
        d    <- NULL
        dtrn <- NULL
        dtst <- NULL
        modelFit <- NULL
        predictions <-NULL
        cm <- NULL

        try ({
            # apply mu and rho -------------------------------------

            # append train and test set so can murho (and get same levels)
            # but then split again
            kyt$train <- 1
            kyv$train <- 0

            # recombine the train and test sets
            kya <- rbind(kyt, kyv)
            # apply murho to both train and test
            da <- murho(data = kya, mu = mu, rho = rho)
            # resplit into train and test
            d  <- da[da$train==1, ]
            dv <- da[da$train==0, ]

            # count factors
            factor_i = 0
            for (var in names(d)) {
                if (is.factor(d[,var])) {
                    factor_i <- factor_i + 1
                }
            }
            # do I sample within the dataset after mu/rho, or use kyv, the validation set?
            do_OOS_valid = TRUE
            if (!(do_OOS_valid)) { # initially, I used within-sample validation, but don't
                                          # do this, this is wrong. Use the out of sample VALIDATION
                                          # set that was set aside initially.
                    # partition the data -----------------------------------
                    train <- createDataPartition(y = d$hiv03, p=0.80, list=FALSE)
                    # subset using partition
                    dtrn <- d[train, ]    # used for fitting the model, and Mu and Rho etc
                    dtst  <- d[-train, ]  # THE VALIDATION SET

                    # build a predictor ------------------------------------
                    modelFit <- randomForest(hiv03 ~ . , data = dtrn, importance = TRUE)  # will auto center and scale

                    # do within-sample validation --------------------------
                    #  ****** need to confirm with out of sample validation
                    predictions <- predict(modelFit, newdata=dtst)
                    # compare to actual
                    cm <- confusionMatrix(predictions, dtst$hiv03) # use confusion matrix to EVALUATE

                    # store meta results
                    sumstats[i,"mu"] <- mu
                    sumstats[i,"rho"] <- rho
                    sumstats[i,"obs"] <- dim(d)[1]
                    sumstats[i,"vars"] <- dim(d)[2]
                    sumstats[i,"cat_vars"] <- factor_i

                    sumstats[i,"tn"] <- cm$table[1,1]
                    sumstats[i,"fn"] <- cm$table[1,2]
                    sumstats[i,"tp"] <- cm$table[2,2]
                    sumstats[i,"fp"] <- cm$table[2,1]

                    sumstats[i,"se"] <- cm$byClass["Sensitivity"]
                    sumstats[i,"sp"] <- cm$byClass["Specificity"]
                    sumstats[i,"ppv"] <- cm$byClass["Pos Pred Value"]
                    sumstats[i,"npv"] <- cm$byClass["Neg Pred Value"]
            }
            # YES, do THIS one! use the validation data that was set aside prior to murho
            if (do_OOS_valid) { # use OOS validation, the CORRECT way
                # don't partition the data -----------------------------------
                # use whole murho's dataset to train
                # train <- createDataPartition(y = d$hiv03, p=0.80, list=FALSE)
                # subset using partition
                dtrn  <- d      # used for fitting the model, and Mu and Rho etc
                dtst  <- dv  # THE VALIDATION SET,  kyv (validation) was murho'd

                # build a predictor ------------------------------------
                message("running randomForest ...")
                modelFit <- randomForest(hiv03 ~ . , data = dtrn, importance = TRUE)  # will auto center and scale

                # do OUT OF sample validation --------------------------
                message("RF complete! predicting from test set ...")
                predictions <- predict(modelFit, newdata=dtst)
                # compare to actual
                cm <- confusionMatrix(predictions, dtst$hiv03) # use confusion matrix to EVALUATE

                # store meta results
                sumstats[i,"mu"] <- mu
                sumstats[i,"rho"] <- rho
                sumstats[i,"obs"] <- dim(d)[1]
                sumstats[i,"vars"] <- dim(d)[2]
                sumstats[i,"cat_vars"] <- factor_i

                sumstats[i,"tn"] <- cm$table[1,1]
                sumstats[i,"fn"] <- cm$table[1,2]
                sumstats[i,"tp"] <- cm$table[2,2]
                sumstats[i,"fp"] <- cm$table[2,1]

                sumstats[i,"se"] <- cm$byClass["Sensitivity"]
                sumstats[i,"sp"] <- cm$byClass["Specificity"]
                sumstats[i,"ppv"] <- cm$byClass["Pos Pred Value"]
                sumstats[i,"npv"] <- cm$byClass["Neg Pred Value"]
                message("DONE!\n")
            }
            i = i + 1
        })
    }
}

# show results
sumstats

# ------------------------------------------------------------------------------
# plot mu vs. # vars
pdf("figs/mu_vars_oos.pdf", height=4, width=4.7)
    ggplot(sumstats, aes(x=mu, y=vars)) + geom_point(aes(size=obs))
dev.off()

# plot mu vs. Sensitivity
pdf("figs/mu_se_oos.pdf", height=4, width=4.7)
    ggplot(sumstats[1:32,], aes(x=mu, y=se)) + geom_jitter(width=0.005, height=0.001) + geom_smooth(method="loess")
dev.off()










#      ~ fin ~
