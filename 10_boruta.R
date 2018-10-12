# ------------------------------------------------------------
# A random forest prediction using Zimbabwe Standard DHS 2015
# Zimbabwe 2015 selected because it has IR, MR, and AR
# survhiv71_flattened.dta
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
# library(plyr)
library(devtools)   # install.packages("devtools")
library(reprtree)   # install_github('araastat/reprtree')
library(glmnet)     # Library for glmnet ie: LASSO
library(doParallel) # Library for parallel ops
library(Boruta)     # feature selection heuristic
library(pROC)       # for the AUC
# setwd("D:/EricChow/DHS_ranforest")
setwd("~/QSU/DHS_boruta")
# load missinginess and murho fucntions
source("0_fn.R")

# get list of surveys

survey_dir <- "/Users/echow8/DHS_live_abstract/DHS_live_abstract/"
srvRDS_dir <- "/Users/echow8/DHS_live_abstract/DHS_RDS/"
mtadat_dir <- "/Users/echow8/DHS_live_abstract/meta_data/"
survey_files <- list.files(survey_dir)
surveys <- substr(survey_files,1,4) # the survey country/wave ie: zw61
# surveys_ <- surveys[48:59]


survey_filepath <- str_c(srvRDS_dir, surveys[8], ".rds", sep="")  # 27 works
mtadat_filepath <- str_c(mtadat_dir, surveys[8], "_metadat.dta", sep="")
sex = 1
# ------------------------------------------------------------------------------
# THIS IS WHERE FUNCTION WOULD START
# ------------------------------------------------------------------------------
do_boruta <- function(survey_filepath, mtadat_filepath, sex, this_survey, bor_ = 1, seed) {

    # RETURN VARIABLES:
    # PERCENT_W_HIV_DATA - the percent of the original data rows that has HIV results (pre-murho)
    # N_HAS_HIV_DATA     - the number of obs w HIV test result  (pre-murho)
    # N_HIV_POS          - the number of obs that are HIV positive  (pre-murho)
    # PREVALENCE         - the prevalence of HIV+  (pre-murho)
    # N_HAS_HIV_DATA_MURHO     - the number of obs w HIV test result  (pre-murho)
    # N_HIV_POS_MURHO          - the number of obs that are HIV positive  (pre-murho)
    # PREVALENCE_MURHO         - the prevalence of HIV+  (pre-murho)

    survey  <- readRDS(survey_filepath) # survey
    meta    <- data.frame(read_dta(mtadat_filepath)) # meta data
    sex     <- sex


    seed    <- 314
    mu      <- 0.15
    rho     <- 0.2
    split   <- 0.80

    # -------------------------------------------
    # Keep ONLY SEX == sex here
    survey <- survey[((as.numeric(survey$hv104) == sex) & !is.na(survey$hv104)),  ]
    SEX <- sex
    # -------------------------------------------

    # -------------------------------------------
    # keep observations with HIV test result data
    survhiv <- survey[!is.na(survey$hiv03), ]
    # keep only if HIV status == negative or positive
    survhiv <- survhiv[as.numeric(survhiv$hiv03) %in% c(1,2), ]
    survhiv$hiv03 <- factor(survhiv$hiv03)
    # dim(survhiv)                # 32,192 obs
    PERCENT_W_HIV_DATA <- nrow(survhiv)/nrow(survey)  # 73.7% of obs had HIV

    # what percent had HIV data?
    # if (nrow(survey) != nrow(survhiv)) {message( nrow(survey) - nrow(survhiv)," (", round(100*(nrow(survey)-nrow(survhiv))/nrow(survey)), "%) observations were dropped." )}

    N_HAS_HIV_DATA     <- nrow(survhiv)               # still 32,192 obs - so they were all 0,1
    N_HIV_POS  <- table(survhiv$hiv03)[2]
    PREVALENCE <- table(survhiv$hiv03)[2]/nrow(survhiv)

    # how many men/women?  HV104 is sex
        #
        # table(survhiv$hv104)/nrow(survhiv)
        # table(survhiv$hiv03, survhiv$hv104)  # 6.8% of men are hiv[+], 10.3% of females are hiv[+]
        # table(survhiv$hiv03, survhiv$hv104)[2,]/table(survhiv$hv104) # proportion of men/women HIV+
        # prop.test(table(survhiv$hv104, survhiv$hiv03))
        #
    # put survhiv into another dataset to preserve survhiv for later
    # survhiv_ <- survhiv

    # -------------------------------------------
    # subset out a validation set
    set.seed(seed)
    train <- createDataPartition(y = survhiv$hiv03, p = split, list=FALSE)

    # Subset to an out-of-sample validation set, and dataset for training
    survhivt <- survhiv[train, ]    # used for fitting the model, and Mu and Rho etc
    dim(survhivt)               # 25,690 obs (~80%)
    survhivv <- survhiv[-train, ]  # THE OOS VALIDATION SET
    dim(survhivv)               # 6,422 obs  (~20%)

    # pre-process with mu and rho
    survhivt$train <- 1
    survhivv$train <- 0
    # recombine the train and test sets
    survhiva <- rbind(survhivt, survhivv)

    # ----------------------------------------------
    # apply murho to both train and test
    # ----------------------------------------------
    da <- murho(data = survhiva, mu = mu, rho = rho, meta = meta)
    # da$hiv03 <- factor(da$hiv03)

    # ----------------------------------------------
    # drop variables that have too many levels (char)
    # ----------------------------------------------
    da <- da[,!(names(da) %in% c("shdist","sdist","_append"))]

    # -----------------------------------------------
    # drop variables that would perfectly predict HIV
    # -----------------------------------------------
    # shiv51 - confirmed HIV status
    # sh279  - result of determine hiv rdt
    # sh279a - result of unigold hiv rdt
    # sh278  - result measurement code of rapid hiv test
    sh279_vars <- names(da)[isin("sh279", names(da))]
    sh278_vars <- names(da)[isin("sh278", names(da))]

    da <- da[,!(names(da) %in% c("shiv51", sh279_vars, sh278_vars))]


    # resplit into train and test
    survhivt_murhoed <- da[da$train==1, ]
    survhivv_murhoed <- da[da$train==0, ]

    N_HAS_HIV_DATA_MURHO <- nrow(survhivt_murhoed)
    N_HIV_POS_MURHO <- table(survhivt_murhoed$hiv03)[2]           # 91% HIV (-), 8.7% HIV (+)
    PREVALENCE_MURHO <- table(survhivt_murhoed$hiv03)[2]/nrow(survhivt_murhoed)
    N_VARS_MURHO <- ncol(survhivt_murhoed)

    #prevalence of HIV in VALIDATION
    # table(survhivv_murhoed$hiv03)
    # table(survhivv$hiv03)/nrow(survhivv)

    # how many men/women?  HV104 is sex
    # table(survhivt_murhoed$hv104)/nrow(survhivt_murhoed)
    # table(survhivt_murhoed$hiv03, survhivt_murhoed$hv104)  # 6.8% of men are hiv[+], 10.3% of females are hiv[+]
    # table(survhivt_murhoed$hiv03, survhivt_murhoed$hv104)[2,]/table(survhivt_murhoed$hv104) # proportion of men/women HIV+
    # prop.test(table(survhivt_murhoed$hv104, survhivt_murhoed$hiv03))



        # -----------------------------------------------------
        # 1. tweak the max # of vars - sqrt, % of vars, any number, mtry
        # 2. minium leaf size - avoid overfitting ie: node size
        # 3. number of trees, ntree
         # will auto center and scale
        nvars    <- ncol(survhivt_murhoed) # sample about half of the variables (200/435)
        nobs     <- round(nrow(survhivt_murhoed),0)
        nodesize <- 12 # min number in each leaf

        nvars # number of vars
        nobs # number of obseravations
        nodesize

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


    cat("random foresting...  ")
    # Run a random forest with all features
    rf_fit <- NULL
    rf_fit <- randomForest(hiv03 ~ . , data = survhivt_murhoed, importance = TRUE,
        ntree=ntree, mtry=mtry, nodesize=nodesize) #, strata = hv104)
    rf_fit
    rf_fit$type
    NTREE <- rf_fit$ntree
    MTRY  <- rf_fit$mtry

    # calculate an ROC
    # head(rf_fit$votes[,2]) # the HIV+ votes
    # rf_fit.roc <- roc(survhivt_murhoed$hiv03, rf_fit$votes[,2])
    # plot(rf_fit.roc)
    # auc(rf_fit.roc)

    # list important vars sorted
    imp_rf <- rf_fit$importance[,3]; imp_rf <- data.frame(imp_rf[order(-imp_rf)])
    head(imp_rf,15)

    get_se(rf_fit)

    cat("predicting...  ")
    # do validation
    predictions <- predict(rf_fit, newdata = survhivv_murhoed)              # classification
    predi_probs <- predict(rf_fit, newdata = survhivv_murhoed, type="prob") # predicted probabilities
    length(predictions)

    # make an roc dataset from validation set
    roc.data <- data.frame(ref=survhivv_murhoed$hiv03, prob=predi_probs[,2])
    roc.data$survey <- this_survey
    roc.data$sex    <- sex
    write.dta(roc.data, str_c("data/roc_",this_survey,"_",sex,".dta"))
    # head(roc.data)
    predicted.roc <- roc(roc.data$ref, roc.data$prob)
    # plot(predicted.roc)

    # compare to actual
    cm <- confusionMatrix(data = predictions, reference = survhivv_murhoed$hiv03, positive = as.character(levels(survhivv_murhoed$hiv03)[2])) # use confusion matrix to EVALUATE

    # the outcomes
    SE <- cm$byClass["Sensitivity"]
    SP <- cm$byClass["Specificity"]
    PPV <- cm$byClass["Pos Pred Value"]
    NPV <- cm$byClass["Neg Pred Value"]
    AUC <- auc(predicted.roc)

    tn <- cm$table[1,1]
    tp <- cm$table[2,2]
    fn <- cm$table[1,2]
    fp <- cm$table[2,1]

    # ------------------------------------------------------------------------------
    # Run Boruta to get important features
    imprtnt_vars <- NA
    if (bor_ == 1) {
        cat("boruting... \n")
        br_fit <- Boruta(hiv03 ~ . , data = survhivt_murhoed) # , doTrace=2) #, strata=hv104)
        summary(br_fit)
        br_fit
        br_fit$pValue
        br_fit$impSource
        br_fit$mcAdj

        # print out BOruta plot
        pdf(file = str_c("figs/",this_survey,"_",sex,".pdf"), width=14)
            plot(br_fit)
        dev.off()

        # get important vars
        imprtnt_vars <- c(getSelectedAttributes(br_fit, withTentative=FALSE),"hiv03")
    }
                    # --------------------------------------------------------------------------
                    if (1==0) {
                        # Get the selected attributes
                        imprtnt_vars <- c(getSelectedAttributes(br_fit, withTentative=FALSE),"hiv03")
                        survhivt_murhoed_b <- survhivt_murhoed[ ,imprtnt_vars]

                        nvars_b <- ncol(survhivt_murhoed_b)
                        nobs_b  <- round(nrow(survhivt_murhoed_b),0)
                        mtry = nvars_b # round(nvars_b/2,0)
                        ntree = round(nobs_b * nvars_b * 0.000131,0)

                        # RANDOM FOREST WITH ONLY BORUTA-SELECTED VARIABLES
                        modelFit <- NULL
                        modelFit <- randomForest(hiv03 ~ . , data = survhivt_murhoed_b, importance = TRUE,
                                    ntree = ntree, mtry = mtry, nodesize=12)

                        get_se(modelFit)
                        modelFit

                        # do validation
                        predictions <- predict(modelFit, newdata = survhivv_murhoed)
                        length(predictions)
                        # compare to actual
                        cm <- confusionMatrix(data = predictions, reference = survhivv_murhoed$hiv03, positive = "[1]hiv  positive") # use confusion matrix to EVALUATE
                    } # ------------------------------------------------------------------------



    return(
        data.frame(
        PERCENT_W_HIV_DATA = PERCENT_W_HIV_DATA ,
        N_HAS_HIV_DATA = N_HAS_HIV_DATA ,
        N_HIV_POS= N_HIV_POS,
        PREVALENCE = PREVALENCE ,
        N_HAS_HIV_DATA_MURHO = N_HAS_HIV_DATA_MURHO ,
        N_VARS_MURHO = N_VARS_MURHO,
        N_HIV_POS_MURHO= N_HIV_POS_MURHO,
        PREVALENCE_MURHO = PREVALENCE_MURHO ,
        SEX = SEX,
        SE = SE,
        SP = SP,
        PPV = PPV,
        NPV = NPV,
        AUC = AUC,
        tn  = tn,
        tp  = tp,
        fn  = fn,
        fp  = fp,
        NTREE = NTREE,
        MTRY = MTRY,
        BORUTA = imprtnt_vars
        )
    )
} # end do_boruta









# ------------------------------------------------------------------------------
# for each survey, run boruta
# ------------------------------------------------------------------------------
all_results <- NULL
surveys
surveys_se1 <- c("ao71", "mw7a", "et71", "zm63", "ls72")
file <- surveys_se1[4]
# ------------------------------------------------------------------------------
# FOR THE ALL THE MALES/FEMALES ...

cat(date())
for (file in surveys) {
    # get the survey filepath and it's corresponding metadata filepath
    # survey_filepath <- str_c(survey_dir, survey, "_flattened.dta", sep="")
    survey_filepath <- str_c(srvRDS_dir, file, ".rds", sep="")
    mtadat_filepath <- str_c(mtadat_dir, file, "_metadat.dta", sep="")

    # print it out
    message("\nANALYZING ", file, " ------------------------------------------------------------")

    # try BORUTA
    this_result_m <- NULL
    this_result_f <- NULL

    # 1 male, 2 female
    try(   this_result_m <- cbind(do_boruta(survey_filepath, mtadat_filepath, sex=1, this_survey = file, bor_ = 0, seed = 314),data.frame(SURVEY = file)) ) # MALES
    try(   this_result_f <- cbind(do_boruta(survey_filepath, mtadat_filepath, sex=2, this_survey = file, bor_ = 0, seed = 314),data.frame(SURVEY = file)) ) # FEMALES
    cat(date()) # output the date finished
    all_results <- rbind(all_results, this_result_m)
    all_results <- rbind(all_results, this_result_f)
}
all_results_AUC <- all_results
# write to a file -----------------------
write.dta(all_results, "all_results_AUC.dta")

# show all results
head(all_results)
nrow(all_results)
date()


# ------------------------------------------------------------------------------
# post analysis of boruta'd survey important variables
# what is the median (IQR) of the number of important variables?
# ------------------------------------------------------------------------------
#
#
# if (1==0) {
#     # how many features were selected
#     summary(as.numeric(table(all_results$SURVEY)))
#     # table(all_results$BORUTA) # count of variable frequency
#     features <- data.frame(table(all_results$BORUTA)/length(table(all_results$SURVEY, all_results$SEX))) # frequency percent
#     features <- features[order(-1*features$Freq), ] #sorted
#
#     # list important features
#     head(features, 15)
#     features
# }


all_results <- read.dta("output_9_18/all_results_9_18.dta")
# analysis of prediction: Se, Sp, PPV, NPV
surv_results <- all_results[ , names(all_results)[!(names(all_results) %in% c("BORUTA"))] ] # drop BORUTA
surv_results <- unique(surv_results) # and get rid of repeated rows (used to be 1 per important var)
surv_results <- surv_results[order(-1*surv_results$PREVALENCE_MURHO), ] # ordered by prevalence

# take a look at results
head(surv_results)
summary(surv_results)
# summary of only those with meaningful predictions (ie: Sensitivty not 0, or 1)
meaningful_results <- surv_results[!(surv_results$SE %in% c(1,0)),]
meaningful_results$SEX <- factor(meaningful_results$SEX)
summary(meaningful_results)
nrow(meaningful_results)
# how many surveys represented?
length(table(meaningful_results$SURVEY)[table(meaningful_results$SURVEY)!=0])
# CreateTableOne(vars = names(meaningful_results), strata = c("SEX"), data = meaningful_results)

# total population possible
sum(meaningful_results$N_HAS_HIV_DATA)
# total population analyzed
sum(meaningful_results$N_HAS_HIV_DATA_MURHO)
# total population HIV+
sum(meaningful_results$N_HIV_POS_MURHO)

# feature list from Boruta ---------------------------


# read in labels for variables (from zm71)
labels <- read.csv("/Users/echow/QSU/DHS_boruta/data/dictionary.csv", sep=",")
head(labels)
labels$var <- as.character(labels$var)
labels$label <- as.character(labels$label)

# how many features were selected
meaningful_features <- all_results[!(all_results$SE %in% c(0,1)), ]
meaningful_features$SURVEY <- factor(meaningful_features$SURVEY)
# how many features selected per survey?
summary(as.numeric(table(meaningful_features$SURVEY)))

# frequency of features
num_meaningful_surveys <- nrow(unique(meaningful_features[,c("SURVEY", "SEX")]))
num_meaningful_surveys

# FEATURES OVERALL
features <- data.frame(table(meaningful_features$BORUTA)/num_meaningful_surveys) # frequency percent
features <- features[order(-1*features$Freq), ] #sorted
features$var <- as.character(features$Var1)
features <- join(features, labels, by=c("var"), type="left", match="first")
head(features, 15) # list important features


# FEATURES FOR MALES
meaningful_features_m <- meaningful_features[meaningful_features$SEX == 1, ]
num_meaningful_surveys_m <- nrow(unique(meaningful_features_m[,c("SURVEY", "SEX")]))
features_m <- data.frame(table(meaningful_features_m$BORUTA)/num_meaningful_surveys_m) # frequency percent
features_m <- features_m[order(-1*features_m$Freq), ] #sorted
features_m$var <- as.character(features_m$Var1)
features_m <- join(features_m, labels, by=c("var"), type="left", match="first")
head(features_m, 30)
features_m

# FEATURES FOR FEMALES
meaningful_features_f <- meaningful_features[meaningful_features$SEX == 2, ]
num_feaningful_surveys_f <- nrow(unique(meaningful_features_f[,c("SURVEY", "SEX")]))
features_f <- data.frame(table(meaningful_features_f$BORUTA)/num_feaningful_surveys_f) # frequency percent
features_f <- features_f[order(-1*features_f$Freq), ] #sorted
features_f$var <- as.character(features_f$Var1)
features_f <- join(features_f, labels, by=c("var"), type="left", match="first")
head(features_f, 60)
features_f



# ------------------------------------------------------------------------------
# PPV and Se, Sp analysis
# ROC curves
roc_ <- read.dta("data/roc_ALL.dta")
meaning_ <- meaningful_results[,c("SURVEY", "SEX")]
meaning_$survey <- meaning_$SURVEY
meaning_$sex <- meaning_$SEX
meaning_$flag <- 1
meaning_ <- meaning_[,c("survey", "sex", "flag")]

roc <- join(roc_, meaning_, by=c("survey", "sex"), type="left", match="first")
roc <- roc[!is.na(roc$flag), ]
head(roc)
summary(roc)
unique(roc$survey)

# ROC plot
roc_p <- roc(roc$ref, roc$prob)
plot(roc_p, col=0, lty=c(2), main="", grid=c(0.2, 0.2))
legend(x=0.1, y=0.85,unique(roc$survey),pch=15, col=1:length(unique(roc$survey)), bty="n", cex=0.6)
legend("bottomright",c("male","female"),lty=1:2, col="black", bty="n", cex=0.6)
text(0.85,0.80, srt=45, "max AUC=0.87")
text(0.65,0.5, srt=45, "min AUC=0.65")
i = 1
for (survey in unique(roc$survey)) {
    roc_t <- roc[roc$survey == survey & roc$sex == 1, ]
    if (nrow(roc_t) != 0) {
        roc_p <- roc(roc_t$ref, roc_t$prob)
        plot(roc_p, col=i, lty=1, lwd=0.8, add=TRUE)
    }
    i = i+1
}
i = 1
for (survey in unique(roc$survey)) {
    roc_t <- roc[roc$survey == survey & roc$sex == 2, ]
    if (nrow(roc_t) != 0) {
        roc_p <- roc(roc_t$ref, roc_t$prob)
        plot(roc_p, col=i, lty=2, lwd=0.8, add=TRUE)
    }
    i = i+1
}


# get survey list for year
survey_list <- read.delim("/Users/echow/QSU/DHS_boruta/survey_list/survey_list.csv", sep=",")
head(survey_list)
survey_list$SURVEY <- substr(survey_list$filename, 1, 4)

mr <- join(meaningful_results, survey_list[,c("SURVEY","year","country")], by=c("SURVEY"), type="left", match="first")
ar <- join(all_results, survey_list[,c("SURVEY","year")], by=c("SURVEY"), type="left", match="first")
head(mr)
ggplot(mr, aes(year, AUC)) + geom_smooth(method="lm") + geom_point(color="blue")
auc_time <- lm(AUC ~ year, data=mr)
summary(auc_time)
cor(mr$year, mr$AUC)

# get Se/Sp that maximizes the PPV ...

# survey <- unique(roc$survey)[1]

result <- data.frame(survey = NA, sex=NA, se=NA, sp=NA, ppv=NA, thresh=NA, youden=NA)
i <- 1
for (survey in unique(roc$survey)) {
    for (sex in 1:2) {
        roc_t <- roc[roc$survey == survey & roc$sex == sex, ]
        if (nrow(roc_t) != 0) {
            J <- 0          # reset Youden's
            J_max <- 0
            opt_thresh <- 0

            # SOLVE FOR BEST YOUDEN'S
            for (j in 1:500) {
                # thresh <- i/100
                # subset predictions to this survey
                # calculate and apply threshold to probability
                thresh <- j/500 # table(roc_t$ref)[2]/nrow(roc_t) # set the threshold to the prevalence}
                roc_t$class <- factor(as.numeric(roc_t$prob >= thresh))
                # make the same factor levels
                levels(roc_t$class) <- levels(roc_t$ref)
                # get Se, Sp, PPV, NPV
                cm <- confusionMatrix(data = roc_t$class, reference = roc_t$ref, positive = as.character(levels(roc_t$ref)[2])) # use confusion matrix to EVALUATE
                se <- cm$byClass["Sensitivity"]
                sp <- cm$byClass["Specificity"]
                ppv <- cm$byClass["Pos Pred Value"]
                npv <- cm$byClass["Neg Pred Value"]

                # Youden's index
                J <- se + sp - 1
                if (J > J_max) {
                    J_max <- J
                    opt_thresh <- thresh
                }
            } # finish solving for best Youden's

            # --------------------------------------------
            # calcuate the Se/Sp with best Youden's
            thresh <- opt_thresh # table(roc_t$ref)[2]/nrow(roc_t) # set the threshold to the prevalence}
            roc_t$class <- factor(as.numeric(roc_t$prob >= thresh))
            # make the same factor levels
            levels(roc_t$class) <- levels(roc_t$ref)
            # get Se, Sp, PPV, NPV
            cm <- confusionMatrix(data = roc_t$class, reference = roc_t$ref, positive = as.character(levels(roc_t$ref)[2])) # use confusion matrix to EVALUATE
            se <- cm$byClass["Sensitivity"]
            sp <- cm$byClass["Specificity"]
            ppv <- cm$byClass["Pos Pred Value"]
            npv <- cm$byClass["Neg Pred Value"]

            # ---------------------------------
            # set results se; sp; ppv; npv
            result[i, "survey"] <- survey
            result[i, "se"] <- se
            result[i, "sp"] <- sp
            result[i, "ppv"] <- ppv
            result[i, "thresh"] <- thresh
            result[i, "youden"] <- J_max
            result[i, "sex"] <- sex
            i = i + 1
        } # end if nrow
    } # end sex
} # end survey loop

result
abstract_surveys <- meaningful_results[,c("SURVEY", "SEX")]
abstract_surveys$survey <- abstract_surveys$SURVEY
abstract_surveys$sex <- abstract_surveys$SEX

result_abstract <- join(abstract_surveys, result, by=c("survey", "sex"), type="left", match="first")
result_abstract
nrow(result_abstract)
summary(result_abstract)
#
#
#
#
# # *****************************************************
# # cycle through a bunch of Mu and Rhos, and for each
# # apply the restrictions and see how many are lost
#
# # an array to store results
# sumstats <- data.frame(mu=NA, rho=NA, vars=NA, cat_vars=NA, obs=NA, tn=NA, fn=NA, tp=NA, fp=NA, se=NA, sp=NA, ppv=NA, npv=NA)
#
# # index var for storing meta results
# i = 1
# m = 5; r = 2
# set.seed(314)
# # takes about ... um ... some time to run! 5:26 -
# # loop through Mu and Rho's
# for (m in 1:9) {
#     for (r in 1:4) {
#         mu = m/20; rho = r/10
#         d    <- NULL
#         dtrn <- NULL
#         dtst <- NULL
#         modelFit <- NULL
#         predictions <-NULL
#         cm <- NULL
#
#         try ({
#             # apply mu and rho -------------------------------------
#
#             # append train and test set so can murho (and get same levels)
#             # but then split again
#             survhivt$train <- 1
#             survhivv$train <- 0
#
#             # recombine the train and test sets
#             survhiva <- rbind(survhivt, survhivv)
#             # apply murho to both train and test
#             da <- murho(data = survhiva, mu = mu, rho = rho)
#             # resplit into train and test
#             d  <- da[da$train==1, ]
#             dv <- da[da$train==0, ]
#
#             # count factors
#             factor_i = 0
#             for (var in names(d)) {
#                 if (is.factor(d[,var])) {
#                     factor_i <- factor_i + 1
#                 }
#             }
#             # do I sample within the dataset after mu/rho, or use survhivv, the validation set?
#
#             # YES, do THIS one! use the validation data that was set aside prior to murho
#             if (do_OOS_valid) { # use OOS validation, the CORRECT way
#                 # don't partition the data -----------------------------------
#                 # use whole murho's dataset to train
#                 # train <- createDataPartition(y = d$hiv03, p=0.80, list=FALSE)
#                 # subset using partition
#                 dtrn  <- d      # used for fitting the model, and Mu and Rho etc
#                 dtst  <- dv  # THE VALIDATION SET,  survhivv (validation) was murho'd
#
#                 # build a predictor ------------------------------------
#                 message("running randomForest ...")
#                 modelFit <- randomForest(hiv03 ~ . , data = dtrn, importance = TRUE)  # will auto center and scale
#
#                 # do OUT OF sample validation --------------------------
#                 message("RF complete! predicting from test set ...")
#                 predictions <- predict(modelFit, newdata=dtst)
#                 # compare to actual
#                 cm <- confusionMatrix(predictions, dtst$hiv03) # use confusion matrix to EVALUATE
#
#                 # store meta results
#                 sumstats[i,"mu"] <- mu
#                 sumstats[i,"rho"] <- rho
#                 sumstats[i,"obs"] <- dim(d)[1]
#                 sumstats[i,"vars"] <- dim(d)[2]
#                 sumstats[i,"cat_vars"] <- factor_i
#
#                 sumstats[i,"tn"] <- cm$table[1,1]
#                 sumstats[i,"fn"] <- cm$table[1,2]
#                 sumstats[i,"tp"] <- cm$table[2,2]
#                 sumstats[i,"fp"] <- cm$table[2,1]
#
#                 sumstats[i,"se"] <- cm$byClass["Sensitivity"]
#                 sumstats[i,"sp"] <- cm$byClass["Specificity"]
#                 sumstats[i,"ppv"] <- cm$byClass["Pos Pred Value"]
#                 sumstats[i,"npv"] <- cm$byClass["Neg Pred Value"]
#                 message("DONE!\n")
#             }
#             i = i + 1
#         })
#     }
# }
#
# # show results
# sumstats
#
# # ------------------------------------------------------------------------------
# # plot mu vs. # vars
# pdf("figs/mu_vars_oos.pdf", height=4, width=4.7)
#     ggplot(sumstats, aes(x=mu, y=vars)) + geom_point(aes(size=obs))
# dev.off()
#
# # plot mu vs. Sensitivity
# pdf("figs/mu_se_oos.pdf", height=4, width=4.7)
#     ggplot(sumstats[1:32,], aes(x=mu, y=se)) + geom_jitter(width=0.005, height=0.001) + geom_smooth(method="loess")
# dev.off()
#
#
#
#
#
#
#
#
#
#
# #      ~ fin ~
