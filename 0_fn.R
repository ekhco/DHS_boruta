
# *****************************************************
# A function that identifies the degree of missinginess
# in a variable.
missingness <- function(var) {
    if (is.character(var)) { # for characters
        count_miss <- sum(as.numeric(var == ""))
    }
    else { # for numerics
        count_miss <- sum(as.numeric(is.na(var)))
    }
    # count_blnk <- sum(as.numeric(zw[,var] == ""))
    count_ttl <- length(var)
    return(count_miss/count_ttl)
}

# returns a boolean vector of wether string is in each element of da
isin <- function(string, da){
  !is.na(str_locate(da, string)[,1])
}




# ------------------------------------------------------------------------------
# DO RANDOM FOREST PREDICTION w MU & RHO
# 2018-05-21 discussion w Eran:
# don't do imputation, stratify by gender, use hyper-parameters, Mu (missingness
# of variable) and Rho (for label). Keep only complete cases of vars selected
# for Mu. Select an optimal Mu, Rho. Report a variable importance list.
#
# Mu:  variable kept if missingness < Mu
# Rho: variable is cont if lbls/lvls < Rho

# *****************************************************
# given missinginess, and cat/cont threshold Rho, will
# create a subsetted dataset and return it.
murho <- function(data, mu = 0.2, rho = 0.2, meta) {
    # ****************************************
    # Set metaparameters:                    *
    Mu  <- mu  # missingness of var          *
    Rho <- rho  # whether var is cont/cat    *
    # ****************************************

    # -------------------------------------------
    # stratify random forest by gender
    table(data$hv104)
    # split by gender
    # if (1==0) {zw1 <- data[data$hv104 == 1, ]}   # 48% are sex == 1 (ie: MALE)
    zw1 <- data[!is.na(data$hv104), ]   # acutally, don't subset my gender
    # zw2 <- data[data$hv104 == 2, ]   # 52% are sex == 2 (ie: FEMALE)

    # -------------------------------------------
    # get a list of vars w missingness < Mu
    names.zw1 <- names(zw1)
    var.out     <- names.zw1
    length(var.out)
    for (var in names.zw1) {
        # print(var)
        if (missingness(zw1[,var]) >= Mu) {
            var.out <- var.out[!(var.out %in% c(var))] # remove variable if has missingness >= Mu
        }
    }
    # subset data to include only vars with Missinessness < Mu
    zw1m <- zw1[,var.out]

    # what percent of variables were dropped?
    if (ncol(zw1) != ncol(zw1m)) {cat( ncol(zw1) - ncol(zw1m)," (", round(100*(ncol(zw1)-ncol(zw1m))/ncol(zw1)), "%) of variables had >",Mu," missingness and were dropped.\n", sep="")}

    # -------------------------------------------
    # drop any vars that are 1 value
    nrow(zw1m)
    var.out <- names(zw1m)
    for (var in names(zw1m)) {
        # if >95% of levels are uniquely identifying, drop
        if ((length(unique(zw1m[,var]))) == 1) {
            var.out <- var.out[!(var.out %in% c(var))]
        }
    }
    # subset data to exclude vars with 1 level
    zw1m_ <- zw1m[, var.out]
    # what percent of variables were dropped for having 1 level)
    if (ncol(zw1m) != ncol(zw1m_)) {cat( ncol(zw1m) - ncol(zw1m_)," (", round(100*(ncol(zw1m)-ncol(zw1m_))/ncol(zw1m)), "%) variables had 1 level and were dropped.\n", sep="")}
    # summary(zw1m_)



    # -------------------------------------------
    # drop variables that are character
    var.out <- names(zw1m_)
    for (var in names(zw1m_)) {
        # if >95% of levels are uniquely identifying, drop
        if (is.character(zw1m_[,var])) {
            var.out <- var.out[!(var.out %in% c(var))]
        }
    }
    zw1m_c <- zw1m_[, var.out]
    # what percent of vars dropped because character?
    if (ncol(zw1m_c) != ncol(zw1m_)) {cat( ncol(zw1m_) - ncol(zw1m_c)," (", round(100*(ncol(zw1m_)-ncol(zw1m_c))/ncol(zw1m_)), "%) variables were character and dropped.\n", sep="")}



    # -------------------------------------------
    # determine if variable is cat or continuous
    zw1m_cc <- zw1m_c
    # keep vars from metadat that are in zw1m_
    insubset <- meta$var_name %in% names(zw1m_cc)
    # subset metadat to these vars
    meta_ <- meta[insubset, ]
    # get name of vars with Rho < Rho
    meta_cat  <- meta_[meta_$rho >= Rho, "var_name" ]  # categorical
    meta_cont <- meta_[meta_$rho <  Rho, "var_name" ]  # continuous

    # turn vars into factors and recode missingness
    # as a new category
    for (var in meta_cat) {
        zw1m_cc[,var] <- as.factor(zw1m_cc[,var]) # to factors
        # if there are missing categories...
        if (sum(as.numeric(is.na(zw1m_cc[,var]))) > 0) { # there are missing levels
            # recode NA to new level
            levels(zw1m_cc[,var])  <- c(levels(zw1m_cc[,var]), "missing")
            zw1m_cc[is.na(zw1m_cc[,var]), var] <- "missing"
        }
        # drop it if there are more than 50 levels because ran forest can't even
        if (length(levels(zw1m_cc[,var])) > 50 ) {
            # zw1m_cc <- subset(zw1m_cc, select = -c(var))
            zw1m_cc <- zw1m_cc[ , names(zw1m_cc)[!(names(zw1m_cc) %in% c(var))] ]
        }
    }
    # turn vars into numeric
    for (var in meta_cont) {
        zw1m_cc[,var] <- as.numeric(zw1m_cc[,var])
    }
    # summary(zw1m_cc)

    # -------------------------------------------
    # keep only complete cases within these vars
    # obs.keep <- 1:nrow(zw1m_cc)
    # for (i in obs.keep) {
    #     if ( missingness(zw1m_cc[i,]) > 0 ) { # not complete case
    #         obs.keep <- obs.keep[!(obs.keep %in% i)]
    #     }
    # }
    obs.keep <- complete.cases(zw1m_cc)

    # subset data to only complete case
    zw1m_ccc <- zw1m_cc[obs.keep,]
    # output
    if (nrow(zw1m_ccc) != nrow(zw1m_cc)) {cat( nrow(zw1m_cc) - nrow(zw1m_ccc)," (", round(100*(nrow(zw1m_cc)-nrow(zw1m_ccc))/nrow(zw1m_cc)), "%) observations were incomplete cases and were dropped.\n", sep="")}
    d <- zw1m_ccc # put it into a short var name

    cat(dim(d)[1], " (", round(100*dim(d)[1]/dim(data)[1]), "%) obs and ", dim(d)[2], " (", round(100*dim(d)[2]/dim(data)[2]), "%) vars remain after mu=", Mu, " and Rho=", Rho, " applied.\n", sep="")
    return(d)
} # end murho
