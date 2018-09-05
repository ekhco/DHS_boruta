# look at random forest tuning PARAMETERS
# how many variables should I consider
# and how many trees should I grow?
rm(list=ls())
gc()
library(ggplot2)

setwd("~/QSU/DHS_ranforest/dhs_hiv")
tuning <- read.csv("~/QSU/DHS_ranforest/dhs_hiv/data/tuning.csv")

# scaled tuning
tuning$se_s <- (tuning$se - mean(tuning$se))/sd(tuning$se)
head(tuning)
tuning <- tuning[,c(1,2,4,5)]


# tuning_heat <- reshape(tuning, idvar="mtry", timevar="ntrees", v.names="vars", direction = "wide", sep="_")
# tuning_heat <- as.matrix(tuning_heat)
# head(tuning_heat)
#
# if (1==0){
#     library(reshape)
#     tuning_h <- cast(tuning, mtry ~ ntrees)
#     tuning_hm <- as.matrix(tuning_h)
#     heatmap(tuning_hm) # but does dendogram
#
#     # using plotly (makes a webpage?)
#     library(plotly)
#     p <- plot_ly(z=tuning_hm, type="heatmap")
# }



# use ggplot
head(tuning)
p <- ggplot(tuning, aes(mtry, ntree)) +
    geom_tile(aes(fill=se)) +
    scale_fill_gradient(low="white", high="tomato") +
    xlim(0,250) +
    labs(x="number of variables") +
    labs(y="number of trees in forest")
p

# as a scatter plot (N VARS)
s <- ggplot(tuning, aes(mtry,se, col=ntree)) +
    geom_point() +
    geom_smooth() +
    xlim(0,430) +
    ylim(0,0.3) +
    labs(x="number of variables") +
    labs(y="sensitivity")
s

# as a scatter plot (NTREES)
nt <- ggplot(tuning, aes(ntree,se)) +
    geom_point() +
    geom_smooth() +
    xlim(0,250) +
    labs(x="number of trees") +
    labs(y="sensitivity")
nt
