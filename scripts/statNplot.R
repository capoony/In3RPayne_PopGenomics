#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

### Input/Output

O <- args[1]

library(ggplot2)
library(plyr)
library(emmeans)

D <- read.table(O, header = T)

Dpairs <- subset(D, D$Origin != "Victoria")
Dpairs <- subset(Dpairs, Dpairs$Origin != "Sweden")
Dpairs <- subset(Dpairs, Dpairs$Origin != "Maine")
Dpairs <- subset(Dpairs, Dpairs$InvStatus != "All")
Dpairs$ST <- Dpairs[, ncol(Dpairs)]


## remove sites which are not found in all datasets
Dpairs <- data.frame(Dpairs, "CP" = paste(Dpairs$C, Dpairs$P, sep = ""))
Dpos <- ddply(Dpairs, c("CP", "location"), summarise,
  N = length(ST)
)
missing <- subset(Dpos, Dpos$N != 8)$CP
Dpairs <- Dpairs[!(Dpairs$CP %in% missing), ]

## randomly subsample sites outside
inside <- subset(Dpairs, Dpairs$location == "inside")
outside <- unique(subset(Dpairs, Dpairs$location == "outside")$CP)
I <- length(inside$CP) / 8
Keep <- sample(outside, I)
outside <- Dpairs[Dpairs$CP %in% Keep, ]
Dpairs <- rbind(inside, outside)


sink(paste(O, ".anova.txt", sep = ""))
res <- aov(ST ~ Origin * InvStatus * location, data = Dpairs)
summary(res)
pairs(emmeans(res, c("Origin", "InvStatus", "location")))
TukeyHSD(res)
sink()


Dsum <- ddply(Dpairs, c("Origin", "InvStatus", "location"), summarise,
  N    = length(ST),
  mean = mean(ST),
  sd   = sd(ST),
  se   = sd / sqrt(N)
)

pdf(paste(O, ".pdf", sep = ""), width = 8, height = 5)

ggplot(data = Dsum, aes(x = Origin, y = mean, fill = InvStatus)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_manual(values = c("grey80", "grey20")) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), position = position_dodge(.9), width = .1) +
  facet_wrap(~location) +
  theme_classic()

dev.off()

sink(paste(O, ".Wilcox.txt", sep = ""))

wilcoxtest <- by(Dpairs, list(Dpairs$Origin, Dpairs$location), function(x) {
  wilcox.test(ST ~ InvStatus, data = x, alternative = "less")
})

wilcoxtest

cat("Now test for Differences between regions inside and outside of inversion")

wilcoxtest <- by(Dpairs, list(Dpairs$Origin), function(x) {
  wilcox.test(ST ~ location, data = x, alternative = "less")
})

wilcoxtest


sink()
