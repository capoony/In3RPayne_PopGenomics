#!/usr/bin/env Rscript

library(ggplot2)
library(plyr)
library(gtools)

args <- commandArgs(trailingOnly = TRUE)

### Input/Output

O <- args[1]

D <- read.table(O, header = T)

# convert negative FST values to zero
r2z <- function(x) {
  ifelse(x < 0, 0.0, x)
}
## remove Africa
D <- subset(D, D$Origin != "Africa")
## remove sites which are not found in all datasets
D <- data.frame(D, "CP" = paste(D$C, D$P, sep = ""))
Dpos <- ddply(D, c("CP", "location"), summarise,
  N    = length(fst),
  M = mean(r2z(fst))
)
missing <- subset(Dpos, Dpos$N != 9 & M == 0.0)$CP
D <- D[!(D$CP %in% missing), ]

## randomly subsample sites outside
inside <- subset(D, D$location == "inside")
outside <- unique(subset(D, D$location == "outside")$CP)
I <- length(inside$CP) / 9
Keep <- sample(outside, I)
Drop <- outside[!(outside %in% Keep)]
D <- D[!(D$CP %in% Drop), ]


Dsum <- ddply(D, c("Origin", "Comparison", "location"), summarise,
  N    = length(fst),
  mean = mean(r2z(fst)),
  sd   = sd(r2z(fst)),
  se   = sd / sqrt(N)
)

ggplot(data = Dsum, aes(x = Comparison, y = mean, fill = location)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_manual(values = c("grey80", "grey20")) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), position = position_dodge(.9), width = .1) +
  facet_wrap(~Origin) +
  theme_classic()

OUT <- paste0(O, ".pdf")
ggsave(OUT)

sink(paste0(O, ".anova.txt"))

wilcoxtest <- by(D, list(D$Origin, D$location), function(x) {
  res <- aov(r2z(fst) ~ Comparison, data = x)
  anova(res)
})

wilcoxtest

wilcoxtest <- by(D, list(D$Origin, D$location), function(x) {
  res <- aov(r2z(fst) ~ Comparison, data = x)
  summary(res)
  TukeyHSD(res)
}, simplify = F)

wilcoxtest

sink()

sink(paste0(O, ".ANOVA-2way.txt"))

wilcoxtest <- by(D, list(D$Origin), function(x) {
  res <- aov(r2z(fst) ~ Comparison + location, data = x)
  anova(res)
})

wilcoxtest

wilcoxtest <- by(D, list(D$Origin), function(x) {
  res <- aov(r2z(fst) ~ Comparison + location, data = x)
  summary(res)
  TukeyHSD(res)
}, simplify = F)

wilcoxtest

sink()


sink(paste0(O, ".ANOVA-3way.txt"))
res <- aov(r2z(fst) ~ Comparison + Origin + location, data = D)
summary(res)
TukeyHSD(res)
sink()

## compare correlation among windows

D <- read.table(O, header = T)
D <- data.frame(D, "CP" = paste(D$C, D$P, sep = ""))
X <- subset(D, D$C == "3R" & D$P > 14432209 & D$P < 26744010)[, 4:ncol(D) - 1]
Y <- na.omit(subset(D, D$C == "3L" & D$P > 4432209 & D$P < 16744010)[, 4:ncol(D) - 1])

pdf(paste(O, ".comparison.pdf", sep = ""), width = 8, heigh = 5)

XY <- rbind(X, Y)

pairs(XY,
  lower.panel = function(x, y, ...) {
    Xx <- x[seq_len(nrow(X))] # corresponds to X subset
    Xy <- y[seq_len(nrow(X))] # corresponds to X subset
    usr <- par("usr")
    on.exit(par(usr))
    par(usr = c(range(X[, -ncol(X)]), range(X[, -1]))) # set up limits
    points(Xx, Xy, col = rgb(0, 0, 1, 0.3), pch = 16, cex = 2)
    reg <- lm(Xy ~ Xx)
    abline(reg, lwd = 3, lty = 2)
    r2 <- round(summary(reg)$r.squared, 2)
    pval <- stars.pval(summary(reg)$coefficients[[8]])
    legend("topleft", legend = substitute(paste(italic(R)^2, "= ", r2, pval, sep = "")), bg = "white")
    if (par("mfg")[2] == 1) axis(2) # if left plot, add left axis
    if (par("mfg")[1] == ncol(X)) axis(1)
  }, # if bottom plot add bottom axis
  upper.panel = function(x, y, ...) {
    Yx <- x[(nrow(X) + 1):length(x)] # Y subset
    Yy <- y[(nrow(X) + 1):length(y)] # Y subset
    usr <- par("usr")
    on.exit(par(usr))
    par(usr = c(range(Y[, -1]), range(Y[, -ncol(Y)]))) # set up limits
    points(Yx, Yy, col = rgb(1, 0, 0, 0.3), pch = 16, cex = 2)
    reg <- lm(Yy ~ Yx)
    abline(reg, lwd = 3, lty = 2)
    r2 <- round(summary(reg)$r.squared, 2)
    pval <- stars.pval(summary(reg)$coefficients[[8]])
    legend("topleft", legend = substitute(paste(italic(R)^2, "= ", r2, pval, sep = "")), bg = "white")
    if (par("mfg")[2] == ncol(Y)) axis(4) # if right plot, add right axis
    if (par("mfg")[1] == 1) axis(3) # if top plot, add top axis
  },
  tick = FALSE, # suppress the default tick marks
  line = 10
)

dev.off()
