#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

library(nlme)
library(nlshelper)

I <- args[1]

## read data
D <- read.table(paste0(I, ".txt"), header = T)
D$Group <- as.factor(D$Group)

## set parameter
HW.st <- c(C = 0.00001)

sink(file = paste0(I, "-stats.txt"))

## calculate with grouping factor
HW.nonlinear_full <- nlsList(LD.data ~ ((10 + C * distance) / ((2 + C * distance) * (11 + C * distance))) * (1 + ((3 + C * distance) * (12 + 12 * C * distance + (C * distance)^2)) / (n * (2 + C * distance) * (11 + C * distance))) | Group, control = nls.control(maxiter = 1000), start = HW.st, data = D)

cat("**********full model**********")
summary(HW.nonlinear_full)

## calculate without grouping factor
HW.nonlinear_red <- nls(LD.data ~ ((10 + C * distance) / ((2 + C * distance) * (11 + C * distance))) * (1 + ((3 + C * distance) * (12 + 12 * C * distance + (C * distance)^2)) / (n * (2 + C * distance) * (11 + C * distance))), start = HW.st, control = nls.control(maxiter = 1000, warnOnly = TRUE), data = D)

cat("*********reduced model**********")
summary(HW.nonlinear_red)

## compare models

cat("**********ANALYSIS OF DEVIANCE**************")
anova_nlslist(HW.nonlinear_full, HW.nonlinear_red)

## comparsion of coefficients
coef.vmax <- coef(summary(HW.nonlinear_full))[, , 1][1:2]
coef.vmax.se <- coef(summary(HW.nonlinear_full))[, , 1][3:4]

t.wert.IS <- as.numeric((coef.vmax[1] - coef.vmax[2]) / sqrt((coef.vmax.se[1] + coef.vmax.se[2])^2))
p.IS <- pt(t.wert.IS, nrow(D))

cat("********* comparsison of coefficients **********")
print(p.IS)

sink()
