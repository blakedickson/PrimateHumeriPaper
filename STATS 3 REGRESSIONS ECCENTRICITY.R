

rm(list = ls())

library(lmodel2)
load("boneVars.Rdata")
eco <- as.data.frame(read.csv("./DATA/specimenEco2.csv"))


# MAX THICKNESS -----------------------------------------------------------
# using reduced major axis

dat <- data.frame(boneVars$ecoDF, boneVars$maxDF, boneVars$meanDF, mass = eco$MASS..KG.)



# ESSENTRICITY ----------------------------------------------------------------

dat <- cbind(boneVars$ecoDF, boneVars$abs50DF, mass = eco$MASS..KG.)

sink("./stats/essentricity regression stats.txt")
# ### var.CURVATURE

LM_massCurvature <- lmodel2((abs.ecc50.Curvature)~log(mass), "interval", "interval",data = dat, 99)
LM_massCurvature

LM_LengthCurvature <- lmodel2((abs.ecc50.Curvature)~log(ZLENGTH), "interval", "interval",data = dat, 99)
LM_LengthCurvature

LM_CircCurvature <- lmodel2((abs.ecc50.Curvature)~log(MIDSHAFTCIRC), "interval", "interval",data = dat, 99)
LM_CircCurvature

# ### abs.ecc50.CORTICAL AREA

LM_massCA <- lmodel2((abs.ecc50.CA)~log(mass), "interval", "interval",data = dat, 99)
LM_massCA
LM_LengthCA <- lmodel2((abs.ecc50.CA)~log(ZLENGTH), "interval", "interval",data = dat, 99)
LM_LengthCA

LM_CircCA <- lmodel2((abs.ecc50.CA)~log(MIDSHAFTCIRC), "interval", "interval",data = dat, 99)
LM_CircCA

# ### abs.ecc50.THICKNESS

LM_massthickness <- lmodel2((abs.ecc50.thickness)~log(mass), "interval", "interval",data = dat, 99)
LM_massthickness

LM_Lengththickness <- lmodel2((abs.ecc50.thickness)~log(ZLENGTH), "interval", "interval",data = dat, 99)
LM_Lengththickness

LM_Circthickness <- lmodel2((abs.ecc50.thickness)~log(MIDSHAFTCIRC), "interval", "interval",data = dat, 99)
LM_Circthickness

# ### abs.ecc50.Ix

LM_massIx <- lmodel2((abs.ecc50.Ix)~log(mass), "interval", "interval",data = dat, 99)
LM_massIx

LM_LengthIx <- lmodel2((abs.ecc50.Ix)~log(ZLENGTH), "interval", "interval",data = dat, 99)
LM_LengthIx

LM_CircIx <- lmodel2((abs.ecc50.Ix)~log(MIDSHAFTCIRC), "interval", "interval",data = dat, 99)
LM_CircIx

# ### abs.ecc50.Iy

LM_massIy <- lmodel2((abs.ecc50.Iy)~log(mass), "interval", "interval",data = dat, 99)
LM_massIy

LM_LengthIy <- lmodel2((abs.ecc50.Iy)~log(ZLENGTH), "interval", "interval",data = dat, 99)
LM_LengthIy

LM_CircIy <- lmodel2((abs.ecc50.Iy)~log(MIDSHAFTCIRC), "interval", "interval",data = dat, 99)
LM_CircIy

# ### abs.ecc50.I.ratio

LM_massI.ratio <- lmodel2((abs.ecc50.I.ratio)~log(mass), "interval", "interval",data = dat, 99)
LM_massI.ratio

LM_LengthI.ratio <- lmodel2((abs.ecc50.I.ratio)~log(ZLENGTH), "interval", "interval",data = dat, 99)
LM_LengthI.ratio

LM_CircI.ratio <- lmodel2((abs.ecc50.I.ratio)~log(MIDSHAFTCIRC), "interval", "interval",data = dat, 99)
LM_CircI.ratio

# ### abs.ecc50.Zpol

LM_massZpol <- lmodel2((abs.ecc50.Zpol)~log(mass), "interval", "interval",data = dat, 99)
LM_massZpol

LM_LengthZpol <- lmodel2((abs.ecc50.Zpol)~log(ZLENGTH), "interval", "interval",data = dat, 99)
LM_LengthZpol

LM_CircZpol <- lmodel2((abs.ecc50.Zpol)~log(MIDSHAFTCIRC), "interval", "interval",data = dat, 99)
LM_CircZpol

# ### abs.ecc50.J



LM_massJ <- lmodel2((abs.ecc50.J)~log(mass), "interval", "interval",data = dat, 99)
LM_massJ


LM_LengthJ <- lmodel2((abs.ecc50.J)~log(ZLENGTH), "interval", "interval",data = dat, 99)
LM_LengthJ

LM_CircJ <- lmodel2((abs.ecc50.J)~log(MIDSHAFTCIRC), "interval", "interval",data = dat, 99)
LM_CircJ

sink()

######################
# plot all together
#####################



# plot together -----------------------------------------------------------

par(mfrow = c(6,4), mar = c(5,4,1,2))

# using reduced major axis


plot(LM_massCurvature, "RMA", ylab = "Curvature essentricity", xlab = "body mass (log)", main = "")
legend("bottomleft", 
       legend = c(
         paste("intercept =", round(LM_massCurvature$regression.results[4,2], digits = 2)),
         paste("slope =", round(LM_massCurvature$regression.results[4,3], digits = 2),
               "±", round(LM_massCurvature$confidence.intervals[4,5] -LM_massCurvature$confidence.intervals[4,4], 2)),
         paste("p =", LM_massCurvature$regression.results$`P-perm (1-tailed)`[4]),
         paste("Rsq =", round(LM_massCurvature$rsquare, 2))
       ),
       bty= "n"
)

# ### CORTICAL AREA


plot(LM_massCA, "RMA", ylab = "Cortical Area essentricity", xlab = "body mass (log)", main = "")
legend("bottomleft", 
       legend = c(
         paste("intercept =", round(LM_massCA$regression.results[4,2], digits = 2)),
         paste("slope =", round(LM_massCA$regression.results[4,3], digits = 2),
               "±", round(LM_massCA$confidence.intervals[4,5] -LM_massCA$confidence.intervals[4,4], 2)),
         paste("p =", LM_massCA$regression.results$`P-perm (1-tailed)`[4]),
         paste("Rsq =", round(LM_massCA$rsquare, 2))
       ),
       bty= "n"
)


# ### THICKNESS

plot(LM_massthickness, "RMA", ylab = "mean Cortical Thickness essentricity", xlab = "body mass (log)", main = "")
legend("bottomleft", 
       legend = c(
         paste("intercept =", round(LM_massthickness$regression.results[4,2], digits = 2)),
         paste("slope =", round(LM_massthickness$regression.results[4,3], digits = 2),
               "±", round(LM_massthickness$confidence.intervals[4,5] -LM_massthickness$confidence.intervals[4,4], 2)),
         paste("p =", LM_massthickness$regression.results$`P-perm (1-tailed)`[4]),
         paste("Rsq =", round(LM_massthickness$rsquare, 2))
       ),
       bty= "n"
)

# ### Ix

plot(LM_massIx, "RMA", ylab = "Ix essentricity", xlab = "body mass (log)", main = "")
legend("bottomleft", 
       legend = c(
         paste("intercept =", round(LM_massIx$regression.results[4,2], digits = 2)),
         paste("slope =", round(LM_massIx$regression.results[4,3], digits = 2),
               "±", round(LM_massIx$confidence.intervals[4,5] -LM_massIx$confidence.intervals[4,4], 2)),
         paste("p =", LM_massIx$regression.results$`P-perm (1-tailed)`[4]),
         paste("Rsq =", round(LM_massIx$rsquare, 2))
       ),
       bty= "n"
)

# ### Iy

plot(LM_massIy, "RMA", ylab = "Iy essentricity", xlab = "body mass (log)", main = "")
legend("bottomleft", 
       legend = c(
         paste("intercept =", round(LM_massIy$regression.results[4,2], digits = 2)),
         paste("slope =", round(LM_massIy$regression.results[4,3], digits = 2),
               "±", round(LM_massIy$confidence.intervals[4,5] -LM_massIy$confidence.intervals[4,4], 2)),
         paste("p =", LM_massIy$regression.results$`P-perm (1-tailed)`[4]),
         paste("Rsq =", round(LM_massIy$rsquare, 2))
       ),
       bty= "n"
)

# ### I.ratio

plot(LM_massI.ratio, "RMA", ylab = "I.ratio essentricity", xlab = "body mass (log)", main = "")
legend("bottomleft", 
       legend = c(
         paste("intercept =", round(LM_massI.ratio$regression.results[4,2], digits = 2)),
         paste("slope =", round(LM_massI.ratio$regression.results[4,3], digits = 2),
               "±", round(LM_massI.ratio$confidence.intervals[4,5] -LM_massI.ratio$confidence.intervals[4,4], 2)),
         paste("p =", LM_massI.ratio$regression.results$`P-perm (1-tailed)`[4]),
         paste("Rsq =", round(LM_massI.ratio$rsquare, 2))
       ),
       bty= "n"
)


# ### Zpol

plot(LM_massZpol, "RMA", ylab = "Zpol essentricity", xlab = "body mass (log)", main = "")
legend("bottomleft", 
       legend = c(
         paste("intercept =", round(LM_massZpol$regression.results[4,2], digits = 2)),
         paste("slope =", round(LM_massZpol$regression.results[4,3], digits = 2),
               "±", round(LM_massZpol$confidence.intervals[4,5] -LM_massZpol$confidence.intervals[4,4], 2)),
         paste("p =", LM_massZpol$regression.results$`P-perm (1-tailed)`[4]),
         paste("Rsq =", round(LM_massZpol$rsquare, 2))
       ),
       bty= "n"
)


# ### abs.ecc50.J

plot(LM_massJ, "RMA", ylab = "J essentricity", xlab = "body mass (log)", main = "")
legend("bottomleft", 
       legend = c(
         paste("intercept =", round(LM_massJ$regression.results[4,2], digits = 2)),
         paste("slope =", round(LM_massJ$regression.results[4,3], digits = 2),
               "±", round(LM_massJ$confidence.intervals[4,5] -LM_massJ$confidence.intervals[4,4], 2)),
         paste("p =", LM_massJ$regression.results$`P-perm (1-tailed)`[4]),
         paste("Rsq =", round(LM_massJ$rsquare, 2))
       ),
       bty= "n"
)



# plot together -----------------------------------------------------------

# par(mfrow = c(4,4), mar = c(5,4,1,2))

# using reduced major axis


plot(LM_LengthCurvature, "RMA", ylab = "Curvature essentricity", xlab = "bone length (log)", main = "")
legend("bottomleft", 
       legend = c(
         paste("intercept =", round(LM_LengthCurvature$regression.results[4,2], digits = 2)),
         paste("slope =", round(LM_LengthCurvature$regression.results[4,3], digits = 2),
               "±", round(LM_LengthCurvature$confidence.intervals[4,5] -LM_LengthCurvature$confidence.intervals[4,4], 2)),
         paste("p =", LM_LengthCurvature$regression.results$`P-perm (1-tailed)`[4]),
         paste("Rsq =", round(LM_LengthCurvature$rsquare, 2))
       ),
       bty= "n"
)

# ### CORTICAL AREA


plot(LM_LengthCA, "RMA", ylab = "Cortical Area essentricity", xlab = "bone length (log)", main = "")
legend("bottomleft", 
       legend = c(
         paste("intercept =", round(LM_LengthCA$regression.results[4,2], digits = 2)),
         paste("slope =", round(LM_LengthCA$regression.results[4,3], digits = 2),
               "±", round(LM_LengthCA$confidence.intervals[4,5] -LM_LengthCA$confidence.intervals[4,4], 2)),
         paste("p =", LM_LengthCA$regression.results$`P-perm (1-tailed)`[4]),
         paste("Rsq =", round(LM_LengthCA$rsquare, 2))
       ),
       bty= "n"
)


# ### THICKNESS

plot(LM_Lengththickness, "RMA", ylab = "mean Cortical Thickness essentricity", xlab = "bone length (log)", main = "")
legend("bottomleft", 
       legend = c(
         paste("intercept =", round(LM_Lengththickness$regression.results[4,2], digits = 2)),
         paste("slope =", round(LM_Lengththickness$regression.results[4,3], digits = 2),
               "±", round(LM_Lengththickness$confidence.intervals[4,5] -LM_Lengththickness$confidence.intervals[4,4], 2)),
         paste("p =", LM_Lengththickness$regression.results$`P-perm (1-tailed)`[4]),
         paste("Rsq =", round(LM_Lengththickness$rsquare, 2))
       ),
       bty= "n"
)

# ### Ix

plot(LM_LengthIx, "RMA", ylab = "Ix essentricity", xlab = "bone length (log)", main = "")
legend("bottomleft", 
       legend = c(
         paste("intercept =", round(LM_LengthIx$regression.results[4,2], digits = 2)),
         paste("slope =", round(LM_LengthIx$regression.results[4,3], digits = 2),
               "±", round(LM_LengthIx$confidence.intervals[4,5] -LM_LengthIx$confidence.intervals[4,4], 2)),
         paste("p =", LM_LengthIx$regression.results$`P-perm (1-tailed)`[4]),
         paste("Rsq =", round(LM_LengthIx$rsquare, 2))
       ),
       bty= "n"
)

# ### Iy

plot(LM_LengthIy, "RMA", ylab = "Iy essentricity", xlab = "bone length (log)", main = "")
legend("bottomleft", 
       legend = c(
         paste("intercept =", round(LM_LengthIy$regression.results[4,2], digits = 2)),
         paste("slope =", round(LM_LengthIy$regression.results[4,3], digits = 2),
               "±", round(LM_LengthIy$confidence.intervals[4,5] -LM_LengthIy$confidence.intervals[4,4], 2)),
         paste("p =", LM_LengthIy$regression.results$`P-perm (1-tailed)`[4]),
         paste("Rsq =", round(LM_LengthIy$rsquare, 2))
       ),
       bty= "n"
)

# ### I.ratio

plot(LM_LengthI.ratio, "RMA", ylab = "I.ratio essentricity", xlab = "bone length (log)", main = "")
legend("bottomleft", 
       legend = c(
         paste("intercept =", round(LM_LengthI.ratio$regression.results[4,2], digits = 2)),
         paste("slope =", round(LM_LengthI.ratio$regression.results[4,3], digits = 2),
               "±", round(LM_LengthI.ratio$confidence.intervals[4,5] -LM_LengthI.ratio$confidence.intervals[4,4], 2)),
         paste("p =", LM_LengthI.ratio$regression.results$`P-perm (1-tailed)`[4]),
         paste("Rsq =", round(LM_LengthI.ratio$rsquare, 2))
       ),
       bty= "n"
)


# ### Zpol

plot(LM_LengthZpol, "RMA", ylab = "Zpol essentricity", xlab = "bone length (log)", main = "")
legend("bottomleft", 
       legend = c(
         paste("intercept =", round(LM_LengthZpol$regression.results[4,2], digits = 2)),
         paste("slope =", round(LM_LengthZpol$regression.results[4,3], digits = 2),
               "±", round(LM_LengthZpol$confidence.intervals[4,5] -LM_LengthZpol$confidence.intervals[4,4], 2)),
         paste("p =", LM_LengthZpol$regression.results$`P-perm (1-tailed)`[4]),
         paste("Rsq =", round(LM_LengthZpol$rsquare, 2))
       ),
       bty= "n"
)


# ### abs.ecc50.J

plot(LM_LengthJ, "RMA", ylab = "J essentricity", xlab = "bone length (log)", main = "")
legend("bottomleft", 
       legend = c(
         paste("intercept =", round(LM_LengthJ$regression.results[4,2], digits = 2)),
         paste("slope =", round(LM_LengthJ$regression.results[4,3], digits = 2),
               "±", round(LM_LengthJ$confidence.intervals[4,5] -LM_LengthJ$confidence.intervals[4,4], 2)),
         paste("p =", LM_LengthJ$regression.results$`P-perm (1-tailed)`[4]),
         paste("Rsq =", round(LM_LengthJ$rsquare, 2))
       ),
       bty= "n"
)




# plot together -----------------------------------------------------------

# par(mfrow = c(4,4), mar = c(5,4,1,2))

# using reduced major axis

plot(LM_CircCurvature, "RMA", ylab = "Curvature essentricity", xlab = "midshaft circumferance (log)", main = "")
legend("bottomleft", 
       legend = c(
         paste("intercept =", round(LM_CircCurvature$regression.results[4,2], digits = 2)),
         paste("slope =", round(LM_CircCurvature$regression.results[4,3], digits = 2),
               "±", round(LM_CircCurvature$confidence.intervals[4,5] -LM_CircCurvature$confidence.intervals[4,4], 2)),
         paste("p =", LM_CircCurvature$regression.results$`P-perm (1-tailed)`[4]),
         paste("Rsq =", round(LM_CircCurvature$rsquare, 2))
       ),
       bty= "n"
)


# ### CORTICAL AREA


plot(LM_CircCA, "RMA", ylab = "Cortical Area essentricity", xlab = "midshaft circumferance (log)", main = "")
legend("bottomleft", 
       legend = c(
         paste("intercept =", round(LM_CircCA$regression.results[4,2], digits = 2)),
         paste("slope =", round(LM_CircCA$regression.results[4,3], digits = 2),
               "±", round(LM_CircCA$confidence.intervals[4,5] -LM_CircCA$confidence.intervals[4,4], 2)),
         paste("p =", LM_CircCA$regression.results$`P-perm (1-tailed)`[4]),
         paste("Rsq =", round(LM_CircCA$rsquare, 2))
       ),
       bty= "n"
)

plot(LM_Circthickness, "RMA", ylab = "mean Cortical Thickness essentricity", xlab = "midshaft circumferance (log)", main = "")
legend("bottomleft", 
       legend = c(
         paste("intercept =", round(LM_Circthickness$regression.results[4,2], digits = 2)),
         paste("slope =", round(LM_Circthickness$regression.results[4,3], digits = 2),
               "±", round(LM_Circthickness$confidence.intervals[4,5] -LM_Circthickness$confidence.intervals[4,4], 2)),
         paste("p =", LM_Circthickness$regression.results$`P-perm (1-tailed)`[4]),
         paste("Rsq =", round(LM_Circthickness$rsquare, 2))
       ),
       bty= "n"
)


# ### Ix

plot(LM_CircIx, "RMA", ylab = "Ix essentricity", xlab = "midshaft circumferance (log)", main = "")
legend("bottomleft", 
       legend = c(
         paste("intercept =", round(LM_CircIx$regression.results[4,2], digits = 2)),
         paste("slope =", round(LM_CircIx$regression.results[4,3], digits = 2),
               "±", round(LM_CircIx$confidence.intervals[4,5] -LM_CircIx$confidence.intervals[4,4], 2)),
         paste("p =", LM_CircIx$regression.results$`P-perm (1-tailed)`[4]),
         paste("Rsq =", round(LM_CircIx$rsquare, 2))),
       bty= "n")


# ### Iy

plot(LM_CircIy, "RMA", ylab = "Iy essentricity", xlab = "midshaft circumferance (log)", main = "")
legend("bottomleft", 
       legend = c(
         paste("intercept =", round(LM_CircIy$regression.results[4,2], digits = 2)),
         paste("slope =", round(LM_CircIy$regression.results[4,3], digits = 2),
               "±", round(LM_CircIy$confidence.intervals[4,5] -LM_CircIy$confidence.intervals[4,4], 2)),
         paste("p =", LM_CircIy$regression.results$`P-perm (1-tailed)`[4]),
         paste("Rsq =", round(LM_CircIy$rsquare, 2))),
       bty= "n")


# ### I.ratio

plot(LM_CircI.ratio, "RMA", ylab = "I.ratio essentricity", xlab = "midshaft circumferance (log)", main = "")
legend("bottomleft", 
       legend = c(
         paste("intercept =", round(LM_CircI.ratio$regression.results[4,2], digits = 2)),
         paste("slope =", round(LM_CircI.ratio$regression.results[4,3], digits = 2),
               "±", round(LM_CircI.ratio$confidence.intervals[4,5] -LM_CircI.ratio$confidence.intervals[4,4], 2)),
         paste("p =", LM_CircI.ratio$regression.results$`P-perm (1-tailed)`[4]),
         paste("Rsq =", round(LM_CircI.ratio$rsquare, 2))),
       bty= "n")



# ### Zpol

plot(LM_CircZpol, "RMA", ylab = "Zpol essentricity", xlab = "midshaft circumferance (log)", main = "")
legend("bottomleft", 
       legend = c(
         paste("intercept =", round(LM_CircZpol$regression.results[4,2], digits = 2)),
         paste("slope =", round(LM_CircZpol$regression.results[4,3], digits = 2),
               "±", round(LM_CircZpol$confidence.intervals[4,5] -LM_CircZpol$confidence.intervals[4,4], 2)),
         paste("p =", LM_CircZpol$regression.results$`P-perm (1-tailed)`[4]),
         paste("Rsq =", round(LM_CircZpol$rsquare, 2))),
       bty= "n")




# ### abs.ecc50.J


plot(LM_CircJ, "RMA", ylab = "J essentricity", xlab = "midshaft circumferance (log)", main = "")
legend("bottomleft", 
       legend = c(
         paste("intercept =", round(LM_CircJ$regression.results[4,2], digits = 2)),
         paste("slope =", round(LM_CircJ$regression.results[4,3], digits = 2),
               "±", round(LM_CircJ$confidence.intervals[4,5] -LM_CircJ$confidence.intervals[4,4], 2)),
         paste("p =", LM_CircJ$regression.results$`P-perm (1-tailed)`[4]),
         paste("Rsq =", round(LM_CircJ$rsquare, 2))),
       bty= "n")

