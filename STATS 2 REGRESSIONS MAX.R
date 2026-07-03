
rm(list = ls())

library(lmodel2)
load("boneVars.Rdata")
eco <- as.data.frame(read.csv("./DATA/specimenEco2.csv"))


# MAX THICKNESS -----------------------------------------------------------
# using reduced major axis

dat <- data.frame(boneVars$ecoDF, boneVars$maxDF, boneVars$meanDF, mass = eco$MASS..KG.)



  
# plot --------------------------------------------------------------------




###########################################
# CALCULATE SCALING USING REDUCED MAJOR AXIS
############################################

# install.packages("lmodel2")

dat <- data.frame(boneVars$ecoDF, boneVars$maxDF, boneVars$meanDF, mass = eco$MASS..KG.)


# ### CURVATURE

# sink("./stats/max regression stats.txt.")

LM_massCurvature <- lmodel2((max.Curvature)~(mass), "interval", "interval",data = dat, 99)
LM_massCurvature

LM_LengthCurvature <- lmodel2((max.Curvature)~(ZLENGTH), "interval", "interval",data = dat, 99)
LM_LengthCurvature

LM_CircCurvature <- lmodel2((max.Curvature)~(MIDSHAFTCIRC), "interval", "interval",data = dat, 99)
LM_CircCurvature

# ### CORTICAL AREA

LM_massCA <- lmodel2((max.CA)~(mass), "interval", "interval",data = dat, 99)
LM_massCA

LM_LengthCA <- lmodel2((max.CA)~(ZLENGTH), "interval", "interval",data = dat, 99)
LM_LengthCA

LM_CircCA <- lmodel2((max.CA)~(MIDSHAFTCIRC), "interval", "interval",data = dat, 99)
LM_CircCA

# ### THICKNESS

LM_massthickness <- lmodel2((max.thickness)~(mass), "interval", "interval",data = dat, 99)
LM_massthickness


LM_Lengththickness <- lmodel2((max.thickness)~(ZLENGTH), "interval", "interval",data = dat, 99)
LM_Lengththickness

LM_Circthickness <- lmodel2((max.thickness)~(MIDSHAFTCIRC), "interval", "interval",data = dat, 99)
LM_Circthickness

# ### I.X

LM_massIx <- lmodel2((max.Ix)~(mass), "interval", "interval",data = dat, 99)
LM_massIx

LM_LengthIx <- lmodel2((max.Ix)~(ZLENGTH), "interval", "interval",data = dat, 99)
LM_LengthIx

LM_CircIx <- lmodel2((max.Ix)~(MIDSHAFTCIRC), "interval", "interval",data = dat, 99)
LM_CircIx

# ### I.Y

LM_massIy <- lmodel2((max.Iy)~(mass), "interval", "interval",data = dat, 99)
LM_massIy

LM_LengthIy <- lmodel2((max.Iy)~(ZLENGTH), "interval", "interval",data = dat, 99)
LM_LengthIy

LM_CircIy <- lmodel2((max.Iy)~(MIDSHAFTCIRC), "interval", "interval",data = dat, 99)
LM_CircIy

# ### I.ratio

LM_massI.ratio <- lmodel2((max.I.ratio)~(mass), "interval", "interval",data = dat, 99)
LM_massI.ratio


LM_LengthI.ratio <- lmodel2((max.I.ratio)~(ZLENGTH), "interval", "interval",data = dat, 99)
LM_LengthI.ratio

LM_CircI.ratio <- lmodel2((max.I.ratio)~(MIDSHAFTCIRC), "interval", "interval",data = dat, 99)
LM_CircI.ratio

# ### Zpol


LM_massZpol <- lmodel2((max.Zpol)~(mass), "interval", "interval",data = dat, 99)
LM_massZpol


LM_LengthZpol <- lmodel2((max.Zpol)~(ZLENGTH), "interval", "interval",data = dat, 99)
LM_LengthZpol

LM_LengthZpol$regression.results

LM_CircZpol <- lmodel2((max.Zpol)~(MIDSHAFTCIRC), "interval", "interval",data = dat, 99)
LM_CircZpol

# ### max.J
LM_massJ <- lmodel2((max.J)~(mass), "interval", "interval",data = dat, 99)
LM_massJ


LM_LengthJ <- lmodel2((max.J)~(ZLENGTH), "interval", "interval",data = dat, 99)
LM_LengthJ

LM_CircJ <- lmodel2((max.J)~(MIDSHAFTCIRC), "interval", "interval",data = dat, 99)
LM_CircJ

# sink()


# -------------------------------------------------------------------------


## SAVE RESIDUALS

# res_massCurvature <- residuals(lm(max.Curvature~mass, data = dat))
# res_massCA <- residuals(lm(max.CA~mass, data = dat))
# res_massthickness <- residuals(lm(max.thickness~mass, data = dat))
# res_massIx <- residuals(lm(max.Ix~mass, data = dat))
# res_massIy <- residuals(lm(max.Iy~mass, data = dat))
# res_massI.ratio <- residuals(lm(max.I.ratio~mass, data = dat))
# res_massZpol <- residuals(lm(max.Zpol~mass, data = dat))
# res_massJ <- residuals(lm(max.J~mass, data = dat))

res_massCurvature <- residuals(lm(log(max.Curvature)~log(mass), data = dat))
res_massCA <- residuals(lm(log(max.CA)~log(mass), data = dat))
res_massthickness <- residuals(lm(log(max.thickness)~log(mass), data = dat))
res_massIx <- residuals(lm(log(max.Ix)~log(mass), data = dat))
res_massIy <- residuals(lm(log(max.Iy)~log(mass), data = dat))
res_massI.ratio <- residuals(lm(log(max.I.ratio)~log(mass), data = dat))
res_massZpol <- residuals(lm(log(max.Zpol)~log(mass), data = dat))
res_massJ <- residuals(lm(log(max.J)~log(mass), data = dat))



res_data <- as.data.frame(cbind(res_massCurvature,
                                res_massCA,
                                res_massthickness, 
                                res_massIx, 
                                res_massIy, 
                                res_massI.ratio, 
                                res_massZpol, 
                                res_massJ))



# -------------------------------------------------------------------------





res_LengthCurvature <- residuals(lm(max.Curvature~ZLENGTH, data = dat))
res_CircCurvature   <- residuals(lm(max.Curvature~MIDSHAFTCIRC, data = dat))

res_LengthCA <- residuals(lm(max.CA~ZLENGTH, data = dat))
res_CircCA   <- residuals(lm(max.CA~MIDSHAFTCIRC, data = dat))

res_Lengththickness <- residuals(lm(max.thickness~ZLENGTH, data = dat))
res_Circthickness   <- residuals(lm(max.thickness~MIDSHAFTCIRC, data = dat))

res_LengthIx <- residuals(lm(max.Ix~ZLENGTH, data = dat))
res_CircIx   <- residuals(lm(max.Ix~MIDSHAFTCIRC, data = dat))

res_LengthIy <- residuals(lm(max.Iy~ZLENGTH, data = dat))
res_CircIy   <- residuals(lm(max.Iy~MIDSHAFTCIRC, data = dat))

res_LengthI.ratio <- residuals(lm(max.I.ratio~ZLENGTH, data = dat))
res_CircI.ratio   <- residuals(lm(max.I.ratio~MIDSHAFTCIRC, data = dat))

res_LengthZpol <- residuals(lm(max.Zpol~ZLENGTH, data = dat))
res_CircZpol   <- residuals(lm(max.Zpol~MIDSHAFTCIRC, data = dat))

res_LengthJ <- residuals(lm(max.J~ZLENGTH, data = dat))
res_CircJ   <- residuals(lm(max.J~MIDSHAFTCIRC, data = dat))

res_data <- as.data.frame(cbind(res_massCurvature,res_LengthCurvature, res_CircCurvature, 
                                res_massCA,res_LengthCA, res_CircCA,
                                res_massthickness, res_Lengththickness, res_Circthickness,
                                res_massIx, res_LengthIx, res_CircIx,
                                res_massIy, res_LengthIy, res_CircIy,
                                res_massI.ratio, res_LengthI.ratio, res_CircI.ratio,
                                res_massZpol, res_LengthZpol, res_CircZpol,
                                res_massJ, res_LengthJ, res_CircJ))





par(mfrow = c(8,3))

for(i in 1:ncol(res_data)){
  
  plot(res_data[,i], xlab = colnames(res_data)[i])
  
}

save(res_data, file= "res_data.Rdata")



# log-log regressions -----------------------------------------------------


dat <- data.frame(boneVars$ecoDF, boneVars$maxDF, boneVars$meanDF, mass = eco$MASS..KG.)


# ### CURVATURE

# sink("./stats/max regression stats.txt.")

LM_massCurvature <- lmodel2(log(max.Curvature)~log(mass), "interval", "interval",data = dat, 99)
LM_massCurvature

LM_LengthCurvature <- lmodel2(log(max.Curvature)~log(ZLENGTH), "interval", "interval",data = dat, 99)
LM_LengthCurvature

LM_CircCurvature <- lmodel2(log(max.Curvature)~log(MIDSHAFTCIRC), "interval", "interval",data = dat, 99)
LM_CircCurvature

# ### CORTICAL AREA

LM_massCA <- lmodel2(log(max.CA)~log(mass), "interval", "interval",data = dat, 99)
LM_massCA

LM_LengthCA <- lmodel2(log(max.CA)~log(ZLENGTH), "interval", "interval",data = dat, 99)
LM_LengthCA

LM_CircCA <- lmodel2(log(max.CA)~log(MIDSHAFTCIRC), "interval", "interval",data = dat, 99)
LM_CircCA

# ### THICKNESS

LM_massthickness <- lmodel2(log(max.thickness)~log(mass), "interval", "interval",data = dat, 99)
LM_massthickness


LM_Lengththickness <- lmodel2(log(max.thickness)~log(ZLENGTH), "interval", "interval",data = dat, 99)
LM_Lengththickness

LM_Circthickness <- lmodel2(log(max.thickness)~log(MIDSHAFTCIRC), "interval", "interval",data = dat, 99)
LM_Circthickness

# ### Ix

LM_massIx <- lmodel2(log(max.Ix)~log(mass), "interval", "interval",data = dat, 99)
LM_massIx

LM_LengthIx <- lmodel2(log(max.Ix)~log(ZLENGTH), "interval", "interval",data = dat, 99)
LM_LengthIx

LM_CircIx <- lmodel2(log(max.Ix)~log(MIDSHAFTCIRC), "interval", "interval",data = dat, 99)
LM_CircIx
# ### Iy
LM_massIy <- lmodel2(log(max.Iy)~log(mass), "interval", "interval",data = dat, 99)
LM_massIy

LM_LengthIy <- lmodel2(log(max.Iy)~log(ZLENGTH), "interval", "interval",data = dat, 99)
LM_LengthIy

LM_CircIy <- lmodel2(log(max.Iy)~log(MIDSHAFTCIRC), "interval", "interval",data = dat, 99)
LM_CircIy
# ### I.ratio

LM_massI.ratio <- lmodel2((max.I.ratio)~log(mass), "interval", "interval",data = dat, 99)
LM_massI.ratio


LM_LengthI.ratio <- lmodel2((max.I.ratio)~log(ZLENGTH), "interval", "interval",data = dat, 99)
LM_LengthI.ratio

LM_CircI.ratio <- lmodel2((max.I.ratio)~log(MIDSHAFTCIRC), "interval", "interval",data = dat, 99)
LM_CircI.ratio

# ### Zpol

LM_massZpol <- lmodel2(log(max.Zpol)~log(mass), "interval", "interval",data = dat, 99)
LM_massZpol


LM_LengthZpol <- lmodel2(log(max.Zpol)~log(ZLENGTH), "interval", "interval",data = dat, 99)
LM_LengthZpol


LM_CircZpol <- lmodel2(log(max.Zpol)~log(MIDSHAFTCIRC), "interval", "interval",data = dat, 99)
LM_CircZpol

# ### max.J

LM_massJ <- lmodel2(log(max.J)~log(mass), "interval", "interval",data = dat, 99)
LM_massJ

LM_LengthJ <- lmodel2(log(max.J)~log(ZLENGTH), "interval", "interval",data = dat, 99)
LM_LengthJ

LM_CircJ <- lmodel2(log(max.J)~log(MIDSHAFTCIRC), "interval", "interval",data = dat, 99)
LM_CircJ

# sink()



# make p val list ---------------------------------------------------------
pvals <- c(round(LM_massCurvature$regression.results$Slope, 3),
           round(LM_massCA$P.param, 3),
           round(LM_massthickness$P.param, 3),
           round(LM_massIx$P.param, 3),
           round(LM_massIy$P.param, 3),
           round(LM_massI.ratio$P.param, 3),
           round(LM_CircZpol$P.param, 3),
           round(LM_massJ$P.param, 3))


pvals <- c(round(LM_massCurvature$P.param, 3),
  round(LM_massCA$P.param, 3),
  round(LM_massthickness$P.param, 3),
  round(LM_massIx$P.param, 3),
  round(LM_massIy$P.param, 3),
  round(LM_massI.ratio$P.param, 3),
  round(LM_CircZpol$P.param, 3),
  round(LM_massJ$P.param, 3))

pvals <- c(round(LM_massCurvature$rsquare, 3),
           round(LM_massCA$rsquare, 3),
           round(LM_massthickness$rsquare, 3),
           round(LM_massIx$rsquare, 3),
           round(LM_massIy$rsquare, 3),
           round(LM_massI.ratio$rsquare, 3),
           round(LM_CircZpol$rsquare, 3),
           round(LM_massJ$rsquare, 3))





# plot body mass regs -----------------------------------------------------
# LM_massCurvature$
  

par(mfrow = c(2,4), mar = c(5,4,1,2))


plot(LM_massCurvature, "RMA", ylab = "Curvature (log)", xlab = "body mass (log)", main = "")
legend("bottomright", 
       legend = c(
         paste("intercept =", round(LM_massCurvature$regression.results[4,2], digits = 2)),
         paste("slope =", round(LM_massCurvature$regression.results[4,3], digits = 2),
               "±", round(LM_massCurvature$confidence.intervals[4,5] -LM_massCurvature$confidence.intervals[4,4], 2)),
         paste("p < 0.001"),
         paste("Rsq =", round(LM_massCurvature$rsquare, 2))
       ),
       bty= "n"
)


# ### CORTICAL AREA


plot(LM_massCA, "RMA", ylab = "Cortical Area (log)", xlab = "body mass (log)", main = "")
legend("bottomright", 
       legend = c(
         paste("intercept =", round(LM_massCA$regression.results[4,2], digits = 2)),
         paste("slope =", round(LM_massCA$regression.results[4,3], digits = 2),
               "±", round(LM_massCA$confidence.intervals[4,5] -LM_massCA$confidence.intervals[4,4], 2)),
         paste("p < 0.001"),
         paste("Rsq =", round(LM_massCA$rsquare, 2))
       ),
       bty= "n"
)



# ### THICKNESS

plot(LM_massthickness, "RMA", ylab = "mean Cortical Thickness (log)", xlab = "body mass (log)", main = "")
legend("bottomright", 
       legend = c(
         paste("intercept =", round(LM_massthickness$regression.results[4,2], digits = 2)),
         paste("slope =", round(LM_massthickness$regression.results[4,3], digits = 2),
               "±", round(LM_massthickness$confidence.intervals[4,5] -LM_massthickness$confidence.intervals[4,4], 2)),
         paste("p < 0.001"),
         paste("Rsq =", round(LM_massthickness$rsquare, 2))
       ),
       bty= "n"
)


# ### I.x
# LM_massIx$P.param
plot(LM_massIx, "RMA", ylab = "I.x (log)", xlab = "body mass (log)", main = "")
legend("bottomright", 
       legend = c(
         paste("intercept =", round(LM_massIx$regression.results[4,2], digits = 2)),
         paste("slope =", round(LM_massIx$regression.results[4,3], digits = 2),
               "±", round(LM_massIx$confidence.intervals[4,5] -LM_massIx$confidence.intervals[4,4], 2)),
         paste("p < 0.001"),
         paste("Rsq =", round(LM_massIx$rsquare, 2))
       ),
       bty= "n"
)


# ### I.y
# LM_massIy$P.param
plot(LM_massIy, "RMA", ylab = "I.y (log)", xlab = "body mass (log)", main = "")
legend("bottomright", 
       legend = c(
         paste("intercept =", round(LM_massIy$regression.results[4,2], digits = 2)),
         paste("slope =", round(LM_massIy$regression.results[4,3], digits = 2),
               "±", round(LM_massIy$confidence.intervals[4,5] -LM_massIy$confidence.intervals[4,4], 2)),
         paste("p < 0.001"),
         paste("Rsq =", round(LM_massIy$rsquare, 2))
       ),
       bty= "n"
)


# ### I.ratio
# LM_massI.ratio$P.param
plot(LM_massI.ratio, "RMA", ylab = "I.ratio ", xlab = "body mass (log)", main = "")
legend("bottomright", 
       legend = c(
         paste("intercept =", round(LM_massI.ratio$regression.results[4,2], digits = 2)),
         paste("slope =", round(LM_massI.ratio$regression.results[4,3], digits = 2),
               "±", round(LM_massI.ratio$confidence.intervals[4,5] -LM_massI.ratio$confidence.intervals[4,4], 2)),
         paste("p =",  round(LM_massI.ratio$P.param, 3)),
         paste("Rsq =", round(LM_massI.ratio$rsquare, 2))
       ),
       bty= "n"
)



# ### Zpol

plot(LM_massZpol, "RMA", ylab = "Zpol (log)", xlab = "body mass (log)", main = "")
legend("bottomright", 
       legend = c(
         paste("intercept =", round(LM_massZpol$regression.results[4,2], digits = 2)),
         paste("slope =", round(LM_massZpol$regression.results[4,3], digits = 2),
               "±", round(LM_massZpol$confidence.intervals[4,5] -LM_massZpol$confidence.intervals[4,4], 2)),
         paste("p < 0.001"),
         paste("Rsq =", round(LM_massZpol$rsquare, 2))
       ),
       bty= "n"
)



# ### max.J

plot(LM_massJ, "RMA", ylab = "J (log)", xlab = "body mass (log)", main = "")
legend("bottomright", 
       legend = c(
         paste("intercept =", round(LM_massJ$regression.results[4,2], digits = 2)),
         paste("slope =", round(LM_massJ$regression.results[4,3], digits = 2),
               "±", round(LM_massJ$confidence.intervals[4,5] -LM_massJ$confidence.intervals[4,4], 2)),
         paste("p < 0.001"),
         paste("Rsq =", round(LM_massJ$rsquare, 2))
       ),
       bty= "n"
)




# -------------------------------------------------------------------------




# plot together -----------------------------------------------------------



# length and circ ---------------------------------------------------------


par(mfrow = c(4,4), mar = c(5,4,1,2))


plot(LM_LengthCurvature, "RMA", ylab = "Curvature (log)", xlab = "bone length (log)", main = "")
legend("bottomright", 
       legend = c(
         paste("intercept =", round(LM_LengthCurvature$regression.results[4,2], digits = 2)),
         paste("slope =", round(LM_LengthCurvature$regression.results[4,3], digits = 2),
               "±", round(LM_LengthCurvature$confidence.intervals[4,5] -LM_LengthCurvature$confidence.intervals[4,4], 2)),
         paste("p < 0.001"),
         paste("Rsq =", round(LM_LengthCurvature$rsquare, 2))
       ),
       bty= "n"
)


# ### CORTICAL AREA


plot(LM_LengthCA, "RMA", ylab = "Cortical Area (log)", xlab = "bone length (log)", main = "")
legend("bottomright", 
       legend = c(
         paste("intercept =", round(LM_LengthCA$regression.results[4,2], digits = 2)),
         paste("slope =", round(LM_LengthCA$regression.results[4,3], digits = 2),
               "±", round(LM_LengthCA$confidence.intervals[4,5] -LM_LengthCA$confidence.intervals[4,4], 2)),
         paste("p < 0.001"),
         paste("Rsq =", round(LM_LengthCA$rsquare, 2))
       ),
       bty= "n"
)


# ### THICKNESS

plot(LM_Lengththickness, "RMA", ylab = "mean Cortical Thickness (log)", xlab = "bone length (log)", main = "")
legend("bottomright", 
       legend = c(
         paste("intercept =", round(LM_Lengththickness$regression.results[4,2], digits = 2)),
         paste("slope =", round(LM_Lengththickness$regression.results[4,3], digits = 2),
               "±", round(LM_Lengththickness$confidence.intervals[4,5] -LM_Lengththickness$confidence.intervals[4,4], 2)),
         paste("p < 0.001"),
         paste("Rsq =", round(LM_Lengththickness$rsquare, 2))
       ),
       bty= "n"
)


# ### I.x
# LM_LengthIx$P.param
plot(LM_LengthIx, "RMA", ylab = "I.x ", xlab = "bone length(log)", main = "")
legend("bottomright", 
       legend = c(
         paste("intercept =", round(LM_LengthIx$regression.results[4,2], digits = 2)),
         paste("slope =", round(LM_LengthIx$regression.results[4,3], digits = 2),
               "±", round(LM_LengthIx$confidence.intervals[4,5] -LM_LengthIx$confidence.intervals[4,4], 2)),
         paste("p =",  round(LM_LengthIx$P.param, 3)),
         paste("Rsq =", round(LM_LengthIx$rsquare, 2))
       ),
       bty= "n"
)

# ### I.y
# LM_LengthIy$P.param
plot(LM_LengthIy, "RMA", ylab = "I.y ", xlab = "bone length(log)", main = "")
legend("bottomright", 
       legend = c(
         paste("intercept =", round(LM_LengthIy$regression.results[4,2], digits = 2)),
         paste("slope =", round(LM_LengthIy$regression.results[4,3], digits = 2),
               "±", round(LM_LengthIy$confidence.intervals[4,5] -LM_LengthIy$confidence.intervals[4,4], 2)),
         paste("p =",  round(LM_LengthIy$P.param, 3)),
         paste("Rsq =", round(LM_LengthIy$rsquare, 2))
       ),
       bty= "n"
)

# ### I.ratio
# LM_LengthI.ratio$P.param
plot(LM_LengthI.ratio, "RMA", ylab = "I.ratio ", xlab = "bone length(log)", main = "")
legend("bottomleft", 
       legend = c(
         paste("intercept =", round(LM_LengthI.ratio$regression.results[4,2], digits = 2)),
         paste("slope =", round(LM_LengthI.ratio$regression.results[4,3], digits = 2),
               "±", round(LM_LengthI.ratio$confidence.intervals[4,5] -LM_LengthI.ratio$confidence.intervals[4,4], 2)),
         paste("p =",  round(LM_LengthI.ratio$P.param, 3)),
         paste("Rsq =", round(LM_LengthI.ratio$rsquare, 2))
       ),
       bty= "n"
)


# ### Zpol

plot(LM_LengthZpol, "RMA", ylab = "Zpol ", xlab = "bone length(log)", main = "")
legend("bottomright", 
       legend = c(
         paste("intercept =", round(LM_LengthZpol$regression.results[4,2], digits = 2)),
         paste("slope =", round(LM_LengthZpol$regression.results[4,3], digits = 2),
               "±", round(LM_LengthZpol$confidence.intervals[4,5] -LM_LengthZpol$confidence.intervals[4,4], 2)),
         paste("p < 0.001"),
         paste("Rsq =", round(LM_LengthZpol$rsquare, 2))
       ),
       bty= "n"
)


# ### max.J

plot(LM_LengthJ, "RMA", ylab = "J ", xlab = "bone length (log)", main = "")
legend("bottomright", 
       legend = c(
         paste("intercept =", round(LM_LengthJ$regression.results[4,2], digits = 2)),
         paste("slope =", round(LM_LengthJ$regression.results[4,3], digits = 2),
               "±", round(LM_LengthJ$confidence.intervals[4,5] -LM_LengthJ$confidence.intervals[4,4], 2)),
         paste("p < 0.001"),
         paste("Rsq =", round(LM_LengthJ$rsquare, 2))
       ),
       bty= "n"
)


# -------------------------------------------------------------------------


# next --------------------------------------------------------------------



# length and circ ---------------------------------------------------------


# par(mfrow = c(4,4), mar = c(5,4,1,2))



plot(LM_CircCurvature, "RMA", ylab = "Curvature (log)", xlab = "midshaft circumferance (log)", main = "")
legend("bottomright", 
       legend = c(
         paste("intercept =", round(LM_CircCurvature$regression.results[4,2], digits = 2)),
         paste("slope =", round(LM_CircCurvature$regression.results[4,3], digits = 2),
               "±", round(LM_CircCurvature$confidence.intervals[4,5] -LM_CircCurvature$confidence.intervals[4,4], 2)),
         paste("p < 0.001"),
         paste("Rsq =", round(LM_CircCurvature$rsquare, 2))
       ),
       bty= "n"
)


# ### CORTICAL AREA



plot(LM_CircCA, "RMA", ylab = "Cortical Area (log)", xlab = "midshaft circumferance (log)", main = "")
legend("bottomright", 
       legend = c(
         paste("intercept =", round(LM_CircCA$regression.results[4,2], digits = 2)),
         paste("slope =", round(LM_CircCA$regression.results[4,3], digits = 2),
               "±", round(LM_CircCA$confidence.intervals[4,5] -LM_CircCA$confidence.intervals[4,4], 2)),
         paste("p < 0.001"),
         paste("Rsq =", round(LM_CircCA$rsquare, 2))
       ),
       bty= "n"
)


# ### THICKNESS


plot(LM_Circthickness, "RMA", ylab = "mean Cortical Thickness (log)", xlab = "midshaft circumferance(log)", main = "")
legend("bottomright", 
       legend = c(
         paste("intercept =", round(LM_Circthickness$regression.results[4,2], digits = 2)),
         paste("slope =", round(LM_Circthickness$regression.results[4,3], digits = 2),
               "±", round(LM_Circthickness$confidence.intervals[4,5] -LM_Circthickness$confidence.intervals[4,4], 2)),
         paste("p < 0.001"),
         paste("Rsq =", round(LM_Circthickness$rsquare, 2))
       ),
       bty= "n"
)


# ### I.x
# LM_LengthIx$P.param


plot(LM_CircIx, "RMA", ylab = "I.x ", xlab = "midshaft circumferance (log)", main = "")
legend("bottomright", 
       legend = c(
         paste("intercept =", round(LM_CircIx$regression.results[4,2], digits = 2)),
         paste("slope =", round(LM_CircIx$regression.results[4,3], digits = 2),
               "±", round(LM_CircIx$confidence.intervals[4,5] -LM_CircIx$confidence.intervals[4,4], 2)),
         paste("p =",  round(LM_CircIx$P.param, 3)),
         paste("Rsq =", round(LM_CircIx$rsquare, 2))),
       bty= "n")

# ### I.y
# LM_LengthIy$P.param

plot(LM_CircIy, "RMA", ylab = "I.y ", xlab = "midshaft circumferance (log)", main = "")
legend("bottomright", 
       legend = c(
         paste("intercept =", round(LM_CircIy$regression.results[4,2], digits = 2)),
         paste("slope =", round(LM_CircIy$regression.results[4,3], digits = 2),
               "±", round(LM_CircIy$confidence.intervals[4,5] -LM_CircIy$confidence.intervals[4,4], 2)),
         paste("p =",  round(LM_CircIy$P.param, 3)),
         paste("Rsq =", round(LM_CircIy$rsquare, 2))),
       bty= "n")

# ### I.ratio
# LM_LengthI.ratio$P.param

plot(LM_CircI.ratio, "RMA", ylab = "I.ratio ", xlab = "midshaft circumferance (log)", main = "")
legend("bottomleft", 
       legend = c(
         paste("intercept =", round(LM_CircI.ratio$regression.results[4,2], digits = 2)),
         paste("slope =", round(LM_CircI.ratio$regression.results[4,3], digits = 2),
               "±", round(LM_CircI.ratio$confidence.intervals[4,5] -LM_CircI.ratio$confidence.intervals[4,4], 2)),
         paste("p =",  round(LM_CircI.ratio$P.param, 3)),
         paste("Rsq =", round(LM_CircI.ratio$rsquare, 2))),
       bty= "n")



# ### Zpol



plot(LM_CircZpol, "RMA", ylab = "Zpol ", xlab = "midshaft circumferance (log)", main = "")
legend("bottomright", 
       legend = c(
         paste("intercept =", round(LM_CircZpol$regression.results[4,2], digits = 2)),
         paste("slope =", round(LM_CircZpol$regression.results[4,3], digits = 2),
               "±", round(LM_CircZpol$confidence.intervals[4,5] -LM_CircZpol$confidence.intervals[4,4], 2)),
         paste("p < 0.001"),
         paste("Rsq =", round(LM_CircZpol$rsquare, 2))),
       bty= "n")




# ### max.J



plot(LM_CircJ, "RMA", ylab = "J ", xlab = "midshaft circumferance (log)", main = "")
legend("bottomright", 
       legend = c(
         paste("intercept =", round(LM_CircJ$regression.results[4,2], digits = 2)),
         paste("slope =", round(LM_CircJ$regression.results[4,3], digits = 2),
               "±", round(LM_CircJ$confidence.intervals[4,5] -LM_CircJ$confidence.intervals[4,4], 2)),
         paste("p < 0.001"),
         paste("Rsq =", round(LM_CircJ$rsquare, 2))),
       bty= "n")


# -------------------------------------------------------------------------


# next --------------------------------------------------------------------

