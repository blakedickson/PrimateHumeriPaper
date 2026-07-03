rm(list = ls())


library(lmodel2)
load("boneVars.Rdata")
eco <- as.data.frame(read.csv("./DATA/specimenEco2.csv"))


# -------------------------------------------------------------------------


dat <- data.frame(mass = eco$MASS..KG., boneVars$abs50DF)

dat <- data.frame(mass = eco$MASS..KG., boneVars$abs50DF$abs.ecc50.Curvature, boneVars$maxDF[,-1])

dat <- data.frame(mass = eco$MASS..KG., boneVars$maxDF$max.Curvature, boneVars$abs.ecc50.Curvature[,-1])

dat[,2] <- boneVars$maxDF$max.Curvature

colnames(dat) <- c("mass", "Curvature", "Thickness", "CA", "Ix", "Iy", "I.ratio", "Zpol", "J")

# -------------------------------------------------------------------------


# res max curve v max csg -----------------------------------------------------


# ### CORTICAL AREA

# LM_CurveCA <- lmodel2(CA~Curvature, "interval", "interval",data = dat, 99)
LM_CurveCA <- lmodel2(CA~Curvature, "interval", "interval",data = dat, 99)
LM_CurveCA

# ### THICKNESS

LM_Curvethickness <- lmodel2(Thickness~Curvature, "interval", "interval",data = dat, 99)
LM_Curvethickness

# ### I.ratio

LM_CurveIx <- lmodel2(Ix~Curvature, "interval", "interval",data = dat, 99)
LM_CurveIx
# ### I.ratio

LM_CurveIy <- lmodel2(Iy~Curvature, "interval", "interval",data = dat, 99)
LM_CurveIy
# ### I.ratio

LM_CurveI.ratio <- lmodel2(I.ratio~Curvature, "interval", "interval",data = dat, 99)
LM_CurveI.ratio

# ### Zpol

LM_CurveZpol <- lmodel2(Zpol~Curvature, "interval", "interval",data = dat, 99)
LM_CurveZpol

# ### max.J

LM_CurveJ <- lmodel2(J~Curvature, "interval", "interval",data = dat, 99)
LM_CurveJ



eccCurveCSG <- list(LM_CurveCA = LM_CurveCA, 
                    LM_Curvethickness= LM_Curvethickness, 
                    LM_CurveIx = LM_CurveIx, 
                    LM_CurveIy = LM_CurveIy, 
                    LM_CurveI.ratio = LM_CurveI.ratio, 
                    LM_CurveZpol = LM_CurveZpol, 
                    LM_CurveJ = LM_CurveJ)

# resCurveCSG

sink("./stats/LMs residuals Curvature vs max CSG.txt")
eccCurveCSG
sink()


eccCurveCSG


# plot --------------------------------------------------------------------





par(mfrow = c(2,4), mar = c(5,4,1,2))


plot(eccCurveCSG$LM_CurveCA
     , "RMA", ylab = "ecc CA ", xlab = "Curvature ", main = "")
legend("bottomright", 
       legend = c(
         paste("intercept =", round(eccCurveCSG$LM_CurveCA$regression.results[4,2], digits = 2)),
         paste("slope =", round(eccCurveCSG$LM_CurveCA$regression.results[4,3], digits = 2),
               "±", round(eccCurveCSG$LM_CurveCA$confidence.intervals[4,5]-eccCurveCSG$LM_CurveCA$confidence.intervals[4,4], 2)),
         paste("p =",  round(eccCurveCSG$LM_CurveCA$P.param, 2)),
         paste("Rsq =", round(eccCurveCSG$LM_CurveCA$rsquare, 2))
       ),
       bty= "n"
)

plot(eccCurveCSG$LM_Curvethickness
     , "RMA", ylab = "ecc Thickness ", xlab = "Curvature ", main = "")
legend("bottomright", 
       legend = c(
         paste("intercept =", round(eccCurveCSG$LM_Curvethickness$regression.results[4,2], digits = 2)),
         paste("slope =", round(eccCurveCSG$LM_Curvethickness$regression.results[4,3], digits = 2),
               "±", round(eccCurveCSG$LM_Curvethickness$confidence.intervals[4,5]-eccCurveCSG$LM_Curvethickness$confidence.intervals[4,4], 2)),
         paste("p =",  round(eccCurveCSG$LM_Curvethickness$P.param, 2)),
         paste("Rsq =", round(eccCurveCSG$LM_Curvethickness$rsquare, 2))
       ),
       bty= "n"
)


plot(eccCurveCSG$LM_CurveIx
     , "RMA", ylab = "ecc Ix ", xlab = "Curvature ", main = "")
legend("bottomright", 
       legend = c(
         paste("intercept =", round(eccCurveCSG$LM_CurveIx$regression.results[4,2], digits = 2)),
         paste("slope =", round(eccCurveCSG$LM_CurveIx$regression.results[4,3], digits = 2),
               "±", round(eccCurveCSG$LM_CurveIx$confidence.intervals[4,5]-eccCurveCSG$LM_CurveIx$confidence.intervals[4,4], 2)),
         paste("p =",  round(eccCurveCSG$LM_CurveIx$P.param, 2)),
         paste("Rsq =", round(eccCurveCSG$LM_CurveIx$rsquare, 2))
       ),
       bty= "n"
)



plot(eccCurveCSG$LM_CurveIy
     , "RMA", ylab = "ecc Iy ", xlab = "Curvature ", main = "")
legend("bottomright", 
       legend = c(
         paste("intercept =", round(eccCurveCSG$LM_CurveIy$regression.results[4,2], digits = 2)),
         paste("slope =", round(eccCurveCSG$LM_CurveIy$regression.results[4,3], digits = 2),
               "±", round(eccCurveCSG$LM_CurveIy$confidence.intervals[4,5]-eccCurveCSG$LM_CurveIy$confidence.intervals[4,4], 2)),
         paste("p =",  round(eccCurveCSG$LM_CurveIy$P.param, 2)),
         paste("Rsq =", round(eccCurveCSG$LM_CurveIy$rsquare, 2))
       ),
       bty= "n"
)



plot(eccCurveCSG$LM_CurveI.ratio
     , "RMA", ylab = "ecc I-Ratio ", xlab = "Curvature ", main = "")
legend("bottomright", 
       legend = c(
         paste("intercept =", round(eccCurveCSG$LM_CurveI.ratio$regression.results[4,2], digits = 2)),
         paste("slope =", round(eccCurveCSG$LM_CurveI.ratio$regression.results[4,3], digits = 2),
               "±", round(eccCurveCSG$LM_CurveI.ratio$confidence.intervals[4,5]-eccCurveCSG$LM_CurveI.ratio$confidence.intervals[4,4], 2)),
         paste("p =",  round(eccCurveCSG$LM_CurveI.ratio$P.param, 2)),
         paste("Rsq =", round(eccCurveCSG$LM_CurveI.ratio$rsquare, 2))
       ),
       bty= "n"
)


plot(eccCurveCSG$LM_CurveZpol
     , "RMA", ylab = "ecc Zpol ", xlab = "Curvature ", main = "")
legend("bottomright", 
       legend = c(
         paste("intercept =", round(eccCurveCSG$LM_CurveZpol$regression.results[4,2], digits = 2)),
         paste("slope =", round(eccCurveCSG$LM_CurveZpol$regression.results[4,3], digits = 2),
               "±", round(eccCurveCSG$LM_CurveZpol$confidence.intervals[4,5]-eccCurveCSG$LM_CurveZpol$confidence.intervals[4,4], 2)),
         paste("p =",  round(eccCurveCSG$LM_CurveZpol$P.param, 2)),
         paste("Rsq =", round(eccCurveCSG$LM_CurveZpol$rsquare, 2))
       ),
       bty= "n"
)


plot(eccCurveCSG$LM_CurveJ
     , "RMA", ylab = "ecc J ", xlab = "Curvature ", main = "")
legend("bottomright", 
       legend = c(
         paste("intercept =", round(eccCurveCSG$LM_CurveJ$regression.results[4,2], digits = 2)),
         paste("slope =", round(eccCurveCSG$LM_CurveJ$regression.results[4,3], digits = 2),
               "±", round(eccCurveCSG$LM_CurveJ$confidence.intervals[4,5]-eccCurveCSG$LM_CurveJ$confidence.intervals[4,4], 2)),
         paste("p =",  round(eccCurveCSG$LM_CurveJ$P.param, 2)),
         paste("Rsq =", round(eccCurveCSG$LM_CurveJ$rsquare, 2))
       ),
       bty= "n"
)




##### CURVE VS ecc VARIANCE

###### CURVE VS ESSENTRICITY/POSITION




