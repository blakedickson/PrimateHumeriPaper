rm(list = ls())


library(lmodel2)
load("boneVars.Rdata")
eco <- as.data.frame(read.csv("./DATA/specimenEco2.csv"))

# -------------------------------------------------------------------------


dat <- data.frame(boneVars$ecoDF, boneVars$maxDF, boneVars$meanDF, mass = eco$MASS..KG.)


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



save(res_data, file= "res_data.Rdata")
# -------------------------------------------------------------------------




###########################################
# CALCULATE Curvature vs CSG USING REDUCED MAJOR AXIS
############################################


load("res_data.Rdata")

# 
# resDAT<- data.frame(Curvature = res_data$res_massCurvature,
#            CA = res_data$res_massCA,
#            Thickness = res_data$res_massthickness, 
#            Ix = res_data$res_massIx,
#            Iy = res_data$res_massIy,
#            Iratio = res_data$res_massI.ratio,
#            Zpol = res_data$res_massZpol,
#            J = res_data$res_massJ)

resDAT<- data.frame(Curvature = dat$max.Curvature,
           CA = res_data$res_massCA,
           Thickness = res_data$res_massthickness, 
           Ix = res_data$res_massIx,
           Iy = res_data$res_massIy,
           Iratio = res_data$res_massI.ratio,
           Zpol = res_data$res_massZpol,
           J = res_data$res_massJ)

# res max curve v max csg -----------------------------------------------------


# ### CORTICAL AREA

# LM_CurveCA <- lmodel2(CA~Curvature, "interval", "interval",data = resDAT, 99)
LM_CurveCA <- lmodel2(CA~Curvature, "interval", "interval",data = resDAT, 99)
LM_CurveCA

# ### THICKNESS

LM_Curvethickness <- lmodel2(Thickness~Curvature, "interval", "interval",data = resDAT, 99)
LM_Curvethickness

# ### I.ratio

LM_CurveIx <- lmodel2(Ix~Curvature, "interval", "interval",data = resDAT, 99)
LM_CurveIx
# ### I.ratio

LM_CurveIy <- lmodel2(Iy~Curvature, "interval", "interval",data = resDAT, 99)
LM_CurveIy
# ### I.ratio

LM_CurveI.ratio <- lmodel2(Iratio~Curvature, "interval", "interval",data = resDAT, 99)
LM_CurveI.ratio

# ### Zpol

LM_CurveZpol <- lmodel2(Zpol~Curvature, "interval", "interval",data = resDAT, 99)
LM_CurveZpol

# ### max.J

LM_CurveJ <- lmodel2(J~Curvature, "interval", "interval",data = resDAT, 99)
LM_CurveJ



resCurveCSG <- list(LM_CurveCA = LM_CurveCA, 
                    LM_Curvethickness= LM_Curvethickness, 
                    LM_CurveIx = LM_CurveIx, 
                    LM_CurveIy = LM_CurveIy, 
                    LM_CurveI.ratio = LM_CurveI.ratio, 
                    LM_CurveZpol = LM_CurveZpol, 
                    LM_CurveJ = LM_CurveJ)

# resCurveCSG

sink("./stats/LMs residuals Curvature vs max CSG.txt")
resCurveCSG
sink()


resCurveCSG


# plot --------------------------------------------------------------------





par(mfrow = c(4,2), mar = c(5,4,1,2))


plot(resCurveCSG$LM_CurveCA
     , "RMA", ylab = "res CA ", xlab = "Curvature ", main = "")
legend("bottomright", 
       legend = c(
         paste("intercept =", round(resCurveCSG$LM_CurveCA$regression.results[4,2], digits = 2)),
         paste("slope =", round(resCurveCSG$LM_CurveCA$regression.results[4,3], digits = 2),
               "±", round(resCurveCSG$LM_CurveCA$confidence.intervals[4,5]-resCurveCSG$LM_CurveCA$confidence.intervals[4,4], 2)),
         paste("p =",  round(resCurveCSG$LM_CurveCA$P.param, 2)),
         paste("Rsq =", round(resCurveCSG$LM_CurveCA$rsquare, 2))
       ),
       bty= "n"
)

plot(resCurveCSG$LM_Curvethickness
     , "RMA", ylab = "res Thickness ", xlab = "Curvature ", main = "")
legend("bottomright", 
       legend = c(
         paste("intercept =", round(resCurveCSG$LM_Curvethickness$regression.results[4,2], digits = 2)),
         paste("slope =", round(resCurveCSG$LM_Curvethickness$regression.results[4,3], digits = 2),
               "±", round(resCurveCSG$LM_Curvethickness$confidence.intervals[4,5]-resCurveCSG$LM_Curvethickness$confidence.intervals[4,4], 2)),
         paste("p =",  round(resCurveCSG$LM_Curvethickness$P.param, 2)),
         paste("Rsq =", round(resCurveCSG$LM_Curvethickness$rsquare, 2))
       ),
       bty= "n"
)


plot(resCurveCSG$LM_CurveIx
     , "RMA", ylab = "res Ix ", xlab = "Curvature ", main = "")
legend("bottomright", 
       legend = c(
         paste("intercept =", round(resCurveCSG$LM_CurveIx$regression.results[4,2], digits = 2)),
         paste("slope =", round(resCurveCSG$LM_CurveIx$regression.results[4,3], digits = 2),
               "±", round(resCurveCSG$LM_CurveIx$confidence.intervals[4,5]-resCurveCSG$LM_CurveIx$confidence.intervals[4,4], 2)),
         paste("p =",  round(resCurveCSG$LM_CurveIx$P.param, 2)),
         paste("Rsq =", round(resCurveCSG$LM_CurveIx$rsquare, 2))
       ),
       bty= "n"
)



plot(resCurveCSG$LM_CurveIy
     , "RMA", ylab = "res Iy ", xlab = "Curvature ", main = "")
legend("bottomright", 
       legend = c(
         paste("intercept =", round(resCurveCSG$LM_CurveIy$regression.results[4,2], digits = 2)),
         paste("slope =", round(resCurveCSG$LM_CurveIy$regression.results[4,3], digits = 2),
               "±", round(resCurveCSG$LM_CurveIy$confidence.intervals[4,5]-resCurveCSG$LM_CurveIy$confidence.intervals[4,4], 2)),
         paste("p =",  round(resCurveCSG$LM_CurveIy$P.param, 2)),
         paste("Rsq =", round(resCurveCSG$LM_CurveIy$rsquare, 2))
       ),
       bty= "n"
)



plot(resCurveCSG$LM_CurveI.ratio
     , "RMA", ylab = "res I-Ratio ", xlab = "Curvature ", main = "")
legend("bottomright", 
       legend = c(
         paste("intercept =", round(resCurveCSG$LM_CurveI.ratio$regression.results[4,2], digits = 2)),
         paste("slope =", round(resCurveCSG$LM_CurveI.ratio$regression.results[4,3], digits = 2),
               "±", round(resCurveCSG$LM_CurveI.ratio$confidence.intervals[4,5]-resCurveCSG$LM_CurveI.ratio$confidence.intervals[4,4], 2)),
         paste("p =",  round(resCurveCSG$LM_CurveI.ratio$P.param, 2)),
         paste("Rsq =", round(resCurveCSG$LM_CurveI.ratio$rsquare, 2))
       ),
       bty= "n"
)


plot(resCurveCSG$LM_CurveZpol
     , "RMA", ylab = "res Zpol ", xlab = "Curvature ", main = "")
legend("bottomright", 
       legend = c(
         paste("intercept =", round(resCurveCSG$LM_CurveZpol$regression.results[4,2], digits = 2)),
         paste("slope =", round(resCurveCSG$LM_CurveZpol$regression.results[4,3], digits = 2),
               "±", round(resCurveCSG$LM_CurveZpol$confidence.intervals[4,5]-resCurveCSG$LM_CurveZpol$confidence.intervals[4,4], 2)),
         paste("p =",  round(resCurveCSG$LM_CurveZpol$P.param, 2)),
         paste("Rsq =", round(resCurveCSG$LM_CurveZpol$rsquare, 2))
       ),
       bty= "n"
)


plot(resCurveCSG$LM_CurveJ
     , "RMA", ylab = "res J ", xlab = "Curvature ", main = "")
legend("bottomright", 
       legend = c(
         paste("intercept =", round(resCurveCSG$LM_CurveJ$regression.results[4,2], digits = 2)),
         paste("slope =", round(resCurveCSG$LM_CurveJ$regression.results[4,3], digits = 2),
               "±", round(resCurveCSG$LM_CurveJ$confidence.intervals[4,5]-resCurveCSG$LM_CurveJ$confidence.intervals[4,4], 2)),
         paste("p =",  round(resCurveCSG$LM_CurveJ$P.param, 2)),
         paste("Rsq =", round(resCurveCSG$LM_CurveJ$rsquare, 2))
       ),
       bty= "n"
)




##### CURVE VS RES VARIANCE

###### CURVE VS ESSENTRICITY/POSITION






