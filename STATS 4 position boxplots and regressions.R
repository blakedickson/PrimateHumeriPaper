eco <- as.data.frame(read.csv("./DATA/specimenEco2.csv"))

load("./morphomapSnapshots85Procrustes/allCSG85Procrustes_centroid.Rdata")
load("./morphomapSnapshots85Procrustes/allLands85Procrustes.Rdata")


load("boneVars.Rdata")


data <- data.frame(LOCSUB = eco$LOCSUB, boneVars)

LOCSUB <- eco$LOCSUB

boneVars$posDF

data$posDF.pos.CA

# -------------------------------------------------------------------------

require(ggplot2)


par(mfrow = c(2,3))


# data <- data.frame(LOCSUB = eco$LOCSUB, boneVars$posDF)

data <- data.frame(LOCSUB = eco$LOCSUB, boneVars$maxDF, boneVars$posDF, boneVars$varDF)
par(mfcol = c(8,2))

which(eco$LOCSUB == levels(as.factor(eco$LOCSUB))[1])
boneVars$ecc50DF

mean(as.matrix(boneVars$abs50DF))

colMeans(boneVars$abs50DF[which(eco$LOCSUB == levels(as.factor(eco$LOCSUB))[1]),])
colMeans(boneVars$abs50DF[which(eco$LOCSUB == levels(as.factor(eco$LOCSUB))[1]),])
colMeans(boneVars$abs50DF[which(eco$LOCSUB == levels(as.factor(eco$LOCSUB))[1]),])
colMeans(boneVars$abs50DF[which(eco$LOCSUB == levels(as.factor(eco$LOCSUB))[1]),])

boneVars$posDF
colMeans(boneVars$posDF[,])



colMeans(boneVars$posDF[which(eco$LOCSUB == levels(as.factor(eco$LOCSUB))[1]),])
colMeans(boneVars$posDF[which(eco$LOCSUB == levels(as.factor(eco$LOCSUB))[2]),])
colMeans(boneVars$posDF[which(eco$LOCSUB == levels(as.factor(eco$LOCSUB))[3]),])
colMeans(boneVars$posDF[which(eco$LOCSUB == levels(as.factor(eco$LOCSUB))[4]),])




# position ----------------------------------------------------------------


par(mfrow= c(4,2),mar= c(2, 10, 2, 2) + 0.1)


## POSTITION
maxCurveAtLoc1 <- aov(pos.Curvature~LOCSUB, data = data)
s <- summary(maxCurveAtLoc1)
boxplot(pos.Curvature~LOCSUB, data = data, horizontal =T, ylab = "",
        ylim = c(0.15,0.85), las = 1, col = palette.colors()[2:5],
        main = paste("Curve:", paste("p = ", round(s[[1]]$`Pr(>F)`[1], 3))))


legend("bottomright", bty = "n",
       legend = paste("p = ", round(s[[1]]$`Pr(>F)`[1], 3)))
# legend("topleft", bty = "n", legend = paste("max Curve Position:", paste("p = ", round(s[[1]]$`Pr(>F)`[1], 3))))
summary(maxCurveAtLoc1)


data$pos.CA

pos.CA <- aov(pos.CA~LOCSUB, data = data)
s <- summary(pos.CA)
boxplot(pos.CA~LOCSUB, data = data, horizontal =T, ylab = "",
        ylim = c(0.15,0.85), las = 1, col = palette.colors()[2:5],
        main = paste("CA:", paste("p = ", round(s[[1]]$`Pr(>F)`[1], 3))))


# legend("topleft", bty = "n", legend = "max CA Position")
# legend("topleft", bty = "n", legend = paste("max CA Position:", paste("p = ", round(s[[1]]$`Pr(>F)`[1], 3))))

# legend("bottomright", bty = "n",
#        legend = paste("p = ", round(s[[1]]$`Pr(>F)`[1], 3)))
summary(maxCurveAtLoc1)

data$pos.thickness

pos.thickness <- aov(pos.thickness~LOCSUB, data = data)
s <- summary(pos.thickness)
boxplot(pos.thickness~LOCSUB, data = data, horizontal =T, ylab = "",
        ylim = c(0.15,0.85), las = 1, col = palette.colors()[2:5],
        main = paste("Thick.:", paste("p = ", round(s[[1]]$`Pr(>F)`[1], 3))))


# legend("topleft", bty = "n", legend = "max Thick. Position")
# legend("topleft", bty = "n", legend = paste("max Thick. Position:", paste("p = ", round(s[[1]]$`Pr(>F)`[1], 3))))


# legend("bottomright", bty = "n",
#        legend = paste("p = ", round(s[[1]]$`Pr(>F)`[1], 3)))
summary(maxCurveAtLoc1)


data$pos.Ix

pos.Ix <- aov(pos.Ix~LOCSUB, data = data)
s <- summary(pos.Ix)
boxplot(pos.Ix~LOCSUB, data = data, horizontal =T, ylab = "",
        ylim = c(0.15,0.85), las = 1, col = palette.colors()[2:5],
        main = paste("Ix:", paste("p = ", round(s[[1]]$`Pr(>F)`[1], 3))))

# legend("topleft", bty = "n", legend = paste("max Ix Position:", paste("p = ", round(s[[1]]$`Pr(>F)`[1], 3))))

# 
# legend("bottomright", bty = "n",
#        legend = paste("p = ", round(s[[1]]$`Pr(>F)`[1], 3)))

data$pos.Iy



pos.Iy <- aov(pos.Iy~LOCSUB, data = data)
s <- summary(pos.Iy)
boxplot(pos.Iy~LOCSUB, data = data, horizontal =T, ylab = "",
        ylim = c(0.15,0.85), las = 1, col = palette.colors()[2:5],
        main = paste("Iy:", paste("p = ", round(s[[1]]$`Pr(>F)`[1], 3))))



# legend("topleft", bty = "n", legend = "max Iy Position")
# legend("topleft", bty = "n", legend = paste("max Iy Position:", paste("p = ", round(s[[1]]$`Pr(>F)`[1], 3))))

# legend("bottomright", bty = "n",
#        legend = paste("p = ", round(s[[1]]$`Pr(>F)`[1], 3)))



data$pos.I.ratio
pos.I.ratio <- aov(pos.I.ratio~LOCSUB, data = data)
s <- summary(pos.I.ratio)
boxplot(pos.I.ratio~LOCSUB, data = data, horizontal =T, ylab = "",
        ylim = c(0.15,0.85), las = 1, col = palette.colors()[2:5],
        main = paste("I.ratio:", paste("p = ", round(s[[1]]$`Pr(>F)`[1], 3))))



# legend("topleft", bty = "n", legend = "max I.ratio Position")
# legend("topleft", bty = "n", legend = paste("max I.ratio Position:", paste("p = ", round(s[[1]]$`Pr(>F)`[1], 3))))

# legend("bottomright", bty = "n",
#        legend = paste("p = ", round(s[[1]]$`Pr(>F)`[1], 3)))

data$pos.Zpol
pos.Zpol <- aov(pos.Zpol~LOCSUB, data = data)
s <- summary(pos.Zpol)
boxplot(pos.Zpol~LOCSUB, data = data, horizontal =T, ylab = "",
        ylim = c(0.15,0.85), las = 1, col = palette.colors()[2:5],
        main = paste("Zpol:", paste("p = ", round(s[[1]]$`Pr(>F)`[1], 3))))



# legend("topleft", bty = "n", legend = "max Zpol Position")
# legend("topleft", bty = "n", legend = paste("max Zpol Position:", paste("p = ", round(s[[1]]$`Pr(>F)`[1], 3))))

# legend("bottomright", bty = "n",
#        legend = paste("p = ", round(s[[1]]$`Pr(>F)`[1], 3)))





data$pos.J
pos.J <- aov(pos.J~LOCSUB, data = data)
s <- summary(pos.J)
boxplot(pos.J~LOCSUB, data = data, horizontal =T, ylab = "",
        ylim = c(0.15,0.85), las = 1, col = palette.colors()[2:5],
        main = paste("J:", paste("p = ", round(s[[1]]$`Pr(>F)`[1], 3))))



# legend("topleft", bty = "n", legend = "max J Position")
# legend("topleft", bty = "n", legend = paste("max J Position:", paste("p = ", round(s[[1]]$`Pr(>F)`[1], 3))))

# legend("bottomright", bty = "n",
#        legend = paste("p = ", round(s[[1]]$`Pr(>F)`[1], 3)))



# pairwise comps ----------------------------------------------------------

# ??posthoc

LOCSUB <- as.factor(LOCSUB)
data <- data.frame(boneVars$posDF)

comp <- expand.grid(levels(LOCSUB), levels(LOCSUB))



compMAT <- matrix(NA,4,4)
dimnames(compMAT) <- list(levels(LOCSUB),levels(LOCSUB))

d=1
i=1

j=2

compLIST<- list()

for(d in 1:ncol(data)){
  
  compMAT <- matrix(NA,4,4)
  dimnames(compMAT) <- list(levels(LOCSUB),levels(LOCSUB))
  
  for(i in 1:length(levels(LOCSUB))){
    for(j in 1:length(levels(LOCSUB))){
      
      X <- data[which(eco$LOCSUB == levels(as.factor(eco$LOCSUB))[i]),d]
      Y <- data[which(eco$LOCSUB == levels(as.factor(eco$LOCSUB))[j]),d]
      compMAT[i,j] <- wilcox.test(X, Y)$p.value
      
      
    }
  }
  
  which(compMAT== 0)
  
  compLIST[[d]] <- compMAT
  
  
}


names(compLIST) <- colnames(data)

sink("./stats/pairwiseMaxPos.txt")
lapply(compLIST, round, 3)
sink()

compLIST
# AOV and pairwise comps --------------------------------------------------

pos.C <- aov(pos.Curvature~LOCSUB, data = data)
summary(pos.C)

pos.CA <- aov(pos.CA~LOCSUB, data = data)
summary(pos.CA)

pos.thickness <- aov(pos.thickness~LOCSUB, data = data)
summary(pos.thickness)

pos.Ix <- aov(pos.Ix~LOCSUB, data = data)
summary(pos.Ix)

pos.Iy <- aov(pos.Iy~LOCSUB, data = data)
summary(pos.Iy)

pos.I.ratio <- aov(pos.I.ratio~LOCSUB, data = data)
summary(pos.I.ratio)

pos.Zpol <- aov(pos.Zpol~LOCSUB, data = data)
summary(pos.Zpol)

pos.J <- aov(pos.J~LOCSUB, data = data)
summary(pos.J)




# positions v postiions ---------------------------------------------------------------

library(lmodel2)



dat <- boneVars$posDF

dat <- cbind(boneVars$ecoDF, boneVars$ecc50DF, boneVars$abs50DF ,boneVars$maxDF, boneVars$varDF, boneVars$meanDF, boneVars$posDF)

dat$pos.Curvature



# ### CORTICAL AREA

LM_posCurvePosCA <- lmodel2((pos.CA)~(pos.Curvature), "interval", "interval",data = dat, 99)
LM_posCurvePosCA

# ### THICKNESS

LM_posCurvePosthickness <- lmodel2((pos.thickness)~(pos.Curvature), "interval", "interval",data = dat, 99)
LM_posCurvePosthickness

# ### Ix

LM_posCurvePosIx <- lmodel2((pos.Ix)~(pos.Curvature), "interval", "interval",data = dat, 99)
LM_posCurvePosIx

# ### Iy

LM_posCurvePosIy <- lmodel2((pos.Iy)~(pos.Curvature), "interval", "interval",data = dat, 99)
LM_posCurvePosIy

# ### I.ratio

LM_posCurvePosI.ratio <- lmodel2((pos.I.ratio)~(pos.Curvature), "interval", "interval",data = dat, 99)
LM_posCurvePosI.ratio

# ### Zpol

LM_posCurvePosZpol <- lmodel2((pos.Zpol)~(pos.Curvature), "interval", "interval",data = dat, 99)
LM_posCurvePosZpol

# ### pos.J

LM_posCurvePosJ <- lmodel2((pos.J)~(pos.Curvature), "interval", "interval",data = dat, 99)
LM_posCurvePosJ

CurveposCSGpos <- list(LM_CurveCA = LM_posCurvePosCA, 
                       LM_Curvethickness= LM_posCurvePosthickness, 
                       LM_CurveIx = LM_posCurvePosIx, 
                       LM_CurveIy = LM_posCurvePosIy, 
                       LM_CurveI.ratio = LM_posCurvePosI.ratio, 
                       LM_CurveZpol = LM_posCurvePosZpol, 
                       LM_CurveJ = LM_posCurvePosJ)

CurveposCSGpos



sink("./stats/LMs pos Curvature vs pos CSG.txt")
CurveposCSGpos
sink()




# plot --------------------------------------------------------------------



par(mfrow = c(2,4),mar= c(1, 4, 1, 2) + 0.1)


# plot.new()


plot(CurveposCSGpos$LM_CurveCA
     , "RMA", ylab = "max CA pos", xlab = "", main = "")
legend("bottomright", 
       legend = c(
         paste("intercept =", round(CurveposCSGpos$LM_CurveCA$regression.results[4,2], digits = 2)),
         paste("slope =", round(CurveposCSGpos$LM_CurveCA$regression.results[4,3], digits = 2),
               "±", round(CurveposCSGpos$LM_CurveCA$confidence.intervals[4,5]-CurveposCSGpos$LM_CurveCA$confidence.intervals[4,4], 2)),
         paste("p =",  round(CurveposCSGpos$LM_CurveCA$P.param, 2)),
         paste("Rsq =", round(CurveposCSGpos$LM_CurveCA$rsquare, 2))
       ),
       bty= "n"
)

plot(CurveposCSGpos$LM_Curvethickness
     , "RMA", ylab = "max Thickness pos", xlab = "", main = "")
legend("bottomright", 
       legend = c(
         paste("intercept =", round(CurveposCSGpos$LM_Curvethickness$regression.results[4,2], digits = 2)),
         paste("slope =", round(CurveposCSGpos$LM_Curvethickness$regression.results[4,3], digits = 2),
               "±", round(CurveposCSGpos$LM_Curvethickness$confidence.intervals[4,5]-CurveposCSGpos$LM_Curvethickness$confidence.intervals[4,4], 2)),
         paste("p =",  round(CurveposCSGpos$LM_Curvethickness$P.param, 2)),
         paste("Rsq =", round(CurveposCSGpos$LM_Curvethickness$rsquare, 2))
       ),
       bty= "n"
)


plot(CurveposCSGpos$LM_CurveIx
     , "RMA", ylab = "max Ix pos", xlab = "", main = "")
legend("bottomright", 
       legend = c(
         paste("intercept =", round(CurveposCSGpos$LM_CurveIx$regression.results[4,2], digits = 2)),
         paste("slope =", round(CurveposCSGpos$LM_CurveIx$regression.results[4,3], digits = 2),
               "±", round(CurveposCSGpos$LM_CurveIx$confidence.intervals[4,5]-CurveposCSGpos$LM_CurveIx$confidence.intervals[4,4], 2)),
         paste("p =",  round(CurveposCSGpos$LM_CurveIx$P.param, 2)),
         paste("Rsq =", round(CurveposCSGpos$LM_CurveIx$rsquare, 2))
       ),
       bty= "n"
)


plot(CurveposCSGpos$LM_CurveIy
     , "RMA", ylab = "max Iy pos", xlab = "", main = "")
legend("bottomright", 
       legend = c(
         paste("intercept =", round(CurveposCSGpos$LM_CurveIy$regression.results[4,2], digits = 2)),
         paste("slope =", round(CurveposCSGpos$LM_CurveIy$regression.results[4,3], digits = 2),
               "±", round(CurveposCSGpos$LM_CurveIy$confidence.intervals[4,5]-CurveposCSGpos$LM_CurveIy$confidence.intervals[4,4], 2)),
         paste("p =",  round(CurveposCSGpos$LM_CurveIy$P.param, 2)),
         paste("Rsq =", round(CurveposCSGpos$LM_CurveIy$rsquare, 2))
       ),
       bty= "n"
)
plot(CurveposCSGpos$LM_CurveI.ratio
     , "RMA", ylab = "max I-Ratio pos", xlab = "", main = "")
legend("bottomright", 
       legend = c(
         paste("intercept =", round(CurveposCSGpos$LM_CurveI.ratio$regression.results[4,2], digits = 2)),
         paste("slope =", round(CurveposCSGpos$LM_CurveI.ratio$regression.results[4,3], digits = 2),
               "±", round(CurveposCSGpos$LM_CurveI.ratio$confidence.intervals[4,5]-CurveposCSGpos$LM_CurveI.ratio$confidence.intervals[4,4], 2)),
         paste("p =",  round(CurveposCSGpos$LM_CurveI.ratio$P.param, 2)),
         paste("Rsq =", round(CurveposCSGpos$LM_CurveI.ratio$rsquare, 2))
       ),
       bty= "n"
)


plot(CurveposCSGpos$LM_CurveZpol
     , "RMA", ylab = "max Zpol pos", xlab = "", main = "")
legend("bottomright", 
       legend = c(
         paste("intercept =", round(CurveposCSGpos$LM_CurveZpol$regression.results[4,2], digits = 2)),
         paste("slope =", round(CurveposCSGpos$LM_CurveZpol$regression.results[4,3], digits = 2),
               "±", round(CurveposCSGpos$LM_CurveZpol$confidence.intervals[4,5]-CurveposCSGpos$LM_CurveZpol$confidence.intervals[4,4], 2)),
         paste("p =",  round(CurveposCSGpos$LM_CurveZpol$P.param, 2)),
         paste("Rsq =", round(CurveposCSGpos$LM_CurveZpol$rsquare, 2))
       ),
       bty= "n"
)


plot(CurveposCSGpos$LM_CurveJ
     , "RMA", ylab = "max J pos", xlab = "", main = "")
legend("bottomright", 
       legend = c(
         paste("intercept =", round(CurveposCSGpos$LM_CurveJ$regression.results[4,2], digits = 2)),
         paste("slope =", round(CurveposCSGpos$LM_CurveJ$regression.results[4,3], digits = 2),
               "±", round(CurveposCSGpos$LM_CurveJ$confidence.intervals[4,5]-CurveposCSGpos$LM_CurveJ$confidence.intervals[4,4], 2)),
         paste("p =",  round(CurveposCSGpos$LM_CurveJ$P.param, 2)),
         paste("Rsq =", round(CurveposCSGpos$LM_CurveJ$rsquare, 2))
       ),
       bty= "n"
)













