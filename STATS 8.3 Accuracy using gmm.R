rm(list = (ls()))

require(Morpho)


eco <- as.data.frame(read.csv("./DATA/specimenEco2.csv"))

load(file = "./morphomapSnapshots85Procrustes/allCSG85Procrustes_centroid.Rdata")
load(file = "./morphomapSnapshots85Procrustes/allLands85Procrustes.Rdata")
load(file = "./morphomapSnapshots85Procrustes/ShapeSections.Rdata")


allLands



diaphysisLands <- allShapeSections

# coerce into kxpxn array --------------------------------------------------
XX <- diaphysisLands[[1]]$'3D_out'

nslice = dim(XX)[1]
n <- dim(XX)[1] * dim(XX)[3]

l <- length(XX)



inner <- list()
outer <- list()
comb <- list()

ii=1
i=1

for(ii in 1:length(diaphysisLands)){
  X <- diaphysisLands[[ii]]
  
  Xout <- X$'3D_out'
  Xin <- X$'3D_in'
  
  tmpOUT <- NULL
  tmpIN <- NULL
  
  for(i in 1:dim(XX)[3]){
    
    tmpOUT <- rbind(tmpOUT, Xout[,,i])
    tmpIN <- rbind(tmpIN, Xin[,,i])
    
  }
  
  inner[[ii]] <- tmpIN
  outer[[ii]] <- tmpOUT
  comb[[ii]] <- rbind(tmpOUT, tmpIN)
  
  
}

names(inner) <- names(diaphysisLands)
names(outer) <- names(diaphysisLands)
names(comb) <- names(diaphysisLands)


inner <- list2array(inner)
outer <- list2array(outer)
comb <- list2array(comb)

# procrustes --------------------------------------------------------------

procOUTERscaled <- procSym(outer, scale = T, reflect = F)
procCOMBscaled <- procSym(comb, scale = T, reflect = F)
procOUTERunscaled <- procSym(outer, scale = F, reflect = F)
procCOMBunscaled <- procSym(comb, scale = F, reflect = F)


# procOUTERscaled$rotated
# procCOMBscaled$rotated
# procCOMBunscaled$size
# procCOMBscaled$size

gc()


# CVA ---------------------------------------------------------------------

groups <- eco$LOCSUB
# 
# CVAOUTERscaled <- CVA(procOUTERscaled$rotated, cv = F, groups = groups)
# CVACOMBscaled <- CVA(procCOMBscaled$rotated, cv = F, groups = groups)
# CVAOUTERunscaled <- CVA(procOUTERunscaled$rotated, cv = F, groups = groups)
# CVACOMBunscaled <- CVA(procCOMBunscaled$rotated, cv = F, groups = groups)
# 
# save(CVAOUTERscaled,
#      CVACOMBscaled,
#      CVAOUTERunscaled,
#      CVACOMBunscaled,
#      file = "CVAs.Rdata")
# 
# 


# CVAOUTERscaled <- CVA(procOUTERunscaled$rotated, cv = T, groups = groups)
# CVACOMBscaled <- CVA(procCOMBunscaled$rotated, cv = T, groups = groups)
# CVAOUTERunscaled <- CVA(procOUTERunscaled$rotated, cv = T, groups = groups)
# CVACOMBunscaled <- CVA(procCOMBunscaled$rotated, cv = T, groups = groups)
# 
# 
# save(CVAOUTERscaled,
#      CVACOMBscaled,
#      CVAOUTERunscaled,
#      CVACOMBunscaled,
#      file = "CVAscv.Rdata")

# load("CVAscv.Rdata")

# -------------------------------------------------------------------------



load("CVAs.Rdata")

CVAOUTERscaled
CVACOMBscaled
CVAOUTERunscaled
CVACOMBunscaled


# test slices -------------------------------------------------------------

ii=1

names(diaphysisLands) <- names(allCSG)

lengthAlong <- seq(diaphysisLands[[1]]$start,
                   diaphysisLands[[1]]$end,
    length.out = dim(diaphysisLands[[1]]$"3D_out")[3])




CVAslicesOUTER <- list()
CVAslicesCOMB <- list()


i=20
ii=1

for(i in 1:dim(XX)[3]){
  
  sliceinner <- list()
  sliceouter <- list()
  slicecomb <- list()
  
  for(ii in 1:length(diaphysisLands)){
    X <- diaphysisLands[[ii]]
    
    Xout <- X$'2D_out'
    Xin <- X$'2D_in'
    
      
      tmpOUT <- Xout[,,i]
      tmpIN <- Xin[,,i]
      
  
    
    sliceinner[[ii]] <- tmpIN
    sliceouter[[ii]] <- tmpOUT
    slicecomb[[ii]] <- rbind(tmpOUT, tmpIN)
    
    
    
    
  }
  
  sliceinner <- list2array(sliceinner)
  sliceouter <- list2array(sliceouter)
  slicecomb <- list2array(slicecomb)
  
  sliceprocOUTERscaled <- procSym(sliceouter, scale = T, reflect = F)
  sliceprocCOMBscaled <- procSym(slicecomb, scale = T, reflect = F)
  
  sliceCVAOUTERscaled <- CVA(sliceprocOUTERscaled$rotated, cv = F, weighting = F, rounds = 9999, groups = as.factor(groups))
  
  classify(sliceCVAOUTERscaled, cv = F)

  sliceCVACOMBscaled <- CVA(sliceprocCOMBscaled$rotated, cv = F, weighting = F, rounds = 9999, groups = as.factor(groups))
  
  print(classify(sliceCVACOMBscaled, cv = F))$accuracy
  
  CVAslicesOUTER[[i]] <- sliceCVAOUTERscaled
  CVAslicesCOMB[[i]] <- sliceCVACOMBscaled
  
  
}

names(CVAslicesOUTER) <- lengthAlong
names(CVAslicesCOMB) <- lengthAlong


sink("./stats/CVA gmm slices OUTER.txt")
print(CVAslicesOUTER)

sink()
sink("./stats/CVA gmm slices COMB.txt")
print(CVAslicesCOMB)

sink()



SLICESprobsOUTER <- unlist(lapply(CVAslicesOUTER, FUN = function(X) print(classify(X, cv = F))$accuracy))
SLICESprobsCOMB <- unlist(lapply(CVAslicesCOMB, FUN = function(X) print(classify(X, cv = F))$accuracy))

# write.csv()

write.csv(data.frame(outer = SLICESprobsOUTER, comb = SLICESprobsCOMB), file = "cvagmmprobs.csv")


# write.csv(SLICESprobsOUTER, file = "cvagmmprobs.csv")



sliceAcc <- read.csv(file = "./CVA CSG unscaled res/cvaprobs.csv", row.names = 1)

nrow(sliceAcc)

CVAsliceaccuracy <- data.frame(CSG = sliceAcc, GMMOUTER = SLICESprobsOUTER, GMMCOMB = SLICESprobsCOMB)
colMeans(CVAsliceaccuracy)

colnames(CVAsliceaccuracy)[1] <- "CSG"


# -------------------------------------------------------------------------



par(mfrow = c(1,1), mar = c(5,4,4,2)+0.1)

plot(unique(lengthAlong), CVAsliceaccuracy$CSG, type = "n", ylim = range(CVAsliceaccuracy),
     main = "Classification accuracy of cortical CSG along the diaphysis",
     xlab = "% diaphysis", ylab = "classification accuracy (%)")

palette.colors()

lines(unique(lengthAlong), CVAsliceaccuracy$CSG, pch = 21)
points(unique(lengthAlong), CVAsliceaccuracy$CSG, pch = 21, bg = "black")
lines(unique(lengthAlong), CVAsliceaccuracy$GMMOUTER)
points(unique(lengthAlong), CVAsliceaccuracy$GMMOUTER, pch = 24, bg = "white")
lines(unique(lengthAlong), CVAsliceaccuracy$GMMCOMB, pch = 25)
points(unique(lengthAlong), CVAsliceaccuracy$GMMCOMB, pch = 25, bg = "grey")

legend("right", pch = c(25,24,21), cex = 0.8, pt.bg = c( "grey", "white", "black"), legend = c("GMM combined", "GMM outer", "CSG unscl"))


# load CSG data -----------------------------------------------------------







# visualize ---------------------------------------------------------------


sink("./stats/CVA gmm.txt")

print("CVAOUTERscaled")
CVAOUTERscaled
print("CVAOUTERscaled")
CVACOMBscaled
print("CVAOUTERunscaled")
CVAOUTERunscaled
print("CVACOMBunscaled")
CVACOMBunscaled
sink()



# -------------------------------------------------------------------------



palette(palette.colors()[2:5])

par(mfrow = c(2,2))

plot(CVAOUTERscaled$CVscores, col=as.factor(groups), pch=as.numeric(as.factor(groups))+14, typ="p",asp=1,
     xlab=paste("1st canonical axis", paste(round(CVAOUTERscaled$Var[1,2],1),"%")),
     ylab=paste("2nd canonical axis", paste(round(CVAOUTERscaled$Var[2,2],1),"%")),
     main= "CVA using OUTER landmarks")
# text(CVAOUTERscaled$CVscores, as.character(data[,1]), col=as.numeric(data[,1]), cex=.7)

for(i in 1:length(levels(CVAOUTERscaled$groups))){
  X <- CVAOUTERscaled$CVscores[CVAOUTERscaled$groups == levels(CVAOUTERscaled$groups)[i],]
  hpts <- chull(X)
  hpts <- c(hpts, hpts[1])
  lines(X[hpts, ],col = palette()[i])
  
}

a <- print(CVAOUTERscaled)

a$accuracy


legend("topright", legend = levels(CVAOUTERscaled$groups), cex = 0.8, pch = (1:4)+14, col = palette()[1:4])
legend("bottomright", legend = paste("Accuracy =", a$accuracy, "%"), bty = "n")



plot(CVAOUTERscaled$CVscores[,3:2], col=as.factor(groups), pch=as.numeric(as.factor(groups))+14, typ="p",asp=1,
     xlab=paste("3rd canonical axis", paste(round(CVAOUTERscaled$Var[3,2],1),"%")),
     ylab=paste("2nd canonical axis", paste(round(CVAOUTERscaled$Var[2,2],1),"%")),
     main= "CVA using OUTER landmarks")
# text(CVAOUTERscaled$CVscores, as.character(data[,1]), col=as.numeric(data[,1]), cex=.7)

for(i in 1:length(levels(CVAOUTERscaled$groups))){
  X <- CVAOUTERscaled$CVscores[CVAOUTERscaled$groups == levels(CVAOUTERscaled$groups)[i],3:2]
  hpts <- chull(X)
  hpts <- c(hpts, hpts[1])
  lines(X[hpts, ],col = palette()[i])
  
}



plot(CVACOMBscaled$CVscores, col=as.factor(groups), pch=as.numeric(as.factor(groups))+14, typ="p",asp=1,
     xlab=paste("1st canonical axis", paste(round(CVACOMBscaled$Var[1,2],1),"%")),
     ylab=paste("2nd canonical axis", paste(round(CVACOMBscaled$Var[2,2],1),"%")),
     main= "CVA using COMB landmarks")
# text(CVACOMBscaled$CVscores, as.character(data[,1]), col=as.numeric(data[,1]), cex=.7)

for(i in 1:length(levels(CVACOMBscaled$groups))){
  X <- CVACOMBscaled$CVscores[CVACOMBscaled$groups == levels(CVACOMBscaled$groups)[i],]
  hpts <- chull(X)
  hpts <- c(hpts, hpts[1])
  lines(X[hpts, ],col = palette()[i])
  
}

a <- print(CVACOMBscaled)

a$accuracy


# legend("topright", legend = levels(CVACOMBscaled$groups), cex = 0.8, pch = (1:4)+14, col = palette()[1:4])
legend("bottomright", legend = paste("Accuracy =", a$accuracy, "%"), bty = "n")



plot(CVACOMBscaled$CVscores[,3:2], col=as.factor(groups), pch=as.numeric(as.factor(groups))+14, typ="p",asp=1,
     xlab=paste("3rd canonical axis", paste(round(CVACOMBscaled$Var[3,2],1),"%")),
     ylab=paste("2nd canonical axis", paste(round(CVACOMBscaled$Var[2,2],1),"%")),
     main= "CVA using COMB landmarks")
# text(CVACOMBscaled$CVscores, as.character(data[,1]), col=as.numeric(data[,1]), cex=.7)

for(i in 1:length(levels(CVACOMBscaled$groups))){
  X <- CVACOMBscaled$CVscores[CVACOMBscaled$groups == levels(CVACOMBscaled$groups)[i],3:2]
  hpts <- chull(X)
  hpts <- c(hpts, hpts[1])
  lines(X[hpts, ],col = palette()[i])
  
}




# -------------------------------------------------------------------------

layout(rbind(c(1,2)), 
       widths = c(1,dist(range(CVACOMBscaled$CVscores[,3]))/
                    dist(range(CVACOMBscaled$CVscores[,1]))
))

par(mar = c(5, 4, 4, 0) + 0.1)



plot(CVACOMBscaled$CVscores[,1:2], col=as.factor(groups), pch=as.numeric(as.factor(groups))+14, 
     xlab=paste("1st canonical axis", paste(round(CVACOMBscaled$Var[1,2],1),"%")),
     ylab=paste("2nd canonical axis", paste(round(CVACOMBscaled$Var[2,2],1),"%")),
     ylim = range(CVACOMBscaled$CVscores[,2]), 
     axes = F, 
     type = "n", 
     main= "CVA using COMB landmarks")


box()
axis(1, at = seq(-8, 8, by = 2))
axis(2, at = seq(-10, 10, by = 2))

abline(h = 0, lty = 5, col = "grey")
abline(v = 0, lty = 5, col = "grey")



# text(CVACOMBscaled$CVscores, as.character(data[,1]), col=as.numeric(data[,1]), cex=.7)


for(i in 1:length(levels(CVACOMBscaled$groups))){
  X <- CVACOMBscaled$CVscores[CVACOMBscaled$groups == levels(CVACOMBscaled$groups)[i],]
  hpts <- chull(X)
  hpts <- c(hpts, hpts[1])
  lines(X[hpts, ],col = palette()[i])
  
}

points(CVACOMBscaled$CVscores[,1:2], col=as.factor(groups), pch=as.numeric(as.factor(groups))+14,
       xlab=paste("1st canonical axis", paste(round(CVACOMBscaled$Var[1,2],1),"%")),
       ylab=paste("2nd canonical axis", paste(round(CVACOMBscaled$Var[2,2],1),"%")),
       ylim = range(CVACOMBscaled$CVscores[,2]), 
       # axes = F, 
       # type = "n", 
       main= "CVA using COMB landmarks")

# a <- print(CVACOMBscaled)

# a$accuracy

legend("topright", legend = levels(CVAOUTERscaled$groups), cex = 0.8, pch = (1:4)+14, col = palette()[1:4])
# legend("bottomright", legend = paste("Accuracy =", a$accuracy, "%"), bty = "n")

par(mar = c(5, 0, 4, 1) + 0.1)


plot(CVACOMBscaled$CVscores[,3:2], col=as.factor(groups), pch=as.numeric(as.factor(groups))+14, 
     type="n",
     xlab=paste("3rd canonical axis", paste(round(CVACOMBscaled$Var[3,2],1),"%")),
     ylab="",
     ylim = range(CVACOMBscaled$CVscores[,2]), 
     axes = F)

box()
axis(1, at = seq(-8, 8, by = 2))
# axis(2, at = seq(-10, 10, by = 2))

abline(h = 0, lty = 5, col = "grey")
abline(v = 0, lty = 5, col = "grey")



# axis(2)
# text(CVACOMBscaled$CVscores, as.character(data[,1]), col=as.numeric(data[,1]), cex=.7)

for(i in 1:length(levels(CVACOMBscaled$groups))){
  X <- CVACOMBscaled$CVscores[CVACOMBscaled$groups == levels(CVACOMBscaled$groups)[i],3:2]
  hpts <- chull(X)
  hpts <- c(hpts, hpts[1])
  lines(X[hpts, ],col = palette()[i])
  
}

points(CVACOMBscaled$CVscores[,3:2], col=as.factor(groups), pch=as.numeric(as.factor(groups))+14, 
       type="p",
       xlab=paste("3rd canonical axis", paste(round(CVACOMBscaled$Var[3,2],1),"%")),
       ylab="",
       ylim = range(CVACOMBscaled$CVscores[,2]), 
)


# -------------------------------------------------------------------------


