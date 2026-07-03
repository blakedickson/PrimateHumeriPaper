rm(list = (ls()))

eco <- as.data.frame(read.csv("./DATA/specimenEco2.csv"))

groups <- eco$LOCSUB

# dir.create("./CVA CSG scaled body res/")

load("allDATA.Rdata")
load("boneVars.Rdata")

load("res_data.Rdata")

##########################################################################
# ECOLOGY PREDICTIONS -----------------------------------------------------
###########################################################################

boneVars$ecoDF
lengthAlong = allDATA$lengthAlong

boneVars$ecoDF

lengthAlong = allDATA$lengthAlong

midCirc <- unlist(allDATA$MidShaft.Ext_perim)

ALLdat <- data.frame(Curvature = allDATA$DistCentroid,
                     thickness = allDATA$Mean_thick,
                     CA = allDATA$CA,
                     Ix = (allDATA$Ix),
                     Iy = (allDATA$Iy),
                     I.ratio = (allDATA$Imax/allDATA$Imin),
                     Zpol = allDATA$Zpol.Zpol,
                     J = allDATA$J.Ix
)

DF <- ALLdat[lengthAlong == 0.5,]

ECO <- boneVars$ecoDF[,c(3:5)]


length(unique(lengthAlong))

MASS <- unlist(lapply(eco$MASS..KG., rep, length(unique(lengthAlong))))

SIZES <- boneVars$ecoDF[,c(6:7)]

data <- data.frame(LOCSUB = as.factor(paste(ECO$LOCOMOTION, ECO$SUBSTRATE)), DF, mass = eco$MASS..KG.)
dat <- data.frame(LOCSUB = as.factor(paste(ECO$LOCOMOTION, ECO$SUBSTRATE)), DF, mass = eco$MASS..KG.)

require(Morpho)

require(MASS)

par(mfrow= c(1,1), mar = c(5.1,4.1, 4.1, 2.1))
b <- barplot(summary(data$LOCSUB), col = palette.colors())
text(b, summary(data$LOCSUB) - 2 , summary(data$LOCSUB), font=2, col=c("grey", "black", "black", "black"))


# export pub data file ----------------------------------------------------


OUT <- data.frame(ECO = LOCSUB,
           MIDCIRC = midCirc,
           lengthAlong = lengthAlong,
           "Body Mass" = MASS,
           Curvature = allDATA$DistCentroid,
           thickness = allDATA$Mean_thick,
           CA = allDATA$CA,
           Ix = (allDATA$Ix),
           Iy = (allDATA$Iy),
           I.ratio = (allDATA$Imax/allDATA$Imin),
           Zpol = allDATA$Zpol.Zpol,
           J = allDATA$J.Ix
)

rownames(OUT) <- rownames(boneVars[[1]]$ID)



write.csv(OUT, file = "../0MANUSCRIPT/CSGdata.csv")
write.csv(eco, file = "../0MANUSCRIPT/ECOdata.csv")





# CVA USING JUST MIDSHAFT -------------------------------------------------

lengthAlong = allDATA$lengthAlong

dat <- data.frame(Curvature = allDATA$DistCentroid,
                  thickness = allDATA$Mean_thick,
                  CA = allDATA$CA,
                  I.x = (allDATA$Ix),
                  I.y = (allDATA$Iy),
                  I.ratio = (allDATA$Imax/allDATA$Imin),
                  Zpol = allDATA$Zpol.Zpol,
                  J = allDATA$J.Ix
)





rawdat <- data.frame(Curvature = allDATA$DistCentroid/(MASS^(1/3)),
                     thickness = allDATA$Mean_thick/(MASS^(1/3)),
                     CA = allDATA$CA/(MASS^(2/3)),
                     Ix = (allDATA$Ix)/(MASS^(4/3)),
                     Iy = (allDATA$Iy)/(MASS^(4/3)),
                     I.ratio = (allDATA$Imax/allDATA$Imin),
                     Zpol = allDATA$Zpol.Zpol/(MASS^(3/3)),
                     J = allDATA$J.Ix/(MASS^(4/3))
)



# calc and export eco group means -----------------------------------------

rawdat

LOCSUB <- as.factor(paste(allDATA$LOCOMOTION, allDATA$SUBSTRATE))


write
outMAT <- matrix(NA, ncol = 8, nrow = 4)

colnames(outMAT) <- colnames(rawdat)
rownames(outMAT) <- levels(LOCSUB)

pos = c(0.15, 0.5, 0.85)

i=1
ii=1

outLIST <- list()


for(ii in 1:3){
  
  MEANs <- outMAT
  SDs <- outMAT
  
  for(i in 1:4){
  
    INDD <- (LOCSUB == levels(LOCSUB)[i]) & (lengthAlong == pos[ii])
    
    which(INDD)
    
    tmpdat<- rawdat[which(INDD),]
    
    MEANs[i,] <- apply(tmpdat, 2, mean)
    SDs[i,] <- apply(tmpdat, 2, sd)
    
  }  
  
  outLIST[[ii]] <- list(MEANs = MEANs, SDs = SDs)
  
  names(outLIST)[ii] <- pos[ii]
}

outLIST

sink("./scaledCSGProxMidDist.txt")
outLIST
sink()

rownames(outMAT)<- levels(LOCSUB)


 
LOCSUB


# Calculate at 50% -------------------------------------------------------------------------

ALLdat <- rawdat
ALLdat <- scale(ALLdat, center = F)

REVdat = t(apply(ALLdat, 1, function(r)r*attr(ALLdat,'scaled:scale')))

round(REVdat - rawdat, 3)


identical(ALLdat*attr(ALLdat, "scaled:scale"), rawdat)

DF <- ALLdat[lengthAlong == 0.5,]

ECO <- boneVars$ecoDF[,c(3:5)]

LOCSUB = as.factor(paste(ECO$LOCOMOTION, ECO$SUBSTRATE))
SIZES <- boneVars$ecoDF[,c(6:7)]


# data <- data.frame(LOCSUB = as.factor(paste(ECO$LOCOMOTION, ECO$SUBSTRATE)), DF, SIZES$MIDSHAFTCIRC )

require(Morpho)


groups <- eco$LOCSUB

CVAmidshaft50 <- Morpho::CVA(DF, groups = LOCSUB, weighting = F, cv = F, rounds = 9999)
CVAmidshaft50$CV
CVAmidshaft50$groupmeans
CVAmidshaft50$CVvis


par(mfrow = c(1,2))


sink("./CVA CSG scaled body res/CVAmidshaft50.txt")
CVAmidshaft50
sink()


par(mfrow = c(2,2))
plot(CVAmidshaft50$CVscores, col=data[,1], pch=as.numeric(data[,1])+14, typ="p",asp=1,
     xlab=paste("1st canonical axis", paste(round(CVAmidshaft50$Var[1,2],1),"%")),
     ylab=paste("2nd canonical axis", paste(round(CVAmidshaft50$Var[2,2],1),"%")),
     main = "CVA using midshaft50 values")

for(i in 1:length(levels(CVAmidshaft50$groups))){
  X <- CVAmidshaft50$CVscores[CVAmidshaft50$groups == levels(CVAmidshaft50$groups)[i],]
  # X <- CVAmidshaft50$Scores[CVAmidshaft50$groups == levels(CVAmidshaft50$groups)[i],]
  hpts <- chull(X)
  hpts <- c(hpts, hpts[1])
  lines(X[hpts, ],col = palette()[i])
  
}


acc <- print(CVAmidshaft50)$accuracy


legend("topleft", legend = levels(CVAmidshaft50$groups), pch = (1:4)+14, cex = 0.7,col = palette()[1:4])
legend("bottomright", legend = paste("Accuracy = ", round(acc,2), "%", sep = ""), bty = "n")



plot(CVAmidshaft50$CVscores[,c(1,3)], col=data[,1], pch=as.numeric(data[,1])+14, typ="p",asp=1,
     xlab=paste("1st canonical axis", paste(round(CVAmidshaft50$Var[1,2],1),"%")),
     ylab=paste("3rd canonical axis", paste(round(CVAmidshaft50$Var[3,2],1),"%")))


for(i in 1:length(levels(CVAmidshaft50$groups))){
  X <- CVAmidshaft50$CVscores[CVAmidshaft50$groups == levels(CVAmidshaft50$groups)[i],c(1,3)]
  # X <- CVAmidshaft50$Scores[CVAmidshaft50$groups == levels(CVAmidshaft50$groups)[i],c(1,3)]
  hpts <- chull(X)
  hpts <- c(hpts, hpts[1])
  lines(X[hpts, ],col = palette()[i])
  
}


# Calculate at 40% -------------------------------------------------------------------------

lengthAlong = allDATA$lengthAlong


seq(0,1, by = 0.005) == 0.4

lengthAlong <- as.numeric(as.character(lengthAlong)) 



str(lengthAlong)

DF <- ALLdat[lengthAlong == 0.400,]

ECO <- boneVars$ecoDF[,c(3:5)]
LOCSUB
SIZES <- boneVars$ecoDF[,c(6:7)]



CVAmidshaft40 <- Morpho::CVA(DF, groups = LOCSUB, weighting = F, cv = F, rounds = 9999)
# CVAmidshaft40 <- Morpho::groupPCA(data[,-1], groups = data[,1], weighting = F,cv = T, rounds = 9999)
CVAmidshaft40
CVAmidshaft40$CV
CVAmidshaft40$groupmeans
CVAmidshaft40$CVvis


CVAmidshaft40
CVAmidshaft50


sink("./CVA CSG scaled body res/CVAmidshaft40.txt")
CVAmidshaft40
sink()


palette(palette.colors())

# par(mfrow = c(2,1))
plot(CVAmidshaft40$CVscores, col=data[,1], pch=as.numeric(data[,1])+14, typ="p",asp=1,
     xlab=paste("1st canonical axis", paste(round(CVAmidshaft40$Var[1,2],1),"%")),
     ylab=paste("2nd canonical axis", paste(round(CVAmidshaft40$Var[2,2],1),"%")),
     main = "CVA using midshaft40 values")


for(i in 1:length(levels(CVAmidshaft40$groups))){
  X <- CVAmidshaft40$CVscores[CVAmidshaft40$groups == levels(CVAmidshaft40$groups)[i],]
  # X <- CVAmidshaft40$Scores[CVAmidshaft40$groups == levels(CVAmidshaft40$groups)[i],]
  hpts <- chull(X)
  hpts <- c(hpts, hpts[1])
  lines(X[hpts, ],col = palette()[i])
  
}



acc <- print(CVAmidshaft40)$accuracy

# legend("topright", legend = levels(CVAmidshaft40$groups), pch = (1:4)+14, col = palette()[1:4])
legend("bottomright", legend = paste("Accuracy = ", round(acc,2), "%", sep = ""), bty = "n")


plot(CVAmidshaft40$CVscores[,c(1,3)], col=data[,1], pch=as.numeric(data[,1])+14, typ="p",asp=1,
     xlab=paste("1st canonical axis", paste(round(CVAmidshaft40$Var[1,2],1),"%")),
     ylab=paste("3rd canonical axis", paste(round(CVAmidshaft40$Var[3,2],1),"%")))


for(i in 1:length(levels(CVAmidshaft40$groups))){
  X <- CVAmidshaft40$CVscores[CVAmidshaft40$groups == levels(CVAmidshaft40$groups)[i],c(1,3)]
  # X <- CVAmidshaft40$Scores[CVAmidshaft40$groups == levels(CVAmidshaft40$groups)[i],c(1,3)]
  hpts <- chull(X)
  hpts <- c(hpts, hpts[1])
  lines(X[hpts, ],col = palette()[i])
  
}


# which position predicts ecology best? -----------------------------------

ID <- allDATA$X

LOCSUB <- as.factor(paste(allDATA$LOCOMOTION, allDATA$SUBSTRATE))

sliceAcc <- rep(NA, length(unique(lengthAlong)))
names(sliceAcc) <- unique(lengthAlong)

CVlist <- list()

which(unique(lengthAlong) == 0.5)

for (i in 1:length(unique(lengthAlong))){
  
  lengthAlong[i]
  
  INDEX <- which(lengthAlong == unique(lengthAlong)[i])
  
  datIND <- ALLdat[INDEX,]
  IDIND <- allDATA$X[INDEX]
  LOCSUBIND <- LOCSUB[INDEX]
  
  levels(LOCSUBIND)
  
  resCV <- Morpho::CVA(datIND, groups = LOCSUBIND, plot= T, weighting = F, cv = F, rounds = 9999)

  
  cl <- print(classify(resCV))
  
  
  sliceAcc[i] <- print(resCV)$accuracy
  
  CVlist[[i]] <- resCV
  # SLICESprobsOUTER <- unlist(lapply(CVAslicesOUTER, FUN = function(X) print(X)$accuracy))
  # SLICESprobsCOMB <- unlist(lapply(CVAslicesCOMB, FUN = function(X) print(X)$accuracy))
  # 
  # sliceAcc[i] <- cl$accuracy
  
}

sliceAcc

write.csv(sliceAcc, file = "./CVA CSG scaled body res/cvaprobs.csv")

# plot

par(mfrow = c(1,1))
# pdf("./CVA CSG scaled body res/Classification accuracy of cortical CSG along the diaphysis.pdf",
#     width = 11.83, height = 5.74)


plot(unique(lengthAlong), sliceAcc, type = "l", 
     main = "Classification accuracy of cortical CSG along the diaphysis",
     xlab = "% diaphysis", ylab = "classification accuracy (%)")
points(unique(lengthAlong), sliceAcc, pch = 16)


# dev.off()

names(which.max(sliceAcc))


max(sliceAcc)

INDEX <- which(lengthAlong == unique(lengthAlong)[which.max(sliceAcc)])

datIND <- ALLdat[INDEX,]
IDIND <- allDATA$X[INDEX]
LOCSUBIND <- LOCSUB[INDEX]

levels(LOCSUBIND)


cvaMAX <- Morpho::CVA(datIND, groups = LOCSUBIND, plot= T, weighting = F ,cv = F,rounds = 9999)
# cvaMAX <- Morpho::groupPCA(datIND, groups = LOCSUBIND, weighting = F,cv = T,rounds = 9999)

sink("./CVA CSG scaled body res/CVAmaxclassification.txt")
cvaMAX
sink()

cl <- print(classify(cvaMAX))


write.csv(cl$probtable, file = "cvaMAXprobs.csv")
palette(palette.colors())



# PLOT MAX ACCURACY CVA ---------------------------------------------------
# dat <- data.frame(Curvature = allDATA$DistCentroid/SIZES$MIDSHAFTCIRC,
#                   # Circumferance = allDATA$Ext_perim,
#                   thickness = allDATA$Mean_thick/SIZES$MIDSHAFTCIRC,
#                   CA = allDATA$CA/(SIZES$MIDSHAFTCIRC^2),
#                   I.x = (allDATA$Ix/(SIZES$MIDSHAFTCIRC^4)),
#                   I.y = (allDATA$Iy/(SIZES$MIDSHAFTCIRC^4)),
#                   I.ratio = (allDATA$Imax/allDATA$Imin),
#                   Zpol = allDATA$Zpol.Zpol/(SIZES$MIDSHAFTCIRC^3),
#                   J = allDATA$J.Ix/(SIZES$MIDSHAFTCIRC^4))

# dat<- scale(dat, center = T, scale = T)

par(mfrow = c(1,2))
plot(cvaMAX$CVscores, col=data[,1], pch=as.numeric(data[,1])+14, typ="p",asp=1,
     xlab=paste("1st canonical axis", paste(round(cvaMAX$Var[1,2],1),"%")),
     ylab=paste("2nd canonical axis", paste(round(cvaMAX$Var[2,2],1),"%")),
     main = paste("CVA using diaphysis", names(which.max(sliceAcc)), "values"))

# text(CVAmaxvarpos$CVscores, as.character(data[,1]), col=as.numeric(data[,1]), cex=.7)

for(i in 1:length(levels(cvaMAX$groups))){
  X <- cvaMAX$CVscores[cvaMAX$groups == levels(cvaMAX$groups)[i],]
  hpts <- chull(X)
  hpts <- c(hpts, hpts[1])
  lines(X[hpts, ],col = palette()[i])
  
}

legend("topright", legend = levels(CVAmidshaft40$groups),
       pch = (1:4)+14, col = palette()[1:4],
       cex = 0.7)
legend("bottomright",
       legend = paste("Accuracy = ", round(max(sliceAcc),2), "%"),
       bty = "n")


plot(cvaMAX$CVscores[,c(3,2)], col=data[,1], pch=as.numeric(data[,1])+14, typ="p",asp=1,
     xlab=paste("3rd canonical axis", paste(round(cvaMAX$Var[3,2],1),"%")),
     ylab=paste("2nd canonical axis", paste(round(cvaMAX$Var[2,2],1),"%")))
# text(CVAmaxvarpos$CVscores, as.character(data[,1]), col=as.numeric(data[,1]), cex=.7)

for(i in 1:length(levels(cvaMAX$groups))){
  X <- cvaMAX$CVscores[cvaMAX$groups == levels(cvaMAX$groups)[i],c(3,2)]
  hpts <- chull(X)
  hpts <- c(hpts, hpts[1])
  lines(X[hpts, ],col = palette()[i])
  
}

for(i in 1:length(levels(cvaMAX$groups))){
  # X <- cvaMAX$CVscores[cvaMAX$groups == levels(cvaMAX$groups)[i],c(1,3)]
  X <- cvaMAX$Scores[cvaMAX$groups == levels(cvaMAX$groups)[i],c(1,3)]
  hpts <- chull(X)
  hpts <- c(hpts, hpts[1])
  lines(X[hpts, ],col = palette()[i])
  
}





# plot groupmeans

cvaMAX$groupmeans
CVAmidshaft50$groupmeans

write.csv(CVAmidshaft50$groupmeans/cvaMAX$groupmeans, 
          "./CVA CSG scaled body res/CSGmmidshaft-CSGdistal.csv")



# calc group means and SD -------------------------------------------------


### calc group means and standard deviation

maxIND<- which.max(sliceAcc)

INDEX <- which(lengthAlong == names(maxIND))

# ALLdat <- scale(ALLdat, center = F)


datIND <- rawdat[INDEX,]
IDIND <- allDATA$X[INDEX]
LOCSUBIND <- LOCSUB[INDEX]

grpMeansMAX <- matrix(ncol = ncol(datIND), nrow = 4)

colnames(grpMeansMAX) <- colnames(datIND)
rownames(grpMeansMAX) <- levels(LOCSUBIND)


grpSDMAX <- grpMeansMAX

for(i in 1:length(levels(LOCSUBIND))){
  
  grpMeansMAX[i,] <- colMeans(datIND[LOCSUBIND==levels(LOCSUBIND)[i],])
  grpSDMAX[i,] <- apply(datIND[LOCSUBIND==levels(LOCSUBIND)[i],], 2, sd)
  
}

grpMeansMAX
grpSDMAX

##

INDEX <- which(lengthAlong == 0.5)

datIND <- rawdat[INDEX,]
IDIND <- allDATA$X[INDEX]
LOCSUBIND <- LOCSUB[INDEX]

grpMeans50 <- matrix(ncol = ncol(datIND), nrow = 4)

colnames(grpMeans50) <- colnames(datIND)
rownames(grpMeans50) <- levels(LOCSUBIND)


grpSD50  <- grpMeans50 

for(i in 1:length(levels(LOCSUBIND))){
  
  grpMeans50[i,] <- colMeans(datIND[LOCSUBIND==levels(LOCSUBIND)[i],])
  grpSD50[i,] <- apply(datIND[LOCSUBIND==levels(LOCSUBIND)[i],], 2, sd)
  
}


# transform CVA data ------------------------------------------------------
#transform groupmeans into original scale


grpMeans50
colnames(grpSD50) <- rep("sd", 8)

grpMeansMAX
colnames(grpSDMAX) <- rep("sd", 8)

grpMSD50 <- cbind(grpMeans50, grpSD50)[,c(1,1+8,2,2+8,3,3+8, 4,4+8,5,5+8,6,6+8,7,7+8,8,+8)]
grpMSDMAX <- cbind(grpMeansMAX, grpSDMAX)[,c(1,1+8,2,2+8,3,3+8, 4,4+8,5,5+8,6,6+8,7,7+8,8,+8)]

write.csv((grpMSD50), "./CVA CSG scaled body res/groupCSGMeansSD50.csv")
write.csv((grpMSDMAX), "./CVA CSG scaled body res/groupCSGMeansSDMAX.csv")
round(grpMSDMAX,4)

par(mfrow = c(4,2), mar = c(5,1,4,1) + 0.1)

labs = colnames(CVAmidshaft50$groupmeans)

for(i in 1:nrow(CVAmidshaft50$groupmeans)){
  
  tmp <- cbind(CVAmidshaft50$groupmeans[i,]-colMeans(CVAmidshaft50$groupmeans),
               cvaMAX$groupmeans[i,]-colMeans(cvaMAX$groupmeans))
  
  colnames(tmp) <- c("midshaft","distal")
  
  
  bp1 <- barplot(tmp[,1],beside = T, col = palette.colors()[1:8], 
                 main = paste(rownames(CVAmidshaft50$groupmeans)[i], "(midshaft)"), 
                 space = 0, axes = F, las = 2)
  
  # text(cex=1, x=bp1-.25, y=-1.25, labs, xpd=TRUE, srt=45)
  
  bp2 <- barplot(tmp[,2],beside = T, col = palette.colors()[1:8], 
                 main = paste(rownames(CVAmidshaft50$groupmeans)[i], "(distal)"), 
                 space = 0, axes = F, las = 2)  
  
  
  
}




###########################################################################

# run permutation model to optimize site combos for best model fit --------

# test all combinations of two sites to find bestdiscrimination
# 
# lengthAlong = allDATA$lengthAlong
# 
# 
# ID <- allDATA$X
# 
# LOCSUB <- as.factor(paste(allDATA$LOCOMOTION, allDATA$SUBSTRATE))
# 
# 
# ALLdat
# 
# pairs <- expand.grid(unique(lengthAlong), unique(lengthAlong))
# 
# rmInd<- apply(pairs, 1, FUN = function(X){
#   
#   X[1] == X[2]
#   
# })
# 
# # pairs <-pairs[-rmInd,]
# 
# 
# 
# sliceAccPAIRS <- rep(NA, nrow(pairs))
# 
# sliceAccPAIRS <- data.frame(pairs, sliceAccPAIRS)
# 
# 
# i =1 
# 
# for (i in 1:nrow(sliceAccPAIRS)){
#   
#   
#   if(sliceAccPAIRS[i,1] == sliceAccPAIRS[i,2]){
#     
#     INDEX <- grep(sliceAccPAIRS[i,1], lengthAlong)
#     
#     
#   } else{
#     INDEX <- grep(sliceAccPAIRS[i,1], lengthAlong)
#     INDEX <- c(INDEX,grep(sliceAccPAIRS[i,2], lengthAlong))
#     
#   }
#   
#   lengthAlong[INDEX]
#   datIND <- ALLdat[INDEX,]
#   IDIND <- allDATA$X[INDEX]
#   LOCSUBIND <- LOCSUB[INDEX]
#   
#   levels(LOCSUBIND)
#   
#   resCV <- Morpho::CVA(datIND, groups = LOCSUBIND, weighting = F, cv = F, plot= T,rounds = 9999)
#   # resCV <- Morpho::groupPCA(datIND, groups = LOCSUBIND, weighting = F, cv = F,rounds = 9999)
#   # res <- Morpho::CVA(datIND, groups = LOCSUBIND, rounds = 9999)
#   
#   # 
#   # cl <- print(classify(resCV))
#   # 
#   # sliceAccPAIRS[i,3] <- cl$accuracy
#   sliceAccPAIRS[i,3] <- print(resCV)$accuracy
#   
#   
# }
# 
# 
# sliceAccPAIRS <- sliceAccPAIRS[order(sliceAccPAIRS[,3], decreasing = T),]
# 
# save(sliceAccPAIRS, file = "./CVA CSG scaled body res/CVAsliceAccPAIRS.Rdata")
# write.csv(sliceAccPAIRS, "./CVA CSG scaled body res/CVAsliceAccPAIRS.csv")
# 
# 
# sliceAccPAIRS

###########################################################################
# CVA USING FULL DIAPHYSIS CSG --------------------------------------------
###########################################################################


# make wide dataframe -----------------------------------------------------

ID <- allDATA$X

LOCSUB <- as.factor(paste(allDATA$LOCOMOTION, allDATA$SUBSTRATE))

# dat <- data.frame(Curvature = allDATA$DistCentroid/max(allDATA$DistCentroid),
#                   Thickness = allDATA$Mean_thick/max(allDATA$Mean_thick),
#                   Circumferance = allDATA$Ext_perim/max(allDATA$Ext_perim),
#                   CA = allDATA$CA/max(allDATA$CA),
#                   Ix = (allDATA$Ix)/max(allDATA$Ix),
#                   Iy = (allDATA$Iy)/max(allDATA$Iy),
#                   I.ratio = (allDATA$Imax/allDATA$Imin),
#                   Zpol = allDATA$Zpol.Zpol/max(allDATA$Zpol.Zpol),
#                   J = allDATA$J.Ix/max(allDATA$J.Ix)
# )


# 
# dat <- data.frame(Curvature = allDATA$DistCentroid,
#                   Circumferance = allDATA$Ext_perim,
#                   thickness = allDATA$Mean_thick,
#                   CA = allDATA$CA,
#                   I.x = (allDATA$Ix),
#                   I.y = (allDATA$Iy),
#                   I.ratio = (allDATA$Imax/allDATA$Imin),
#                   Zpol = allDATA$Zpol.Zpol,
#                   J = allDATA$J.Ix
# )


lengthAlong

SIZES <- boneVars$ecoDF[,c(6:7)]
ECO <- boneVars$ecoDF[,c(3:5)]
LOCSUB = as.factor(paste(ECO$LOCOMOTION, ECO$SUBSTRATE))
MIDSHAFTCIRC = unlist(allDATA$MidShaft.Ext_perim)
# ZLENGTH = ZlengthDF

unique(lengthAlong)

WIDE <- NULL

i=1 # set ID value

for(i in 1:length(levels(ID))){
  
  II <- levels(ID)[i]
  
  XX <- ALLdat[ID == II,]
  combNames <- expand.grid(unique(lengthAlong), colnames(XX))[,c(2,1)]
  
  CC <- c((as.matrix(XX)))
  names(CC) <- paste(combNames[,1],combNames[,2])
  
  WIDE <- rbind(WIDE, CC)
  # rownames(WIDE)[i] <- ID 
  
}

allDATAwide <- data.frame(ID = levels(ID), ECO = LOCSUB, SIZES, WIDE)

rownames(allDATAwide) <- levels(ID)

allDATAwide

save(allDATAwide, file= "./CVA CSG scaled body res/allDATAwide.Rdata")



# calc CVA ----------------------------------------------------------------

load("./CVA CSG scaled body res/allDATAwide.Rdata")


require(Morpho)


allDATAwide

LOCSUB = as.factor(paste(ECO$LOCOMOTION, ECO$SUBSTRATE))


PCA <- prcomp(allDATAwide[,-c(1,2,3,4)])

par(mfrow = c(1,1))
plot(PCA$x, col = allDATAwide$ECO)
PCA$x

CVAallCSG <- Morpho::CVA(allDATAwide[,-c(1,2,3,4)], groups = allDATAwide$ECO, weighting= F,cv = F, rounds = 9999)


# CVAallCSGcv <- Morpho::CVA(allDATAwide[,-c(1,2,3,4)], cv = F, groups = allDATAwide$ECO, rounds = 9999)
# bgPCAallCSG <- Morpho::groupPCA(allDATAwide[,-c(1,2,3,4)], groups = allDATAwide$ECO, rounds = 9999)





# -------------------------------------------------------------------------



CVAallCSG
# bgPCAallCSG

CVAallCSG$groupmeans

# CVAallCSG$CV
# CVAallCSG$groupmeans
# CVAallCSG$CVvis

save(CVAallCSG, file = "./CVA CSG scaled body res/CVAallCSG.Rdata")


sink("./CVA CSG scaled body res/CVAallCSG.txt")
CVAallCSG
sink()


CVAallCSG$CV


cl <- print(classify(CVAallCSG))


write.csv(cl$probtable, file = "./CVA CSG scaled body res/cvaALLCSGprobs.csv")



# plot --------------------------------------------------------------------


palette(palette.colors())

par(mfrow = c(1,2), mar = c(5,4,4,2))
plot(CVAallCSG$CVscores, col=data[,1], pch=as.numeric(data[,1])+14, typ="p",
     xlab=paste("1st canonical axis", paste(round(CVAallCSG$Var[1,2],1),"%")),
     ylab=paste("2nd canonical axis", paste(round(CVAallCSG$Var[2,2],1),"%")),
     main = "CVA using all CSG values")
# text(CVAmaxvarpos$CVscores, as.character(data[,1]), col=as.numeric(data[,1]), cex=.7)

for(i in 1:length(levels(CVAallCSG$groups))){
  X <- CVAallCSG$CVscores[CVAallCSG$groups == levels(CVAallCSG$groups)[i],]
  hpts <- chull(X)
  hpts <- c(hpts, hpts[1])
  lines(X[hpts, ],col = palette()[i])
  
}
# 
# if (require(car)) {
#   for(ii in 1:length(levels(data[,1]))){
#     dataEllipse(CVAmaxvarpos$CVscores[data[,1]==levels(data[,1])[ii],1],
#                 CVAmaxvarpos$CVscores[data[,1]==levels(data[,1])[ii],2], 
#                 add=TRUE,levels=.80, col=c(1:7)[ii])}
# }



legend("topleft", legend = "accuracy = 100%", bty = "n")



plot(CVAallCSG$CVscores[,c(3,2)], col=data[,1], pch=as.numeric(data[,1])+14, typ="p",asp=1,
     xlab=paste("1st canonical axis", paste(round(CVAallCSG$Var[3,2],1),"%")),
     ylab=paste("3rd canonical axis", paste(round(CVAallCSG$Var[2,2],1),"%")))
# text(CVAmaxvarpos$CVscores, as.character(data[,1]), col=as.numeric(data[,1]), cex=.7)

for(i in 1:length(levels(CVAallCSG$groups))){
  X <- CVAallCSG$CVscores[CVAallCSG$groups == levels(CVAallCSG$groups)[i],c(3,2)]
  hpts <- chull(X)
  hpts <- c(hpts, hpts[1])
  lines(X[hpts, ],col = palette()[i])
  
}
legend("topright", legend = levels(CVAallCSG$groups),
       cex = 0.7,
       pch = (1:4)+14, col = palette()[1:4])







# plot all CVA together ---------------------------------------------------


# par(mfrow = c(4,2))
par(mfrow = c(4,2), mar = c(4,4,3,1))

### 50 CVA


plot(CVAmidshaft50$CVscores, col=data[,1], pch=as.numeric(data[,1])+14, typ="p",asp=1,
     xlab=paste("1st canonical axis", paste(round(CVAmidshaft50$Var[1,2],1),"%")),
     ylab=paste("2nd canonical axis", paste(round(CVAmidshaft50$Var[2,2],1),"%")),
     main = "CVA using diaphysis 50 values")
# text(CVAmaxvarpos$CVscores, as.character(data[,1]), col=as.numeric(data[,1]), cex=.7)

for(i in 1:length(levels(CVAmidshaft50$groups))){
  X <- CVAmidshaft50$CVscores[CVAmidshaft50$groups == levels(CVAmidshaft50$groups)[i],]
  hpts <- chull(X)
  hpts <- c(hpts, hpts[1])
  lines(X[hpts, ],col = palette()[i])
  
}
# 
# if (require(car)) {
#   for(ii in 1:length(levels(data[,1]))){
#     dataEllipse(CVAmaxvarpos$CVscores[data[,1]==levels(data[,1])[ii],1],
#                 CVAmaxvarpos$CVscores[data[,1]==levels(data[,1])[ii],2], 
#                 add=TRUE,levels=.80, col=c(1:7)[ii])}
# }


acc <- print(CVAmidshaft50)$accuracy


# legend("bottomright", legend = levels(CVAmidshaft50$groups), pch = (1:4)+14, cex = 0.7,col = palette()[1:4])
legend("bottomright", legend = paste("Accuracy = ", round(acc,2), "%", sep = ""), bty = "n")
# legend("bottomright", legend = levels(CVAallCSG$groups),
#        cex = 1,
#        pch = (1:4)+14, col = palette()[1:4])


plot(CVAmidshaft50$CVscores[,c(1,3)], col=data[,1], pch=as.numeric(data[,1])+14, typ="p",asp=1,
     xlab=paste("1st canonical axis", paste(round(CVAmidshaft50$Var[1,2],1),"%")),
     ylab=paste("3rd canonical axis", paste(round(CVAmidshaft50$Var[3,2],1),"%")))
# text(CVAmaxvarpos$CVscores, as.character(data[,1]), col=as.numeric(data[,1]), cex=.7)

for(i in 1:length(levels(CVAmidshaft50$groups))){
  X <- CVAmidshaft50$CVscores[CVAmidshaft50$groups == levels(CVAmidshaft50$groups)[i],c(1,3)]
  hpts <- chull(X)
  hpts <- c(hpts, hpts[1])
  lines(X[hpts, ],col = palette()[i])
  
}



### 40 CVA

plot(CVAmidshaft40$CVscores, col=data[,1], pch=as.numeric(data[,1])+14, typ="p",asp=1,
     xlab=paste("1st canonical axis", paste(round(CVAmidshaft40$Var[1,2],1),"%")),
     ylab=paste("2nd canonical axis", paste(round(CVAmidshaft40$Var[2,2],1),"%")),
     main = "CVA using diaphysis 40 values")
# text(CVAmaxvarpos$CVscores, as.character(data[,1]), col=as.numeric(data[,1]), cex=.7)

for(i in 1:length(levels(CVAmidshaft40$groups))){
  X <- CVAmidshaft40$CVscores[CVAmidshaft40$groups == levels(CVAmidshaft40$groups)[i],]
  hpts <- chull(X)
  hpts <- c(hpts, hpts[1])
  lines(X[hpts, ],col = palette()[i])
  
}


#
# if (require(car)) {
#   for(ii in 1:length(levels(data[,1]))){
#     dataEllipse(CVAmaxvarpos$CVscores[data[,1]==levels(data[,1])[ii],1],
#                 CVAmaxvarpos$CVscores[data[,1]==levels(data[,1])[ii],2],
#                 add=TRUE,levels=.80, col=c(1:7)[ii])}
# }


acc <- print(CVAmidshaft40)$accuracy

# legend("topright", legend = levels(CVAmidshaft40$groups), pch = (1:4)+14, col = palette()[1:4])
legend("bottomright", legend = paste("Accuracy = ", round(acc,2), "%", sep = ""), bty = "n")


plot(CVAmidshaft40$CVscores[,c(1,3)], col=data[,1], pch=as.numeric(data[,1])+14, typ="p",asp=1,
     xlab=paste("1st canonical axis", paste(round(CVAmidshaft40$Var[1,2],1),"%")),
     ylab=paste("3rd canonical axis", paste(round(CVAmidshaft40$Var[3,2],1),"%")))
# text(CVAmaxvarpos$CVscores, as.character(data[,1]), col=as.numeric(data[,1]), cex=.7)

for(i in 1:length(levels(CVAmidshaft40$groups))){
  X <- CVAmidshaft40$CVscores[CVAmidshaft40$groups == levels(CVAmidshaft40$groups)[i],c(1,3)]
  hpts <- chull(X)
  hpts <- c(hpts, hpts[1])
  lines(X[hpts, ],col = palette()[i])
  
}


### max CVA

plot(cvaMAX$CVscores, col=data[,1], pch=as.numeric(data[,1])+14, typ="p",asp=1,
     xlab=paste("1st canonical axis", paste(round(cvaMAX$Var[1,2],1),"%")),
     ylab=paste("2nd canonical axis", paste(round(cvaMAX$Var[2,2],1),"%")),
     main = paste("CVA using diaphysis", names(which.max(sliceAcc)), "values"))
# text(CVAmaxvarpos$CVscores, as.character(data[,1]), col=as.numeric(data[,1]), cex=.7)

for(i in 1:length(levels(cvaMAX$groups))){
  X <- cvaMAX$CVscores[cvaMAX$groups == levels(cvaMAX$groups)[i],]
  hpts <- chull(X)
  hpts <- c(hpts, hpts[1])
  lines(X[hpts, ],col = palette()[i])
  
}
# 
# if (require(car)) {
#   for(ii in 1:length(levels(data[,1]))){
#     dataEllipse(CVAmaxvarpos$CVscores[data[,1]==levels(data[,1])[ii],1],
#                 CVAmaxvarpos$CVscores[data[,1]==levels(data[,1])[ii],2], 
#                 add=TRUE,levels=.80, col=c(1:7)[ii])}
# }




# legend("topright", legend = levels(CVAmidshaft40$groups), 
#        pch = (1:4)+14, col = palette()[1:4],
#        cex = 0.7)
legend("bottomright", 
       legend = paste("Accuracy = ", round(max(sliceAcc),2), "%"), 
       bty = "n")


plot(cvaMAX$CVscores[,c(3,2)], col=data[,1], pch=as.numeric(data[,1])+14, typ="p",asp=1,
     xlab=paste("3rd canonical axis", paste(round(cvaMAX$Var[3,2],1),"%")),
     ylab=paste("2nd canonical axis", paste(round(cvaMAX$Var[2,2],1),"%")))
# text(CVAmaxvarpos$CVscores, as.character(data[,1]), col=as.numeric(data[,1]), cex=.7)

for(i in 1:length(levels(cvaMAX$groups))){
  X <- cvaMAX$CVscores[cvaMAX$groups == levels(cvaMAX$groups)[i],c(3,2)]
  hpts <- chull(X)
  hpts <- c(hpts, hpts[1])
  lines(X[hpts, ],col = palette()[i])
  
}



### all CVA

palette(palette.colors())

plot(CVAallCSG$CVscores, col=data[,1], pch=as.numeric(data[,1])+14, typ="p",
     xlab=paste("1st canonical axis", paste(round(CVAallCSG$Var[1,2],1),"%")),
     ylab=paste("2nd canonical axis", paste(round(CVAallCSG$Var[2,2],1),"%")),
     main = "CVA using all diaphysis CSG values")
# text(CVAmaxvarpos$CVscores, as.character(data[,1]), col=as.numeric(data[,1]), cex=.7)

for(i in 1:length(levels(CVAallCSG$groups))){
  X <- CVAallCSG$CVscores[CVAallCSG$groups == levels(CVAallCSG$groups)[i],]
  hpts <- chull(X)
  hpts <- c(hpts, hpts[1])
  lines(X[hpts, ],col = palette()[i])
  
}
# 
# if (require(car)) {
#   for(ii in 1:length(levels(data[,1]))){
#     dataEllipse(CVAmaxvarpos$CVscores[data[,1]==levels(data[,1])[ii],1],
#                 CVAmaxvarpos$CVscores[data[,1]==levels(data[,1])[ii],2], 
#                 add=TRUE,levels=.80, col=c(1:7)[ii])}
# }


legend("topright", legend = levels(CVAallCSG$groups),
       cex = 0.9,
       pch = (1:4)+14, col = palette()[1:4])

legend("bottomright", legend = "Accuracy = 100%", bty = "n")
# legend(x=-16,y = 17, legend = "Accuracy = 100%", bty = "n")



plot(CVAallCSG$CVscores[,c(3,2)], col=data[,1], pch=as.numeric(data[,1])+14, typ="p",asp=1,
     xlab=paste("1st canonical axis", paste(round(CVAallCSG$Var[3,2],1),"%")),
     ylab=paste("3rd canonical axis", paste(round(CVAallCSG$Var[2,2],1),"%")))
# text(CVAmaxvarpos$CVscores, as.character(data[,1]), col=as.numeric(data[,1]), cex=.7)

for(i in 1:length(levels(CVAallCSG$groups))){
  X <- CVAallCSG$CVscores[CVAallCSG$groups == levels(CVAallCSG$groups)[i],c(3,2)]
  hpts <- chull(X)
  hpts <- c(hpts, hpts[1])
  lines(X[hpts, ],col = palette()[i])
  
}



# plot all CVA together (pdf) ---------------------------------------------------


pdf("./CVA CSG scaled body res/allCSGCVAs.pdf",
    width = 7.78,
    height = 12.82)

# par(mfrow = c(4,2))
# par(mfrow = c(4,2))
par(mfrow = c(4,2), mar = c(4,4,3,1))

### 50 CVA


plot(CVAmidshaft50$CVscores, col=data[,1], pch=as.numeric(data[,1])+14, typ="p",asp=1,
     xlab=paste("1st canonical axis", paste(round(CVAmidshaft50$Var[1,2],1),"%")),
     ylab=paste("2nd canonical axis", paste(round(CVAmidshaft50$Var[2,2],1),"%")),
     main = "CVA using diaphysis 50 values")
# text(CVAmaxvarpos$CVscores, as.character(data[,1]), col=as.numeric(data[,1]), cex=.7)

for(i in 1:length(levels(CVAmidshaft50$groups))){
  X <- CVAmidshaft50$CVscores[CVAmidshaft50$groups == levels(CVAmidshaft50$groups)[i],]
  hpts <- chull(X)
  hpts <- c(hpts, hpts[1])
  lines(X[hpts, ],col = palette()[i])
  
}
# 
# if (require(car)) {
#   for(ii in 1:length(levels(data[,1]))){
#     dataEllipse(CVAmaxvarpos$CVscores[data[,1]==levels(data[,1])[ii],1],
#                 CVAmaxvarpos$CVscores[data[,1]==levels(data[,1])[ii],2], 
#                 add=TRUE,levels=.80, col=c(1:7)[ii])}
# }


acc <- print(CVAmidshaft50)$accuracy


# legend("bottomright", legend = levels(CVAmidshaft50$groups), pch = (1:4)+14, cex = 0.7,col = palette()[1:4])
legend("bottomright", legend = paste("Accuracy = ", round(acc,2), "%", sep = ""), bty = "n")
# legend("bottomright", legend = levels(CVAallCSG$groups),
#        cex = 1,
#        pch = (1:4)+14, col = palette()[1:4])


plot(CVAmidshaft50$CVscores[,c(1,3)], col=data[,1], pch=as.numeric(data[,1])+14, typ="p",asp=1,
     xlab=paste("1st canonical axis", paste(round(CVAmidshaft50$Var[1,2],1),"%")),
     ylab=paste("3rd canonical axis", paste(round(CVAmidshaft50$Var[3,2],1),"%")))
# text(CVAmaxvarpos$CVscores, as.character(data[,1]), col=as.numeric(data[,1]), cex=.7)

for(i in 1:length(levels(CVAmidshaft50$groups))){
  X <- CVAmidshaft50$CVscores[CVAmidshaft50$groups == levels(CVAmidshaft50$groups)[i],c(1,3)]
  hpts <- chull(X)
  hpts <- c(hpts, hpts[1])
  lines(X[hpts, ],col = palette()[i])
  
}



### 40 CVA

plot(CVAmidshaft40$CVscores, col=data[,1], pch=as.numeric(data[,1])+14, typ="p",asp=1,
     xlab=paste("1st canonical axis", paste(round(CVAmidshaft40$Var[1,2],1),"%")),
     ylab=paste("2nd canonical axis", paste(round(CVAmidshaft40$Var[2,2],1),"%")),
     main = "CVA using diaphysis 40 values")
# text(CVAmaxvarpos$CVscores, as.character(data[,1]), col=as.numeric(data[,1]), cex=.7)

for(i in 1:length(levels(CVAmidshaft40$groups))){
  X <- CVAmidshaft40$CVscores[CVAmidshaft40$groups == levels(CVAmidshaft40$groups)[i],]
  hpts <- chull(X)
  hpts <- c(hpts, hpts[1])
  lines(X[hpts, ],col = palette()[i])
  
}


#
# if (require(car)) {
#   for(ii in 1:length(levels(data[,1]))){
#     dataEllipse(CVAmaxvarpos$CVscores[data[,1]==levels(data[,1])[ii],1],
#                 CVAmaxvarpos$CVscores[data[,1]==levels(data[,1])[ii],2],
#                 add=TRUE,levels=.80, col=c(1:7)[ii])}
# }


acc <- print(CVAmidshaft40)$accuracy

# legend("topright", legend = levels(CVAmidshaft40$groups), pch = (1:4)+14, col = palette()[1:4])
legend("bottomright", legend = paste("Accuracy = ", round(acc,2), "%", sep = ""), bty = "n")


plot(CVAmidshaft40$CVscores[,c(1,3)], col=data[,1], pch=as.numeric(data[,1])+14, typ="p",asp=1,
     xlab=paste("1st canonical axis", paste(round(CVAmidshaft40$Var[1,2],1),"%")),
     ylab=paste("3rd canonical axis", paste(round(CVAmidshaft40$Var[3,2],1),"%")))
# text(CVAmaxvarpos$CVscores, as.character(data[,1]), col=as.numeric(data[,1]), cex=.7)

for(i in 1:length(levels(CVAmidshaft40$groups))){
  X <- CVAmidshaft40$CVscores[CVAmidshaft40$groups == levels(CVAmidshaft40$groups)[i],c(1,3)]
  hpts <- chull(X)
  hpts <- c(hpts, hpts[1])
  lines(X[hpts, ],col = palette()[i])
  
}


### max CVA

plot(cvaMAX$CVscores, col=data[,1], pch=as.numeric(data[,1])+14, typ="p",asp=1,
     xlab=paste("1st canonical axis", paste(round(cvaMAX$Var[1,2],1),"%")),
     ylab=paste("2nd canonical axis", paste(round(cvaMAX$Var[2,2],1),"%")),
     main = paste("CVA using diaphysis", names(which.max(sliceAcc)), "values"))
# text(CVAmaxvarpos$CVscores, as.character(data[,1]), col=as.numeric(data[,1]), cex=.7)

for(i in 1:length(levels(cvaMAX$groups))){
  X <- cvaMAX$CVscores[cvaMAX$groups == levels(cvaMAX$groups)[i],]
  hpts <- chull(X)
  hpts <- c(hpts, hpts[1])
  lines(X[hpts, ],col = palette()[i])
  
}
# 
# if (require(car)) {
#   for(ii in 1:length(levels(data[,1]))){
#     dataEllipse(CVAmaxvarpos$CVscores[data[,1]==levels(data[,1])[ii],1],
#                 CVAmaxvarpos$CVscores[data[,1]==levels(data[,1])[ii],2], 
#                 add=TRUE,levels=.80, col=c(1:7)[ii])}
# }




# legend("topright", legend = levels(CVAmidshaft40$groups), 
#        pch = (1:4)+14, col = palette()[1:4],
#        cex = 0.7)
legend("bottomright", 
       legend = paste("Accuracy = ", round(max(sliceAcc),2), "%"), 
       bty = "n")


plot(cvaMAX$CVscores[,c(3,2)], col=data[,1], pch=as.numeric(data[,1])+14, typ="p",asp=1,
     xlab=paste("3rd canonical axis", paste(round(cvaMAX$Var[3,2],1),"%")),
     ylab=paste("2nd canonical axis", paste(round(cvaMAX$Var[2,2],1),"%")))
# text(CVAmaxvarpos$CVscores, as.character(data[,1]), col=as.numeric(data[,1]), cex=.7)

for(i in 1:length(levels(cvaMAX$groups))){
  X <- cvaMAX$CVscores[cvaMAX$groups == levels(cvaMAX$groups)[i],c(3,2)]
  hpts <- chull(X)
  hpts <- c(hpts, hpts[1])
  lines(X[hpts, ],col = palette()[i])
  
}



### all CVA

palette(palette.colors())

plot(CVAallCSG$CVscores, col=data[,1], pch=as.numeric(data[,1])+14, typ="p",
     xlab=paste("1st canonical axis", paste(round(CVAallCSG$Var[1,2],1),"%")),
     ylab=paste("2nd canonical axis", paste(round(CVAallCSG$Var[2,2],1),"%")),
     main = "CVA using all diaphysis CSG values")
# text(CVAmaxvarpos$CVscores, as.character(data[,1]), col=as.numeric(data[,1]), cex=.7)

for(i in 1:length(levels(CVAallCSG$groups))){
  X <- CVAallCSG$CVscores[CVAallCSG$groups == levels(CVAallCSG$groups)[i],]
  hpts <- chull(X)
  hpts <- c(hpts, hpts[1])
  lines(X[hpts, ],col = palette()[i])
  
}
# 
# if (require(car)) {
#   for(ii in 1:length(levels(data[,1]))){
#     dataEllipse(CVAmaxvarpos$CVscores[data[,1]==levels(data[,1])[ii],1],
#                 CVAmaxvarpos$CVscores[data[,1]==levels(data[,1])[ii],2], 
#                 add=TRUE,levels=.80, col=c(1:7)[ii])}
# }


legend("topright", legend = levels(CVAallCSG$groups),
       cex = 0.9,
       pch = (1:4)+14, col = palette()[1:4])

legend("bottomright", legend = "Accuracy = 100%", bty = "n")
# legend(x=-16,y = 17, legend = "Accuracy = 100%", bty = "n")



plot(CVAallCSG$CVscores[,c(3,2)], col=data[,1], pch=as.numeric(data[,1])+14, typ="p",asp=1,
     xlab=paste("1st canonical axis", paste(round(CVAallCSG$Var[3,2],1),"%")),
     ylab=paste("3rd canonical axis", paste(round(CVAallCSG$Var[2,2],1),"%")))
# text(CVAmaxvarpos$CVscores, as.character(data[,1]), col=as.numeric(data[,1]), cex=.7)

for(i in 1:length(levels(CVAallCSG$groups))){
  X <- CVAallCSG$CVscores[CVAallCSG$groups == levels(CVAallCSG$groups)[i],c(3,2)]
  hpts <- chull(X)
  hpts <- c(hpts, hpts[1])
  lines(X[hpts, ],col = palette()[i])
  
}
dev.off()
# 
# 
# # plot vectors ------------------------------------------------------------
# 
# 
# 
# 
# parName <- as.factor(unlist(lapply(colnames(dat),1 , FUN = rep, 29)))
# 
# unique(parName) 
# 
# parL <- rep(unique(lengthAlong), 9)
# 
# CVscores <- data.frame(parName, parL,CVAallCSG$CV)
# 
# rownames(CVscores) <- row.names(CVAallCSG$CVvis)
# 
# 
# 
# 
# par(mfrow= c(1,2), mar = c(5,4,4,2), las = 2)
# 
# plot(CVAallCSG$CV[,c(1,2)], type = "n")
# lines(x = range(CVAallCSG$CV[,1])*1.5, 
#       y = c(0,0),
#       lty = 2)
# 
# lines(x = c(0,0),
#       y = range(CVAallCSG$CV[,2])*1.5,
#       lty = 2)
# 
# 
# arrows(0,0, CVAallCSG$CV[,1], CVAallCSG$CV[,2], length = 0.1, col = parName)
# 
# text(CVAallCSG$CV[,c(1:2)]*1.15, labels = rownames(CVAallCSG$CVvis), col = as.numeric(parName))
# 
# plot(pred.all$x[,c(3,2)], col = cols, pch = 16, cex = 1.5, type = "n",
#      xlab = "LDA3 (0.1634%)",asp = 1,
#      ylab = "LDA2 (0.3290%)")
# lines(x = range(pred.all$x[,3])*1.5, 
#       y = c(0,0),
#       lty = 2)
# 
# lines(x = c(0,0),
#       y = range(pred.all$x[,2])*1.5,
#       lty = 2)
# 
# 
# 
# 
# arrows(0,0, scaling[,3], scaling[,2], length = 0.1, col = "red")
# text(scaling[,c(3:2)]*1.15, labels = rownames(scaling), col = "red")
# 
# #### EVEN WORSE
# 
# 
# # plot Bars ---------------------------------------------------------------
# 
# ### looks like shit
# 
# par(mfcol= c(7,3), mar = c(2,4,1,2), las = 2)
# # barplot(CV.1~parName+parL, data = CVscores, beside=T)
# 
# ii=1
# 
# iii = 3
# 
# CVscores < 1
# 
# 
# CVscores <- CVAallCSG$CV
# 
# 
# 
# apply(abs(CVAallCSG$CV), 1,sum)
# 
# 
# 
# 
# 
# 
# for( iii in 1:ncol(CVscores)){
#   
#   # ylim <- range(CVscores[,iii])
#   ylim <- c(0,max(CVscores[,iii]))
#   
#   
#   tmpCV <- CVscores[parName==levels(parName)[ii],iii]
#   
#   cols <- tmpCV
#   
#   cols[which(tmpCV >0)] <- "black"
#   cols[which(tmpCV <0)] <- "grey"
#   
#   tmpCV <- abs(tmpCV)
#   
#   
#   par(mar = c(2,4,3,1), las = 2)
#   
#   plot(unique(lengthAlong),tmpCV, pch = 16, main = paste("CV", iii,  "loadings"),
#        xlab = "", ylab = levels(parName)[1], ylim = ylim, col = cols)
#   
#   if(iii == 3) legend("topleft", fill = c("black", "grey"), legend = c("pos", "neg"))
#   
#   
#   
#   for(i in 1:length(unique(lengthAlong))){
#     
#     y= c(0 ,tmpCV[i])
#     x= c(unique(lengthAlong)[i], unique(lengthAlong)[i])
#     
#     lines(x, y, lwd= 5, col = cols[i])
#     
#   }
#   par(mar = c(2,4,1,1), las = 2)
#   
#   
#   for(ii in 2:length(levels(parName))){
#     tmpCV <- CVscores[parName==levels(parName)[ii],iii]
#     
#     cols <- tmpCV
#     
#     cols[which(tmpCV >0)] <- "black"
#     cols[which(tmpCV <0)] <- "grey"
#     
#     
#     tmpCV <- abs(tmpCV)
#     
#     
#     names(tmpCV) <- unique(lengthAlong)
#     
#     plot(unique(lengthAlong),tmpCV, pch = 16,
#          xlab = "", axes = T, ylab = levels(parName)[ii], ylim = ylim, col = cols)
#     
# 
#     for(i in 1:length(unique(lengthAlong))){
#       
#       y= c(0 ,tmpCV[i])
#       x= c(unique(lengthAlong)[i], unique(lengthAlong)[i])
#       
#       lines(x, y, lwd= 5,col = cols[i])
#       
#     }
#   }
# }
# 
# 
# 
# 
# # calc total CV loadings for each trait -----------------------------------
# 
# i=1
# 
# mm <- matrix(NA, nrow = 9, ncol = 29)
# 
# CVsums<- cbind(rep(NA, 9),rep(NA, 9),rep(NA, 9))
# rownames(CVsums) <- unique(parName)
# colnames(CVsums) <- colnames(CVAallCSG$CV)
# 
# 
# for(i in 1:9){
#   
#   ii <- (29*(i-1) + 1:29)
#   CVsums[i,]<- apply(abs(CVAallCSG$CV[ii,]), 2, sum)
#   
#   
# }
# 
# CVsums
# 
# 
# par(mfrow = c(1,1), mar = c(7,4,4,1))
# barplot(CVsums, beside = T, names.arg = rep(unique(parName),3), col = palette.colors(n=9))
# 
# 
# barplot(t(CVsums), beside = T)
# legend("topleft", fill = grey.colors(3), legend = c("CV1", "CV2", "CV3"))
# 
# 
# # plot ranked CV scores ---------------------------------------------------
# 
# 
# nrow(CVAallCSG$CV)
# nrow(CVAallCSG$CVvis)
# 
# 
# 
# 
# parName <- as.factor(unlist(lapply(colnames(dat),1 , FUN = rep, 29)))
# 
# parL <- rep(unique(lengthAlong), 7)
# 
# 
# par(mfrow= c(1,3), mar = c(4,10,4,2), las = 2)
# 
# i=1
# 
# for( i in 1: ncol(CVAallCSG$CV)){
#   
#   CVscores <- CVAallCSG$CV[,i]
#   names(CVscores) <- row.names(CVAallCSG$CVvis)
#   o <- order(abs(CVscores), decreasing = F)
#   # o <- order((CVscores), decreasing = F)
#   
#   
#   tmpparName <- parName[o]
#   tmpparL <- parL[o]
#   tmpCV <- CVscores[o]
#   names(tmpCV) <- names(CVscores)[o]
#   
#   SS <- quantile(tmpCV, probs = c(0.10, 0.90))
#   
#   oo <- which(((tmpCV < SS[1]) | (tmpCV > SS[2])))
# 
#   
#   # cols<- palette.colors(length(unique(tmpparName)),"R4")[as.numeric( as.factor(tmpparName[oo]))]
#   
#   
#   cols <- tmpCV
#   
#   cols[which(tmpCV >0)] <- "black"
#   cols[which(tmpCV <0)] <- "grey"
#   
#   
#   # cols <- col2rgb(cols, alpha = tmpparL[oo])
#   
#   colours()
#   
#   barplot(abs(tmpCV[oo]), horiz = T, main = paste("CV", i, " Loadings", sep = ""), col = cols[oo], cex.names=1)
#   
#   
# }
# 
# 
# ###########################################################################
# # PLOT CURVES FOR GROUPMEANS FOR ALLCSG CVA -------------------------------
# 
# 
# # COERCE INTO LONG df AGAIN
# CVAallCSG$groupmeans
# 
# 

