rm(list = (ls()))


gc()

library(morphomap)
load("./morphomapSnapshots85Procrustes/AlignedMeshes.Rdata")


# load("./InertialAligned.Rdata")
load("./morphomapSnapshots85Procrustes/allCSG85Procrustes_centroid.Rdata")
load("./morphomapSnapshots85Procrustes/allLands85Procrustes.Rdata")

# eco <- as.data.frame(read.csv("./DATA/specimenEco.csv"))
# eco <- eco[index,]
# eco


# eco2 <- data.frame(eco, LOCSUB = as.factor(paste(eco$SUBSTRATE, eco$LOCOMOTION)))

# write.csv(eco2, "./DATA/specimenEco2.csv")


eco <- read.csv("./DATA/specimenEco2.csv")

str(eco)
as.factor(eco[,2])

eco <- data.frame(lapply(eco, as.factor))

str(eco)

# convert to long table ---------------------------------------------------




COLNAMES <- colnames(allCSG[[1]])

COLNAMES
COLNAMES[1] <- "lengthAlong" 
COLNAMES[26] <- "SegCenX" 
COLNAMES[27] <- "SegCenY" 
COLNAMES[28] <- "SegCenZ"
COLNAMES[30] <- "DistSegment"
COLNAMES[31] <- "DistCentroid"
COLNAMES

lengthAlong <- allCSG[[1]][[1]]

allDATA <- list()
i=1


YY <- NULL


for(i in 1:length(allCSG)){
  
  colnames(allCSG[[i]]) <-COLNAMES
  
  
  tmp <- allCSG[[i]]
  
  midShaft <- matrix(tmp[15,], nrow = nrow(tmp), ncol = length(tmp[15,]), byrow=T)
  colnames(midShaft) <- paste("MidShaft." ,colnames(tmp), sep="")
  
  Max <- c(apply(tmp, 2, max))
  MaxMAT<- matrix(rep(Max, nrow(tmp)), nrow = 29, byrow =T )
  colnames(MaxMAT) <- paste("Max." ,names(Max), sep="")
  
  maxIND <- which.max(tmp$DistCentroid)
  maxCpos <- unlist(tmp[maxIND,])
  maxCposMAT <- matrix(rep(maxCpos, nrow(tmp)), nrow = 29, byrow = T )
  colnames(maxCposMAT) <- paste("MaxCpos." ,colnames(tmp), sep="")
  
  
  
  # ScaleMidCirc <- median(tmp$Ext_perim)
  
  XX <- cbind(eco[i,], tmp, midShaft)
  
  XX <- cbind(XX, 
              MaxMAT,
              maxCposMAT)
  
  YY <-   rbind(YY, XX)
  
  
  
  
}

colnames(YY)


allDATA <- YY

str(YY)

save(allDATA, file = "allDATA.Rdata")



##########################
########################
#  rewrite plots and stats for long table ---------------------------------


load("allDATA.Rdata")

load("./morphomapSnapshots85Procrustes/AlignedMeshes.Rdata")

# load("./newAlign.Rdata")

A <- newAlign$LANDMESH$MESHES

# plot lengthAlong data using latice

allDATA$Ix.Ix
allDATA$Iy.Iy
i=1


ZLength = NULL

i=1
for(i in 1:length(A)){
  ZLength = c(ZLength, as.numeric(dist(range(A[[i]]$mesh[[1]]$vb[3,]))))
}

ZLength

save(ZLength, file = "ZLength.RData")

i=1
length(lengthAlong)

ZlengthDF <- NULL

for(i in 1:length(ZLength)){
  
  ZlengthDF <- c(ZlengthDF ,rep(ZLength[i], length(lengthAlong)))
  
  
}


CCC <- unlist(allDATA$MidShaft.Ext_perim)

dat <- data.frame(X = allDATA$X,
                  TAXA = allDATA$TAXA,
                  GENUS = allDATA$GENUS,
                  LOCOMOTION = allDATA$LOCOMOTION,
                  SUBSTRATE = allDATA$SUBSTRATE,
                  GRP = allDATA$GRP,
                  MIDSHAFTCIRC = unlist(allDATA$MidShaft.Ext_perim),
                  ZLENGTH = ZlengthDF,
                  lengthAlong = allDATA$lengthAlong,
                  Curvature = allDATA$DistCentroid,
                  thickness = allDATA$Mean_thick,
                  CA = allDATA$CA,
                  Ix = allDATA$Ix.Ix,
                  Iy = allDATA$Iy.Iy,
                  I.ratio = (allDATA$Imax/allDATA$Imin),
                  Zpol = allDATA$Zpol.Zpol,
                  J = allDATA$J.Ix
)



dat


maxCol <- paste("max" ,colnames(dat[10:15]))
varCol <- paste("var" ,colnames(dat[10:15]))
eccCol <- paste("ecc" ,colnames(dat[10:15]))


ecoDF <-NULL
meanDF <- NULL
maxDF <- NULL
minDF <- NULL
varDF <- NULL
posDF <- NULL
ecc50DF <- NULL
ecc40DF <- NULL
abs50DF <- NULL
abs40DF <- NULL
mid50DF <- NULL
mid40DF <- NULL


i=1
for(i in 1:length(levels(dat$X))){
  
  ID  <- levels(dat$X)[i]
  tmp <- dat[dat$X == ID,] 
  
  TAXA <- tmp$TAXA[1]
  LOCOMOTION <- tmp$LOCOMOTION[1]
  SUBSTRATE <- tmp$SUBSTRATE[1]
  GRP <- tmp$GRP[1]
  MIDSHAFTCIRC <- tmp$MIDSHAFTCIRC[1]
  ZLENGTH <- ZLength[i]
  
  V <- tmp[,10:ncol(tmp)]
  
  ecoDF <- rbind(ecoDF, data.frame(ID, TAXA, LOCOMOTION, SUBSTRATE, GRP,MIDSHAFTCIRC,ZLENGTH))
  
  # V[which(tmp$lengthAlong == 0.5),]
  
  
  mid50VEC <- V[which(tmp$lengthAlong == 0.5),]
  names(mid50VEC) <- paste("MID50", names(mid50VEC))
  mid50DF <- rbind(mid50DF, mid50VEC)
  rownames(mid50DF)[i] <- ID
  
  
  mid40VEC <- V[11,]
  names(mid40VEC) <- paste("MID40", names(mid40VEC))
  mid40DF <- rbind(mid40DF, mid40VEC)
  rownames(mid40DF)[i] <- ID
  
  
  
  meanVec <- apply(V,2, mean)
  names(meanVec) <- paste("mean", names(meanVec))
  meanDF <- rbind(meanDF, meanVec)
  rownames(meanDF)[i] <- ID
  
  
  maxVec <- apply(V,2, max)
  names(maxVec) <- paste("max", names(maxVec))
  maxDF <- rbind(maxDF, maxVec)
  rownames(maxDF)[i] <- ID
  
  
  minVec <- apply(V,2, max)
  names(minVec) <- paste("min", names(minVec))
  minDF <- rbind(minDF, minVec)
  rownames(minDF)[i] <- ID
  
  
  varVec <- apply(V,2, var)
  names(varVec) <- paste("var", names(varVec))
  varDF <- rbind(varDF, varVec)
  rownames(varDF)[i] <- ID
  
  
  posVec <- apply(V,2, function(X){
    as.numeric(lengthAlong[which.max(X)])
  }
  )
  names(posVec) <- paste("pos", names(posVec))
  posDF <- rbind(posDF, posVec)
  rownames(posDF)[i] <- ID
  
  
  
  ecc40Vec <- apply(V,2, function(X){
    as.numeric(lengthAlong[which.max(X)] - 0.4)
  }
  )
  names(ecc40Vec) <- paste("ecc40", names(ecc40Vec))
  ecc40DF <- rbind(ecc40DF, ecc40Vec)
  rownames(ecc40DF)[i] <- ID
  
  
  ecc50Vec <- apply(V,2, function(X){
    as.numeric(lengthAlong[which.max(X)] - 0.5)
  }
  )
  names(ecc50Vec) <- paste("ecc50", names(ecc50Vec))
  ecc50DF <- rbind(ecc50DF, ecc50Vec)
  rownames(ecc50DF)[i] <- ID
  
  
  abs50Vec <- abs(ecc50Vec)
  names(abs50Vec) <- paste("abs", names(abs50Vec))
  abs50DF <- rbind(abs50DF, abs50Vec)
  rownames(abs50DF)[i] <- ID
  
  abs40Vec <- abs(ecc40Vec)
  names(abs40Vec) <- paste("abs", names(abs40Vec))
  abs40DF <- rbind(abs40DF, abs40Vec)
  rownames(abs40DF)[i] <- ID
  
  
  
}

boneVars <- list( ecoDF ,  data.frame(mid50DF), data.frame(mid40DF), data.frame(meanDF) ,  data.frame(maxDF) ,  data.frame(minDF) , data.frame(varDF), data.frame(posDF),  data.frame(ecc40DF) ,  data.frame(ecc50DF) , data.frame(abs40DF) ,  data.frame(abs50DF))
names(boneVars) <- c("ecoDF", "mid50DF", "mid40DF","meanDF", "maxDF", "minDF", "varDF", "posDF", "ecc40DF", "ecc50DF","abs40DF" ,"abs50DF")

save(boneVars, file = "boneVars.Rdata")

