
# calc group means and SD -------------------------------------------------




### calc group means and standard deviation



INDEX <- which(lengthAlong == 0.2)

IDIND <- allDATA$X[INDEX]
colnames(grpMeans20) <- colnames(datIND)
rownames(grpMeans20) <- levels(LOCSUBIND)


grpSD20 <- grpMeans20

for(i in 1:length(levels(LOCSUBIND))){
  
  grpMeans20[i,] <- colMeans(datIND[LOCSUBIND==levels(LOCSUBIND)[i],])
  grpSD20[i,] <- apply(datIND[LOCSUBIND==levels(LOCSUBIND)[i],], 2, sd)
  
}

grpMeans20
grpSD20

INDEX <- which(lengthAlong == 0.8)

datIND <- rawdat[INDEX,]
IDIND <- allDATA$X[INDEX]
LOCSUBIND <- LOCSUB[INDEX]

grpMeans80 <- matrix(ncol = ncol(datIND), nrow = 4)

colnames(grpMeans80) <- colnames(datIND)
rownames(grpMeans80) <- levels(LOCSUBIND)


grpSD80 <- grpMeans80

for(i in 1:length(levels(LOCSUBIND))){
  
  grpMeans80[i,] <- colMeans(datIND[LOCSUBIND==levels(LOCSUBIND)[i],])
  grpSD80[i,] <- apply(datIND[LOCSUBIND==levels(LOCSUBIND)[i],], 2, sd)
  
}

grpMeans80
grpSD80

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

# -------------------------------------------------------------------------

# CAlC GROUP MEANS for TABLE

allDATAwide







# -------------------------------------------------------------------------

load(file = "./CVA CSG scaled midhaft res/CVAallCSG.Rdata")


1:length(unique(lengthAlong))

GroupMeansARR <- array(NA, dim = c(length(unique(lengthAlong)), nrow(dim), nrow(CVAallCSG$groupmeans)))


dim <-matrix(1:ncol(CVAallCSG$groupmeans[,]), ncol = length(unique(lengthAlong)), byrow= T)

for(ii in 1:nrow(CVAallCSG$groupmeans)){
  
  for(i in 1:nrow(dim)){
    
    GroupMeansARR[,i,ii] <- CVAallCSG$groupmeans[ii,dim[i,]]
    
    
  }
  
  
}

dimnames(GroupMeansARR) <- list(unique(lengthAlong), 
                                gsub(colnames(CVAallCSG$groupmeans[,dim[,1]]), pattern = ".0.15", replacement = ""),
                                rownames(CVAallCSG$groupmeans))


which(unique(lengthAlong)==0.2)

GroupMeansARR[which(unique(lengthAlong)==0.2),,]
GroupMeansARR[which(unique(lengthAlong)==0.5),,]
GroupMeansARR[which(unique(lengthAlong)==0.8),,]

ylim <- range(GroupMeansARR/colMeans(GroupMeansARR[which(unique(lengthAlong)==0.5),,]))


par(mfrow = c(4,3), mar = c(5,1,4,1) + 0.1)



for(i in 1:nrow(CVAallCSG$groupmeans)){
  
  # tmp <- cbind(GroupMeansARR[which(unique(lengthAlong)==0.2),,i]-GroupMeansARR[which(unique(lengthAlong)==0.5),,i],
  #              GroupMeansARR[which(unique(lengthAlong)==0.5),,i]-colMeans(GroupMeansARR[which(unique(lengthAlong)==0.5),,]),
  #              GroupMeansARR[which(unique(lengthAlong)==0.8),,i]-GroupMeansARR[which(unique(lengthAlong)==0.5),,i])
  # 
  
  tmp <- cbind(GroupMeansARR[which(unique(lengthAlong)==0.15),,i],
               GroupMeansARR[which(unique(lengthAlong)==0.5),,i],
               GroupMeansARR[which(unique(lengthAlong)==0.85),,i])
  
  # tmp <- tmp - min(tmp)
  
  # 
  # tmp <- cbind(GroupMeansARR[which(unique(lengthAlong)==0.2),,i]-colMeans(GroupMeansARR[which(unique(lengthAlong)==0.5),,]),
  #              GroupMeansARR[which(unique(lengthAlong)==0.5),,i]-colMeans(GroupMeansARR[which(unique(lengthAlong)==0.5),,]),
  #              GroupMeansARR[which(unique(lengthAlong)==0.8),,i]-colMeans(GroupMeansARR[which(unique(lengthAlong)==0.5),,]))
  
  # tmp <- cbind(CVAmidshaft20$groupmeans[i,]-(CVAmidshaft50$groupmeans[i,]),
  #              CVAmidshaft50$groupmeans[i,]-colMeans(CVAmidshaft50$groupmeans),
  #              CVAmidshaft80$groupmeans[i,]-(CVAmidshaft50$groupmeans[i,]))
  
  
  par(mar = c(5,3,4,1) + 0.1)
  
  
  plot(x = unique(lengthAlong), y = GroupMeansARR[,1,i], 
       ylim = range(GroupMeansARR), type = "n", 
       main = dimnames(GroupMeansARR)[[3]][i],
       ylab = "", xlab = "Length Along Diaphysis")
  
  abline(h = 0.5, col = "grey", lty = 1, lwd = 1)
  abline(h = 1, col = "grey", lty = 1, lwd = 1)
  abline(h = 1.5, col = "grey", lty = 1, lwd = 1)
  abline(h = 2, col = "grey", lty = 1, lwd = 1)
  abline(h = 2.5, col = "grey", lty = 1, lwd = 1)
  
  
  
  for(ii in 1:8){
    lines(x = unique(lengthAlong), y = GroupMeansARR[,ii,i], lwd = 2, col = palette.colors(8)[ii])
    points(x = unique(lengthAlong), y = GroupMeansARR[,ii,i], cex = 0.5,col = palette.colors(8)[ii])
    
  }
  
  
  colnames(tmp) <- c("proximal","midshaft","distal")
  
  par(mar = c(5,0,4,0.5) + 0.1)
  
  
  
  barplot(tmp[,1],beside = T, col = palette.colors()[1:8], xpd=F,
                 main = "proximal", plot = T, 
                 space = 0, axes = F, las = 2, ylim = range(GroupMeansARR))  
  
  abline(h = 0.5, col = "grey", lty = 1, lwd = 1)
  abline(h = 1, col = "grey", lty = 1, lwd = 1)
  abline(h = 1.5, col = "grey", lty = 1, lwd = 1)
  abline(h = 2, col = "grey", lty = 1, lwd = 1)
  abline(h = 2.5, col = "grey", lty = 1, lwd = 1)
  
  barplot(tmp[,1],beside = T, col = palette.colors()[1:8], xpd=F,
                main = "proximal", add = T,
                space = 0, axes = F, las = 2, ylim = range(GroupMeansARR))  
  
  
  bp2 <- barplot(tmp[,2],beside = T, col = palette.colors()[1:8], xpd=F,
                 main = "midshaft", 
                 space = 0, axes = F, las = 2, ylim = range(GroupMeansARR))
  
  abline(h = 0.5, col = "grey", lty = 1, lwd = 1)
  abline(h = 1, col = "grey", lty = 1, lwd = 1)
  abline(h = 1.5, col = "grey", lty = 1, lwd = 1)
  abline(h = 2, col = "grey", lty = 1, lwd = 1)
  abline(h = 2.5, col = "grey", lty = 1, lwd = 1)
  # text(cex=1, x=bp1-.25, y=-1.25, labs, xpd=TRUE, srt=45)
  
  bp2 <- barplot(tmp[,2],beside = T, col = palette.colors()[1:8], xpd=F,
                 main = "midshaft", add = T,
                 space = 0, axes = F, las = 2, ylim = range(GroupMeansARR))
  
  bp3 <- barplot(tmp[,3],beside = T, col = palette.colors()[1:8], xpd=F,
                 main = "distal", 
                 space = 0, axes = F, las = 2, ylim = range(GroupMeansARR))  
  
  abline(h = 0.5, col = "grey", lty = 1, lwd = 1)
  abline(h = 1, col = "grey", lty = 1, lwd = 1)
  abline(h = 1.5, col = "grey", lty = 1, lwd = 1)
  abline(h = 2, col = "grey", lty = 1, lwd = 1)
  abline(h = 2.5, col = "grey", lty = 1, lwd = 1)
  # bp1
  
  bp3 <- barplot(tmp[,3],beside = T, col = palette.colors()[1:8], xpd=F,
                 main = "distal", add = T,
                 space = 0, axes = F, las = 2, ylim = range(GroupMeansARR))  
  
  
}


# plot curve --------------------------------------------------------------

par(mfcol = c(4,1), mar = c(5,3,4,1))

i=1

for(i in 1:4){
  plot(x = unique(lengthAlong), y = GroupMeansARR[,1,i], 
       ylim = range(GroupMeansARR), type = "n", 
       main = dimnames(GroupMeansARR)[[3]][i],
       ylab = "", xlab = "Length Along Diaphysis")
  
  for(ii in 1:8){
    lines(x = unique(lengthAlong), y = GroupMeansARR[,ii,i], lwd = 2, col = palette.colors(8)[ii])
    points(x = unique(lengthAlong), y = GroupMeansARR[,ii,i], cex = 0.5, col = palette.colors(8)[ii])
    
  }
  
  
}






# -------------------------------------------------------------------------

which(unique(lengthAlong)==0.5)


match(c(0.2,0.5,0.8), unique(lengthAlong))
ylim <- range(GroupMeansARR[match(c(0.2,0.5,0.8), unique(lengthAlong)),,]-rowMeans(GroupMeansARR[which(unique(lengthAlong)==0.5),,]))


# pdf("./CVA CSG scaled midhaft res/CSG group means at midshaft and distal.pdf", width = 8, height = 12)

par(mfrow = c(4,3), mar = c(5,1,4,1) + 0.1)

i=1
for(i in 1:nrow(CVAallCSG$groupmeans)){
  
  # tmp <- cbind(GroupMeansARR[which(unique(lengthAlong)==0.2),,i]-GroupMeansARR[which(unique(lengthAlong)==0.5),,i],
  #              GroupMeansARR[which(unique(lengthAlong)==0.5),,i]-colMeans(GroupMeansARR[which(unique(lengthAlong)==0.5),,]),
  #              GroupMeansARR[which(unique(lengthAlong)==0.8),,i]-GroupMeansARR[which(unique(lengthAlong)==0.5),,i])
  # 
  
  # tmp <- cbind(GroupMeansARR[which(unique(lengthAlong)==0.2),,i],
  #              GroupMeansARR[which(unique(lengthAlong)==0.5),,i],
  #              GroupMeansARR[which(unique(lengthAlong)==0.8),,i])
  # # 
  tmp <- cbind(GroupMeansARR[which(unique(lengthAlong)==0.2),,i]-rowMeans(GroupMeansARR[which(unique(lengthAlong)==0.5),,]),
               GroupMeansARR[which(unique(lengthAlong)==0.5),,i]-rowMeans(GroupMeansARR[which(unique(lengthAlong)==0.5),,]),
               GroupMeansARR[which(unique(lengthAlong)==0.8),,i]-rowMeans(GroupMeansARR[which(unique(lengthAlong)==0.5),,]))
  
  
  # tmp <- cbind(CVAmidshaft20$groupmeans[i,]-(CVAmidshaft50$groupmeans[i,]),
  #              CVAmidshaft50$groupmeans[i,]-colMeans(CVAmidshaft50$groupmeans),
  #              CVAmidshaft80$groupmeans[i,]-(CVAmidshaft50$groupmeans[i,]))
  
  colnames(tmp) <- c("proximal","midshaft","distal")
  
  
  # bp1<- barplot(apply(M, 2, max), plot = F,
  #               beside = T, col = palette.colors()[1:8], 
  #                main = paste(rownames(CVAallCSG$groupmeans)[i], "(proximal)"), 
  #                space = 0, axes = t, las = 2)  

    bp1<- barplot(tmp[,1], xpd = F,
                  beside = T, col = palette.colors()[1:8], 
                 main = paste(rownames(CVAallCSG$groupmeans)[i], "(proximal)"), 
                 space = 0, axes = F, las = 2, ylim = ylim)  
  
  
  
  
  bp2 <- barplot(tmp[,2],beside = T, col = palette.colors()[1:8], xpd = F,
                 main = paste(rownames(CVAallCSG$groupmeans)[i], "(midshaft)"), 
                 space = 0, axes = F, las = 2, ylim = ylim)
  
  # text(cex=1, x=bp1-.25, y=-1.25, labs, xpd=TRUE, srt=45)
  
  bp3 <- barplot(tmp[,3],beside = T, col = palette.colors()[1:8], xpd = F,
                 main = paste(rownames(CVAallCSG$groupmeans)[i], "(distal)"), 
                 space = 0, axes = F, las = 2, ylim = ylim)  
  
  # bp1
  
}

# dev.off()


# -------------------------------------------------------------------------

# pdf("./CVA CSG scaled midhaft res/CSG group means at midshaft and distal_CHANGE.pdf", width = 8, height = 12)

par(mfrow = c(4,3), mar = c(5,1,4,1) + 0.1)

i=1
for(i in 1:nrow(CVAallCSG$groupmeans)){
  
  tmp <- cbind(GroupMeansARR[which(unique(lengthAlong)==0.2),,i]-GroupMeansARR[which(unique(lengthAlong)==0.5),,i],
               GroupMeansARR[which(unique(lengthAlong)==0.5),,i]-rowMeans(GroupMeansARR[which(unique(lengthAlong)==0.5),,]),
               GroupMeansARR[which(unique(lengthAlong)==0.8),,i]-GroupMeansARR[which(unique(lengthAlong)==0.5),,i])
  # 
  
  # tmp <- cbind(GroupMeansARR[which(unique(lengthAlong)==0.2),,i],
  #              GroupMeansARR[which(unique(lengthAlong)==0.5),,i],
  #              GroupMeansARR[which(unique(lengthAlong)==0.8),,i])
  # # 
  # tmp <- cbind(GroupMeansARR[which(unique(lengthAlong)==0.2),,i]-rowMeans(GroupMeansARR[which(unique(lengthAlong)==0.5),,]),
  #              GroupMeansARR[which(unique(lengthAlong)==0.5),,i]-rowMeans(GroupMeansARR[which(unique(lengthAlong)==0.5),,]),
  #              GroupMeansARR[which(unique(lengthAlong)==0.8),,i]-rowMeans(GroupMeansARR[which(unique(lengthAlong)==0.5),,]))
  
  
  # tmp <- cbind(CVAmidshaft20$groupmeans[i,]-(CVAmidshaft50$groupmeans[i,]),
  #              CVAmidshaft50$groupmeans[i,]-colMeans(CVAmidshaft50$groupmeans),
  #              CVAmidshaft80$groupmeans[i,]-(CVAmidshaft50$groupmeans[i,]))
  
  colnames(tmp) <- c("proximal","midshaft","distal")
  
  
  # bp1<- barplot(apply(M, 2, max), plot = F,
  #               beside = T, col = palette.colors()[1:8], 
  #                main = paste(rownames(CVAallCSG$groupmeans)[i], "(proximal)"), 
  #                space = 0, axes = t, las = 2)  

    bp1<- barplot(tmp[,1], xpd = F,
                  beside = T, col = palette.colors()[1:8], 
                 main = paste(rownames(CVAallCSG$groupmeans)[i], "(proximal)"), 
                 space = 0, axes = F, las = 2, ylim = ylim)  
  
  
  
  
  bp2 <- barplot(tmp[,2],beside = T, col = palette.colors()[1:8], xpd = F,
                 main = paste(rownames(CVAallCSG$groupmeans)[i], "(midshaft)"), 
                 space = 0, axes = F, las = 2, ylim = ylim)
  
  # text(cex=1, x=bp1-.25, y=-1.25, labs, xpd=TRUE, srt=45)
  
  bp3 <- barplot(tmp[,3],beside = T, col = palette.colors()[1:8], xpd = F,
                 main = paste(rownames(CVAallCSG$groupmeans)[i], "(distal)"), 
                 space = 0, axes = F, las = 2, ylim = ylim)  
  
  # bp1
  
}

# dev.off()


# -------------------------------------------------------------------------


load(file = "./CVA CSG unscaled res/CVAallCSG.Rdata")



1:length(unique(lengthAlong))

GroupMeansARR <- array(NA, dim = c(length(unique(lengthAlong)), nrow(dim), nrow(CVAallCSG$groupmeans)))


dim <-matrix(1:ncol(CVAallCSG$groupmeans[,]), ncol = length(unique(lengthAlong)), byrow= T)

for(ii in 1:nrow(CVAallCSG$groupmeans)){
  
  for(i in 1:nrow(dim)){
    
    GroupMeansARR[,i,ii] <- CVAallCSG$groupmeans[ii,dim[i,]]
    
    
  }
  
  
}

dimnames(GroupMeansARR) <- list(unique(lengthAlong), 
                                gsub(colnames(CVAallCSG$groupmeans[,dim[,1]]), pattern = ".0.15", replacement = ""),
                                rownames(CVAallCSG$groupmeans))


which(unique(lengthAlong)==0.2)

GroupMeansARR[which(unique(lengthAlong)==0.2),,]
GroupMeansARR[which(unique(lengthAlong)==0.5),,]
GroupMeansARR[which(unique(lengthAlong)==0.8),,]


GroupMeansARR

ylim <- range(GroupMeansARR[match(c(0.2,0.5,0.8), unique(lengthAlong)),,]-rowMeans(GroupMeansARR[which(unique(lengthAlong)==0.5),,]))
# ylim <- c(0, max(GroupMeansARR[match(c(0.2,0.5,0.8), unique(lengthAlong)),,]rowMeans(GroupMeansARR[which(unique(lengthAlong)==0.5),,])))


# pdf("./CVA CSG scaled midhaft res/CSG group means at midshaft and distal",

# pdf("./CVA CSG unscaled res/CSG group means at midshaft and distal.pdf",  width = 8, height = 12)

par(mfrow = c(4,3), mar = c(5,1,4,1) + 0.1)


i=1
for(i in 1:nrow(CVAallCSG$groupmeans)){
  

    
  tmp <- cbind(GroupMeansARR[which(unique(lengthAlong)==0.2),,i]-rowMeans(GroupMeansARR[which(unique(lengthAlong)==0.5),,]),
               GroupMeansARR[which(unique(lengthAlong)==0.5),,i]-rowMeans(GroupMeansARR[which(unique(lengthAlong)==0.5),,]),
               GroupMeansARR[which(unique(lengthAlong)==0.8),,i]-rowMeans(GroupMeansARR[which(unique(lengthAlong)==0.5),,]))
  
  

  
  colnames(tmp) <- c("proximal","midshaft","distal")
  
  
  bp1<- barplot(tmp[,1],beside = T, col = palette.colors()[1:8], 
                main = paste(rownames(CVAallCSG$groupmeans)[i], "(proximal)"), 
                space = 0, axes = F, las = 2, ylim = ylim)  
  
  bp2 <- barplot(tmp[,2],beside = T, col = palette.colors()[1:8], 
                 main = paste(rownames(CVAallCSG$groupmeans)[i], "(midshaft)"), 
                 space = 0, axes = F, las = 2, ylim = ylim)
  
  # text(cex=1, x=bp1-.25, y=-1.25, labs, xpd=TRUE, srt=45)
  
  bp3 <- barplot(tmp[,3],beside = T, col = palette.colors()[1:8], 
                 main = paste(rownames(CVAallCSG$groupmeans)[i], "(distal)"), 
                 space = 0, axes = F, las = 2, ylim = ylim)  
  

  
}

# dev.off()



# -------------------------------------------------------------------------




# pdf("./CVA CSG unscaled res/CSG group means at midshaft and distal_min.pdf", width = 8, height = 12)

par(mfrow = c(4,3), mar = c(5,1,4,1) + 0.1)


i=1
for(i in 1:nrow(CVAallCSG$groupmeans)){

  tmp <- cbind(GroupMeansARR[which(unique(lengthAlong)==0.2),,i],
               GroupMeansARR[which(unique(lengthAlong)==0.5),,i],
               GroupMeansARR[which(unique(lengthAlong)==0.8),,i])
  
  tmp <- tmp - min(tmp)
 
  
  colnames(tmp) <- c("proximal","midshaft","distal")
  
  
  bp1<- barplot(tmp[,1],beside = T, col = palette.colors()[1:8], 
                main = paste(rownames(CVAallCSG$groupmeans)[i], "(proximal)"), 
                space = 0, axes = F, las = 2, ylim = range(tmp))  
  
  bp2 <- barplot(tmp[,2],beside = T, col = palette.colors()[1:8], 
                 main = paste(rownames(CVAallCSG$groupmeans)[i], "(midshaft)"), 
                 space = 0, axes = F, las = 2, ylim = range(tmp))
  
  # text(cex=1, x=bp1-.25, y=-1.25, labs, xpd=TRUE, srt=45)
  
  bp3 <- barplot(tmp[,3],beside = T, col = palette.colors()[1:8], 
                 main = paste(rownames(CVAallCSG$groupmeans)[i], "(distal)"), 
                 space = 0, axes = F, las = 2, ylim = range(tmp))  
  
  # bp1
  
}

# dev.off()


# -------------------------------------------------------------------------



# pdf("./CVA CSG unscaled res/CSG group means at midshaft and distal_RAW.pdf",
#     # pdf("./CVA CSG scaled midhaft res/CSG group means at midshaft and distal",
#     width = 8, height = 12)

par(mfrow = c(4,3), mar = c(5,1,4,1) + 0.1)


i=1
for(i in 1:nrow(CVAallCSG$groupmeans)){

  
  max(GroupMeansARR[,,])
  
  
  tmp <- cbind(GroupMeansARR[which(unique(lengthAlong)==0.2),,i],
               GroupMeansARR[which(unique(lengthAlong)==0.5),,i],
               GroupMeansARR[which(unique(lengthAlong)==0.8),,i])

  
  colnames(tmp) <- c("proximal","midshaft","distal")
  
  
  bp1<- barplot(tmp[,1],beside = T, col = palette.colors()[1:8], 
                main = paste(rownames(CVAallCSG$groupmeans)[i], "(proximal)"), 
                space = 0, axes = F, las = 2, ylim = c(0, max(GroupMeansARR[,,])))  
  
  bp2 <- barplot(tmp[,2],beside = T, col = palette.colors()[1:8], 
                 main = paste(rownames(CVAallCSG$groupmeans)[i], "(midshaft)"), 
                 space = 0, axes = F, las = 2, ylim = c(0, max(GroupMeansARR[,,])))
  

  bp3 <- barplot(tmp[,3],beside = T, col = palette.colors()[1:8], 
                 main = paste(rownames(CVAallCSG$groupmeans)[i], "(distal)"), 
                 space = 0, axes = F, las = 2, ylim = c(0, max(GroupMeansARR[,,])))  
  

}

# dev.off()

