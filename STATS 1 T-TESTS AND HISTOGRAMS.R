
###########################################################################
# T-TESTS -----------------------------------------------------------------
###########################################################################

rm(list = ls())

load("boneVars.Rdata")

eco <- as.data.frame(read.csv("./DATA/specimenEco2.csv"))


# PLOT HISTOGRAMS ---------------------------------------------------------

par(mfcol = c(8,5), mar = c(5,4,1,2))

# MIDSHAFT GSC VS MAX CSG

boneVars$mid40DF
boneVars$mid50DF
boneVars$maxDF


diff40 <- (boneVars$maxDF/boneVars$ecoDF$MIDSHAFTCIRC) - (boneVars$mid40DF/boneVars$ecoDF$MIDSHAFTCIRC)
diff50 <- (boneVars$maxDF/boneVars$ecoDF$MIDSHAFTCIRC) - (boneVars$mid50DF/boneVars$ecoDF$MIDSHAFTCIRC)
diff4050 <- (boneVars$mid40DF/boneVars$ecoDF$MIDSHAFTCIRC) - (boneVars$mid50DF/boneVars$ecoDF$MIDSHAFTCIRC)



varnames <- c("Curvature", "m. Thickness", "Cort. Area", "Ix","Iy", "I.ratio", "Zpol", "J")

colnames(diff40) <- varnames
colnames(diff50) <- varnames
colnames(diff4050) <- varnames



i=1
for(i in 1:ncol(boneVars$mid40DF)){

  # hist(maxDF[,i], col=rgb(0.5,0.5,0.5,1/4), breaks = 10, main = varnames[i])
  # hist(mid40DF[,i], col = rgb(0,0,1,1/4), add =T)
  # hist(mid50DF[,i], col=rgb(1,0,0,1/4), add = T)
  
  hist(diff4050[,i], breaks = 20, xlab = "dif. 40-50", main = varnames[i])
  # hist(diff40[,i], breaks = 20, xlab = "dif. max-40", main = varnames[i])
  # hist(diff40[,i], breaks = 20, xlab = "dif. max-50", main = varnames[i])



  
}

for(i in 1:ncol(boneVars$mid40DF)){

  # hist(maxDF[,i], col=rgb(0.5,0.5,0.5,1/4), breaks = 10, main = varnames[i])
  # hist(mid40DF[,i], col = rgb(0,0,1,1/4), add =T)
  # hist(mid50DF[,i], col=rgb(1,0,0,1/4), add = T)
  
  # hist(diff4050[,i], breaks = 20, xlab = "dif. 40-50", main = varnames[i])
  hist(diff40[,i], breaks = 20, xlab = "dif. max-40", main = varnames[i])
  # hist(diff40[,i], breaks = 20, xlab = "dif. max-50", main = varnames[i])



  
}

for(i in 1:ncol(boneVars$mid40DF)){

  # hist(maxDF[,i], col=rgb(0.5,0.5,0.5,1/4), breaks = 10, main = varnames[i])
  # hist(mid40DF[,i], col = rgb(0,0,1,1/4), add =T)
  # hist(mid50DF[,i], col=rgb(1,0,0,1/4), add = T)
  
  # hist(diff4050[,i], breaks = 20, xlab = "dif. 40-50", main = varnames[i])
  # hist(diff40[,i], breaks = 20, xlab = "dif. max-40", main = varnames[i])
  hist(diff50[,i], breaks = 20, xlab = "dif. max-50", main = varnames[i])



  
}



diff4050T.tests <- lapply(X= diff4050, FUN = wilcox.test, alternative = "g" , conf.int = T, mu= 0)
diff40T.tests <- lapply(X= diff40, FUN = wilcox.test, alternative = "g" , conf.int = T, mu= 0)
diff50T.tests <- lapply(X= diff50, FUN = wilcox.test, alternative = "g" , conf.int = T, mu= 0)




sink("./stats/T-TESTS 40-50-max.txt")
print("difference between positions 40-50")
print(diff4050T.tests)

print("difference between position 40 and MAX")
print(diff40T.tests)

print("difference between position 50 and MAX")
print(diff50T.tests)

sink()


p.valDF <- NULL

for(i in 1:length(diff4050T.tests)){
  p.vals <- c(diff4050T.tests[[i]]$p.value,
  diff40T.tests[[i]]$p.value,
  diff50T.tests[[i]]$p.value)
  
  p.valDF <- cbind(p.valDF, p.vals)
  
}

colnames(p.valDF) <- names(diff4050T.tests)
rownames(p.valDF) <- c("diff4050", "diff40max", "diff50max")

t(p.valDF)

p.valDF <- round(p.valDF, 3) 

p.valDF[which(p.valDF < 0.001)] <- "< 0.001"




write.csv(cbind(round(colMeans(boneVars$posDF[,]), 3), t(p.valDF)), "./stats/wilcoxon pvals for difference between positions.csv")




# -------------------------------------------------------------------------


# eccentricity

n <- nrow(boneVars$ecoDF)


D <- boneVars$posDF

# par(mfcol = c(8,1), mar = c(5,4,1,2))


# for(i in 1: ncol(D)){
#   
#   hist(D[,i], main = "", xlab = colnames(D)[i], breaks = 20)
#   
# }
# 
# # D <- boneVars$ecc50DF
# # 
# # for(i in 1: ncol(D)){
# # 
# #   hist(D[,i], main = "", xlab = colnames(D)[i], breaks = 20)
# #   
# # }

D <- boneVars$abs50DF

for(i in 1: ncol(D)){
  
  hist(D[,i], main = varnames[i], xlab = "", las = 2 ,breaks = 20)
  
}

# 

abs40T.tests <- lapply(X= boneVars$abs40DF, FUN = wilcox.test, alternative = "g" , conf.int = T, mu= 0)
abs50T.tests <- lapply(X= boneVars$abs50DF, FUN = wilcox.test, alternative = "g", mu= 0)

eccT.tests<- c(abs40T.tests, abs50T.tests)

sink("./stats/wilcoxon tests for essentricity.txt")
eccT.tests
sink()



p.valDF <- NULL

for(i in 1:length(diff4050T.tests)){
  p.vals <- c(abs40T.tests[[i]]$p.value,
              abs50T.tests[[i]]$p.value)
  
  p.valDF <- cbind(p.valDF, p.vals)
  
}

colnames(p.valDF) <- names(diff4050T.tests)
rownames(p.valDF) <- c("abs40T.tests", "abs50T.tests")

t(p.valDF)

p.valDF <- round(p.valDF, 3) 

p.valDF[which(p.valDF < 0.001)] <- "< 0.001"


colMeans(boneVars$posDF[which(eco$LOCSUB == levels(as.factor(eco$LOCSUB))[1]),])
colMeans(boneVars$posDF[which(eco$LOCSUB == levels(as.factor(eco$LOCSUB))[1]),])
colMeans(boneVars$posDF[which(eco$LOCSUB == levels(as.factor(eco$LOCSUB))[1]),])
colMeans(boneVars$posDF[which(eco$LOCSUB == levels(as.factor(eco$LOCSUB))[1]),])




write.csv(cbind(round(colMeans(boneVars$abs50DF[,]), 3), t(p.valDF)), "./stats/wilcoxon pvals for eccentricity.csv")



# variance ----------------------------------------------------------------



# variance

n <- nrow(boneVars$ecoDF)


D <- boneVars$varDF

# par(mfrow = c(2,4), mar = c(5,4,1,2))


# for(i in 1: ncol(D)){
#   
#   hist(D[,i], main = "", xlab = colnames(D)[i], breaks = 20)
#   
# }
# 
# # D <- boneVars$ecc50DF
# # 
# # for(i in 1: ncol(D)){
# # 
# #   hist(D[,i], main = "", xlab = colnames(D)[i], breaks = 20)
# #   
# # }


for(i in 1: ncol(D)){
  
  hist(D[,i], xlab = "", main = varnames[i], las = 2 ,breaks = 20)
  
}

# 

var40T.tests <- lapply(X= boneVars$abs40DF, FUN = wilcox.test, alternative = "g" , conf.int = T, mu= 0)
var50T.tests <- lapply(X= boneVars$abs50DF, FUN = wilcox.test, alternative = "g", mu= 0)

varT.tests<- c(var40T.tests, var50T.tests)

sink("./stats/wilcoxon tests for essentricity.txt")
varT.tests
sink()



p.valDF <- NULL

for(i in 1:length(diff4050T.tests)){
  p.vals <- c(var40T.tests[[i]]$p.value,
              var50T.tests[[i]]$p.value)
  
  p.valDF <- cbind(p.valDF, p.vals)
  
}

colnames(p.valDF) <- names(diff4050T.tests)
rownames(p.valDF) <- c("var40T.tests", "var50T.tests")

t(p.valDF)

p.valDF <- round(p.valDF, 3) 

p.valDF[which(p.valDF < 0.001)] <- "< 0.001"

write.csv(t(p.valDF), "./stats/wilcoxon pvals for eccentricity.csv")




