
load("./allDATA.rdata")
lengthAlong = allDATA$lengthAlong


###############################################################
# EFFECTS OF ECOLOGICAL FACTORS -------------------------------------------
######################################################

###########################################################################
# PREDICTING ECOLOGY FROM TRAITS ------------------------------------------
LOCSUB <- as.factor(paste(allDATA$LOCOMOTION, allDATA$SUBSTRATE))
dat <- data.frame(Curvature = allDATA$DistCentroid,
                  Thickness = allDATA$Mean_thick,
                  Circumferance = allDATA$Ext_perim,
                  CA = allDATA$CA,
                  I.ratio = (allDATA$Imax/allDATA$Imin),
                  Zpol = allDATA$Zpol.Zpol,
                  J = allDATA$J.Ix
)

lengthAlong <- allDATA$lengthAlong


# SIZES <- boneVars$ecoDF[,c(6:7)]

# data <- data.frame(, DF, SIZES$MIDSHAFTCIRC )




MANOVA <- adonis2(dat~LOCSUB+lengthAlong)

MANOVA

summary(MANOVA)

# q: How well do traits, and positions predict ecology?

#CVA DOES THIS










# locomotion --------------------------------------------------------------
palette(c("firebrick", "greenyellow", "skyblue"))


res_data



boneVars$residuals <-res_data



# MAX ANOVAS --------------------------------------------------------------

ECO <- boneVars$ecoDF[,c(3:5)]
DF <- boneVars$maxDF
SIZES <- boneVars$ecoDF[,c(6:7)]




sink("./stats/MAX anovaMODELS.txt")


anovaMODELS <-list()

for(i in 1:ncol(DF)){
  
  anovaMODELS[[i]] <- list()
  names(anovaMODELS)[i] <- colnames(DF)[i]
  ECOLIST <-list()
  
  
  
  X <- DF[,i]
  
  
  data <- data.frame(X, ECO, SIZE, as.factor(paste(data$LOCOMOTION, data$SUBSTRATE)))
  pars <- c(colnames(DF)[i],colnames(ECO),  colnames(SIZES)[1], "LOCSUB")
  colnames(data) <- pars
  
  model1 <- paste(pars[1], "~", pars[2], sep = "")
  print(paste("MODEL =", model1))
  ANOVA1 <- lm(log(X)~LOCOMOTION, data = data)
  summ1 <- summary(ANOVA1)
  print(summ1)
  tukey1 <- TukeyHSD(aov(X~LOCOMOTION, data = data))  
  print(tukey1)
  
  
  model2 <- paste(pars[1], "~", pars[3], sep = "")
  print(paste("MODEL =", model2))
  ANOVA2 <- lm(log(X)~SUBSTRATE, data = data)
  summ2 <- summary(ANOVA2)
  print(summ2)
  tukey2 <-  TukeyHSD(aov(X~SUBSTRATE, data = data))  
  print(tukey2)
  
  model3 <- paste(pars[1], "~", pars[2], "*", pars[3], sep = "")
  print(paste("MODEL =", model3))
  ANOVA3 <- lm(log(X)~LOCOMOTION*SUBSTRATE, data = data)
  summ3 <- summary(ANOVA3)
  print(summ3)
  tukey3 <- TukeyHSD(aov(X~LOCSUB, data = data))  
  print(tukey3)
  
  model4 <- paste(pars[1], "~", pars[4], sep = "")
  print(paste("MODEL =", model4))
  ANOVA4 <- lm(log(X)~GRP, data = data)
  summ4 <- summary(ANOVA4)
  print(summ4)
  tukey4 <- TukeyHSD(aov(X~GRP, data = data))  
  print(tukey4)
  
  model5 <- paste(pars[1], "~", pars[2], "*", pars[5], sep = "")
  print(paste("MODEL =", model5))
  ANOVA5 <- lm(log(X)~LOCOMOTION*log(MIDSHAFTCIRC), data = data)
  summ5 <- summary(ANOVA5)
  print(summ5)
  
  model6 <- paste(pars[1], "~", pars[3], "*", pars[5], sep = "")
  print(paste("MODEL =", model6))
  ANOVA6 <- lm(log(X)~SUBSTRATE*log(MIDSHAFTCIRC), data = data)
  summ6 <- summary(ANOVA6)
  print(summ6)
  
  model7 <- paste(pars[1], "~", pars[4], "*", pars[5], sep = "")
  print(paste("MODEL =", model7))
  ANOVA7 <- lm(log(X)~GRP*log(MIDSHAFTCIRC), data = data)
  summ7 <- summary(ANOVA7)
  print(summ7)
  
  model8 <- paste(pars[1], "~", pars[2], "*",pars[3], "*", pars[5], sep = "")
  print(paste("MODEL =", model8))
  ANOVA8 <- lm(log(X)~LOCOMOTION*SUBSTRATE*log(MIDSHAFTCIRC), data = data)
  summ8 <- summary(ANOVA8)
  print(summ8)
  
  model9 <- paste(pars[1], "~", pars[2], "*",pars[3], "*", pars[4], "*", pars[5], sep = "")
  print(paste("MODEL =", model9))
  ANOVA9 <- lm(log(X)~LOCOMOTION*SUBSTRATE*GRP*log(MIDSHAFTCIRC), data = data)
  summ9 <- summary(ANOVA9)
  print(summ9)
  
  MODELS<- data.frame("model" = c(model1, model2, model3, 
                                  model4, model5, model6, 
                                  model7, model8, model9),
                      "AdjRsq" = c(summ1$adj.r.squared, 
                                   summ2$adj.r.squared, 
                                   summ3$adj.r.squared, 
                                   summ4$adj.r.squared, 
                                   summ5$adj.r.squared, 
                                   summ6$adj.r.squared, 
                                   summ7$adj.r.squared, 
                                   summ8$adj.r.squared, 
                                   summ9$adj.r.squared),
                      "AIC" = c(AIC(ANOVA1), AIC(ANOVA2), AIC(ANOVA3),
                                AIC(ANOVA4), AIC(ANOVA5), AIC(ANOVA6),
                                AIC(ANOVA7), AIC(ANOVA8), AIC(ANOVA9))
  )
  
  anovaMODELS[[i]] <- MODELS
}


sink()

anovaMODELS <- lapply(anovaMODELS, FUN = function(X){
  
  ret <- X[order(X$AIC),]
  rownames(ret) <- c(1:nrow(X))
  return = ret
  
} 
)

save(anovaMODELS,file = "./stats/MAX anovaMODELS.Rdata")

sink("./stats/MAX AIC anovaMODELS.txt")
print(anovaMODELS)
sink()

###### ADD ANOVA LOOPS FOR OTHER VARIABLES (VAR, ECC)



# VAR ANOVAS --------------------------------------------------------------



ECO <- boneVars$ecoDF[,c(3:5)]
DF <- boneVars$varDF
SIZES <- boneVars$ecoDF[,c(6:7)]




sink("./stats/VAR anovaMODELS.txt")


anovaMODELS <-list()

for(i in 1:ncol(DF)){
  
  anovaMODELS[[i]] <- list()
  names(anovaMODELS)[i] <- colnames(DF)[i]
  ECOLIST <-list()
  
  
  
  X <- DF[,i]
  
  
  data <- data.frame(X, ECO, SIZE, as.factor(paste(data$LOCOMOTION, data$SUBSTRATE)))
  pars <- c(colnames(DF)[i],colnames(ECO),  colnames(SIZES)[1], "LOCSUB")
  colnames(data) <- pars
  
  model1 <- paste(pars[1], "~", pars[2], sep = "")
  print(paste("MODEL =", model1))
  ANOVA1 <- lm(log(X)~LOCOMOTION, data = data)
  summ1 <- summary(ANOVA1)
  print(summ1)
  tukey1 <- TukeyHSD(aov(X~LOCOMOTION, data = data))  
  print(tukey1)
  
  
  model2 <- paste(pars[1], "~", pars[3], sep = "")
  print(paste("MODEL =", model2))
  ANOVA2 <- lm(log(X)~SUBSTRATE, data = data)
  summ2 <- summary(ANOVA2)
  print(summ2)
  tukey2 <-  TukeyHSD(aov(X~SUBSTRATE, data = data))  
  print(tukey2)
  
  model3 <- paste(pars[1], "~", pars[2], "*", pars[3], sep = "")
  print(paste("MODEL =", model3))
  ANOVA3 <- lm(log(X)~LOCOMOTION*SUBSTRATE, data = data)
  summ3 <- summary(ANOVA3)
  print(summ3)
  tukey3 <- TukeyHSD(aov(X~LOCSUB, data = data))  
  print(tukey3)
  
  model4 <- paste(pars[1], "~", pars[4], sep = "")
  print(paste("MODEL =", model4))
  ANOVA4 <- lm(log(X)~GRP, data = data)
  summ4 <- summary(ANOVA4)
  print(summ4)
  tukey4 <- TukeyHSD(aov(X~GRP, data = data))  
  print(tukey4)
  
  model5 <- paste(pars[1], "~", pars[2], "*", pars[5], sep = "")
  print(paste("MODEL =", model5))
  ANOVA5 <- lm(log(X)~LOCOMOTION*log(MIDSHAFTCIRC), data = data)
  summ5 <- summary(ANOVA5)
  print(summ5)
  
  model6 <- paste(pars[1], "~", pars[3], "*", pars[5], sep = "")
  print(paste("MODEL =", model6))
  ANOVA6 <- lm(log(X)~SUBSTRATE*log(MIDSHAFTCIRC), data = data)
  summ6 <- summary(ANOVA6)
  print(summ6)
  
  model7 <- paste(pars[1], "~", pars[4], "*", pars[5], sep = "")
  print(paste("MODEL =", model7))
  ANOVA7 <- lm(log(X)~GRP*log(MIDSHAFTCIRC), data = data)
  summ7 <- summary(ANOVA7)
  print(summ7)
  
  model8 <- paste(pars[1], "~", pars[2], "*",pars[3], "*", pars[5], sep = "")
  print(paste("MODEL =", model8))
  ANOVA8 <- lm(log(X)~LOCOMOTION*SUBSTRATE*log(MIDSHAFTCIRC), data = data)
  summ8 <- summary(ANOVA8)
  print(summ8)
  
  model9 <- paste(pars[1], "~", pars[2], "*",pars[3], "*", pars[4], "*", pars[5], sep = "")
  print(paste("MODEL =", model9))
  ANOVA9 <- lm(log(X)~LOCOMOTION*SUBSTRATE*GRP*log(MIDSHAFTCIRC), data = data)
  summ9 <- summary(ANOVA9)
  print(summ9)
  
  MODELS<- data.frame("model" = c(model1, model2, model3, 
                                  model4, model5, model6, 
                                  model7, model8, model9),
                      "AdjRsq" = c(summ1$adj.r.squared, 
                                   summ2$adj.r.squared, 
                                   summ3$adj.r.squared, 
                                   summ4$adj.r.squared, 
                                   summ5$adj.r.squared, 
                                   summ6$adj.r.squared, 
                                   summ7$adj.r.squared, 
                                   summ8$adj.r.squared, 
                                   summ9$adj.r.squared),
                      "AIC" = c(AIC(ANOVA1), AIC(ANOVA2), AIC(ANOVA3),
                                AIC(ANOVA4), AIC(ANOVA5), AIC(ANOVA6),
                                AIC(ANOVA7), AIC(ANOVA8), AIC(ANOVA9))
  )
  
  anovaMODELS[[i]] <- MODELS
}

sink()

anovaMODELS <- lapply(anovaMODELS, FUN = function(X){
  
  ret <- X[order(X$AIC),]
  rownames(ret) <- c(1:nrow(X))
  return = ret
  
} 
)

save(anovaMODELS,file = "./stats/VAR AIC anovaMODELS.Rdata")

sink("./stats/VAR AIC anovaMODELS.txt")
print(anovaMODELS)
sink()




# ECC POS ANOVAS --------------------------------------------------------------



ECO <- boneVars$ecoDF[,c(3:5)]
DF <- boneVars$posDF
SIZES <- boneVars$ecoDF[,c(6:7)]




sink("./stats/ECC POS anovaMODELS.txt")


anovaMODELS <-list()

for(i in 1:ncol(DF)){
  
  anovaMODELS[[i]] <- list()
  names(anovaMODELS)[i] <- colnames(DF)[i]
  ECOLIST <-list()
  
  
  
  X <- DF[,i]
  
  
  data <- data.frame(X, ECO, SIZE, as.factor(paste(data$LOCOMOTION, data$SUBSTRATE)))
  pars <- c(colnames(DF)[i],colnames(ECO),  colnames(SIZES)[1], "LOCSUB")
  colnames(data) <- pars
  
  model1 <- paste(pars[1], "~", pars[2], sep = "")
  print(paste("MODEL =", model1))
  ANOVA1 <- lm(log(X)~LOCOMOTION, data = data)
  summ1 <- summary(ANOVA1)
  print(summ1)
  tukey1 <- TukeyHSD(aov(X~LOCOMOTION, data = data))  
  print(tukey1)
  
  
  model2 <- paste(pars[1], "~", pars[3], sep = "")
  print(paste("MODEL =", model2))
  ANOVA2 <- lm(log(X)~SUBSTRATE, data = data)
  summ2 <- summary(ANOVA2)
  print(summ2)
  tukey2 <-  TukeyHSD(aov(X~SUBSTRATE, data = data))  
  print(tukey2)
  
  model3 <- paste(pars[1], "~", pars[2], "*", pars[3], sep = "")
  print(paste("MODEL =", model3))
  ANOVA3 <- lm(log(X)~LOCOMOTION*SUBSTRATE, data = data)
  summ3 <- summary(ANOVA3)
  print(summ3)
  tukey3 <- TukeyHSD(aov(X~LOCSUB, data = data))  
  print(tukey3)
  
  model4 <- paste(pars[1], "~", pars[4], sep = "")
  print(paste("MODEL =", model4))
  ANOVA4 <- lm(log(X)~GRP, data = data)
  summ4 <- summary(ANOVA4)
  print(summ4)
  tukey4 <- TukeyHSD(aov(X~GRP, data = data))  
  print(tukey4)
  
  model5 <- paste(pars[1], "~", pars[2], "*", pars[5], sep = "")
  print(paste("MODEL =", model5))
  ANOVA5 <- lm(log(X)~LOCOMOTION*log(MIDSHAFTCIRC), data = data)
  summ5 <- summary(ANOVA5)
  print(summ5)
  
  model6 <- paste(pars[1], "~", pars[3], "*", pars[5], sep = "")
  print(paste("MODEL =", model6))
  ANOVA6 <- lm(log(X)~SUBSTRATE*log(MIDSHAFTCIRC), data = data)
  summ6 <- summary(ANOVA6)
  print(summ6)
  
  model7 <- paste(pars[1], "~", pars[4], "*", pars[5], sep = "")
  print(paste("MODEL =", model7))
  ANOVA7 <- lm(log(X)~GRP*log(MIDSHAFTCIRC), data = data)
  summ7 <- summary(ANOVA7)
  print(summ7)
  
  model8 <- paste(pars[1], "~", pars[2], "*",pars[3], "*", pars[5], sep = "")
  print(paste("MODEL =", model8))
  ANOVA8 <- lm(log(X)~LOCOMOTION*SUBSTRATE*log(MIDSHAFTCIRC), data = data)
  summ8 <- summary(ANOVA8)
  print(summ8)
  
  model9 <- paste(pars[1], "~", pars[2], "*",pars[3], "*", pars[4], "*", pars[5], sep = "")
  print(paste("MODEL =", model9))
  ANOVA9 <- lm(log(X)~LOCOMOTION*SUBSTRATE*GRP*log(MIDSHAFTCIRC), data = data)
  summ9 <- summary(ANOVA9)
  print(summ9)
  
  MODELS<- data.frame("model" = c(model1, model2, model3, 
                                  model4, model5, model6, 
                                  model7, model8, model9),
                      "AdjRsq" = c(summ1$adj.r.squared, 
                                   summ2$adj.r.squared, 
                                   summ3$adj.r.squared, 
                                   summ4$adj.r.squared, 
                                   summ5$adj.r.squared, 
                                   summ6$adj.r.squared, 
                                   summ7$adj.r.squared, 
                                   summ8$adj.r.squared, 
                                   summ9$adj.r.squared),
                      "AIC" = c(AIC(ANOVA1), AIC(ANOVA2), AIC(ANOVA3),
                                AIC(ANOVA4), AIC(ANOVA5), AIC(ANOVA6),
                                AIC(ANOVA7), AIC(ANOVA8), AIC(ANOVA9))
  )
  
  anovaMODELS[[i]] <- MODELS
}

sink()

anovaMODELS <- lapply(anovaMODELS, FUN = function(X){
  
  ret <- X[order(X$AIC),]
  rownames(ret) <- c(1:nrow(X))
  return = ret
  
} 
)

save(anovaMODELS,file = "./stats/ECC POS AIC anovaMODELS.Rdata")

sink("./stats/ECC POS AIC anovaMODELS.txt")
print(anovaMODELS)
sink()

