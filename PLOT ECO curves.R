
load("allDATA.Rdata")


require(ggplot2)
library(gridExtra)
library(scales)


# create data plots -------------------------------------------------------

midCirc <- unlist(allDATA$MidShaft.Ext_perim)

# scale data
ALLdat <- data.frame(Curvature = allDATA$DistCentroid/(midCirc^1),
                     thickness = allDATA$Mean_thick/(midCirc^1),
                     CA = allDATA$CA/(midCirc^2),
                     Ix = (allDATA$Ix)/(midCirc^4),
                     Iy = (allDATA$Iy)/(midCirc^4),
                     I.ratio = (allDATA$Imax/allDATA$Imin),
                     Zpol = allDATA$Zpol.Zpol/(midCirc^3),
                     J = allDATA$J.Ix/(midCirc^4)
                     
)

ALLdat <- scale(ALLdat, center = F)


spID <-  gsub(pattern =   "-humerus" ,
              replacement = "", 
              x = allDATA$X)

ECO <- data.frame(lengthAlong = allDATA$lengthAlong,
           LOCSUB = allDATA$LOCSUB,
           TAXA = allDATA$TAXA,
           spID)


ALLdatLONG <- NULL 


# convert to long

i=1
for(i in 1:ncol(ALLdat)){
  
  ALLdatLONG <- rbind(ALLdatLONG, data.frame(X = ALLdat[,i], ECO, VAR = rep(colnames(ALLdat)[i], 29)))
  
}

ALLdatLONG


range(ALLdatLONG$X)

c(0.4, 2.5)
# -------------------------------------------------------------------------

level = 0.975

span = 0.75
levels(ALLdatLONG$LOCSUB)
ALLdatLONG[ALLdatLONG$LOCSUB==levels(ALLdatLONG$LOCSUB)[1],]

gg_AQ <- ggplot(ALLdatLONG[ALLdatLONG$LOCSUB==levels(ALLdatLONG$LOCSUB)[1],]
, aes(x = lengthAlong , y = X), group = VAR) +
    # geom_point(aes(col = ID)) +
    geom_smooth(se = T, aes(col = VAR) ,method= "loess", level = level, method.args = list(span = span)) +
    xlab("Length Along Bone") +
  coord_cartesian(ylim= c(0.4, 2.4)) +
    # ylab("Curvature mom.") + 
    # scale_y_continuous(labels = label_number(accuracy = 0.1)) +
    # ggtitle(ecoID) + 
    theme_bw()+
  ggtitle(levels(ALLdatLONG$LOCSUB)[1])
  
gg_AQ <-  gg_AQ +
    theme(axis.text.y = element_text(angle = 90, vjust =0.5, hjust = 0.5)) +
  theme(legend.position = "none") + 
  scale_color_manual(values= c("#E69F00", "#000000", "#56B4E9", "#009E73",
                               "#F0E442", "#0072B2", "#D55E00", "#CC79A7"))
  
gg_AQ

gg_AS <- ggplot(ALLdatLONG[ALLdatLONG$LOCSUB==levels(ALLdatLONG$LOCSUB)[2],]
, aes(x = lengthAlong , y = X), group = VAR) +
    # geom_point(aes(col = ID)) +
    geom_smooth(se = T, aes(col = VAR) ,method= "loess", level = level, method.args = list(span = span)) +
    xlab("Length Along Bone") +
  coord_cartesian(ylim= c(0.4, 2.4)) +
  
    # ylab("Curvature mom.") + 
    # scale_y_continuous(labels = label_number(accuracy = 0.1)) +
    # ggtitle(ecoID) + 
    theme_bw()+
  ggtitle(levels(ALLdatLONG$LOCSUB)[2])
  
gg_AS <-  gg_AS +
    theme(axis.text.y = element_text(angle = 90, vjust =0.5, hjust = 0.5)) +
  theme(legend.position = "none") + 
  scale_color_manual(values= c("#E69F00", "#000000", "#56B4E9", "#009E73",
                               "#F0E442", "#0072B2", "#D55E00", "#CC79A7"))

gg_AS

gg_AV <- ggplot(ALLdatLONG[ALLdatLONG$LOCSUB==levels(ALLdatLONG$LOCSUB)[3],]
, aes(x = lengthAlong , y = X), group = VAR) +
    # geom_point(aes(col = ID)) +
    geom_smooth(se = T, aes(col = VAR) ,method= "loess", level = level, method.args = list(span = span)) +
    xlab("Length Along Bone") +
  coord_cartesian(ylim= c(0.4, 2.4)) +
  
    # ylab("Curvature mom.") + 
    # scale_y_continuous(labels = label_number(accuracy = 0.1)) +
    # ggtitle(ecoID) + 
    theme_bw()+
  ggtitle(levels(ALLdatLONG$LOCSUB)[3])
  
gg_AV <-  gg_AV +
    theme(axis.text.y = element_text(angle = 90, vjust =0.5, hjust = 0.5)) +
  theme(legend.position = "none") + 
  scale_color_manual(values= c("#E69F00", "#000000", "#56B4E9", "#009E73",
                               "#F0E442", "#0072B2", "#D55E00", "#CC79A7"))

gg_AV

gg_TQ <- ggplot(ALLdatLONG[ALLdatLONG$LOCSUB==levels(ALLdatLONG$LOCSUB)[4],]
, aes(x = lengthAlong , y = X), group = VAR) +
    # geom_point(aes(col = ID)) +
    geom_smooth(se = T, aes(col = VAR) ,method= "loess", level = level, method.args = list(span = span)) +
    xlab("Length Along Bone") +
  coord_cartesian(ylim= c(0.4, 2.4)) +
  
    # ylab("Curvature mom.") + 
    # scale_y_continuous(labels = label_number(accuracy = 0.1)) +
    # ggtitle(ecoID) + 
    theme_bw()+
  ggtitle(levels(ALLdatLONG$LOCSUB)[4])
  
gg_TQ <-  gg_TQ +
    theme(axis.text.y = element_text(angle = 90, vjust =0.5, hjust = 0.5)) +
    theme(legend.position = "none") + 
  scale_color_manual(values= c("#E69F00", "#000000", "#56B4E9", "#009E73",
                               "#F0E442", "#0072B2", "#D55E00", "#CC79A7"))

gg_TQ


gg_list <- list(gg_AQ,
                gg_AS,
                gg_AV,
                gg_TQ)




# 


require(gridExtra)
# grid.arrange(arrangeGrob(grobs = taxaPlots[d[2,]], nrow = 9, ncol = 1))



pl <- marrangeGrob(grobs = gg_list, nrow = 4, ncol = 1)
pl



ggsave("../0MANUSCRIPT/0FIGURES/Fig X_long.pdf", pl, width = 5, height = 20, units = "cm")


  