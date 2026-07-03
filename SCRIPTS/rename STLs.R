

setwd("./STLs")



extfiles <- list.files(pattern = "ext")
intfiles <- list.files(pattern = "int")
XMLfiles <- list.files(pattern = "xml")


file.access(extfiles, mode =2)



extfiles1 <- gsub(" ", extfiles, replacement = "_")
intfiles1 <- gsub(" ", intfiles, replacement = "_")
XMLfiles1 <- gsub(" ", XMLfiles, replacement = "_")


extfiles2 <- gsub("_", extfiles1, replacement = "-")
intfiles2 <- gsub("_", intfiles1, replacement = "-")
XMLfiles2 <- gsub("_", XMLfiles1, replacement = "-")


extfiles3 <- gsub("humerus", extfiles2, replacement = "")
intfiles3 <- gsub("humerus", intfiles2, replacement = "")
XMLfiles3 <- gsub("humerus", XMLfiles2, replacement = "")


extfiles4 <- gsub("-ext.stl", extfiles3, replacement = "-humerus_ext.stl")
intfiles4 <- gsub("-int.stl", intfiles3, replacement = "-humerus_int.stl")
XMLfiles4 <- gsub(".xml",     XMLfiles3, replacement = "-humerus.xml")



extfiles4 <- gsub("--", extfiles4, replacement = "-")
intfiles4 <- gsub("--", intfiles4, replacement = "-")
XMLfiles4 <- gsub("--", XMLfiles4, replacement = "-")


extfiles4
intfiles4
XMLfiles4


file.rename(extfiles, extfiles4)
file.rename(intfiles, intfiles4)
file.rename(XMLfiles, XMLfiles4)


nameEXP<-gsub("_int.stl", intfiles3, replacement = "")

write.csv(nameEXP, file = "../DATA/specimenECO")

