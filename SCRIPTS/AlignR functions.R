readXMLlands <- function(filelist, dir, save = T, file.out = "coords.Rdata"){
  for(i in 1: length( filelist)){
    XMLTMP <- xmlToDataFrame(file.path(dir, filelist[i]), )
    
    coordsTMP <-  as.data.frame(matrix(as.numeric(unlist(strsplit(XMLTMP$Coordinate, " "))), ncol = 3, byrow = T))
    rownames(coordsTMP) <- XMLTMP$Name
    coords[[i]] <- coordsTMP
  }
  
  names(coords) <- gsub(pattern = ".xml", replacement = "", filelist)
  
  save(file = file.out, coords)
  return(coords)  
  
}

# morphomapCheck <- function(mesh, lands = NULL, col = "grey"){
#   clear3d()
#   shade3d(mesh, col = col)
#   
#   if(!is.null(lands)){
#     points3d(lands, col = rainbow(nrow(lands)))
#   }
#   
#   axis3d("x", pos = c(NA, 0, 0), lwd = 5, col = "red")
#   axis3d("y", pos = c(0, NA, 0), lwd = 5, col = "green")
#   axis3d("z", pos = c(0, 0, NA), lwd = 5, col = "blue")
#   title3d(main = NULL, xlab = "X axis", zlab = "Z (Biomechanical length)")
#   rgl.bbox()
#   view3d(userMatrix = rotationMatrix(pi, 1, 0, 1), fov = 0)
#   U <- par3d("userMatrix")
#   par3d(userMatrix = rotate3d(U, pi, 0, 0, 1))
#   p1 <- c(0, 0, 0)
#   p2 <- c(0, 0, 100)
#   p3 <- c(100, 0, 100)
#   normal <- crossProduct(p2 - p1, p3 - p1)
#   planes3d(normal[1], normal[2], normal[3], d = 0, alpha = 0.5, 
#            col = "violet")
#   p1 <- c(0, 0, 0)
#   p2 <- c(100, 0, 0)
#   p3 <- c(0, 100, 0)
#   normal <- crossProduct(p2 - p1, p3 - p1)
#   planes3d(normal[1], normal[2], normal[3], d = 0, alpha = 0.5, 
#            col = "yellow")
#   p1 <- c(0, 0, 0)
#   p2 <- c(0, 0, 100)
#   p3 <- c(0, 100, 100)
#   normal <- crossProduct(p2 - p1, p3 - p1)
#   planes3d(normal[1], normal[2], normal[3], d = 0, alpha = 0.5, 
#            col = "lightblue")
# }


ProcGPA2 <- function (dat.array, tol = 1e-05, scale = TRUE, CSinit = FALSE, 
          silent = TRUE, weights = NULL, centerweight = FALSE, reflection = TRUE, 
          pcAlign = TRUE) 
{
  if (!is.null(weights)) 
    weights <- weights/sum(weights)
  t0 <- Sys.time()
  x <- dat.array
  p1 <- 1e+10
  p2 <- p1
  n <- dim(dat.array)[3]
  k <- dim(dat.array)[1]
  m <- dim(dat.array)[2]
  x1 <- gdif(dat.array)
  arr.list <- list(0)
  for (i in 1:n) arr.list[[i]] <- list(x[, , i], 1)
  if (CSinit) {
    arr.list <- lapply(arr.list, function(x) {
      x[[1]] <- scale(x[[1]], scale = FALSE)
      x[[1]] <- x[[1]]/cSize(x[[1]])
      return(list(x[[1]], x[[2]]))
    })
  }
  else {
    arr.list <- lapply(arr.list, function(x) {
      x[[1]] <- scale(x[[1]], scale = FALSE)
      return(list(x[[1]], x[[2]]))
    })
  }
  mshape <- x[, , 1]
  mshape <- scale(mshape, scale = FALSE)
  if (pcAlign) 
    mshape <- pcAlign(mshape)
  if (centerweight && !is.null(weights)) {
    mcent <- apply(mshape, 2, weighted.mean, w = weights)
    mshape <- scale(mshape, scale = F, center = mcent)
  }
  while (p1 > tol) {
    arr.list <- lapply(arr.list, function(x) {
      x[[1]] <- rot.proc(x[[1]], x = mshape, scale = F, 
                         weights = weights, centerweight = centerweight, 
                         reflection = reflection)
      return(list(x[[1]], x[[2]]))
    })
    for (i in 1:n) x[, , i] <- arr.list[[i]][[1]]
    x2 <- gdif(x)
    p1 <- abs(x1 - x2)
    x1 <- x2
    if (scale) {
      for (i in 1:n) arr.list[[i]] <- list(x[, , i], 1)
      while (p2 > tol) {
        for (i in 1:n) if (!is.null(weights)) 
          x[, , i] <- arr.list[[i]][[1]] * weights
        else x[, , i] <- arr.list[[i]][[1]]
        eigc <- scaleproc(x)
        for (i in 1:n) arr.list[[i]][[2]] <- eigc[i]
        arr.list <- lapply(arr.list, function(x) {
          x[[1]] <- x[[1]] * x[[2]]
          return(list(x[[1]], x[[2]]))
        })
        arr.list <- lapply(arr.list, function(x) {
          x[[1]] <- rot.proc(x[[1]], x = mshape, scale = F, 
                             weights = weights, centerweight = centerweight, 
                             reflection = reflection)
          return(list(x[[1]], x[[2]]))
        })
        for (i in 1:n) x[, , i] <- arr.list[[i]][[1]]
        x2 <- gdif(x)
        p2 <- abs(x1 - x2)
        x1 <- x2
      }
    }
    mshape <- arrMean3(x)
    if (CSinit) {
      msize <- cSize(mshape)
      mshape <- mshape/msize
      if (scale) 
        x <- x/msize
    }
  }
  t1 <- Sys.time()
  if (!silent) 
    message(paste("in... ", format(t1 - t0)[[1]], "\n"))
  return(list(rotated = x, mshape = mshape))
}


AlignByMesh <- function(X1, X2 = NULL, lands = NULL, alignIND = NULL ,center= T, check.align = F){
  
  if(!is.null(X2)){
    X1dim <- data.frame(lapply(X1, FUN= dim))
    X2dim <- data.frame(lapply(X2, FUN= dim))
    X <- mergeMeshes(X1,X2)
    # data.frame(lapply(Xcomb, FUN= dim))
    
  }else{
    X <- X1
  }
  

    
    Xprime <-pcAlign(X)
    M = rbind(c(0, 0, 1),
              c(0, 1, 0),
              c(1, 0, 0))
    Z <- applyTransform(Xprime, trafo = M)
    
    Zdim <- data.frame(lapply(Z, FUN= dim))
    

  
  
  if(center){
    Z$vb[3,] <- Z$vb[3,] - min(Z$vb[3,])
    
    
  }

  
  if(check.align){
    repeat {
      morphomapCheck(Z)
      if(askYesNo("Flip Z axis? (blue)", default = F)){
        M = rbind(c(1, 0, 0),
                  c(0, 1, 0),
                  c(0, 0, -1))
        Z <- applyTransform(Z, trafo = M)
        Z$vb[3,] <- Z$vb[3,] - min(Z$vb[3,])
      }
      
      morphomapCheck(Z)
      
      if(askYesNo("Flip X axis? (red)", default = F)){
        M = rbind(c(-1, 0, 0),
                  c(0, 1, 0),
                  c(0, 0, 1))
        Z <- applyTransform(Z, trafo = M)
        Z$vb[3,] <- Z$vb[3,] - min(Z$vb[3,])
      }
      
      morphomapCheck(Z)
      
      if(askYesNo("Flip Y axis? (green)", default = F)){
        M = rbind(c(1, 0, 0),
                  c(0, -1, 0),
                  c(0, 0, 1))
        Z <- applyTransform(Z, trafo = M)
        Z$vb[3,] <- Z$vb[3,] - min(Z$vb[3,])
      }
      morphomapCheck(Z)
      
      if(askYesNo("Correct?", default = F)){
        break
      }
      
    }
    
  }
    
    
    if(!is.null(lands)){
      
      trafo <- computeTransform(Z, X, type = "rigid", reflection = T)
      Zlands <- applyTransform(lands, (trafo))

    }
    
    # clear3d()
    # 
    # 
    # plot3d(Z, col = "red", alpha = 0.2, asp = "iso")
    # points3d(Zlands, size = 5, col = "red")
    # 
    # shade3d(X, col = "grey", alpha = "0.2")
    # points3d(lands, size = 5, col = "black")
    # 
  
  if(!is.null(X2)){
    Z1 <- X1
    Z2 <- X2
    
    Z1$vb <- Z$vb[,1:X1dim[2,1]]
    Z1$it <- Z$it[,1:X1dim[2,2]]
    Z1$normals <- Z$normals[,1:X1dim[2,3]]
    
    Z2$vb <- Z$vb[,(X1dim[2,1]+1):Zdim[2,1]]
    Z2$it <- (Z$it[,(X1dim[2,2]+1):Zdim[2,2]])
    Z2$it <- Z2$it-(min(Z2$it)-1)
    Z2$normals <- Z$normals[,(X1dim[2,3]+1):Zdim[2,3]]
    
    # X1$it[,1]
    # Z$it[,1]
    # 
    # X2$it[,1]
    # Z2$it[,1]
    # 
    
    if(!is.null(lands)){
      
      return(list(mesh = list(ext = Z1,int = Z2), lands = Zlands))
      
    } else  return(list(ext = Z1,int = Z2))
    
  }else{
     if(!is.null(lands)){
          
    return(list(mesh = Z, lands = Zlands))
        
    } else  return(Z)
  }
}


AlignbyLands <- function(MESHES, lands = NULL, alignIND = NULL ,center= T, Zalign=F, zero = T, scale = F,  check.align = F){
  
  if(class(MESHES) == "LANDMESH"){
    
    lands <- MESHES$LANDS
    
    MESHES <- lapply(MESHES$MESHES, function(X){
      X$mesh
    })
    
  }
  
  if(!is.array(lands)) error("lands is not an array")
  if( length(MESHES) != dim(lands)[3]) error("number of meshes does not equal number of landmark sets")
  if(!is.null(alignIND)){
    
    alignLands <- lands[alignIND, ,]
    
  } else{
    alignLands <- lands
  }
  
  mshape <- lands[,,1]
  
  
  REF <- pcAlign(mshape)
  
  if(Zalign){
    flipXZ = rbind(c(0, 0, 1),
                   c(0, 1, 0),
                   c(1, 0, 0))
    REF <- applyTransform(REF, trafo = flipXZ)
    
    # REF[,3] <- REF[,3] - min(REF[,3])
    
  }
  
  
  rotEXT <- rotmesh.onto(MESHES[[1]][[1]], mshape, REF)
  rotINT <- rotmesh.onto(MESHES[[1]][[2]], mshape, REF)
  if(Zalign){
    
    if(zero){
      floor <- min(rotEXT$mesh$vb[3,])
      rotEXT$mesh$vb[3,] <- rotEXT$mesh$vb[3,] - floor  
      rotEXT$yrot[,3] <- rotEXT$yrot[,3] - floor  
      
      rotINT$mesh$vb[3,] <- rotINT$mesh$vb[3,] - floor  
      rotINT$yrot[,3] <- rotINT$yrot[,3] - floor  
      
    }
  }
  
  
  
  refMESH <- list(list(mesh = list(EXTmesh = rotEXT$mesh, INTmesh = rotINT$mesh),
                       yrot = rotEXT$yrot, trafo = matrix(0, ncol = 4, nrow = 4)  ))
  
  names(refMESH) = names(MESHES)[1]
  
  # if(plotref){
  #   clear3d()
  #   shade3d(refMESH[[1]]$mesh$EXTmesh, col = "grey", alpha= 0.2)
  #   spheres3d(refMESH[[1]]$yrot, col= "red", radius = 0.7)
  #   axis3d("x", pos = c(NA, 0, 0), lwd = 5, col = "red")
  #   axis3d("y", pos = c(0, NA, 0), lwd = 5, col = "green")
  #   axis3d("z", pos = c(0, 0, NA), lwd = 5, col = "blue")
  # 
  # }
  
  
  ind = alignIND
  alignList <- array2list(alignLands)
  landsList <- array2list(lands)
  
  IND <- 1:length(landsList)
  
  # 
  # #TEST
  # x <- IND[1]
  # names(IND) = names(MESHES)
  # #ENDTEST
  # 
  
  rotMESHES <- list()
  rotLANDS <- list()
  CSIZE <- c()
  # rotALL[[1]] <- refMESH
  
  for( i in 1:length(IND)){
    cat("registering", names(MESHES)[i], "\n")
    rotlands <- rotonmat(landsList[[i]], tarmat = REF[alignIND,], refmat = landsList[[i]][alignIND,], 
                         scale = scale, reflection = F, getTrafo = T)
    
    rotLANDS[[i]] <- rotlands$Xrot
    
    rotmeshEXT <- rotmesh.onto(MESHES[[i]][[1]], refmat = landsList[[i]], scale = scale, tarmat = rotlands$Xrot)
    rotmeshINT <- rotmesh.onto(MESHES[[i]][[2]], refmat = landsList[[i]], scale = scale, tarmat = rotlands$Xrot)
    
    rotmesh <- list(mesh = list(EXTmesh = rotmeshEXT$mesh, INTmesh = rotmeshINT$mesh), 
                    yrot = rotmeshEXT$yrot, trafo = rotmeshEXT$trafo)
    
    rotMESHES[[i]] <- rotmesh
    
    CSIZE[i] <- cSize(rotmeshEXT$mesh$vb)
  }
  
  names(rotMESHES) = names(MESHES)
  names(CSIZE) = names(MESHES)
  
  rotLANDS <- list2array(rotLANDS)
  dimnames(rotLANDS)[[3]] <- names(rotMESHES) 
  
  return(list(LANDS = rotLANDS, CSIZE = CSIZE, MESHES = rotMESHES, REF = refMESH, alignIND = alignIND))  
}



checkAlignment <- function(alignedLandMeshes, refMesh = 1, reAlign = T, overwrite = F, exportMESHES = T, filepath = "./NewAligned"){
  
  n <- length(alignedLandMeshes)
  
  if(!dir.exists(filepath)){
    dir.create(filepath)
    dir.create(file.path(filepath, "meshes"))
    
  }
  reflectFP <- file.path(filepath,"reflect.csv")
  #set transformation matrix

  
  if(file.exists(reflectFP)){
    if(overwrite){
      cat("existing file exists - overwriting\n")
      unlink(reflectFP)
      
      reflect <- rep(0, length(alignedLandMeshes$MESHES) )
      names(reflect) <- names( alignedLandMeshes$MESHES)
      reflect <- as.data.frame(reflect)
      reflect[refMesh,] <- 1
      
      write.csv(reflect, reflectFP)
    } else{
      cat("existing file exists - reading existing alignment file\n")
      reflect <- read.csv(reflectFP, header = T, row.names = 1)
    }
  } else{
    reflect <- rep(0, length(alignedLandMeshes$MESHES) )
    names(reflect) <- names( alignedLandMeshes$MESHES)
    reflect <- as.data.frame(reflect)
    reflect[refMesh,] <- 1
    
  }
  
  seqIND <- which(reflect == 0)
  newMESHES <- alignedLandMeshes$MESHES
  newLANDS <- alignedLandMeshes$LANDS
  REF <- alignedLandMeshes$MESHES[[refMesh]]
  
  open3d()
  
  i=2
  for(i in seqIND){
    
    if(reflect[i,1]!=0){
      cat("Skipping",names(alignedLandMeshes$MESHES)[i], "\n")
      next()
    }
    cat(i, "of", length(seqIND),": " ,names(alignedLandMeshes$MESHES)[i], "\n")
    clear3d()
    shade3d(REF$mesh$EXTmesh, col = "grey", alpha= 0.2)
    spheres3d(REF$yrot, col = "blue", alpha= 0.2, radius = 0.4)
    
    shade3d(alignedLandMeshes$MESHES[[i]]$mesh$EXTmesh, col = "red", alpha= 0.2)
    spheres3d(alignedLandMeshes$MESHES[[i]]$yrot, col = "red", alpha = 0.2,  radius = 0.7)
    
    # axis3d("x", pos = c(NA, 0, 0), lwd = 2, col = "red")
    # axis3d("y", pos = c(0, NA, 0), lwd = 2, col = "green")
    # axis3d("z", pos = c(0, 0, NA), lwd = 2, col = "blue")
    
    flip <- NA
    flip <- readline(prompt = "Flip Z axis? (Y/N y/n T/F): ")
    # flip = "Y"
    
    while(!flip %in% c("Y", "y", "T", "N", "n", "F")){
      flip <- readline(prompt = "Input Error - Try again. Flip Z axis? (Y/N y/n T/F): ")
    }
    if(flip %in% c("Y", "y", "T")){
      reflect[i,1] <- -1
      write.csv(reflect, reflectFP)
      
    } 
    if(flip %in% c("N", "n", "F")){
      reflect[i,1] <- 1
      write.csv(reflect, reflectFP)
      
    } 
  } # end user input loop
  
  
  n <- length(alignedLandMeshes$MESHES)
  
  M = rbind(c(-1, 0, 0),
            c(0, 1, 0),
            c(0, 0, 1))
  
  # apply transforms to meshes and landmarks. Optional mesh export
  for(i in 1:n){
    if(reflect[i,1] == -1){
      newmesh <- list()
      newmesh$EXTmesh <- applyTransform(alignedLandMeshes$MESHES[[i]]$mesh$EXTmesh, trafo = M)
      newmesh$INTmesh <- applyTransform(alignedLandMeshes$MESHES[[i]]$mesh$INTmesh, trafo = M)
      
      newlands <- applyTransform(alignedLandMeshes$MESHES[[i]]$yrot, trafo = M)
      
      
      clear3d()
      
      shade3d(newmesh$EXTmesh, col = "blue", alpha = 0.2)
      points3d(newlands)
      
      #perform new registration
      cat("RE-registering", names(alignedLandMeshes$MESHES)[i], "\n")
      
      # LANDMARKS
      rotlands <- rotonmat(newlands, tarmat = REF$yrot[,], refmat = newlands[,],
                           scale = F, reflection = F, getTrafo = T)
      
      rotmesh <- newmesh
      
      rotmesh$EXTmesh <- applyTransform(newmesh$EXTmesh, trafo = rotlands$trafo)
      rotmesh$INTmesh <- applyTransform(newmesh$INTmesh, trafo = rotlands$trafo)
      
      
      # 
      
      clear3d()
      
      shade3d(REF$mesh$EXTmesh, col = "grey", alpha = 0.2)
      points3d(REF$yrot, col = "black")
      shade3d(rotmesh$EXTmesh, col = "blue", alpha = 0.2)
      points3d(rotlands$Xrot, col = "blue")
      
      
      # names(rotLANDS) <- names(alignedLandMeshes$MESHES)[i]
      
      # MESHES
      
      
      
      newLANDS[,,i] <- rotlands$Xrot
      newMESHES[[i]] <- list(mesh = rotmesh, 
                             yrot = rotlands$Xrot, trafo = rotlands$trafo)
    

      
    }
    if(reflect[i,1] ==  1){
      #perform new registration
      # cat("SKIP", names(alignedLandMeshes$MESHES)[i], "\n")
      
    }
    
    if(exportMESHES){
      cat("writing to", file.path(filepath,"MESHES", paste(names(alignedLandMeshes$MESHES)[i], "_EXT.stl", sep ="")), "\n")
      vcgStlWrite(newMESHES[[i]]$mesh$EXTmesh, file.path(filepath,"MESHES", paste(names(alignedLandMeshes$MESHES)[i], "_EXT.stl", sep ="")))
      cat("writing to", file.path(filepath,"MESHES", paste(names(alignedLandMeshes$MESHES)[i], "_INT.stl", sep ="")), "\n")
      vcgStlWrite(newMESHES[[i]]$mesh$INTmesh, file.path(filepath,"MESHES", paste(names(alignedLandMeshes$MESHES)[i], "_INT.stl", sep ="")))
      
    }
    
  }
  
  geomorph::writeland.tps(newLANDS, file = file.path(filepath, "landmarks"))
  
  return(list(LANDMESH = list(LANDS = newLANDS, 
                              MESHES = newMESHES, 
                              REF = REF))
  )
  
}

rotLANDMESH <-function(X, trafo){
  
  MESH <- X$mesh
  YROT <- X$yrot
  ROT <- X$trafo
  
  ZMESH <- MESH
  
  Zlands <- applyTransform(lands, (trafo))
  ZMESH$EXTmesh <- applyTransform(MESH$EXTmesh, (trafo))
  ZMESH$INTmesh <- applyTransform(MESH$INTmesh, (trafo))
  
  out <- list(mesh = ZMESH, yrot = Zlands, trafo = ROT)
  class(out) <- "LANDMESH"
  return(out)
}


plotLandMesh <- function(LANDMESH, plotINTERNAL = F, plotPOINTS = T, plotTEXT = T ,new = F){
  
  # if(class(LANDMESH) != "LANDMESH") stop("ERROR: LANDMESH is not of class 'LANDMESH'")
  
  if(new) open3d()
  
  shade3d(LANDMESH$mesh[[1]], col = "grey", alpha = 0.2)
  
  if(plotPOINTS) points3d(LANDMESH$yrot, size = 10, col = rainbow(nrow(LANDMESH$yrot)))
  if(plotTEXT) text3d(LANDMESH$yrot, texts = 1:nrow(LANDMESH$yrot),size = 10, col = "red")
  
  if(plotINTERNAL)  shade3d(LANDMESH$mesh[[2]], col = "blue", alpha = 1)
  
}

rotLANDMESH <-function(X, trafo){
  
  MESH <- X$mesh
  YROT <- X$yrot
  ROT <- X$trafo
  
  ZMESH <- MESH
  
  Zlands <- applyTransform(YROT, (trafo))
  ZMESH$EXTmesh <- applyTransform(MESH$EXTmesh, (trafo))
  ZMESH$INTmesh <- applyTransform(MESH$INTmesh, (trafo))

  out <- list(mesh = ZMESH, yrot = Zlands, trafo = ROT)
  class(out) <- "LANDMESH"
  return(out)
}


dist3d <- function(a,b,c) {
  
  v1 <- b - c
  v2 <- a - b      
  v3 <- cross3d_prod(v1,v2)
  area <- sqrt(sum(v3*v3))/2
  d <- 2*area/sqrt(sum(v1*v1))
  return(d)
}


cross3d_prod <- function(v1,v2){
  v3 <- vector()
  v3[1] <- v1[2]*v2[3]-v1[3]*v2[2]
  v3[2] <- v1[3]*v2[1]-v1[1]*v2[3]
  v3[3] <- v1[1]*v2[2]-v1[2]*v2[1]
  return(v3)
}


findMaxC <- function(X, Chordline = NULL,verbose = F, segmentCentroid = F){
  X <- as.matrix(X)
  
  lengthAlong <- seq(from = start, to = end, length.out = nrow(X))
  
  if(segmentCentroid){
    dists <- apply(X, MARGIN = 1, FUN = dist3d, 
                   b = X[1,], c = X[nrow(X),])
    
  } else{
    dists <- apply(X, MARGIN = 1, FUN = function(x){
      
      dists <- dist(rbind(x[c(1,2)], c(0,0)))
      
      })
    
  }
  
  
  
  Xdists <- cbind(X,lengthAlong,dists)
  
  maxI <- which.max(dists)
  maxAlong <- lengthAlong[maxI]
  
  Y <- Xdists[maxI,]
  
  d <- max(dists) 
  Y <- c(Y,d)
  
  
  if(verbose){
    return(list(maxC = Y, distAlong = Xdists))
  }else{
    return(Y)
  }

}

# converts an nx3 matrix of n landmarks to a landmark slice array with dimintions nlandx3xnslice
# X an nx3 matrix
# nland number of landmarks within eac slice
# nslice number of slices
# inner logical. does the matrix contain an inner layer of landmarks? If yes, exports a list with outer and inner landmarks

mat2sliceArr <- function(X, nland, nslice, inner = T){
  
  n = nland*nslice
  SLICEind <- matrix(1:n, ncol = nland, nrow = nslice, byrow= T)
  
  
  SLICEarr <- array(NA, dim = c(nland, p, nslice))
  
  for(i in 1:nslice){
    
    SLICEarr[,,i] <- X[SLICEind[i,] ,]
    
  }
  
  if(inner){
    
    INNERind <- SLICEind + n
    INNERarr <- SLICEarr
    
    for(i in 1:nslice){
      INNERarr[,,i] <- X[INNERind[i,] ,]
    }
    
    return(list(OUTER = SLICEarr, INNER = INNERarr))
  }  else{ return(SLICEarr)}
  
}


# calculates each slice as a mesh using INNER and OUTER landmark arrays
# OUTER, INNER landmark slice arrays in nlandx3xnslice format


calcSliceMeshes <-function(OUTER, INNER){
  
  nland <- dim(OUTER)[[1]]
  nslice <- dim(OUTER)[[3]]
  
  indMAT <- 1:nland
  
  it <- NULL
  vb <- NULL
  # nland <- (dim(set)[1]/n)
  maxit <- 0
  
  # a=1
  # i=1
  
  for (i in 1:(nslice)) {
    
    pro_z <- morphomapTri2sects(OUTER[, ,i] , 
                                INNER[, ,i])
    it <- rbind(it, pro_z$tri + maxit)
    vb <- rbind(vb, pro_z$matrix)
    maxit <- maxit + nland * 2
  }
  mesh_d <- list(vb = t(cbind(vb, 1)), it = t(it))
  class(mesh_d) <- "mesh3d"
  
  mesh <- mesh_d
  
  
  shade3d(mesh, alpha = 0.5, add= F, col = palette.colors()[a+1])
  if(a != 4) next3d()
  
  return(mesh)
  
  
  
}

# Calculates a mesh from a set of slice landmarks in array format.
# POINTS - an array with dimensions nlandsx3xnslices 
# plot, logical, plot triangles during meshing
# ... arguments to pass to as.mesh3d

calcMesh <- function(POINTS, plot = F, ...){
  nland <- dim(POINTS)[[1]]
  nslice <- dim(POINTS)[[3]]
  
  m <- apply(POINTS, 2,c)
  
  slicePairs <- cbind(1:(nslice-1),
                      2:(nslice))
  
  LMIND <- matrix(1:n, ncol = nland, nrow = nslice, byrow= T)
  LMIND <- cbind(LMIND,LMIND[,1])
  
  if(plot){
    
    open3d()
    points3d(m)
    
  }
  
  Tri <- lapply(1:nrow(slicePairs), #MARGIN = 1, 
                FUN = function(X, LMIND, m, slicePairs){
                  
                  
                  sIND <- LMIND[slicePairs[X,],]
                  S <- lapply(X = (1:(nland)), FUN = function(X, sIND, m){
                    
                    ss <- m[c(c(sIND[1,X],
                                sIND[1,X+1],
                                sIND[2,X+1]),
                              c(sIND[1,X],
                                sIND[2,X],
                                sIND[2,X+1])),]
                    
                    if(plot){
                      as.mesh3d(ss,type = "triangles", merge= F)
                      wire3d(as.mesh3d(ss,type = "triangles"), col = c("red"))
                    }
                    
                    return(ss)
                  }, sIND, m)
                  
                  S <- do.call("rbind", S)
                  
                  return(S)
                  
                }
                , LMIND, m, slicePairs)
  
  Tri <- do.call("rbind", Tri)
  
  return(as.mesh3d(Tri,type = "triangles", ...))
  
  
}




# findXYZ  <- function(chordLine,){
#   
#   p <- prcomp(chordLine, center = chordLine[1,])
#   
#   p$x
#   
#   
# }

