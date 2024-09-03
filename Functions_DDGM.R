gram_matrix_bspline <- function(nt=100,rangeval=c(0,1),nbasis=5){
  #return a nbasis*nbasis matrix
  phi.basis=create.bspline.basis(rangeval=rangeval, nbasis=nbasis)
  t = seq(rangeval[1], rangeval[2], length.out = nt)
  phi = predict(phi.basis, t)
  gram_matrix <- matrix(0,nbasis,nbasis)
  for(i in 1:nbasis){
    for(j in 1:i){
      gram_matrix[i,j]=gram_matrix[j,i]=sum(phi[,i]*phi[,j])/nt
    }
  }
  gram_matrix
}
fit_func <- function(ytrain,xtrain,xpre,basis){
  f <- Data2fd(xtrain,ytrain,basis)
  ypre <- predict(f, newdata=xpre, Lfdobj=0, returnMatrix=FALSE)
  ypre
}
select_nbasis_cv <- function(funclist, x, nbasislist, Kfold = 5) {
  #nbasis>3
  index <- round(seq(1, length(x), length.out = Kfold + 1))
  err <- c()
  for (b in 1:length(nbasislist)) {
    
    phi.basis=create.bspline.basis(rangeval=range(x), nbasis=nbasislist[b])
    for (k in 1:Kfold) {
      testindex <- index[k]:index[k + 1]
      sub_err <- 0
      for (f in 1:length(funclist)) {
        func <- funclist[[f]]
        trainx <- x[-testindex]
        trainy <- func[, -testindex]
        testx <- x[testindex]
        testy <- func[, testindex]
        prey <- t(do.call(cbind,lapply(1:nrow(func), function(xx) fit_func(trainy[xx, ], trainx, testx,phi.basis))))
        sub_err <- sub_err+sum((prey-testy)^2)
      }
    }
    err <- c(err,sub_err)
  }
  best_nbasis <- nbasislist[which(err==min(err))]
  best_nbasis
}
func_covmatrix <- function(funclist,x,nbasis=5,gram_matrix = NULL){
  
  if(is.null(gram_matrix)){
    gram_matrix <- gram_matrix_bspline(nt=500,rangeval=range(x),nbasis=nbasis)
  }
  phi.basis=create.bspline.basis(rangeval=range(x), nbasis=nbasis)
  bspline_fit <- function(func,x,basis){
    res <- lapply(1:nrow(func),function(xx) Data2fd(x,func[xx,],basis)$coef)
    return(t(do.call(cbind,res)))
  }
  veclist <- lapply(1:length(funclist),function(xx) bspline_fit(funclist[[xx]],x,phi.basis))
  
  gram_inner <- function(vec1,vec2,gram_matrix){
    sum(unlist(lapply(1:nrow(vec1),function(xx) (vec1[xx,]%*%gram_matrix)%*%vec2[xx,])))/(nrow(vec1))
  }
  covmatrix <- matrix(0,length(veclist),length(veclist))
  for(i in 1:length(veclist)){
    for(j in 1:i){
      covmatrix[i,j]=covmatrix[j,i]=gram_inner(veclist[[i]],veclist[[j]],gram_matrix)
    }
  }
  
  covmatrix
}

lqdtrans <- function(points,xrange=NULL,breaknum=50){
  if(length(unique(points))==1){
    points[1] <- points[1]+1
  }
  prob <- seq(0,1,length.out = breaknum)
  points <- as.numeric(na.omit(points))
  #points <- points/mean(points,na.rm =T)
  xrange <- range(points)
  h1=1/(breaknum-1)
  newd=density(points,from = xrange[1],to = xrange[2],n=breaknum,na.rm =T)
  dmatrix <- t(newd$y)
  dSup <- newd$x
  density <- normaliseDensities(dmatrix, dSup)
  quan <- dens2quantile(dens=density, dSup = dSup,qSup = prob)
  qd=gradient(quan,h1=h1)
  lqd <- log(qd)
  list(lqd=lqd,prob=prob,density=density,quan=quan,qd=qd)
}
FuncPCstable <- function (suffStat, indepTest, p, alpha, verbose = FALSE, edgeWeights, 
                          fixedEdges = NULL, NAdelete = FALSE) 
{
  stopifnot((p = as.integer(p)) >= 2)
  cl = match.call()
  fixedGaps = (edgeWeights == 0)
  if (is.null(edgeWeights)) {
    stop("no edgeWeights")
  }
  else {
    if (!identical(dim(fixedGaps), c(p, p))) {
      stop("Dimensions of the dataset and fixedGaps do not agree.")
    }
    else {
      if (!all(fixedGaps == t(fixedGaps))) 
        stop("fixedGaps must be symmetric")
      G = !fixedGaps
    }
  }
  m.max = max(rowSums(G))
  if (is.null(fixedEdges)) {
    fixedEdges = matrix(FALSE, p, p)
  }
  else {
    if (!(identical(dim(fixedEdges), c(p, p)))) 
      stop("Dimensions of the dataset and fixedEdges do not agree.")
    if (fixedEdges != t(fixedEdges)) 
      stop("fixedEdges must be symmetric")
  }
  seq_p = seq_len(p)
  sepset = pl = vector("list", p)
  for (i in seq_p) sepset[[i]] = pl
  if (sum(G) > 0) {
    w.upper = which(upper.tri(diag(p)), arr.ind = T)
    ind = w.upper[G[w.upper], ]
    for (i in 1:dim(ind)[1]) {
      x <- ind[i, 1]
      y <- ind[i, 2]
      pval <- indepTest(x, y, integer(0), suffStat)
      if (verbose) 
        cat("Marginal test", x, "and", y, "pvalue = ", 
            pval, "\n")
      if (pval >= alpha) {
        G[x, y] = G[y, x] = FALSE
        sepset[[x]][[y]] = sepset[[y]][[x]] = integer(0)
      }
      else {
        sepset[[x]][[y]] = setdiff(seq_p, c(x, y))
        isect = (G[x, ] & G[y, ])
        if (sum(isect) > 0) {
          allZ = seq_p[isect]
          for (z in allZ) {
            pval = indepTest(x, y, z, suffStat)
            if (verbose) 
              cat("text", x, "and", y, "with commone vertices", 
                  allZ, "pvalue = ", pval, "\n")
            if (is.na(pval)) 
              pval <- if (NAdelete) {
                1
              }
            else {
              0
            }
            if (pval < alpha) {
              sepset[[x]][[y]] = sepset[[y]][[x]] = setdiff(sepset[[x]][[y]], 
                                                            z)
            }
          }
        }
      }
    }
  }
  pMax = matrix(-Inf, p, p)
  diag(pMax) = 1
  done = FALSE
  ord = 0
  n.edgetests = numeric(1)
  n.edgerms = numeric(1)
  while (any(G) && ord <= m.max) {
    #cat(date(), ": order=", ord, ", # of edges remaining = ", 
    #sum(G)/2, "\n")
    ord1 = ord + 2L
    n.edgetests[ord1] = 0
    n.edgerms[ord1] = 0
    tmpG = G
    w.lower = which(lower.tri(tmpG), arr.ind = T)
    tmpG[w.lower] = 0
    ind = matrix(which(tmpG != 0, arr.ind = T), ncol = 2)
    remainingEdgeTests = nrow(ind)
    aG = G
    for (i in 1:remainingEdgeTests) {
      if (verbose) {
        if (i%%100 == 0) 
          cat("|i=", i, "|iMax=", nrow(ind), "\n")
      }
      x = ind[i, 1]
      y = ind[i, 2]
      if (G[y, x] && !fixedEdges[y, x]) {
        nbrsUnion = (aG[, x] | aG[, y])
        nbrsUnion[x] = FALSE
        nbrsUnion[y] = FALSE
        nbrsIsect = (aG[, x] & aG[, y])
        nbrsIsect[x] = FALSE
        nbrsIsect[y] = FALSE
        nbrsdec = nbrsIsect
        if (sum(nbrsIsect) != 0) {
          Gsub = aG
          Gsub[, c(x, y)] = Gsub[c(x, y), ] = FALSE
          nbrsdec = (connectedComp(Gsub, seq_p[nbrsIsect]) & 
                       nbrsUnion)
        }
        nbrs = seq_p[nbrsUnion]
        nbrs2rm = seq_p[nbrsdec]
        length_nbrs = length(nbrs2rm)
        if (length_nbrs >= ord) {
          if (length_nbrs > ord) 
            done = FALSE
          if (ord == 0) {
            nbrs2test = nbrs
          }
          else {
            S = seq_len(ord)
            nbrs2test = setdiff(nbrs, nbrs2rm[S])
          }
          repeat {
            n.edgetests[ord1] = n.edgetests[ord1] + 
              1
            pval = indepTest(x, y, nbrs2test, suffStat)
            if (verbose) {
              cat("x=", x, " y=", y, " nbrs=", nbrs, 
                  " nbrs2rm=", nbrs2rm, " nbrs2test=", 
                  nbrs2test, ": pval =", pval, "\n")
            }
            if (is.na(pval)) {
              print("pval is NA")
              pval = if (NAdelete) 
                1
              else 0
            }
            if (pval > pMax[x, y]) 
              pMax[x, y] = pval
            if (pval >= alpha) {
              n.edgerms[ord1] = n.edgerms[ord1] + 1
              G[x, y] = G[y, x] = FALSE
              sepset[[x]][[y]] = nbrs2test
              break
            }
            else if (ord <= 0) {
              break
            }
            else {
              nextSet = getNextSet(length_nbrs, ord, 
                                   S)
              if (nextSet$wasLast) {
                break
              }
              S = nextSet$nextSet
            }
          }
        }
      }
    }
    ord = ord + 1
  }
  for (i in 1:(p - 1)) {
    for (j in 2:p) {
      pMax[i, j] = pMax[j, i] = max(pMax[i, j], pMax[j, 
                                                     i])
    }
  }
  nnms = as.character(seq_p)
  Gobject = if (sum(G) == 0) {
    new("graphNEL", nodes = nnms)
  }
  else {
    colnames(G) = rownames(G) = nnms
    as(G, "graphNEL")
  }
  new("pcAlgo", graph = Gobject, call = cl, n = integer(0), 
      max.ord = as.integer(ord - 1), n.edgetests = n.edgetests, 
      sepset = sepset, pMax = pMax, zMin = matrix(NA, 1, 1))
}
network_estimate <- function(covmatrix,n,alpha=0.05,ifr=F){
  p <- nrow(covmatrix)
  #d <- diag(sqrt(diag(covmatrix)))
  #r <- solve(d)%*%covmatrix%*%solve(d)
  if(ifr){
    r <- covmatrix
  }else{
    r <- cov2cor(covmatrix)
  }
  
  edgeWeights <- matrix(1,p,p)
  diag(edgeWeights) <- 0
  skeletom <- FuncPCstable(suffStat = list(C = r, n = n),
                           indepTest = gaussCItest, ## (partial correlations)
                           alpha = alpha,p=p,edgeWeights=edgeWeights)
  adj <- as(skeletom,"matrix")
  #sparsity <- sum(adj)/(p*p-p)
  return(adj)
}