##################################
#   Here is an example of DDGM
##################################
library(fdadensity)
library(fda)
library(pracma)
library(fdapace)
library(corpcor)
library(glasso)
library(MASS)
library(PenPC)
library(pcalg)
library(lubridate)
library(igraph)
library(zoo)
library(RColorBrewer)
library(ggplot2)

source('Functions_DDGM.R')

load("example_data.RData")
#The example_data is a dataset with six variables, each containing 413 samples, and each sample includes several observations.

func_data <- lapply(example_data,function(tt) do.call(rbind,lapply(tt,function(xx) lqdtrans(xx)$lqd)))
#Utilize the lqdtrans() function to estimate the distribution function from the observations and perform feature embedding.


cleanfunc <- function(mat){
  mat[is.na(mat)] <- 0
  mat[mat==-Inf] <- 0
  mat[mat==Inf] <- 0
  mat
}
func_data <- lapply(func_data,cleanfunc)#shape information
#Cleanse outliers

m_data <- lapply(example_data,function(tt) do.call(rbind,lapply(tt,function(xx) median(xx))))
m_data <- do.call(cbind,m_data)#position information

BestK <- select_nbasis_cv(funclist = func_data,x,4:9,Kfold = 5)
# the result of BestK is 7
# BestK <- 7
x <- seq(0,1,length.out=50)
#Calculate the correlation coefficient matrix.
gram_matrix <- gram_matrix_bspline(nt=500,nbasis = BestK)
f_covmatrix <- func_covmatrix(func_data,x,nbasis = BestK,gram_matrix = gram_matrix)
#estimate the results
adjf <- network_estimate(f_covmatrix,nrow(func_data[[1]]),alpha = 0.01)
adjm <- network_estimate(cov(m_data),nrow(m_data),alpha = 0.01)
#Merge results
adj <- adjf+adjm
adj <- adj>0
adj[adj] <- 1
#plot
g <- graph.adjacency(adj,mode="undirected",weighted=T)
plot(g)

