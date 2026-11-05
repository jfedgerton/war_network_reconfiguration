## This script conducts SAOM models for "The Coevolution of Networks of Interstate Support, Interstate Threat and Civil War" 
## The starting point for the script was derived from Kinne & Bunte, "Guns or Money," BJPS 2018

## Since the estimation involves simulations, results may vary slightly from the published results.

library(reshape2)
library(network)
library(sna)
library(ggplot2)
library(countrycode)
library(parallel)
library(RSiena)



#dir <- "[Directory Information Here]"




## Import the network and monadic data in two separate objects
nets.total <- read.csv(file.path(ROOT, "work_saom/X060/relation_network.tab"), header=TRUE, row.names=1)
mons.total <- read.csv(file.path(ROOT, "work_saom/X060/relation_nodal.tab"), header=TRUE, row.names=1)
compchanges <- read.csv(file.path(ROOT, "work_saom/X060/compchanges.tab"), header=TRUE, row.names=1)



## Set some parameters
yr <- sort(unique(mons.total$year))
yr.lo <- yr[-length(yr)]

## List of names against which all data must be sorted. Also use this
## to build an empty matrix, i.e., for merging data
id <- sort(unique(mons.total$ccode1))
n <- length(id)
emat <- matrix(NA, length(id), length(id), dimnames=list(id,id))

## Make composition change list 
write.table(compchanges, "cc.dat",row.names=FALSE, col.names=FALSE,
            quote=FALSE)
changes <- sienaCompositionChangeFromFile("cc.dat")
names(changes) <- id

#########################################################
## Start with making the DV network matrices
#########################################################

##First the support matrices

for (jj in yr) {
  if ( jj == yr[1] ) {
    sm <- 1
    objs <- vector()
  }
  net <- nets.total[nets.total$year==jj,c("ccode1","ccode2","support_bin")]
  net <- acast(net, ccode1~ccode2, value.var="support_bin")
  out <- emat
  out[rownames(out) %in% rownames(net), colnames(out) %in% colnames(net)] <- net
  ## Order rows/columns
  out <- out[order(id),order(id)]
  ## Rename
  assign(paste("S.NET.", jj, sep=""), out)
  rm(net,out)
  objs[sm] <- paste("S.NET.", jj, sep="")
  sm <- sm + 1
}; rm(jj)
assign(paste("Ss.NET", sep=""), lapply(objs,get))
names(Ss.NET) <- objs
rm(objs,sm)
rm(list = ls(pattern = "^S.NET."))

## Now the threat matrices
for (jj in yr) {
  if ( jj == yr[1] ) {
    sm <- 1
    objs <- vector()
  }
  net <- nets.total[nets.total$year==jj,c("ccode1","ccode2","threat_bin_3")]
  net <- acast(net, ccode1~ccode2, value.var="threat_bin_3")
  out <- emat
  out[rownames(out) %in% rownames(net), colnames(out) %in% colnames(net)] <- net
  ## Order rows/columns
  out <- out[order(id),order(id)]
  ## Rename
  assign(paste("T.NET.", jj, sep=""), out)
  rm(net,out)
  objs[sm] <- paste("T.NET.", jj, sep="")
  sm <- sm + 1
}; rm(jj)
assign(paste("Ts.NET", sep=""), lapply(objs,get))
names(Ts.NET) <- objs
rm(objs,sm)
rm(list = ls(pattern = "^T.NET."))


#########################################################
## Then do all the covariates for both networks
#########################################################

## Names of dyadic and monadic covariates, plus behavior variable
d.cons <- names(nets.total)[ !names(nets.total) %in% c("ccode1","ccode2","year","support_bin","threat_bin_3") ]
m.cons <- names(mons.total)[ !names(mons.total) %in% c("ccode1","year","intra_intens") ]
b.cons <- names(mons.total)[ names(mons.total) %in% c("intra_intens") ]

## Set the monadic covariates as NxK matrices, where K is number of years. Remember to
## subtract the last year for RSiena estimation.
for (k in 1:length(m.cons)) {
  out <- matrix(nrow=length(id),ncol=0)
  rownames(out) <- id
  for (i in yr.lo) {
    tmp <- as.matrix(mons.total[mons.total$year==i,c(m.cons[k])])
    rownames(tmp) <- mons.total[mons.total$year==i,]$ccode1
    ## Bind on row names
    out <- cbind(out, tmp[match(rownames(out),rownames(tmp))])
    rm(tmp)
  }; rm(i)
  colnames(out) <- yr.lo
  out <- out[order(id),]
  assign(paste(m.cons[k],sep=""), out)
  rm(out)
}; rm(k)

#Now make behavior variables into matrix form, not subtracting last year
for (k in 1:length(b.cons)) {
  out <- matrix(nrow=length(id),ncol=0)
  rownames(out) <- id
  for (i in yr) {
    tmp <- as.matrix(mons.total[mons.total$year==i,c(b.cons[k])])
    rownames(tmp) <- mons.total[mons.total$year==i,]$ccode1
    ## Bind on row names
    out <- cbind(out, tmp[match(rownames(out),rownames(tmp))])
    rm(tmp)
  }; rm(i)
  colnames(out) <- yr
  out <- out[order(id),]
  assign(paste(b.cons[k],sep=""), out)
  rm(out)
}; rm(k)

## Save dyadic covariates as lists of matrices
for (k in 1:length(d.cons)) {
  for (i in yr.lo) {
    cov <- nets.total[nets.total$year==i, c("ccode1","ccode2",d.cons[k])]
    ## Here's the adjacency matrix
    cov <- acast(cov, ccode1~ccode2, value.var=d.cons[k])
    out <- emat
    ## Add cov data to the full matrix with all actors
    out[rownames(out) %in% rownames(cov), colnames(out) %in% colnames(cov)] <- cov
    diag(out) <- 0
    ## Order rows/columns
    out <- out[order(id),order(id)]
    assign(paste(d.cons[k], i, sep=""), out)
    rm(out,cov)
  }; rm(i)
  out.list <- lapply(paste(d.cons[k], yr.lo, sep=""), get)
  names(out.list) <- paste(d.cons[k], yr.lo, sep="")
  assign(paste(d.cons[k], sep=""), out.list)
  rm(out.list, list=paste(d.cons[k], yr.lo, sep=""))
}; rm(k)

####################################
## Convert objects to RSiena objects
####################################

yrhi <- yr
yrlo <- yr.lo
nlo <- length(yrlo)
nhi <- length(yrhi)
N <- n

## Dependent networkS
assign("support.net", sienaNet(array(do.call("c", Ss.NET), dim=c(N,N,nhi))))
assign("threat.net", sienaNet(array(do.call("c", Ts.NET), dim=c(N,N,nhi)),allowOnly=FALSE))
assign("intra_intens.net", sienaNet(intra_intens, type="behavior"))

## Changing dyadic covariates, omitting the constant ones
dyad.effects <- d.cons[ !d.cons %in% c("contiguous") ]
for (x in dyad.effects) {
  assign(paste(x, ".net", sep=""), varDyadCovar(array(do.call("c", get(x)), dim=c(N,N,nlo))))
}; rm(x)

## Constant dyadic covariates
for (x in d.cons[ d.cons %in% c("contiguous")]) {
  assign(paste(x, ".net", sep=""), coDyadCovar(get(x)[[nlo]]))
}

## Monadic covariates
monad.effects <- m.cons
for (x in monad.effects) {
  assign(paste(x, ".net", sep=""), varCovar(get(x)))
}; rm(x)


## Create an RSiena network object

netdata <- sienaDataCreate(
  support.net,
  threat.net,
  intra_intens.net,
  
  mi1.net,
  energypc.net,
  mmp_change.net,
  dem.net,
  s_t_ipdistance.net,
  support_IGO.net,
  support_ipcloseness.net,
  support_mi.net,
  threat_mi.net,
  support_mi_pot.net,
  threat_density_nbrs.net,
  intra_intens_nbrs.net,
  
  contiguous.net,
  idealpointdistance.net,
  IGOcommon.net,
  support_i_j_ipdistance.net,
  
  changes

)

