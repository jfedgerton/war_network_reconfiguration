## function to correct for an issue in the Weeks model
## for the AMEN continuation from the previous starting points

function(Z,EZ,rho,Y) {

    ## simulates Z under the contraints
    ## (1)  Y[i,j]=1   => Z[i,j]>0
    ## (2)  Y[i,j]=0   => Z[i,j]<0

    sz<-sqrt(1-rho^2)
    ut<-upper.tri(EZ)
    lt<-lower.tri(EZ)

    Y[is.na(Y)]<- -1
    for(y in c((-1):1))
    {
        lb<-c(-Inf,-Inf,0)[y+2] ; ub<-c(Inf,0,Inf)[y+2]

        up<- ut & Y==y
        ez<- EZ[up] + rho*( t(Z)[up]  - t(EZ)[up] )
        lbUnif <- pnorm((lb-ez)/sz)
        ubUnif <- pnorm((ub-ez)/sz)
        unif <- runif(sum(up),lbUnif,ubUnif)
        unif[unif==1]=1-(1e-16)
        Z[up]<-ez+sz*qnorm(unif)

        up<- lt & Y==y
        ez<- EZ[up] + rho*( t(Z)[up]  - t(EZ)[up] )
        lbUnif <- pnorm((lb-ez)/sz)
        ubUnif <- pnorm((ub-ez)/sz)
        unif <- runif(sum(up),lbUnif,ubUnif)
        unif[unif==1]=1-(1e-16)
        Z[up]<-ez+sz*qnorm(unif)
    }

    diag(Z)<-rnorm(nrow(Z),diag(EZ),1)
    Z
}

## when you call this function, use the following code:
#assignInNamespace("rZ_bin_fc", rZ <- bin <- fc2, pos="package:amen")
