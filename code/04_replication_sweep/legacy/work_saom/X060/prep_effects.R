eff <- getEffects(netdata)

###Start with base model, and then add in network effects

## Support network equation first
d.support <- paste(
  c("contiguous",
    "idealpointdistance", "IGOcommon"
  ),
  ".net", sep=""
)
m.support_i <- paste(
  c("mi1", "energypc", "mmp_change", "dem", "intra_intens"),
  ".net", sep=""
)
m.support_j <- paste(
  c("mi1", "energypc", "mmp_change", "dem", "intra_intens"),
  ".net", sep=""
)


for (x in d.support) { # Dyadic effects
  eff <- includeEffects(eff, X, name="support.net", interaction1=x, include=T, type="eval")
}; rm(x)

for (x in m.support_j) { # Monadic effects for supported states
  eff <- includeEffects(eff, altX, name="support.net", interaction1=x, include=T, type="eval")
}; rm(x)

for (x in m.support_i) { # Monadic effects for supporters
  eff <- includeEffects(eff, egoX, name="support.net", interaction1=x, include=T, type="eval")
}; rm(x)


## Do a political similarity measure
eff <- includeEffects(eff, simX, name="support.net", interaction1="dem.net", include=T, type="eval")

## Endogenous support network effects
eff <- includeEffects(eff, recip, name="support.net", include=T, type="eval")

## Now the threat network equation
d.threat <- paste(
  c("contiguous",
    "idealpointdistance", "IGOcommon"),
  ".net", sep=""
)

m.threat_i <- paste(
  c("mi1", "energypc", "mmp_change", "dem", "intra_intens"),
  ".net", sep=""
)

m.threat_j <- paste(
  c("mi1", "energypc", "mmp_change", "dem", "intra_intens"),
  ".net", sep=""
)


for (x in d.threat) { # Dyadic covariates
  eff <- includeEffects(eff, X, name="threat.net", interaction1=x, include=T, type="eval")
}; rm(x)

for (x in m.threat_j) { # Monadic covariates for threatened states
  eff <- includeEffects(eff, altX, name="threat.net", interaction1=x, include=T, type="eval")
}; rm(x)

for (x in m.threat_i) { # Monadic covariates for threatening states
  eff <- includeEffects(eff, egoX, name="threat.net", interaction1=x, include=T, type="eval")
}; rm(x)


## Do a political similarity measure
eff <- includeEffects(eff, simX, name="threat.net", interaction1="dem.net", include=T, type="eval")

## Endogenous threat network effects
eff <- includeEffects(eff, recip, name="threat.net", include=T, type="eval")

#Now with civil war intensity as a dependent behavior variable
m.intra_i <- paste(
  c("mi1", "energypc", "mmp_change", "dem"),
  ".net", sep=""
)

for (x in m.intra_i) { # Monadic covariates for threatened states
  eff <- includeEffects(eff, effFrom, name="intra_intens.net", interaction1=x, include=T, type="eval")
}; rm(x)


## Add in some network variables into all following models
eff2<-eff

## Add the basic power from the support network
eff2 <- includeEffects(eff2, altX, name="threat.net", interaction1="support_mi.net", include=T, type="eval")
eff2 <- includeEffects(eff2, egoX, name="threat.net", interaction1="support_mi.net", include=T, type="eval")
eff2 <- includeEffects(eff2, effFrom, name="intra_intens.net", interaction1="support_mi.net", include=T, type="eval")

## Include the cross-products
eff2 <- includeEffects(eff2, crprodRecip, name="support.net", interaction1="threat.net", include=T, type="eval")
eff2 <- includeEffects(eff2, crprodRecip, name="threat.net", interaction1="support.net", include=T, type="eval")

## Add in some higher order cross-network effects
eff2 <- includeEffects(eff2, inPopIntn, name="support.net",
                       interaction1="threat.net", include=T, type="eval")
eff2 <- includeEffects(eff2, sharedIn, name="support.net",
                      interaction1="threat.net", include=T, type="eval")

## Add in information about the threat network for the behavior equation
eff2 <- includeEffects(eff2, effFrom, name="intra_intens.net", interaction1="threat_mi.net", include=T, type="eval")

## Create two models based on different ways to measure friend of my enemy effects
#eff1a<-eff0
#eff1a<- includeEffects(eff1a, mixedInXW, name="threat.net",
                     # interaction1="support.net", include=T, type="eval")

eff2 <- includeEffects(eff2, cl.XWX2, name="threat.net",
                      interaction1="support.net", include=T, type="eval")

##Add in ideal point distance in the support networks, which will switch out in next models
eff2 <- includeEffects(eff2, altX, name="threat.net", interaction1="support_ipcloseness.net", include=T, type="eval")
eff2 <- includeInteraction(eff2, altX, altX, name="threat.net", interaction1=c("support_mi.net", "support_ipcloseness.net"), include=T, type="eval")
eff2 <- includeEffects(eff2, egoX, name="threat.net", interaction1="support_ipcloseness.net", include=T, type="eval")
eff2 <- includeInteraction(eff2, egoX, egoX, name="threat.net", interaction1=c("support_mi.net", "support_ipcloseness.net"), include=T, type="eval")
eff2 <- includeEffects(eff2, effFrom, name="intra_intens.net", interaction1="support_ipcloseness.net", include=T, type="eval")
eff2 <- includeInteraction(eff2, effFrom, effFrom, name="intra_intens.net", interaction1=c("support_mi.net", "support_ipcloseness.net"), include=T, type="eval")



