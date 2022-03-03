# pop same

set.seed(23091990)
library(stringr)
library(mvtnorm)
library(OpenMx)
library(purrr)
library(testit)
mxRunDave <- possibly(mxRun, otherwise = -999)

groupsize <- 13
nitems <- 1

#define your factor loading parameters (depends, of course, on "nitems")
generalActorLoadings <- 1
generalPartnerLoadings <- 1
generalRelationshipLoadings <- 1

#define lv var-cov parameters (excluding residuals)
grouplevelVarCov <- 0

generalActorVar <- .35
generalPartnerVar <- .15
generalRelationshipVar <- .5
generalActorPartnerCov <-  .04582576 #(implying correlation of .2)
generalRelationshipCov <- .05        #(implying correlation of .1)


meanvector <- c(rep(0, 13*12))



#actors (extended, will be shortened after partner is also defined)
actor <- c()
for (i in 1 : groupsize){
  vekact <- rep(i, groupsize)
  actor <- c(actor, vekact)
}

#partners (extended, will be shortened in next step)
partner <- c()
for (i in 1: groupsize){
  vekpar <- c(1:groupsize)
  partner <- c(partner, vekpar)
}

#make data frame of what you have so far
sofar <- as.data.frame(cbind(actor, partner))
sofar$same <- ifelse(sofar$actor == sofar$partner,1,0)
sofar <-subset(sofar, same ==0)
sofar$same <- NULL

#extend to number of items
sofarlist <- list()
for (i in 1:nitems){
  sofarlist[[i]] <- sofar
}
sofar2 <- data.frame()
for (i in 1:length(sofarlist)){
  sofarnow <- sofarlist[[i]]
  sofar2 <- rbind(sofar2, sofarnow)
}
sofar2$itemnumber <- rep(c(1:nitems), each = groupsize*(groupsize-1))


######

#make latent variable names

#group-level: one for every item
grouplevellv<-c()
for (i in 1:nitems){
  if (i < 10){
    grouplevellv[i]<-paste0("groupi0", i)
  }
  else {grouplevellv[i]<-paste0("groupi", i)}
}

#

#person-level: general actor effects
generalactor <- c()
for (i in 1:groupsize){
  if(i < 10){
    generalactor[i]<- paste0("generalActor_0",i)
  }
  else{
    generalactor[i]<- paste0("generalActor_",i)
  }
}


#

#person-level: general partner effects
generalpartner <- c()
for (i in 1:groupsize){
  if(i < 10){
    generalpartner[i]<- paste0("generalPartner_0",i)
  }
  else{
    generalpartner[i]<- paste0("generalPartner_",i)
  }
}




#dyad-level: general relationship effects
trickysubset <- subset(sofar2, itemnumber==1)
generalRelationship <- c()
for (i in 1:nrow(trickysubset)){
  
  if(trickysubset$actor[i] >= 10 && trickysubset$partner[i] >= 10){
    generalRelationship[i] <- paste0("generalRelationship_ap_",trickysubset$actor[i],trickysubset$partner[i])
  }
  
  if(trickysubset$actor[i] >= 10 && trickysubset$partner[i] < 10){
    generalRelationship[i] <- paste0("generalRelationship_ap_",trickysubset$actor[i],"0",trickysubset$partner[i])
  }
  
  if(trickysubset$actor[i] < 10 && trickysubset$partner[i] >= 10){
    generalRelationship[i] <- paste0("generalRelationship_ap_0",trickysubset$actor[i],trickysubset$partner[i])
  }
  
  if(trickysubset$actor[i] < 10 && trickysubset$partner[i] < 10){
    generalRelationship[i] <- paste0("generalRelationship_ap_0",trickysubset$actor[i],"0",trickysubset$partner[i])
  }
}



#all latent variables
all_lvnames <- c(grouplevellv, generalactor, generalpartner, generalRelationship)
factorloadingmatrix <- as.data.frame(matrix(data=0, nrow = nrow(sofar2),
                                            ncol = length(all_lvnames)))
factorloadingmatrix <- cbind(sofar2, factorloadingmatrix)
names(factorloadingmatrix)<- c(names(sofar2),all_lvnames)
factorloadingmatrix2 <- factorloadingmatrix[,4:ncol(factorloadingmatrix)]



################################################################

#fill in the matrix
for (i in 1: nrow(factorloadingmatrix2)){
  for (j in 1: ncol(factorloadingmatrix2)){
    
    if (
      startsWith(names(factorloadingmatrix2[j]), "group") == T &&
      factorloadingmatrix$itemnumber[i] == as.numeric(str_sub(names(factorloadingmatrix2[j]),7,8))
    )
    {
      factorloadingmatrix2[i,j] <- 1
    }
    
    if (
      startsWith(names(factorloadingmatrix2[j]), "generalActor") == T &&
      factorloadingmatrix$actor[i] == as.numeric(str_sub(names(factorloadingmatrix2[j]),14,15)) 
    )
    {
      factorloadingmatrix2[i,j] <- 1
    }
    
    if (
      startsWith(names(factorloadingmatrix2[j]), "generalPartner") == T &&
      factorloadingmatrix$partner[i] == as.numeric(str_sub(names(factorloadingmatrix2[j]),16,17))
    )
    {
      factorloadingmatrix2[i,j] <- 1
    }
    
    if (
      startsWith(names(factorloadingmatrix2[j]), "generalRelationship_ap_") == T &&
      factorloadingmatrix$actor[i] == as.numeric(str_sub(names(factorloadingmatrix2[j]),24,25)) &&
      factorloadingmatrix$partner[i] == as.numeric(str_sub(names(factorloadingmatrix2[j]),26,27)) 
    )
    {
      factorloadingmatrix2[i,j] <- 1
    }
  }
}

#save the matrix (could save complete matrix,
#but itemspecific Relationship effects are probably better treated as residuals)
pop_factorloadingMatrix <- as.matrix(factorloadingmatrix2)





########################################

#now latent variable covariances (excluding residuals)

lvcovarianceMatrix <- matrix(data = 0, nrow = length(colnames(pop_factorloadingMatrix)), 
                             ncol = length(colnames(pop_factorloadingMatrix)))
colnames(lvcovarianceMatrix) <- colnames(pop_factorloadingMatrix)
rownames(lvcovarianceMatrix) <- colnames(pop_factorloadingMatrix)


#fill in the matrix

for (i in 1:length(rownames(lvcovarianceMatrix))){
  for (j in 1:length(colnames(lvcovarianceMatrix))){
    
    if(i == j){
      if(startsWith(rownames(lvcovarianceMatrix)[i], "group") == T){
        lvcovarianceMatrix[i,j] <- grouplevelVarCov
      }
      if(startsWith(rownames(lvcovarianceMatrix)[i], "generalActor") == T){
        lvcovarianceMatrix[i,j] <- generalActorVar
      }
      if(startsWith(rownames(lvcovarianceMatrix)[i], "generalPartner") == T){
        lvcovarianceMatrix[i,j] <- generalPartnerVar
      }
      if(startsWith(rownames(lvcovarianceMatrix)[i], "generalRelationship") == T){
        lvcovarianceMatrix[i,j] <- generalRelationshipVar
      }
    }
    else {
      
      
      if(startsWith(rownames(lvcovarianceMatrix)[i], "generalActor")==T & 
         startsWith(colnames(lvcovarianceMatrix)[j], "generalPartner")==T &
         str_sub(rownames(lvcovarianceMatrix)[i],14,15) == str_sub(colnames(lvcovarianceMatrix)[j],16,17)
      )
      {
        lvcovarianceMatrix[i,j] <- generalActorPartnerCov
      }
      if(startsWith(rownames(lvcovarianceMatrix)[i], "generalPartner")==T & 
         startsWith(colnames(lvcovarianceMatrix)[j], "generalActor")==T &
         str_sub(rownames(lvcovarianceMatrix)[i],16,17) == str_sub(colnames(lvcovarianceMatrix)[j],14,15)
      )
      {
        lvcovarianceMatrix[i,j] <- generalActorPartnerCov
      }
      if(startsWith(rownames(lvcovarianceMatrix)[i], "generalRelationship")==T & 
         startsWith(colnames(lvcovarianceMatrix)[j], "generalRelationship")==T &
         str_sub(rownames(lvcovarianceMatrix)[i],24,27) == paste0(str_sub(colnames(lvcovarianceMatrix)[j],26,27), str_sub(colnames(lvcovarianceMatrix)[j],24,25))
      )
      {
        lvcovarianceMatrix[i,j] <- generalRelationshipCov
      }
      
      
    }
    
  }
}

pop_lvcovarianceMatrix <- lvcovarianceMatrix



#############################################################


############################################


#make population covariance matrix and population data
pcov <- pop_factorloadingMatrix%*%pop_lvcovarianceMatrix%*%t(pop_factorloadingMatrix)
PopSame  <- rmvnorm(100000, mean = meanvector, sigma = pcov)

PopSame <- as.data.frame(PopSame)

names(PopSame) <- c(
  "xab1",
  "xac1",
  "xad1",
  "xae1",
  "xaf1",
  "xag1",
  "xah1",
  "xai1",
  "xaj1",
  "xak1",
  "xal1",
  "xam1",
  
  "xba1",
  "xbc1",
  "xbd1",
  "xbe1",
  "xbf1",
  "xbg1",
  "xbh1",
  "xbi1",
  "xbj1",
  "xbk1",
  "xbl1",
  "xbm1",
  
  
  "xca1",
  "xcb1",
  "xcd1",
  "xce1",
  "xcf1",
  "xcg1",
  "xch1",
  "xci1",
  "xcj1",
  "xck1",
  "xcl1",
  "xcm1",
  
  "xda1",
  "xdb1",
  "xdc1",
  "xde1",
  "xdf1",
  "xdg1",
  "xdh1",
  "xdi1",
  "xdj1",
  "xdk1",
  "xdl1",
  "xdm1",
  
  "xea1",
  "xeb1",
  "xec1",
  "xed1",
  "xef1",
  "xeg1",
  "xeh1",
  "xei1",
  "xej1",
  "xek1",
  "xel1",
  "xem1",
  
  
  "xfa1",
  "xfb1",
  "xfc1",
  "xfd1",
  "xfe1",
  "xfg1",
  "xfh1",
  "xfi1",
  "xfj1",
  "xfk1",
  "xfl1",
  "xfm1",
  
  "xga1",
  "xgb1",
  "xgc1",
  "xgd1",
  "xge1",
  "xgf1",
  "xgh1",
  "xgi1",
  "xgj1",
  "xgk1",
  "xgl1",
  "xgm1",
  
  
  "xha1",
  "xhb1",
  "xhc1",
  "xhd1",
  "xhe1",
  "xhf1",
  "xhg1",
  "xhi1",
  "xhj1",
  "xhk1",
  "xhl1",
  "xhm1",
  
  
  "xia1",
  "xib1",
  "xic1",
  "xid1",
  "xie1",
  "xif1",
  "xig1",
  "xih1",
  "xij1",
  "xik1",
  "xil1",
  "xim1",
  
  
  "xja1",
  "xjb1",
  "xjc1",
  "xjd1",
  "xje1",
  "xjf1",
  "xjg1",
  "xjh1",
  "xji1",
  "xjk1",
  "xjl1",
  "xjm1",
  
  
  "xka1",
  "xkb1",
  "xkc1",
  "xkd1",
  "xke1",
  "xkf1",
  "xkg1",
  "xkh1",
  "xki1",
  "xkj1",
  "xkl1",
  "xkm1",
  
  
  "xla1",
  "xlb1",
  "xlc1",
  "xld1",
  "xle1",
  "xlf1",
  "xlg1",
  "xlh1",
  "xli1",
  "xlj1",
  "xlk1",
  "xlm1",
  
  
  "xma1",
  "xmb1",
  "xmc1",
  "xmd1",
  "xme1",
  "xmf1",
  "xmg1",
  "xmh1",
  "xmi1",
  "xmj1",
  "xmk1",
  "xml1")

#save population data
save(PopSame, file = "PopSame.rda")





###########################################################################




#pop different



groupsize <- 13
nitems <- 1

#define your factor loading parameters (depends, of course, on "nitems")
generalActorLoadings <- 1
generalPartnerLoadings <- 1
generalRelationshipLoadings <- 1

#define lv var-cov parameters (excluding residuals)
grouplevelVarCov <- 0

generalActorVar <- .15
generalPartnerVar <- .06428571
generalRelationshipVar <- .5
generalActorPartnerCov <- .01963961  #(implying correlation of .2)
generalRelationshipCov <- .05        #(implying correlation of .1)


meanvector <- c(rep(0, 13*12))



#actors (extended, will be shortened after partner is also defined)
actor <- c()
for (i in 1 : groupsize){
  vekact <- rep(i, groupsize)
  actor <- c(actor, vekact)
}

#partners (extended, will be shortened in next step)
partner <- c()
for (i in 1: groupsize){
  vekpar <- c(1:groupsize)
  partner <- c(partner, vekpar)
}

#make data frame of what you have so far
sofar <- as.data.frame(cbind(actor, partner))
sofar$same <- ifelse(sofar$actor == sofar$partner,1,0)
sofar <-subset(sofar, same ==0)
sofar$same <- NULL

#extend to number of items
sofarlist <- list()
for (i in 1:nitems){
  sofarlist[[i]] <- sofar
}
sofar2 <- data.frame()
for (i in 1:length(sofarlist)){
  sofarnow <- sofarlist[[i]]
  sofar2 <- rbind(sofar2, sofarnow)
}
sofar2$itemnumber <- rep(c(1:nitems), each = groupsize*(groupsize-1))


######

#make latent variable names

#group-level: one for every item
grouplevellv<-c()
for (i in 1:nitems){
  if (i < 10){
    grouplevellv[i]<-paste0("groupi0", i)
  }
  else {grouplevellv[i]<-paste0("groupi", i)}
}

#

#person-level: general actor effects
generalactor <- c()
for (i in 1:groupsize){
  if(i < 10){
    generalactor[i]<- paste0("generalActor_0",i)
  }
  else{
    generalactor[i]<- paste0("generalActor_",i)
  }
}


#

#person-level: general partner effects
generalpartner <- c()
for (i in 1:groupsize){
  if(i < 10){
    generalpartner[i]<- paste0("generalPartner_0",i)
  }
  else{
    generalpartner[i]<- paste0("generalPartner_",i)
  }
}




#dyad-level: general relationship effects
trickysubset <- subset(sofar2, itemnumber==1)
generalRelationship <- c()
for (i in 1:nrow(trickysubset)){
  
  if(trickysubset$actor[i] >= 10 && trickysubset$partner[i] >= 10){
    generalRelationship[i] <- paste0("generalRelationship_ap_",trickysubset$actor[i],trickysubset$partner[i])
  }
  
  if(trickysubset$actor[i] >= 10 && trickysubset$partner[i] < 10){
    generalRelationship[i] <- paste0("generalRelationship_ap_",trickysubset$actor[i],"0",trickysubset$partner[i])
  }
  
  if(trickysubset$actor[i] < 10 && trickysubset$partner[i] >= 10){
    generalRelationship[i] <- paste0("generalRelationship_ap_0",trickysubset$actor[i],trickysubset$partner[i])
  }
  
  if(trickysubset$actor[i] < 10 && trickysubset$partner[i] < 10){
    generalRelationship[i] <- paste0("generalRelationship_ap_0",trickysubset$actor[i],"0",trickysubset$partner[i])
  }
}



#all latent variables
all_lvnames <- c(grouplevellv, generalactor, generalpartner, generalRelationship)
factorloadingmatrix <- as.data.frame(matrix(data=0, nrow = nrow(sofar2),
                                            ncol = length(all_lvnames)))
factorloadingmatrix <- cbind(sofar2, factorloadingmatrix)
names(factorloadingmatrix)<- c(names(sofar2),all_lvnames)
factorloadingmatrix2 <- factorloadingmatrix[,4:ncol(factorloadingmatrix)]



################################################################

#fill in the matrix
for (i in 1: nrow(factorloadingmatrix2)){
  for (j in 1: ncol(factorloadingmatrix2)){
    
    if (
      startsWith(names(factorloadingmatrix2[j]), "group") == T &&
      factorloadingmatrix$itemnumber[i] == as.numeric(str_sub(names(factorloadingmatrix2[j]),7,8))
    )
    {
      factorloadingmatrix2[i,j] <- 1
    }
    
    if (
      startsWith(names(factorloadingmatrix2[j]), "generalActor") == T &&
      factorloadingmatrix$actor[i] == as.numeric(str_sub(names(factorloadingmatrix2[j]),14,15)) 
    )
    {
      factorloadingmatrix2[i,j] <- 1
    }
    
    if (
      startsWith(names(factorloadingmatrix2[j]), "generalPartner") == T &&
      factorloadingmatrix$partner[i] == as.numeric(str_sub(names(factorloadingmatrix2[j]),16,17))
    )
    {
      factorloadingmatrix2[i,j] <- 1
    }
    
    if (
      startsWith(names(factorloadingmatrix2[j]), "generalRelationship_ap_") == T &&
      factorloadingmatrix$actor[i] == as.numeric(str_sub(names(factorloadingmatrix2[j]),24,25)) &&
      factorloadingmatrix$partner[i] == as.numeric(str_sub(names(factorloadingmatrix2[j]),26,27)) 
    )
    {
      factorloadingmatrix2[i,j] <- 1
    }
  }
}

#save the matrix (could save complete matrix,
#but itemspecific Relationship effects are probably better treated as residuals)
pop_factorloadingMatrix <- as.matrix(factorloadingmatrix2)





########################################

#now latent variable covariances (excluding residuals)

lvcovarianceMatrix <- matrix(data = 0, nrow = length(colnames(pop_factorloadingMatrix)), 
                             ncol = length(colnames(pop_factorloadingMatrix)))
colnames(lvcovarianceMatrix) <- colnames(pop_factorloadingMatrix)
rownames(lvcovarianceMatrix) <- colnames(pop_factorloadingMatrix)


#fill in the matrix

for (i in 1:length(rownames(lvcovarianceMatrix))){
  for (j in 1:length(colnames(lvcovarianceMatrix))){
    
    if(i == j){
      if(startsWith(rownames(lvcovarianceMatrix)[i], "group") == T){
        lvcovarianceMatrix[i,j] <- grouplevelVarCov
      }
      if(startsWith(rownames(lvcovarianceMatrix)[i], "generalActor") == T){
        lvcovarianceMatrix[i,j] <- generalActorVar
      }
      if(startsWith(rownames(lvcovarianceMatrix)[i], "generalPartner") == T){
        lvcovarianceMatrix[i,j] <- generalPartnerVar
      }
      if(startsWith(rownames(lvcovarianceMatrix)[i], "generalRelationship") == T){
        lvcovarianceMatrix[i,j] <- generalRelationshipVar
      }
    }
    else {
      
      
      if(startsWith(rownames(lvcovarianceMatrix)[i], "generalActor")==T & 
         startsWith(colnames(lvcovarianceMatrix)[j], "generalPartner")==T &
         str_sub(rownames(lvcovarianceMatrix)[i],14,15) == str_sub(colnames(lvcovarianceMatrix)[j],16,17)
      )
      {
        lvcovarianceMatrix[i,j] <- generalActorPartnerCov
      }
      if(startsWith(rownames(lvcovarianceMatrix)[i], "generalPartner")==T & 
         startsWith(colnames(lvcovarianceMatrix)[j], "generalActor")==T &
         str_sub(rownames(lvcovarianceMatrix)[i],16,17) == str_sub(colnames(lvcovarianceMatrix)[j],14,15)
      )
      {
        lvcovarianceMatrix[i,j] <- generalActorPartnerCov
      }
      if(startsWith(rownames(lvcovarianceMatrix)[i], "generalRelationship")==T & 
         startsWith(colnames(lvcovarianceMatrix)[j], "generalRelationship")==T &
         str_sub(rownames(lvcovarianceMatrix)[i],24,27) == paste0(str_sub(colnames(lvcovarianceMatrix)[j],26,27), str_sub(colnames(lvcovarianceMatrix)[j],24,25))
      )
      {
        lvcovarianceMatrix[i,j] <- generalRelationshipCov
      }
      
      
    }
    
  }
}

pop_lvcovarianceMatrix <- lvcovarianceMatrix



#############################################################


############################################


#make population covariance matrix and population data
pcov <- pop_factorloadingMatrix%*%pop_lvcovarianceMatrix%*%t(pop_factorloadingMatrix)
PopDifferent  <- rmvnorm(100000, mean = meanvector, sigma = pcov)

PopDifferent <- as.data.frame(PopDifferent)

names(PopDifferent) <- c(
  "xab1",
  "xac1",
  "xad1",
  "xae1",
  "xaf1",
  "xag1",
  "xah1",
  "xai1",
  "xaj1",
  "xak1",
  "xal1",
  "xam1",
  
  "xba1",
  "xbc1",
  "xbd1",
  "xbe1",
  "xbf1",
  "xbg1",
  "xbh1",
  "xbi1",
  "xbj1",
  "xbk1",
  "xbl1",
  "xbm1",
  
  
  "xca1",
  "xcb1",
  "xcd1",
  "xce1",
  "xcf1",
  "xcg1",
  "xch1",
  "xci1",
  "xcj1",
  "xck1",
  "xcl1",
  "xcm1",
  
  "xda1",
  "xdb1",
  "xdc1",
  "xde1",
  "xdf1",
  "xdg1",
  "xdh1",
  "xdi1",
  "xdj1",
  "xdk1",
  "xdl1",
  "xdm1",
  
  "xea1",
  "xeb1",
  "xec1",
  "xed1",
  "xef1",
  "xeg1",
  "xeh1",
  "xei1",
  "xej1",
  "xek1",
  "xel1",
  "xem1",
  
  
  "xfa1",
  "xfb1",
  "xfc1",
  "xfd1",
  "xfe1",
  "xfg1",
  "xfh1",
  "xfi1",
  "xfj1",
  "xfk1",
  "xfl1",
  "xfm1",
  
  "xga1",
  "xgb1",
  "xgc1",
  "xgd1",
  "xge1",
  "xgf1",
  "xgh1",
  "xgi1",
  "xgj1",
  "xgk1",
  "xgl1",
  "xgm1",
  
  
  "xha1",
  "xhb1",
  "xhc1",
  "xhd1",
  "xhe1",
  "xhf1",
  "xhg1",
  "xhi1",
  "xhj1",
  "xhk1",
  "xhl1",
  "xhm1",
  
  
  "xia1",
  "xib1",
  "xic1",
  "xid1",
  "xie1",
  "xif1",
  "xig1",
  "xih1",
  "xij1",
  "xik1",
  "xil1",
  "xim1",
  
  
  "xja1",
  "xjb1",
  "xjc1",
  "xjd1",
  "xje1",
  "xjf1",
  "xjg1",
  "xjh1",
  "xji1",
  "xjk1",
  "xjl1",
  "xjm1",
  
  
  "xka1",
  "xkb1",
  "xkc1",
  "xkd1",
  "xke1",
  "xkf1",
  "xkg1",
  "xkh1",
  "xki1",
  "xkj1",
  "xkl1",
  "xkm1",
  
  
  "xla1",
  "xlb1",
  "xlc1",
  "xld1",
  "xle1",
  "xlf1",
  "xlg1",
  "xlh1",
  "xli1",
  "xlj1",
  "xlk1",
  "xlm1",
  
  
  "xma1",
  "xmb1",
  "xmc1",
  "xmd1",
  "xme1",
  "xmf1",
  "xmg1",
  "xmh1",
  "xmi1",
  "xmj1",
  "xmk1",
  "xml1")

#save population data
save(PopDifferent, file = "PopDifferent.rda")



############################################################################################




#starting values

generalActorVar <- .35
generalPartnerVar <- .15
generalRelationshipVar <- .5
generalActorPartnerCov <-  .04582576 #(implying correlation of .2)
generalRelationshipCov <-  .05       #(implying correlation of .1)

differentActorVar <- .15
differentPartnerVar <- .06428571
differentRelationshipVar <- .5
differentActorPartnerCov <- .01963961  #(implying correlation of .2)
differentRelationshipCov <- .05        #(implying correlation of .1)



#############################################################################################


#base5

base5 <- mxModel("base5",
                 type="RAM",
                 manifestVars = c(
                   "xab1",
                   "xac1",
                   "xad1",
                   "xae1",
                   
                   "xba1",
                   "xbc1",
                   "xbd1",
                   "xbe1",
                   
                   "xca1",
                   "xcb1",
                   "xcd1",
                   "xce1",
                   
                   "xda1",
                   "xdb1",
                   "xdc1",
                   "xde1",
                   
                   "xea1",
                   "xeb1",
                   "xec1",
                   "xed1"),
                 
                 latentVars = c(
                   "Aon1"     ,      
                   "Bon1"     ,
                   "Con1"     ,
                   "Don1"     ,
                   "Eon1"     ,
                   
                   "onA1"     ,
                   "onB1"     ,
                   "onC1"     ,
                   "onD1"     ,
                   "onE1"),
                 
                 #person level
                 
                 mxPath(from="Aon1", to= "xab1", arrows = 1, free = F, values=1, labels = "la1"),
                 mxPath(from="Aon1", to= "xac1", arrows = 1, free = F, values=1, labels = "la1"),
                 mxPath(from="Aon1", to= "xad1", arrows = 1, free = F, values=1, labels = "la1"),
                 mxPath(from="Aon1", to= "xae1", arrows = 1, free = F, values=1, labels = "la1"),
                 
                 
                 mxPath(from="Bon1", to= "xba1", arrows = 1, free = F, values=1, labels = "la1"),
                 mxPath(from="Bon1", to= "xbc1", arrows = 1, free = F, values=1, labels = "la1"),
                 mxPath(from="Bon1", to= "xbd1", arrows = 1, free = F, values=1, labels = "la1"),
                 mxPath(from="Bon1", to= "xbe1", arrows = 1, free = F, values=1, labels = "la1"),
                 
                 
                 mxPath(from="Con1", to= "xca1", arrows = 1, free = F, values=1, labels = "la1"),
                 mxPath(from="Con1", to= "xcb1", arrows = 1, free = F, values=1, labels = "la1"),
                 mxPath(from="Con1", to= "xcd1", arrows = 1, free = F, values=1, labels = "la1"),
                 mxPath(from="Con1", to= "xce1", arrows = 1, free = F, values=1, labels = "la1"),
                 
                 
                 mxPath(from="Don1", to= "xda1", arrows = 1, free = F, values=1, labels = "la1"),
                 mxPath(from="Don1", to= "xdb1", arrows = 1, free = F, values=1, labels = "la1"),
                 mxPath(from="Don1", to= "xdc1", arrows = 1, free = F, values=1, labels = "la1"),
                 mxPath(from="Don1", to= "xde1", arrows = 1, free = F, values=1, labels = "la1"),
                 
                 
                 mxPath(from="Eon1", to= "xea1", arrows = 1, free = F, values=1, labels = "la1"),
                 mxPath(from="Eon1", to= "xeb1", arrows = 1, free = F, values=1, labels = "la1"),
                 mxPath(from="Eon1", to= "xec1", arrows = 1, free = F, values=1, labels = "la1"),
                 mxPath(from="Eon1", to= "xed1", arrows = 1, free = F, values=1, labels = "la1"),
                 
                 
                 
                 mxPath(from="onA1", to= "xba1", arrows = 1, free = F, values=1, labels = "lp1"),
                 mxPath(from="onA1", to= "xca1", arrows = 1, free = F, values=1, labels = "lp1"),
                 mxPath(from="onA1", to= "xda1", arrows = 1, free = F, values=1, labels = "lp1"),
                 mxPath(from="onA1", to= "xea1", arrows = 1, free = F, values=1, labels = "lp1"),
                 
                 
                 mxPath(from="onB1", to= "xab1", arrows = 1, free = F, values=1, labels = "lp1"),
                 mxPath(from="onB1", to= "xcb1", arrows = 1, free = F, values=1, labels = "lp1"),
                 mxPath(from="onB1", to= "xdb1", arrows = 1, free = F, values=1, labels = "lp1"),
                 mxPath(from="onB1", to= "xeb1", arrows = 1, free = F, values=1, labels = "lp1"),
                 
                 
                 mxPath(from="onC1", to= "xac1", arrows = 1, free = F, values=1, labels = "lp1"),
                 mxPath(from="onC1", to= "xbc1", arrows = 1, free = F, values=1, labels = "lp1"),
                 mxPath(from="onC1", to= "xdc1", arrows = 1, free = F, values=1, labels = "lp1"),
                 mxPath(from="onC1", to= "xec1", arrows = 1, free = F, values=1, labels = "lp1"),
                 
                 
                 mxPath(from="onD1", to= "xad1", arrows = 1, free = F, values=1, labels = "lp1"),
                 mxPath(from="onD1", to= "xbd1", arrows = 1, free = F, values=1, labels = "lp1"),
                 mxPath(from="onD1", to= "xcd1", arrows = 1, free = F, values=1, labels = "lp1"),
                 mxPath(from="onD1", to= "xed1", arrows = 1, free = F, values=1, labels = "lp1"),
                 
                 
                 mxPath(from="onE1", to= "xae1", arrows = 1, free = F, values=1, labels = "lp1"),
                 mxPath(from="onE1", to= "xbe1", arrows = 1, free = F, values=1, labels = "lp1"),
                 mxPath(from="onE1", to= "xce1", arrows = 1, free = F, values=1, labels = "lp1"),
                 mxPath(from="onE1", to= "xde1", arrows = 1, free = F, values=1, labels = "lp1"),
                 
                 
                 
                 mxPath(from="Aon1", arrows=2, free=T,  values=generalActorVar, labels = "av1"),
                 mxPath(from="Bon1", arrows=2, free=T,  values=generalActorVar, labels = "av1"),
                 mxPath(from="Con1", arrows=2, free=T,  values=generalActorVar, labels = "av1"),
                 mxPath(from="Don1", arrows=2, free=T,  values=generalActorVar, labels = "av1"),
                 mxPath(from="Eon1", arrows=2, free=T,  values=generalActorVar, labels = "av1"),
                 
                 
                 mxPath(from="onA1", arrows=2, free=T,  values=generalPartnerVar, labels = "pv1"),
                 mxPath(from="onB1", arrows=2, free=T,  values=generalPartnerVar, labels = "pv1"),
                 mxPath(from="onC1", arrows=2, free=T,  values=generalPartnerVar, labels = "pv1"),
                 mxPath(from="onD1", arrows=2, free=T,  values=generalPartnerVar, labels = "pv1"),
                 mxPath(from="onE1", arrows=2, free=T,  values=generalPartnerVar, labels = "pv1"),
                 
                 
                 mxPath(from="Aon1", to = "onA1", arrows=2, free=T,  values=generalActorPartnerCov, labels = "apc1"),
                 mxPath(from="Bon1", to = "onB1", arrows=2, free=T,  values=generalActorPartnerCov, labels = "apc1"),
                 mxPath(from="Con1", to = "onC1", arrows=2, free=T,  values=generalActorPartnerCov, labels = "apc1"),
                 mxPath(from="Don1", to = "onD1", arrows=2, free=T,  values=generalActorPartnerCov, labels = "apc1"),
                 mxPath(from="Eon1", to = "onE1", arrows=2, free=T,  values=generalActorPartnerCov, labels = "apc1"),
                 
                 mxPath(from="one", to = "onA1", arrows=1, free=F,  values=0, labels = "mean1"),
                 mxPath(from="one", to = "onB1", arrows=1, free=F,  values=0, labels = "mean1"),
                 mxPath(from="one", to = "onC1", arrows=1, free=F,  values=0, labels = "mean1"),
                 mxPath(from="one", to = "onD1", arrows=1, free=F,  values=0, labels = "mean1"),
                 mxPath(from="one", to = "onE1", arrows=1, free=F,  values=0, labels = "mean1"),
                 
                 
                 
                 #dyad level
                 
                 mxPath(from="xab1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                 mxPath(from="xac1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                 mxPath(from="xad1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                 mxPath(from="xae1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                 
                 
                 mxPath(from="xba1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                 mxPath(from="xbc1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                 mxPath(from="xbd1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                 mxPath(from="xbe1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                 
                 
                 mxPath(from="xca1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                 mxPath(from="xcb1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                 mxPath(from="xcd1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                 mxPath(from="xce1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                 
                 
                 mxPath(from="xda1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                 mxPath(from="xdb1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                 mxPath(from="xdc1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                 mxPath(from="xde1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                 
                 
                 mxPath(from="xea1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                 mxPath(from="xeb1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                 mxPath(from="xec1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                 mxPath(from="xed1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                 
                 
                 
                 mxPath(from="xab1", to = "xba1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc1"),
                 mxPath(from="xac1", to = "xca1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc1"),
                 mxPath(from="xad1", to = "xda1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc1"),
                 mxPath(from="xae1", to = "xea1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc1"),
                 
                 
                 mxPath(from="xbc1", to = "xcb1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc1"),
                 mxPath(from="xbd1", to = "xdb1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc1"),
                 mxPath(from="xbe1", to = "xeb1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc1"),
                 
                 
                 mxPath(from="xcd1", to = "xdc1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc1"),
                 mxPath(from="xce1", to = "xec1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc1"),
                 
                 mxPath(from="xde1", to = "xed1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc1")
)


#########################################################################################################################


#base13

base13 <- mxModel("base13",
                  type="RAM",
                  manifestVars = c(
                    "xab1",
                    "xac1",
                    "xad1",
                    "xae1",
                    "xaf1",
                    "xag1",
                    "xah1",
                    "xai1",
                    "xaj1",
                    "xak1",
                    "xal1",
                    "xam1",
                    
                    
                    "xba1",
                    "xbc1",
                    "xbd1",
                    "xbe1",
                    "xbf1",
                    "xbg1",
                    "xbh1",
                    "xbi1",
                    "xbj1",
                    "xbk1",
                    "xbl1",
                    "xbm1",
                    
                    
                    "xca1",
                    "xcb1",
                    "xcd1",
                    "xce1",
                    "xcf1",
                    "xcg1",
                    "xch1",
                    "xci1",
                    "xcj1",
                    "xck1",
                    "xcl1",
                    "xcm1",
                    
                    
                    "xda1",
                    "xdb1",
                    "xdc1",
                    "xde1",
                    "xdf1",
                    "xdg1",
                    "xdh1",
                    "xdi1",
                    "xdj1",
                    "xdk1",
                    "xdl1",
                    "xdm1",
                    
                    
                    "xea1",
                    "xeb1",
                    "xec1",
                    "xed1",
                    "xef1",
                    "xeg1",
                    "xeh1",
                    "xei1",
                    "xej1",
                    "xek1",
                    "xel1",
                    "xem1",
                    
                    
                    "xfa1",
                    "xfb1",
                    "xfc1",
                    "xfd1",
                    "xfe1",
                    "xfg1",
                    "xfh1",
                    "xfi1",
                    "xfj1",
                    "xfk1",
                    "xfl1",
                    "xfm1",
                    
                    
                    "xga1",
                    "xgb1",
                    "xgc1",
                    "xgd1",
                    "xge1",
                    "xgf1",
                    "xgh1",
                    "xgi1",
                    "xgj1",
                    "xgk1",
                    "xgl1",
                    "xgm1",
                    
                    
                    "xha1",
                    "xhb1",
                    "xhc1",
                    "xhd1",
                    "xhe1",
                    "xhf1",
                    "xhg1",
                    "xhi1",
                    "xhj1",
                    "xhk1",
                    "xhl1",
                    "xhm1",
                    
                    
                    "xia1",
                    "xib1",
                    "xic1",
                    "xid1",
                    "xie1",
                    "xif1",
                    "xig1",
                    "xih1",
                    "xij1",
                    "xik1",
                    "xil1",
                    "xim1",
                    
                    
                    "xja1",
                    "xjb1",
                    "xjc1",
                    "xjd1",
                    "xje1",
                    "xjf1",
                    "xjg1",
                    "xjh1",
                    "xji1",
                    "xjk1",
                    "xjl1",
                    "xjm1",
                    
                    
                    "xka1",
                    "xkb1",
                    "xkc1",
                    "xkd1",
                    "xke1",
                    "xkf1",
                    "xkg1",
                    "xkh1",
                    "xki1",
                    "xkj1",
                    "xkl1",
                    "xkm1",
                    
                    
                    "xla1",
                    "xlb1",
                    "xlc1",
                    "xld1",
                    "xle1",
                    "xlf1",
                    "xlg1",
                    "xlh1",
                    "xli1",
                    "xlj1",
                    "xlk1",
                    "xlm1",
                    
                    "xma1",
                    "xmb1",
                    "xmc1",
                    "xmd1",
                    "xme1",
                    "xmf1",
                    "xmg1",
                    "xmh1",
                    "xmi1",
                    "xmj1",
                    "xmk1",
                    "xml1"),
                  
                  latentVars = c(
                    "Aon1"     ,      
                    "Bon1"     ,
                    "Con1"     ,
                    "Don1"     ,
                    "Eon1"     ,
                    "Fon1"     ,
                    "Gon1"     ,
                    "Hon1"     ,
                    "Ion1"     ,
                    "Jon1"     ,
                    "Kon1"     ,
                    "Lon1"     ,
                    "Mon1"     ,
                    
                    "onA1"     ,
                    "onB1"     ,
                    "onC1"     ,
                    "onD1"     ,
                    "onE1"     ,
                    "onF1"     ,
                    "onG1"     ,
                    "onH1"     ,
                    "onI1"     ,
                    "onJ1"     ,
                    "onK1"     ,
                    "onL1"     ,
                    "onM1"     ),
                  
                  
                  #person level
                  
                  mxPath(from="Aon1", to= "xab1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Aon1", to= "xac1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Aon1", to= "xad1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Aon1", to= "xae1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Aon1", to= "xaf1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Aon1", to= "xag1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Aon1", to= "xah1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Aon1", to= "xai1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Aon1", to= "xaj1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Aon1", to= "xak1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Aon1", to= "xal1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Aon1", to= "xam1", arrows = 1, free = F, values=1, labels = "la1"),
                  
                  mxPath(from="Bon1", to= "xba1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Bon1", to= "xbc1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Bon1", to= "xbd1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Bon1", to= "xbe1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Bon1", to= "xbf1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Bon1", to= "xbg1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Bon1", to= "xbh1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Bon1", to= "xbi1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Bon1", to= "xbj1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Bon1", to= "xbk1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Bon1", to= "xbl1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Bon1", to= "xbm1", arrows = 1, free = F, values=1, labels = "la1"),
                  
                  mxPath(from="Con1", to= "xca1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Con1", to= "xcb1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Con1", to= "xcd1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Con1", to= "xce1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Con1", to= "xcf1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Con1", to= "xcg1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Con1", to= "xch1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Con1", to= "xci1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Con1", to= "xcj1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Con1", to= "xck1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Con1", to= "xcl1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Con1", to= "xcm1", arrows = 1, free = F, values=1, labels = "la1"),
                  
                  mxPath(from="Don1", to= "xda1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Don1", to= "xdb1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Don1", to= "xdc1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Don1", to= "xde1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Don1", to= "xdf1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Don1", to= "xdg1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Don1", to= "xdh1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Don1", to= "xdi1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Don1", to= "xdj1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Don1", to= "xdk1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Don1", to= "xdl1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Don1", to= "xdm1", arrows = 1, free = F, values=1, labels = "la1"),
                  
                  mxPath(from="Eon1", to= "xea1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Eon1", to= "xeb1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Eon1", to= "xec1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Eon1", to= "xed1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Eon1", to= "xef1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Eon1", to= "xeg1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Eon1", to= "xeh1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Eon1", to= "xei1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Eon1", to= "xej1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Eon1", to= "xek1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Eon1", to= "xel1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Eon1", to= "xem1", arrows = 1, free = F, values=1, labels = "la1"),
                  
                  mxPath(from="Fon1", to= "xfa1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Fon1", to= "xfb1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Fon1", to= "xfc1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Fon1", to= "xfd1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Fon1", to= "xfe1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Fon1", to= "xfg1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Fon1", to= "xfh1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Fon1", to= "xfi1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Fon1", to= "xfj1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Fon1", to= "xfk1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Fon1", to= "xfl1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Fon1", to= "xfm1", arrows = 1, free = F, values=1, labels = "la1"),
                  
                  mxPath(from="Gon1", to= "xga1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Gon1", to= "xgb1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Gon1", to= "xgc1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Gon1", to= "xgd1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Gon1", to= "xge1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Gon1", to= "xgf1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Gon1", to= "xgh1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Gon1", to= "xgi1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Gon1", to= "xgj1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Gon1", to= "xgk1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Gon1", to= "xgl1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Gon1", to= "xgm1", arrows = 1, free = F, values=1, labels = "la1"),
                  
                  mxPath(from="Hon1", to= "xha1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Hon1", to= "xhb1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Hon1", to= "xhc1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Hon1", to= "xhd1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Hon1", to= "xhe1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Hon1", to= "xhf1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Hon1", to= "xhg1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Hon1", to= "xhi1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Hon1", to= "xhj1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Hon1", to= "xhk1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Hon1", to= "xhl1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Hon1", to= "xhm1", arrows = 1, free = F, values=1, labels = "la1"),
                  
                  mxPath(from="Ion1", to= "xia1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Ion1", to= "xib1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Ion1", to= "xic1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Ion1", to= "xid1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Ion1", to= "xie1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Ion1", to= "xif1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Ion1", to= "xig1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Ion1", to= "xih1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Ion1", to= "xij1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Ion1", to= "xik1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Ion1", to= "xil1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Ion1", to= "xim1", arrows = 1, free = F, values=1, labels = "la1"),
                  
                  mxPath(from="Jon1", to= "xja1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Jon1", to= "xjb1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Jon1", to= "xjc1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Jon1", to= "xjd1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Jon1", to= "xje1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Jon1", to= "xjf1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Jon1", to= "xjg1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Jon1", to= "xjh1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Jon1", to= "xji1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Jon1", to= "xjk1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Jon1", to= "xjl1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Jon1", to= "xjm1", arrows = 1, free = F, values=1, labels = "la1"),
                  
                  mxPath(from="Kon1", to= "xka1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Kon1", to= "xkb1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Kon1", to= "xkc1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Kon1", to= "xkd1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Kon1", to= "xke1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Kon1", to= "xkf1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Kon1", to= "xkg1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Kon1", to= "xkh1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Kon1", to= "xki1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Kon1", to= "xkj1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Kon1", to= "xkl1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Kon1", to= "xkm1", arrows = 1, free = F, values=1, labels = "la1"),
                  
                  mxPath(from="Lon1", to= "xla1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Lon1", to= "xlb1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Lon1", to= "xlc1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Lon1", to= "xld1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Lon1", to= "xle1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Lon1", to= "xlf1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Lon1", to= "xlg1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Lon1", to= "xlh1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Lon1", to= "xli1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Lon1", to= "xlj1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Lon1", to= "xlk1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Lon1", to= "xlm1", arrows = 1, free = F, values=1, labels = "la1"),
                  
                  mxPath(from="Mon1", to= "xma1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Mon1", to= "xmb1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Mon1", to= "xmc1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Mon1", to= "xmd1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Mon1", to= "xme1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Mon1", to= "xmf1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Mon1", to= "xmg1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Mon1", to= "xmh1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Mon1", to= "xmi1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Mon1", to= "xmj1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Mon1", to= "xmk1", arrows = 1, free = F, values=1, labels = "la1"),
                  mxPath(from="Mon1", to= "xml1", arrows = 1, free = F, values=1, labels = "la1"),
                  
                  mxPath(from="onA1", to= "xba1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onA1", to= "xca1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onA1", to= "xda1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onA1", to= "xea1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onA1", to= "xfa1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onA1", to= "xga1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onA1", to= "xha1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onA1", to= "xia1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onA1", to= "xja1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onA1", to= "xka1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onA1", to= "xla1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onA1", to= "xma1", arrows = 1, free = F, values=1, labels = "lp1"),
                  
                  mxPath(from="onB1", to= "xab1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onB1", to= "xcb1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onB1", to= "xdb1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onB1", to= "xeb1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onB1", to= "xfb1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onB1", to= "xgb1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onB1", to= "xhb1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onB1", to= "xib1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onB1", to= "xjb1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onB1", to= "xkb1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onB1", to= "xlb1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onB1", to= "xmb1", arrows = 1, free = F, values=1, labels = "lp1"),
                  
                  mxPath(from="onC1", to= "xac1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onC1", to= "xbc1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onC1", to= "xdc1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onC1", to= "xec1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onC1", to= "xfc1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onC1", to= "xgc1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onC1", to= "xhc1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onC1", to= "xic1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onC1", to= "xjc1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onC1", to= "xkc1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onC1", to= "xlc1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onC1", to= "xmc1", arrows = 1, free = F, values=1, labels = "lp1"),
                  
                  mxPath(from="onD1", to= "xad1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onD1", to= "xbd1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onD1", to= "xcd1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onD1", to= "xed1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onD1", to= "xfd1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onD1", to= "xgd1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onD1", to= "xhd1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onD1", to= "xid1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onD1", to= "xjd1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onD1", to= "xkd1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onD1", to= "xld1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onD1", to= "xmd1", arrows = 1, free = F, values=1, labels = "lp1"),
                  
                  mxPath(from="onE1", to= "xae1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onE1", to= "xbe1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onE1", to= "xce1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onE1", to= "xde1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onE1", to= "xfe1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onE1", to= "xge1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onE1", to= "xhe1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onE1", to= "xie1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onE1", to= "xje1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onE1", to= "xke1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onE1", to= "xle1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onE1", to= "xme1", arrows = 1, free = F, values=1, labels = "lp1"),
                  
                  mxPath(from="onF1", to= "xaf1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onF1", to= "xbf1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onF1", to= "xcf1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onF1", to= "xdf1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onF1", to= "xef1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onF1", to= "xgf1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onF1", to= "xhf1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onF1", to= "xif1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onF1", to= "xjf1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onF1", to= "xkf1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onF1", to= "xlf1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onF1", to= "xmf1", arrows = 1, free = F, values=1, labels = "lp1"),
                  
                  mxPath(from="onG1", to= "xag1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onG1", to= "xbg1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onG1", to= "xcg1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onG1", to= "xdg1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onG1", to= "xeg1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onG1", to= "xfg1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onG1", to= "xhg1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onG1", to= "xig1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onG1", to= "xjg1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onG1", to= "xkg1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onG1", to= "xlg1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onG1", to= "xmg1", arrows = 1, free = F, values=1, labels = "lp1"),
                  
                  mxPath(from="onH1", to= "xah1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onH1", to= "xbh1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onH1", to= "xch1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onH1", to= "xdh1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onH1", to= "xeh1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onH1", to= "xfh1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onH1", to= "xgh1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onH1", to= "xih1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onH1", to= "xjh1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onH1", to= "xkh1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onH1", to= "xlh1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onH1", to= "xmh1", arrows = 1, free = F, values=1, labels = "lp1"),
                  
                  mxPath(from="onI1", to= "xai1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onI1", to= "xbi1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onI1", to= "xci1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onI1", to= "xdi1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onI1", to= "xei1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onI1", to= "xfi1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onI1", to= "xgi1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onI1", to= "xhi1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onI1", to= "xji1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onI1", to= "xki1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onI1", to= "xli1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onI1", to= "xmi1", arrows = 1, free = F, values=1, labels = "lp1"),
                  
                  mxPath(from="onJ1", to= "xaj1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onJ1", to= "xbj1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onJ1", to= "xcj1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onJ1", to= "xdj1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onJ1", to= "xej1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onJ1", to= "xfj1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onJ1", to= "xgj1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onJ1", to= "xhj1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onJ1", to= "xij1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onJ1", to= "xkj1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onJ1", to= "xlj1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onJ1", to= "xmj1", arrows = 1, free = F, values=1, labels = "lp1"),
                  
                  mxPath(from="onK1", to= "xak1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onK1", to= "xbk1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onK1", to= "xck1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onK1", to= "xdk1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onK1", to= "xek1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onK1", to= "xfk1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onK1", to= "xgk1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onK1", to= "xhk1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onK1", to= "xik1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onK1", to= "xjk1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onK1", to= "xlk1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onK1", to= "xmk1", arrows = 1, free = F, values=1, labels = "lp1"),
                  
                  mxPath(from="onL1", to= "xal1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onL1", to= "xbl1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onL1", to= "xcl1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onL1", to= "xdl1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onL1", to= "xel1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onL1", to= "xfl1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onL1", to= "xgl1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onL1", to= "xhl1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onL1", to= "xil1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onL1", to= "xjl1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onL1", to= "xkl1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onL1", to= "xml1", arrows = 1, free = F, values=1, labels = "lp1"),
                  
                  mxPath(from="onM1", to= "xam1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onM1", to= "xbm1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onM1", to= "xcm1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onM1", to= "xdm1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onM1", to= "xem1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onM1", to= "xfm1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onM1", to= "xgm1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onM1", to= "xhm1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onM1", to= "xim1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onM1", to= "xjm1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onM1", to= "xkm1", arrows = 1, free = F, values=1, labels = "lp1"),
                  mxPath(from="onM1", to= "xlm1", arrows = 1, free = F, values=1, labels = "lp1"),
                  
                  
                  mxPath(from="Aon1", arrows=2, free=T,  values=generalActorVar, labels = "av1"),
                  mxPath(from="Bon1", arrows=2, free=T,  values=generalActorVar, labels = "av1"),
                  mxPath(from="Con1", arrows=2, free=T,  values=generalActorVar, labels = "av1"),
                  mxPath(from="Don1", arrows=2, free=T,  values=generalActorVar, labels = "av1"),
                  mxPath(from="Eon1", arrows=2, free=T,  values=generalActorVar, labels = "av1"),
                  mxPath(from="Fon1", arrows=2, free=T,  values=generalActorVar, labels = "av1"),
                  mxPath(from="Gon1", arrows=2, free=T,  values=generalActorVar, labels = "av1"),
                  mxPath(from="Hon1", arrows=2, free=T,  values=generalActorVar, labels = "av1"),
                  mxPath(from="Ion1", arrows=2, free=T,  values=generalActorVar, labels = "av1"),
                  mxPath(from="Jon1", arrows=2, free=T,  values=generalActorVar, labels = "av1"),
                  mxPath(from="Kon1", arrows=2, free=T,  values=generalActorVar, labels = "av1"),
                  mxPath(from="Lon1", arrows=2, free=T,  values=generalActorVar, labels = "av1"),
                  mxPath(from="Mon1", arrows=2, free=T,  values=generalActorVar, labels = "av1"),
                  
                  mxPath(from="onA1", arrows=2, free=T,  values=generalPartnerVar, labels = "pv1"),
                  mxPath(from="onB1", arrows=2, free=T,  values=generalPartnerVar, labels = "pv1"),
                  mxPath(from="onC1", arrows=2, free=T,  values=generalPartnerVar, labels = "pv1"),
                  mxPath(from="onD1", arrows=2, free=T,  values=generalPartnerVar, labels = "pv1"),
                  mxPath(from="onE1", arrows=2, free=T,  values=generalPartnerVar, labels = "pv1"),
                  mxPath(from="onF1", arrows=2, free=T,  values=generalPartnerVar, labels = "pv1"),
                  mxPath(from="onG1", arrows=2, free=T,  values=generalPartnerVar, labels = "pv1"),
                  mxPath(from="onH1", arrows=2, free=T,  values=generalPartnerVar, labels = "pv1"),
                  mxPath(from="onI1", arrows=2, free=T,  values=generalPartnerVar, labels = "pv1"),
                  mxPath(from="onJ1", arrows=2, free=T,  values=generalPartnerVar, labels = "pv1"),
                  mxPath(from="onK1", arrows=2, free=T,  values=generalPartnerVar, labels = "pv1"),
                  mxPath(from="onL1", arrows=2, free=T,  values=generalPartnerVar, labels = "pv1"),
                  mxPath(from="onM1", arrows=2, free=T,  values=generalPartnerVar, labels = "pv1"),
                  
                  mxPath(from="Aon1", to = "onA1", arrows=2, free=T,  values=generalActorPartnerCov, labels = "apc1"),
                  mxPath(from="Bon1", to = "onB1", arrows=2, free=T,  values=generalActorPartnerCov, labels = "apc1"),
                  mxPath(from="Con1", to = "onC1", arrows=2, free=T,  values=generalActorPartnerCov, labels = "apc1"),
                  mxPath(from="Don1", to = "onD1", arrows=2, free=T,  values=generalActorPartnerCov, labels = "apc1"),
                  mxPath(from="Eon1", to = "onE1", arrows=2, free=T,  values=generalActorPartnerCov, labels = "apc1"),
                  mxPath(from="Fon1", to = "onF1", arrows=2, free=T,  values=generalActorPartnerCov, labels = "apc1"),
                  mxPath(from="Gon1", to = "onG1", arrows=2, free=T,  values=generalActorPartnerCov, labels = "apc1"),
                  mxPath(from="Hon1", to = "onH1", arrows=2, free=T,  values=generalActorPartnerCov, labels = "apc1"),
                  mxPath(from="Ion1", to = "onI1", arrows=2, free=T,  values=generalActorPartnerCov, labels = "apc1"),
                  mxPath(from="Jon1", to = "onJ1", arrows=2, free=T,  values=generalActorPartnerCov, labels = "apc1"),
                  mxPath(from="Kon1", to = "onK1", arrows=2, free=T,  values=generalActorPartnerCov, labels = "apc1"),
                  mxPath(from="Lon1", to = "onL1", arrows=2, free=T,  values=generalActorPartnerCov, labels = "apc1"),
                  mxPath(from="Mon1", to = "onM1", arrows=2, free=T,  values=generalActorPartnerCov, labels = "apc1"),
                  
                  mxPath(from="one", to = "onA1", arrows=1, free=F,  values=0, labels = "mean1"),
                  mxPath(from="one", to = "onB1", arrows=1, free=F,  values=0, labels = "mean1"),
                  mxPath(from="one", to = "onC1", arrows=1, free=F,  values=0, labels = "mean1"),
                  mxPath(from="one", to = "onD1", arrows=1, free=F,  values=0, labels = "mean1"),
                  mxPath(from="one", to = "onE1", arrows=1, free=F,  values=0, labels = "mean1"),
                  mxPath(from="one", to = "onF1", arrows=1, free=F,  values=0, labels = "mean1"),
                  mxPath(from="one", to = "onG1", arrows=1, free=F,  values=0, labels = "mean1"),
                  mxPath(from="one", to = "onH1", arrows=1, free=F,  values=0, labels = "mean1"),
                  mxPath(from="one", to = "onI1", arrows=1, free=F,  values=0, labels = "mean1"),
                  mxPath(from="one", to = "onJ1", arrows=1, free=F,  values=0, labels = "mean1"),
                  mxPath(from="one", to = "onK1", arrows=1, free=F,  values=0, labels = "mean1"),
                  mxPath(from="one", to = "onL1", arrows=1, free=F,  values=0, labels = "mean1"),
                  mxPath(from="one", to = "onM1", arrows=1, free=F,  values=0, labels = "mean1"),
                  
                  
                  #dyad level
                  
                  mxPath(from="xab1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xac1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xad1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xae1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xaf1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xag1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xah1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xai1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xaj1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xak1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xal1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xam1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  
                  mxPath(from="xba1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xbc1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xbd1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xbe1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xbf1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xbg1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xbh1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xbi1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xbj1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xbk1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xbl1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xbm1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  
                  mxPath(from="xca1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xcb1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xcd1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xce1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xcf1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xcg1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xch1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xci1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xcj1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xck1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xcl1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xcm1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  
                  mxPath(from="xda1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xdb1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xdc1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xde1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xdf1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xdg1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xdh1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xdi1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xdj1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xdk1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xdl1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xdm1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  
                  mxPath(from="xea1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xeb1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xec1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xed1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xef1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xeg1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xeh1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xei1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xej1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xek1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xel1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xem1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  
                  mxPath(from="xfa1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xfb1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xfc1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xfd1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xfe1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xfg1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xfh1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xfi1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xfj1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xfk1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xfl1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xfm1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  
                  mxPath(from="xga1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xgb1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xgc1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xgd1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xge1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xgf1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xgh1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xgi1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xgj1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xgk1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xgl1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xgm1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  
                  mxPath(from="xha1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xhb1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xhc1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xhd1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xhe1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xhf1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xhg1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xhi1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xhj1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xhk1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xhl1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xhm1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  
                  mxPath(from="xia1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xib1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xic1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xid1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xie1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xif1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xig1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xih1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xij1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xik1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xil1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xim1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  
                  mxPath(from="xja1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xjb1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xjc1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xjd1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xje1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xjf1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xjg1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xjh1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xji1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xjk1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xjl1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xjm1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  
                  mxPath(from="xka1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xkb1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xkc1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xkd1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xke1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xkf1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xkg1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xkh1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xki1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xkj1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xkl1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xkm1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  
                  mxPath(from="xla1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xlb1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xlc1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xld1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xle1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xlf1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xlg1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xlh1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xli1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xlj1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xlk1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xlm1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  
                  mxPath(from="xma1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xmb1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xmc1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xmd1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xme1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xmf1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xmg1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xmh1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xmi1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xmj1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xmk1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  mxPath(from="xml1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv1"),
                  
                  
                  mxPath(from="xab1", to = "xba1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc1"),
                  mxPath(from="xac1", to = "xca1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc1"),
                  mxPath(from="xad1", to = "xda1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc1"),
                  mxPath(from="xae1", to = "xea1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc1"),
                  mxPath(from="xaf1", to = "xfa1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc1"),
                  mxPath(from="xag1", to = "xga1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc1"),
                  mxPath(from="xah1", to = "xha1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc1"),
                  mxPath(from="xai1", to = "xia1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc1"),
                  mxPath(from="xaj1", to = "xja1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc1"),
                  mxPath(from="xak1", to = "xka1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc1"),
                  mxPath(from="xal1", to = "xla1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc1"),
                  mxPath(from="xam1", to = "xma1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc1"),
                  
                  mxPath(from="xbc1", to = "xcb1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc1"),
                  mxPath(from="xbd1", to = "xdb1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc1"),
                  mxPath(from="xbe1", to = "xeb1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc1"),
                  mxPath(from="xbf1", to = "xfb1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc1"),
                  mxPath(from="xbg1", to = "xgb1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc1"),
                  mxPath(from="xbh1", to = "xhb1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc1"),
                  mxPath(from="xbi1", to = "xib1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc1"),
                  mxPath(from="xbj1", to = "xjb1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc1"),
                  mxPath(from="xbk1", to = "xkb1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc1"),
                  mxPath(from="xbl1", to = "xlb1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc1"),
                  mxPath(from="xbm1", to = "xmb1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc1"),
                  
                  mxPath(from="xcd1", to = "xdc1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc1"),
                  mxPath(from="xce1", to = "xec1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc1"),
                  mxPath(from="xcf1", to = "xfc1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc1"),
                  mxPath(from="xcg1", to = "xgc1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc1"),
                  mxPath(from="xch1", to = "xhc1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc1"),
                  mxPath(from="xci1", to = "xic1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc1"),
                  mxPath(from="xcj1", to = "xjc1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc1"),
                  mxPath(from="xck1", to = "xkc1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc1"),
                  mxPath(from="xcl1", to = "xlc1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc1"),
                  mxPath(from="xcm1", to = "xmc1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc1"),
                  
                  mxPath(from="xde1", to = "xed1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc1"),
                  mxPath(from="xdf1", to = "xfd1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc1"),
                  mxPath(from="xdg1", to = "xgd1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc1"),
                  mxPath(from="xdh1", to = "xhd1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc1"),
                  mxPath(from="xdi1", to = "xid1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc1"),
                  mxPath(from="xdj1", to = "xjd1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc1"),
                  mxPath(from="xdk1", to = "xkd1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc1"),
                  mxPath(from="xdl1", to = "xld1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc1"),
                  mxPath(from="xdm1", to = "xmd1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc1"),
                  
                  mxPath(from="xef1", to = "xfe1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc1"),
                  mxPath(from="xeg1", to = "xge1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc1"),
                  mxPath(from="xeh1", to = "xhe1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc1"),
                  mxPath(from="xei1", to = "xie1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc1"),
                  mxPath(from="xej1", to = "xje1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc1"),
                  mxPath(from="xek1", to = "xke1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc1"),
                  mxPath(from="xel1", to = "xle1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc1"),
                  mxPath(from="xem1", to = "xme1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc1"),
                  
                  mxPath(from="xfg1", to = "xgf1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc1"),
                  mxPath(from="xfh1", to = "xhf1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc1"),
                  mxPath(from="xfi1", to = "xif1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc1"),
                  mxPath(from="xfj1", to = "xjf1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc1"),
                  mxPath(from="xfk1", to = "xkf1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc1"),
                  mxPath(from="xfl1", to = "xlf1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc1"),
                  mxPath(from="xfm1", to = "xmf1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc1"),
                  
                  mxPath(from="xgh1", to = "xhg1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc1"),
                  mxPath(from="xgi1", to = "xig1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc1"),
                  mxPath(from="xgj1", to = "xjg1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc1"),
                  mxPath(from="xgk1", to = "xkg1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc1"),
                  mxPath(from="xgl1", to = "xlg1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc1"),
                  mxPath(from="xgm1", to = "xmg1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc1"),
                  
                  mxPath(from="xhi1", to = "xih1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc1"),
                  mxPath(from="xhj1", to = "xjh1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc1"),
                  mxPath(from="xhk1", to = "xkh1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc1"),
                  mxPath(from="xhl1", to = "xlh1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc1"),
                  mxPath(from="xhm1", to = "xmh1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc1"),
                  
                  mxPath(from="xij1", to = "xji1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc1"),
                  mxPath(from="xik1", to = "xki1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc1"),
                  mxPath(from="xil1", to = "xli1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc1"),
                  mxPath(from="xim1", to = "xmi1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc1"),
                  
                  mxPath(from="xjk1", to = "xkj1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc1"),
                  mxPath(from="xjl1", to = "xlj1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc1"),
                  mxPath(from="xjm1", to = "xmj1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc1"),
                  
                  mxPath(from="xkl1", to = "xlk1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc1"),
                  mxPath(from="xkm1", to = "xmk1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc1"),
                  
                  mxPath(from="xlm1", to = "xml1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc1")
)



##############################################################################################################################


#other5same

other5same <- mxModel("other5same",
                      type="RAM",
                      manifestVars = c(
                        "xab1",
                        "xac1",
                        "xad1",
                        "xae1",
                        
                        "xba1",
                        "xbc1",
                        "xbd1",
                        "xbe1",
                        
                        "xca1",
                        "xcb1",
                        "xcd1",
                        "xce1",
                        
                        "xda1",
                        "xdb1",
                        "xdc1",
                        "xde1",
                        
                        "xea1",
                        "xeb1",
                        "xec1",
                        "xed1"),
                      
                      latentVars = c(
                        "Aon1"     ,      
                        "Bon1"     ,
                        "Con1"     ,
                        "Don1"     ,
                        "Eon1"     ,
                        
                        "onA1"     ,
                        "onB1"     ,
                        "onC1"     ,
                        "onD1"     ,
                        "onE1"),
                      
                      #person level
                      
                      mxPath(from="Aon1", to= "xab1", arrows = 1, free = F, values=1, labels = "la2"),
                      mxPath(from="Aon1", to= "xac1", arrows = 1, free = F, values=1, labels = "la2"),
                      mxPath(from="Aon1", to= "xad1", arrows = 1, free = F, values=1, labels = "la2"),
                      mxPath(from="Aon1", to= "xae1", arrows = 1, free = F, values=1, labels = "la2"),
                      
                      
                      mxPath(from="Bon1", to= "xba1", arrows = 1, free = F, values=1, labels = "la2"),
                      mxPath(from="Bon1", to= "xbc1", arrows = 1, free = F, values=1, labels = "la2"),
                      mxPath(from="Bon1", to= "xbd1", arrows = 1, free = F, values=1, labels = "la2"),
                      mxPath(from="Bon1", to= "xbe1", arrows = 1, free = F, values=1, labels = "la2"),
                      
                      
                      mxPath(from="Con1", to= "xca1", arrows = 1, free = F, values=1, labels = "la2"),
                      mxPath(from="Con1", to= "xcb1", arrows = 1, free = F, values=1, labels = "la2"),
                      mxPath(from="Con1", to= "xcd1", arrows = 1, free = F, values=1, labels = "la2"),
                      mxPath(from="Con1", to= "xce1", arrows = 1, free = F, values=1, labels = "la2"),
                      
                      
                      mxPath(from="Don1", to= "xda1", arrows = 1, free = F, values=1, labels = "la2"),
                      mxPath(from="Don1", to= "xdb1", arrows = 1, free = F, values=1, labels = "la2"),
                      mxPath(from="Don1", to= "xdc1", arrows = 1, free = F, values=1, labels = "la2"),
                      mxPath(from="Don1", to= "xde1", arrows = 1, free = F, values=1, labels = "la2"),
                      
                      
                      mxPath(from="Eon1", to= "xea1", arrows = 1, free = F, values=1, labels = "la2"),
                      mxPath(from="Eon1", to= "xeb1", arrows = 1, free = F, values=1, labels = "la2"),
                      mxPath(from="Eon1", to= "xec1", arrows = 1, free = F, values=1, labels = "la2"),
                      mxPath(from="Eon1", to= "xed1", arrows = 1, free = F, values=1, labels = "la2"),
                      
                      
                      
                      mxPath(from="onA1", to= "xba1", arrows = 1, free = F, values=1, labels = "lp2"),
                      mxPath(from="onA1", to= "xca1", arrows = 1, free = F, values=1, labels = "lp2"),
                      mxPath(from="onA1", to= "xda1", arrows = 1, free = F, values=1, labels = "lp2"),
                      mxPath(from="onA1", to= "xea1", arrows = 1, free = F, values=1, labels = "lp2"),
                      
                      
                      mxPath(from="onB1", to= "xab1", arrows = 1, free = F, values=1, labels = "lp2"),
                      mxPath(from="onB1", to= "xcb1", arrows = 1, free = F, values=1, labels = "lp2"),
                      mxPath(from="onB1", to= "xdb1", arrows = 1, free = F, values=1, labels = "lp2"),
                      mxPath(from="onB1", to= "xeb1", arrows = 1, free = F, values=1, labels = "lp2"),
                      
                      
                      mxPath(from="onC1", to= "xac1", arrows = 1, free = F, values=1, labels = "lp2"),
                      mxPath(from="onC1", to= "xbc1", arrows = 1, free = F, values=1, labels = "lp2"),
                      mxPath(from="onC1", to= "xdc1", arrows = 1, free = F, values=1, labels = "lp2"),
                      mxPath(from="onC1", to= "xec1", arrows = 1, free = F, values=1, labels = "lp2"),
                      
                      
                      mxPath(from="onD1", to= "xad1", arrows = 1, free = F, values=1, labels = "lp2"),
                      mxPath(from="onD1", to= "xbd1", arrows = 1, free = F, values=1, labels = "lp2"),
                      mxPath(from="onD1", to= "xcd1", arrows = 1, free = F, values=1, labels = "lp2"),
                      mxPath(from="onD1", to= "xed1", arrows = 1, free = F, values=1, labels = "lp2"),
                      
                      
                      mxPath(from="onE1", to= "xae1", arrows = 1, free = F, values=1, labels = "lp2"),
                      mxPath(from="onE1", to= "xbe1", arrows = 1, free = F, values=1, labels = "lp2"),
                      mxPath(from="onE1", to= "xce1", arrows = 1, free = F, values=1, labels = "lp2"),
                      mxPath(from="onE1", to= "xde1", arrows = 1, free = F, values=1, labels = "lp2"),
                      
                      
                      
                      mxPath(from="Aon1", arrows=2, free=T,  values=generalActorVar, labels = "av2"),
                      mxPath(from="Bon1", arrows=2, free=T,  values=generalActorVar, labels = "av2"),
                      mxPath(from="Con1", arrows=2, free=T,  values=generalActorVar, labels = "av2"),
                      mxPath(from="Don1", arrows=2, free=T,  values=generalActorVar, labels = "av2"),
                      mxPath(from="Eon1", arrows=2, free=T,  values=generalActorVar, labels = "av2"),
                      
                      
                      mxPath(from="onA1", arrows=2, free=T,  values=generalPartnerVar, labels = "pv2"),
                      mxPath(from="onB1", arrows=2, free=T,  values=generalPartnerVar, labels = "pv2"),
                      mxPath(from="onC1", arrows=2, free=T,  values=generalPartnerVar, labels = "pv2"),
                      mxPath(from="onD1", arrows=2, free=T,  values=generalPartnerVar, labels = "pv2"),
                      mxPath(from="onE1", arrows=2, free=T,  values=generalPartnerVar, labels = "pv2"),
                      
                      
                      mxPath(from="Aon1", to = "onA1", arrows=2, free=T,  values=generalActorPartnerCov, labels = "apc2"),
                      mxPath(from="Bon1", to = "onB1", arrows=2, free=T,  values=generalActorPartnerCov, labels = "apc2"),
                      mxPath(from="Con1", to = "onC1", arrows=2, free=T,  values=generalActorPartnerCov, labels = "apc2"),
                      mxPath(from="Don1", to = "onD1", arrows=2, free=T,  values=generalActorPartnerCov, labels = "apc2"),
                      mxPath(from="Eon1", to = "onE1", arrows=2, free=T,  values=generalActorPartnerCov, labels = "apc2"),
                      
                      mxPath(from="one", to = "onA1", arrows=1, free=F,  values=0, labels = "mean2"),
                      mxPath(from="one", to = "onB1", arrows=1, free=F,  values=0, labels = "mean2"),
                      mxPath(from="one", to = "onC1", arrows=1, free=F,  values=0, labels = "mean2"),
                      mxPath(from="one", to = "onD1", arrows=1, free=F,  values=0, labels = "mean2"),
                      mxPath(from="one", to = "onE1", arrows=1, free=F,  values=0, labels = "mean2"),
                      
                      
                      
                      #dyad level
                      
                      mxPath(from="xab1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                      mxPath(from="xac1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                      mxPath(from="xad1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                      mxPath(from="xae1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                      
                      
                      mxPath(from="xba1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                      mxPath(from="xbc1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                      mxPath(from="xbd1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                      mxPath(from="xbe1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                      
                      
                      mxPath(from="xca1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                      mxPath(from="xcb1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                      mxPath(from="xcd1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                      mxPath(from="xce1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                      
                      
                      mxPath(from="xda1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                      mxPath(from="xdb1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                      mxPath(from="xdc1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                      mxPath(from="xde1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                      
                      
                      mxPath(from="xea1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                      mxPath(from="xeb1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                      mxPath(from="xec1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                      mxPath(from="xed1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                      
                      
                      
                      mxPath(from="xab1", to = "xba1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc2"),
                      mxPath(from="xac1", to = "xca1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc2"),
                      mxPath(from="xad1", to = "xda1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc2"),
                      mxPath(from="xae1", to = "xea1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc2"),
                      
                      
                      mxPath(from="xbc1", to = "xcb1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc2"),
                      mxPath(from="xbd1", to = "xdb1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc2"),
                      mxPath(from="xbe1", to = "xeb1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc2"),
                      
                      
                      mxPath(from="xcd1", to = "xdc1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc2"),
                      mxPath(from="xce1", to = "xec1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc2"),
                      
                      mxPath(from="xde1", to = "xed1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc2")
)


##########################################################################################################################


#other5different

other5different <- mxModel("other5different",
                           type="RAM",
                           manifestVars = c(
                             "xab1",
                             "xac1",
                             "xad1",
                             "xae1",
                             
                             "xba1",
                             "xbc1",
                             "xbd1",
                             "xbe1",
                             
                             "xca1",
                             "xcb1",
                             "xcd1",
                             "xce1",
                             
                             "xda1",
                             "xdb1",
                             "xdc1",
                             "xde1",
                             
                             "xea1",
                             "xeb1",
                             "xec1",
                             "xed1"),
                           
                           latentVars = c(
                             "Aon1"     ,      
                             "Bon1"     ,
                             "Con1"     ,
                             "Don1"     ,
                             "Eon1"     ,
                             
                             "onA1"     ,
                             "onB1"     ,
                             "onC1"     ,
                             "onD1"     ,
                             "onE1"),
                           
                           #person level
                           
                           mxPath(from="Aon1", to= "xab1", arrows = 1, free = F, values=1, labels = "la2"),
                           mxPath(from="Aon1", to= "xac1", arrows = 1, free = F, values=1, labels = "la2"),
                           mxPath(from="Aon1", to= "xad1", arrows = 1, free = F, values=1, labels = "la2"),
                           mxPath(from="Aon1", to= "xae1", arrows = 1, free = F, values=1, labels = "la2"),
                           
                           
                           mxPath(from="Bon1", to= "xba1", arrows = 1, free = F, values=1, labels = "la2"),
                           mxPath(from="Bon1", to= "xbc1", arrows = 1, free = F, values=1, labels = "la2"),
                           mxPath(from="Bon1", to= "xbd1", arrows = 1, free = F, values=1, labels = "la2"),
                           mxPath(from="Bon1", to= "xbe1", arrows = 1, free = F, values=1, labels = "la2"),
                           
                           
                           mxPath(from="Con1", to= "xca1", arrows = 1, free = F, values=1, labels = "la2"),
                           mxPath(from="Con1", to= "xcb1", arrows = 1, free = F, values=1, labels = "la2"),
                           mxPath(from="Con1", to= "xcd1", arrows = 1, free = F, values=1, labels = "la2"),
                           mxPath(from="Con1", to= "xce1", arrows = 1, free = F, values=1, labels = "la2"),
                           
                           
                           mxPath(from="Don1", to= "xda1", arrows = 1, free = F, values=1, labels = "la2"),
                           mxPath(from="Don1", to= "xdb1", arrows = 1, free = F, values=1, labels = "la2"),
                           mxPath(from="Don1", to= "xdc1", arrows = 1, free = F, values=1, labels = "la2"),
                           mxPath(from="Don1", to= "xde1", arrows = 1, free = F, values=1, labels = "la2"),
                           
                           
                           mxPath(from="Eon1", to= "xea1", arrows = 1, free = F, values=1, labels = "la2"),
                           mxPath(from="Eon1", to= "xeb1", arrows = 1, free = F, values=1, labels = "la2"),
                           mxPath(from="Eon1", to= "xec1", arrows = 1, free = F, values=1, labels = "la2"),
                           mxPath(from="Eon1", to= "xed1", arrows = 1, free = F, values=1, labels = "la2"),
                           
                           
                           
                           mxPath(from="onA1", to= "xba1", arrows = 1, free = F, values=1, labels = "lp2"),
                           mxPath(from="onA1", to= "xca1", arrows = 1, free = F, values=1, labels = "lp2"),
                           mxPath(from="onA1", to= "xda1", arrows = 1, free = F, values=1, labels = "lp2"),
                           mxPath(from="onA1", to= "xea1", arrows = 1, free = F, values=1, labels = "lp2"),
                           
                           
                           mxPath(from="onB1", to= "xab1", arrows = 1, free = F, values=1, labels = "lp2"),
                           mxPath(from="onB1", to= "xcb1", arrows = 1, free = F, values=1, labels = "lp2"),
                           mxPath(from="onB1", to= "xdb1", arrows = 1, free = F, values=1, labels = "lp2"),
                           mxPath(from="onB1", to= "xeb1", arrows = 1, free = F, values=1, labels = "lp2"),
                           
                           
                           mxPath(from="onC1", to= "xac1", arrows = 1, free = F, values=1, labels = "lp2"),
                           mxPath(from="onC1", to= "xbc1", arrows = 1, free = F, values=1, labels = "lp2"),
                           mxPath(from="onC1", to= "xdc1", arrows = 1, free = F, values=1, labels = "lp2"),
                           mxPath(from="onC1", to= "xec1", arrows = 1, free = F, values=1, labels = "lp2"),
                           
                           
                           mxPath(from="onD1", to= "xad1", arrows = 1, free = F, values=1, labels = "lp2"),
                           mxPath(from="onD1", to= "xbd1", arrows = 1, free = F, values=1, labels = "lp2"),
                           mxPath(from="onD1", to= "xcd1", arrows = 1, free = F, values=1, labels = "lp2"),
                           mxPath(from="onD1", to= "xed1", arrows = 1, free = F, values=1, labels = "lp2"),
                           
                           
                           mxPath(from="onE1", to= "xae1", arrows = 1, free = F, values=1, labels = "lp2"),
                           mxPath(from="onE1", to= "xbe1", arrows = 1, free = F, values=1, labels = "lp2"),
                           mxPath(from="onE1", to= "xce1", arrows = 1, free = F, values=1, labels = "lp2"),
                           mxPath(from="onE1", to= "xde1", arrows = 1, free = F, values=1, labels = "lp2"),
                           
                           
                           
                           mxPath(from="Aon1", arrows=2, free=T,  values=differentActorVar, labels = "av2"),
                           mxPath(from="Bon1", arrows=2, free=T,  values=differentActorVar, labels = "av2"),
                           mxPath(from="Con1", arrows=2, free=T,  values=differentActorVar, labels = "av2"),
                           mxPath(from="Don1", arrows=2, free=T,  values=differentActorVar, labels = "av2"),
                           mxPath(from="Eon1", arrows=2, free=T,  values=differentActorVar, labels = "av2"),
                           
                           
                           mxPath(from="onA1", arrows=2, free=T,  values=differentPartnerVar, labels = "pv2"),
                           mxPath(from="onB1", arrows=2, free=T,  values=differentPartnerVar, labels = "pv2"),
                           mxPath(from="onC1", arrows=2, free=T,  values=differentPartnerVar, labels = "pv2"),
                           mxPath(from="onD1", arrows=2, free=T,  values=differentPartnerVar, labels = "pv2"),
                           mxPath(from="onE1", arrows=2, free=T,  values=differentPartnerVar, labels = "pv2"),
                           
                           
                           mxPath(from="Aon1", to = "onA1", arrows=2, free=T,  values=differentActorPartnerCov, labels = "apc2"),
                           mxPath(from="Bon1", to = "onB1", arrows=2, free=T,  values=differentActorPartnerCov, labels = "apc2"),
                           mxPath(from="Con1", to = "onC1", arrows=2, free=T,  values=differentActorPartnerCov, labels = "apc2"),
                           mxPath(from="Don1", to = "onD1", arrows=2, free=T,  values=differentActorPartnerCov, labels = "apc2"),
                           mxPath(from="Eon1", to = "onE1", arrows=2, free=T,  values=differentActorPartnerCov, labels = "apc2"),
                           
                           mxPath(from="one", to = "onA1", arrows=1, free=F,  values=0, labels = "mean2"),
                           mxPath(from="one", to = "onB1", arrows=1, free=F,  values=0, labels = "mean2"),
                           mxPath(from="one", to = "onC1", arrows=1, free=F,  values=0, labels = "mean2"),
                           mxPath(from="one", to = "onD1", arrows=1, free=F,  values=0, labels = "mean2"),
                           mxPath(from="one", to = "onE1", arrows=1, free=F,  values=0, labels = "mean2"),
                           
                           
                           
                           #dyad level
                           
                           mxPath(from="xab1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                           mxPath(from="xac1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                           mxPath(from="xad1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                           mxPath(from="xae1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                           
                           
                           mxPath(from="xba1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                           mxPath(from="xbc1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                           mxPath(from="xbd1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                           mxPath(from="xbe1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                           
                           
                           mxPath(from="xca1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                           mxPath(from="xcb1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                           mxPath(from="xcd1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                           mxPath(from="xce1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                           
                           
                           mxPath(from="xda1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                           mxPath(from="xdb1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                           mxPath(from="xdc1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                           mxPath(from="xde1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                           
                           
                           mxPath(from="xea1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                           mxPath(from="xeb1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                           mxPath(from="xec1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                           mxPath(from="xed1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                           
                           
                           
                           mxPath(from="xab1", to = "xba1", arrows=2, free=T,  values=differentRelationshipCov, labels = "rc2"),
                           mxPath(from="xac1", to = "xca1", arrows=2, free=T,  values=differentRelationshipCov, labels = "rc2"),
                           mxPath(from="xad1", to = "xda1", arrows=2, free=T,  values=differentRelationshipCov, labels = "rc2"),
                           mxPath(from="xae1", to = "xea1", arrows=2, free=T,  values=differentRelationshipCov, labels = "rc2"),
                           
                           
                           mxPath(from="xbc1", to = "xcb1", arrows=2, free=T,  values=differentRelationshipCov, labels = "rc2"),
                           mxPath(from="xbd1", to = "xdb1", arrows=2, free=T,  values=differentRelationshipCov, labels = "rc2"),
                           mxPath(from="xbe1", to = "xeb1", arrows=2, free=T,  values=differentRelationshipCov, labels = "rc2"),
                           
                           
                           mxPath(from="xcd1", to = "xdc1", arrows=2, free=T,  values=differentRelationshipCov, labels = "rc2"),
                           mxPath(from="xce1", to = "xec1", arrows=2, free=T,  values=differentRelationshipCov, labels = "rc2"),
                           
                           mxPath(from="xde1", to = "xed1", arrows=2, free=T,  values=differentRelationshipCov, labels = "rc2")
)


################################################################################################################################


#other13same

other13same <- mxModel("other13same",
                       type="RAM",
                       manifestVars = c(
                         "xab1",
                         "xac1",
                         "xad1",
                         "xae1",
                         "xaf1",
                         "xag1",
                         "xah1",
                         "xai1",
                         "xaj1",
                         "xak1",
                         "xal1",
                         "xam1",
                         
                         
                         "xba1",
                         "xbc1",
                         "xbd1",
                         "xbe1",
                         "xbf1",
                         "xbg1",
                         "xbh1",
                         "xbi1",
                         "xbj1",
                         "xbk1",
                         "xbl1",
                         "xbm1",
                         
                         
                         "xca1",
                         "xcb1",
                         "xcd1",
                         "xce1",
                         "xcf1",
                         "xcg1",
                         "xch1",
                         "xci1",
                         "xcj1",
                         "xck1",
                         "xcl1",
                         "xcm1",
                         
                         
                         "xda1",
                         "xdb1",
                         "xdc1",
                         "xde1",
                         "xdf1",
                         "xdg1",
                         "xdh1",
                         "xdi1",
                         "xdj1",
                         "xdk1",
                         "xdl1",
                         "xdm1",
                         
                         
                         "xea1",
                         "xeb1",
                         "xec1",
                         "xed1",
                         "xef1",
                         "xeg1",
                         "xeh1",
                         "xei1",
                         "xej1",
                         "xek1",
                         "xel1",
                         "xem1",
                         
                         
                         "xfa1",
                         "xfb1",
                         "xfc1",
                         "xfd1",
                         "xfe1",
                         "xfg1",
                         "xfh1",
                         "xfi1",
                         "xfj1",
                         "xfk1",
                         "xfl1",
                         "xfm1",
                         
                         
                         "xga1",
                         "xgb1",
                         "xgc1",
                         "xgd1",
                         "xge1",
                         "xgf1",
                         "xgh1",
                         "xgi1",
                         "xgj1",
                         "xgk1",
                         "xgl1",
                         "xgm1",
                         
                         
                         "xha1",
                         "xhb1",
                         "xhc1",
                         "xhd1",
                         "xhe1",
                         "xhf1",
                         "xhg1",
                         "xhi1",
                         "xhj1",
                         "xhk1",
                         "xhl1",
                         "xhm1",
                         
                         
                         "xia1",
                         "xib1",
                         "xic1",
                         "xid1",
                         "xie1",
                         "xif1",
                         "xig1",
                         "xih1",
                         "xij1",
                         "xik1",
                         "xil1",
                         "xim1",
                         
                         
                         "xja1",
                         "xjb1",
                         "xjc1",
                         "xjd1",
                         "xje1",
                         "xjf1",
                         "xjg1",
                         "xjh1",
                         "xji1",
                         "xjk1",
                         "xjl1",
                         "xjm1",
                         
                         
                         "xka1",
                         "xkb1",
                         "xkc1",
                         "xkd1",
                         "xke1",
                         "xkf1",
                         "xkg1",
                         "xkh1",
                         "xki1",
                         "xkj1",
                         "xkl1",
                         "xkm1",
                         
                         
                         "xla1",
                         "xlb1",
                         "xlc1",
                         "xld1",
                         "xle1",
                         "xlf1",
                         "xlg1",
                         "xlh1",
                         "xli1",
                         "xlj1",
                         "xlk1",
                         "xlm1",
                         
                         "xma1",
                         "xmb1",
                         "xmc1",
                         "xmd1",
                         "xme1",
                         "xmf1",
                         "xmg1",
                         "xmh1",
                         "xmi1",
                         "xmj1",
                         "xmk1",
                         "xml1"),
                       
                       latentVars = c(
                         "Aon1"     ,      
                         "Bon1"     ,
                         "Con1"     ,
                         "Don1"     ,
                         "Eon1"     ,
                         "Fon1"     ,
                         "Gon1"     ,
                         "Hon1"     ,
                         "Ion1"     ,
                         "Jon1"     ,
                         "Kon1"     ,
                         "Lon1"     ,
                         "Mon1"     ,
                         
                         "onA1"     ,
                         "onB1"     ,
                         "onC1"     ,
                         "onD1"     ,
                         "onE1"     ,
                         "onF1"     ,
                         "onG1"     ,
                         "onH1"     ,
                         "onI1"     ,
                         "onJ1"     ,
                         "onK1"     ,
                         "onL1"     ,
                         "onM1"     ),
                       
                       
                       #person level
                       
                       mxPath(from="Aon1", to= "xab1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Aon1", to= "xac1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Aon1", to= "xad1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Aon1", to= "xae1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Aon1", to= "xaf1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Aon1", to= "xag1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Aon1", to= "xah1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Aon1", to= "xai1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Aon1", to= "xaj1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Aon1", to= "xak1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Aon1", to= "xal1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Aon1", to= "xam1", arrows = 1, free = F, values=1, labels = "la2"),
                       
                       mxPath(from="Bon1", to= "xba1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Bon1", to= "xbc1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Bon1", to= "xbd1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Bon1", to= "xbe1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Bon1", to= "xbf1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Bon1", to= "xbg1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Bon1", to= "xbh1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Bon1", to= "xbi1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Bon1", to= "xbj1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Bon1", to= "xbk1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Bon1", to= "xbl1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Bon1", to= "xbm1", arrows = 1, free = F, values=1, labels = "la2"),
                       
                       mxPath(from="Con1", to= "xca1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Con1", to= "xcb1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Con1", to= "xcd1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Con1", to= "xce1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Con1", to= "xcf1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Con1", to= "xcg1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Con1", to= "xch1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Con1", to= "xci1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Con1", to= "xcj1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Con1", to= "xck1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Con1", to= "xcl1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Con1", to= "xcm1", arrows = 1, free = F, values=1, labels = "la2"),
                       
                       mxPath(from="Don1", to= "xda1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Don1", to= "xdb1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Don1", to= "xdc1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Don1", to= "xde1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Don1", to= "xdf1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Don1", to= "xdg1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Don1", to= "xdh1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Don1", to= "xdi1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Don1", to= "xdj1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Don1", to= "xdk1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Don1", to= "xdl1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Don1", to= "xdm1", arrows = 1, free = F, values=1, labels = "la2"),
                       
                       mxPath(from="Eon1", to= "xea1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Eon1", to= "xeb1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Eon1", to= "xec1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Eon1", to= "xed1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Eon1", to= "xef1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Eon1", to= "xeg1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Eon1", to= "xeh1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Eon1", to= "xei1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Eon1", to= "xej1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Eon1", to= "xek1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Eon1", to= "xel1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Eon1", to= "xem1", arrows = 1, free = F, values=1, labels = "la2"),
                       
                       mxPath(from="Fon1", to= "xfa1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Fon1", to= "xfb1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Fon1", to= "xfc1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Fon1", to= "xfd1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Fon1", to= "xfe1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Fon1", to= "xfg1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Fon1", to= "xfh1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Fon1", to= "xfi1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Fon1", to= "xfj1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Fon1", to= "xfk1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Fon1", to= "xfl1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Fon1", to= "xfm1", arrows = 1, free = F, values=1, labels = "la2"),
                       
                       mxPath(from="Gon1", to= "xga1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Gon1", to= "xgb1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Gon1", to= "xgc1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Gon1", to= "xgd1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Gon1", to= "xge1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Gon1", to= "xgf1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Gon1", to= "xgh1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Gon1", to= "xgi1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Gon1", to= "xgj1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Gon1", to= "xgk1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Gon1", to= "xgl1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Gon1", to= "xgm1", arrows = 1, free = F, values=1, labels = "la2"),
                       
                       mxPath(from="Hon1", to= "xha1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Hon1", to= "xhb1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Hon1", to= "xhc1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Hon1", to= "xhd1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Hon1", to= "xhe1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Hon1", to= "xhf1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Hon1", to= "xhg1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Hon1", to= "xhi1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Hon1", to= "xhj1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Hon1", to= "xhk1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Hon1", to= "xhl1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Hon1", to= "xhm1", arrows = 1, free = F, values=1, labels = "la2"),
                       
                       mxPath(from="Ion1", to= "xia1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Ion1", to= "xib1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Ion1", to= "xic1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Ion1", to= "xid1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Ion1", to= "xie1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Ion1", to= "xif1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Ion1", to= "xig1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Ion1", to= "xih1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Ion1", to= "xij1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Ion1", to= "xik1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Ion1", to= "xil1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Ion1", to= "xim1", arrows = 1, free = F, values=1, labels = "la2"),
                       
                       mxPath(from="Jon1", to= "xja1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Jon1", to= "xjb1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Jon1", to= "xjc1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Jon1", to= "xjd1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Jon1", to= "xje1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Jon1", to= "xjf1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Jon1", to= "xjg1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Jon1", to= "xjh1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Jon1", to= "xji1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Jon1", to= "xjk1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Jon1", to= "xjl1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Jon1", to= "xjm1", arrows = 1, free = F, values=1, labels = "la2"),
                       
                       mxPath(from="Kon1", to= "xka1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Kon1", to= "xkb1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Kon1", to= "xkc1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Kon1", to= "xkd1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Kon1", to= "xke1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Kon1", to= "xkf1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Kon1", to= "xkg1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Kon1", to= "xkh1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Kon1", to= "xki1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Kon1", to= "xkj1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Kon1", to= "xkl1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Kon1", to= "xkm1", arrows = 1, free = F, values=1, labels = "la2"),
                       
                       mxPath(from="Lon1", to= "xla1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Lon1", to= "xlb1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Lon1", to= "xlc1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Lon1", to= "xld1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Lon1", to= "xle1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Lon1", to= "xlf1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Lon1", to= "xlg1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Lon1", to= "xlh1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Lon1", to= "xli1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Lon1", to= "xlj1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Lon1", to= "xlk1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Lon1", to= "xlm1", arrows = 1, free = F, values=1, labels = "la2"),
                       
                       mxPath(from="Mon1", to= "xma1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Mon1", to= "xmb1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Mon1", to= "xmc1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Mon1", to= "xmd1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Mon1", to= "xme1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Mon1", to= "xmf1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Mon1", to= "xmg1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Mon1", to= "xmh1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Mon1", to= "xmi1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Mon1", to= "xmj1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Mon1", to= "xmk1", arrows = 1, free = F, values=1, labels = "la2"),
                       mxPath(from="Mon1", to= "xml1", arrows = 1, free = F, values=1, labels = "la2"),
                       
                       mxPath(from="onA1", to= "xba1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onA1", to= "xca1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onA1", to= "xda1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onA1", to= "xea1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onA1", to= "xfa1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onA1", to= "xga1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onA1", to= "xha1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onA1", to= "xia1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onA1", to= "xja1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onA1", to= "xka1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onA1", to= "xla1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onA1", to= "xma1", arrows = 1, free = F, values=1, labels = "lp2"),
                       
                       mxPath(from="onB1", to= "xab1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onB1", to= "xcb1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onB1", to= "xdb1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onB1", to= "xeb1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onB1", to= "xfb1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onB1", to= "xgb1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onB1", to= "xhb1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onB1", to= "xib1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onB1", to= "xjb1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onB1", to= "xkb1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onB1", to= "xlb1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onB1", to= "xmb1", arrows = 1, free = F, values=1, labels = "lp2"),
                       
                       mxPath(from="onC1", to= "xac1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onC1", to= "xbc1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onC1", to= "xdc1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onC1", to= "xec1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onC1", to= "xfc1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onC1", to= "xgc1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onC1", to= "xhc1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onC1", to= "xic1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onC1", to= "xjc1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onC1", to= "xkc1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onC1", to= "xlc1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onC1", to= "xmc1", arrows = 1, free = F, values=1, labels = "lp2"),
                       
                       mxPath(from="onD1", to= "xad1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onD1", to= "xbd1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onD1", to= "xcd1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onD1", to= "xed1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onD1", to= "xfd1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onD1", to= "xgd1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onD1", to= "xhd1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onD1", to= "xid1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onD1", to= "xjd1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onD1", to= "xkd1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onD1", to= "xld1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onD1", to= "xmd1", arrows = 1, free = F, values=1, labels = "lp2"),
                       
                       mxPath(from="onE1", to= "xae1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onE1", to= "xbe1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onE1", to= "xce1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onE1", to= "xde1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onE1", to= "xfe1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onE1", to= "xge1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onE1", to= "xhe1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onE1", to= "xie1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onE1", to= "xje1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onE1", to= "xke1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onE1", to= "xle1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onE1", to= "xme1", arrows = 1, free = F, values=1, labels = "lp2"),
                       
                       mxPath(from="onF1", to= "xaf1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onF1", to= "xbf1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onF1", to= "xcf1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onF1", to= "xdf1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onF1", to= "xef1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onF1", to= "xgf1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onF1", to= "xhf1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onF1", to= "xif1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onF1", to= "xjf1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onF1", to= "xkf1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onF1", to= "xlf1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onF1", to= "xmf1", arrows = 1, free = F, values=1, labels = "lp2"),
                       
                       mxPath(from="onG1", to= "xag1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onG1", to= "xbg1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onG1", to= "xcg1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onG1", to= "xdg1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onG1", to= "xeg1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onG1", to= "xfg1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onG1", to= "xhg1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onG1", to= "xig1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onG1", to= "xjg1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onG1", to= "xkg1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onG1", to= "xlg1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onG1", to= "xmg1", arrows = 1, free = F, values=1, labels = "lp2"),
                       
                       mxPath(from="onH1", to= "xah1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onH1", to= "xbh1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onH1", to= "xch1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onH1", to= "xdh1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onH1", to= "xeh1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onH1", to= "xfh1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onH1", to= "xgh1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onH1", to= "xih1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onH1", to= "xjh1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onH1", to= "xkh1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onH1", to= "xlh1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onH1", to= "xmh1", arrows = 1, free = F, values=1, labels = "lp2"),
                       
                       mxPath(from="onI1", to= "xai1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onI1", to= "xbi1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onI1", to= "xci1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onI1", to= "xdi1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onI1", to= "xei1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onI1", to= "xfi1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onI1", to= "xgi1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onI1", to= "xhi1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onI1", to= "xji1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onI1", to= "xki1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onI1", to= "xli1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onI1", to= "xmi1", arrows = 1, free = F, values=1, labels = "lp2"),
                       
                       mxPath(from="onJ1", to= "xaj1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onJ1", to= "xbj1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onJ1", to= "xcj1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onJ1", to= "xdj1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onJ1", to= "xej1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onJ1", to= "xfj1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onJ1", to= "xgj1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onJ1", to= "xhj1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onJ1", to= "xij1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onJ1", to= "xkj1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onJ1", to= "xlj1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onJ1", to= "xmj1", arrows = 1, free = F, values=1, labels = "lp2"),
                       
                       mxPath(from="onK1", to= "xak1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onK1", to= "xbk1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onK1", to= "xck1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onK1", to= "xdk1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onK1", to= "xek1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onK1", to= "xfk1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onK1", to= "xgk1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onK1", to= "xhk1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onK1", to= "xik1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onK1", to= "xjk1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onK1", to= "xlk1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onK1", to= "xmk1", arrows = 1, free = F, values=1, labels = "lp2"),
                       
                       mxPath(from="onL1", to= "xal1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onL1", to= "xbl1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onL1", to= "xcl1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onL1", to= "xdl1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onL1", to= "xel1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onL1", to= "xfl1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onL1", to= "xgl1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onL1", to= "xhl1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onL1", to= "xil1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onL1", to= "xjl1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onL1", to= "xkl1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onL1", to= "xml1", arrows = 1, free = F, values=1, labels = "lp2"),
                       
                       mxPath(from="onM1", to= "xam1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onM1", to= "xbm1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onM1", to= "xcm1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onM1", to= "xdm1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onM1", to= "xem1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onM1", to= "xfm1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onM1", to= "xgm1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onM1", to= "xhm1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onM1", to= "xim1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onM1", to= "xjm1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onM1", to= "xkm1", arrows = 1, free = F, values=1, labels = "lp2"),
                       mxPath(from="onM1", to= "xlm1", arrows = 1, free = F, values=1, labels = "lp2"),
                       
                       
                       mxPath(from="Aon1", arrows=2, free=T,  values=generalActorVar, labels = "av2"),
                       mxPath(from="Bon1", arrows=2, free=T,  values=generalActorVar, labels = "av2"),
                       mxPath(from="Con1", arrows=2, free=T,  values=generalActorVar, labels = "av2"),
                       mxPath(from="Don1", arrows=2, free=T,  values=generalActorVar, labels = "av2"),
                       mxPath(from="Eon1", arrows=2, free=T,  values=generalActorVar, labels = "av2"),
                       mxPath(from="Fon1", arrows=2, free=T,  values=generalActorVar, labels = "av2"),
                       mxPath(from="Gon1", arrows=2, free=T,  values=generalActorVar, labels = "av2"),
                       mxPath(from="Hon1", arrows=2, free=T,  values=generalActorVar, labels = "av2"),
                       mxPath(from="Ion1", arrows=2, free=T,  values=generalActorVar, labels = "av2"),
                       mxPath(from="Jon1", arrows=2, free=T,  values=generalActorVar, labels = "av2"),
                       mxPath(from="Kon1", arrows=2, free=T,  values=generalActorVar, labels = "av2"),
                       mxPath(from="Lon1", arrows=2, free=T,  values=generalActorVar, labels = "av2"),
                       mxPath(from="Mon1", arrows=2, free=T,  values=generalActorVar, labels = "av2"),
                       
                       mxPath(from="onA1", arrows=2, free=T,  values=generalPartnerVar, labels = "pv2"),
                       mxPath(from="onB1", arrows=2, free=T,  values=generalPartnerVar, labels = "pv2"),
                       mxPath(from="onC1", arrows=2, free=T,  values=generalPartnerVar, labels = "pv2"),
                       mxPath(from="onD1", arrows=2, free=T,  values=generalPartnerVar, labels = "pv2"),
                       mxPath(from="onE1", arrows=2, free=T,  values=generalPartnerVar, labels = "pv2"),
                       mxPath(from="onF1", arrows=2, free=T,  values=generalPartnerVar, labels = "pv2"),
                       mxPath(from="onG1", arrows=2, free=T,  values=generalPartnerVar, labels = "pv2"),
                       mxPath(from="onH1", arrows=2, free=T,  values=generalPartnerVar, labels = "pv2"),
                       mxPath(from="onI1", arrows=2, free=T,  values=generalPartnerVar, labels = "pv2"),
                       mxPath(from="onJ1", arrows=2, free=T,  values=generalPartnerVar, labels = "pv2"),
                       mxPath(from="onK1", arrows=2, free=T,  values=generalPartnerVar, labels = "pv2"),
                       mxPath(from="onL1", arrows=2, free=T,  values=generalPartnerVar, labels = "pv2"),
                       mxPath(from="onM1", arrows=2, free=T,  values=generalPartnerVar, labels = "pv2"),
                       
                       mxPath(from="Aon1", to = "onA1", arrows=2, free=T,  values=generalActorPartnerCov, labels = "apc2"),
                       mxPath(from="Bon1", to = "onB1", arrows=2, free=T,  values=generalActorPartnerCov, labels = "apc2"),
                       mxPath(from="Con1", to = "onC1", arrows=2, free=T,  values=generalActorPartnerCov, labels = "apc2"),
                       mxPath(from="Don1", to = "onD1", arrows=2, free=T,  values=generalActorPartnerCov, labels = "apc2"),
                       mxPath(from="Eon1", to = "onE1", arrows=2, free=T,  values=generalActorPartnerCov, labels = "apc2"),
                       mxPath(from="Fon1", to = "onF1", arrows=2, free=T,  values=generalActorPartnerCov, labels = "apc2"),
                       mxPath(from="Gon1", to = "onG1", arrows=2, free=T,  values=generalActorPartnerCov, labels = "apc2"),
                       mxPath(from="Hon1", to = "onH1", arrows=2, free=T,  values=generalActorPartnerCov, labels = "apc2"),
                       mxPath(from="Ion1", to = "onI1", arrows=2, free=T,  values=generalActorPartnerCov, labels = "apc2"),
                       mxPath(from="Jon1", to = "onJ1", arrows=2, free=T,  values=generalActorPartnerCov, labels = "apc2"),
                       mxPath(from="Kon1", to = "onK1", arrows=2, free=T,  values=generalActorPartnerCov, labels = "apc2"),
                       mxPath(from="Lon1", to = "onL1", arrows=2, free=T,  values=generalActorPartnerCov, labels = "apc2"),
                       mxPath(from="Mon1", to = "onM1", arrows=2, free=T,  values=generalActorPartnerCov, labels = "apc2"),
                       
                       mxPath(from="one", to = "onA1", arrows=1, free=F,  values=0, labels = "mean2"),
                       mxPath(from="one", to = "onB1", arrows=1, free=F,  values=0, labels = "mean2"),
                       mxPath(from="one", to = "onC1", arrows=1, free=F,  values=0, labels = "mean2"),
                       mxPath(from="one", to = "onD1", arrows=1, free=F,  values=0, labels = "mean2"),
                       mxPath(from="one", to = "onE1", arrows=1, free=F,  values=0, labels = "mean2"),
                       mxPath(from="one", to = "onF1", arrows=1, free=F,  values=0, labels = "mean2"),
                       mxPath(from="one", to = "onG1", arrows=1, free=F,  values=0, labels = "mean2"),
                       mxPath(from="one", to = "onH1", arrows=1, free=F,  values=0, labels = "mean2"),
                       mxPath(from="one", to = "onI1", arrows=1, free=F,  values=0, labels = "mean2"),
                       mxPath(from="one", to = "onJ1", arrows=1, free=F,  values=0, labels = "mean2"),
                       mxPath(from="one", to = "onK1", arrows=1, free=F,  values=0, labels = "mean2"),
                       mxPath(from="one", to = "onL1", arrows=1, free=F,  values=0, labels = "mean2"),
                       mxPath(from="one", to = "onM1", arrows=1, free=F,  values=0, labels = "mean2"),
                       
                       
                       #dyad level
                       
                       mxPath(from="xab1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xac1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xad1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xae1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xaf1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xag1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xah1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xai1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xaj1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xak1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xal1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xam1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       
                       mxPath(from="xba1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xbc1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xbd1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xbe1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xbf1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xbg1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xbh1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xbi1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xbj1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xbk1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xbl1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xbm1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       
                       mxPath(from="xca1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xcb1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xcd1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xce1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xcf1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xcg1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xch1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xci1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xcj1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xck1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xcl1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xcm1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       
                       mxPath(from="xda1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xdb1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xdc1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xde1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xdf1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xdg1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xdh1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xdi1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xdj1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xdk1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xdl1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xdm1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       
                       mxPath(from="xea1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xeb1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xec1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xed1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xef1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xeg1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xeh1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xei1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xej1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xek1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xel1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xem1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       
                       mxPath(from="xfa1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xfb1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xfc1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xfd1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xfe1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xfg1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xfh1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xfi1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xfj1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xfk1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xfl1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xfm1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       
                       mxPath(from="xga1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xgb1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xgc1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xgd1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xge1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xgf1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xgh1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xgi1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xgj1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xgk1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xgl1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xgm1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       
                       mxPath(from="xha1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xhb1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xhc1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xhd1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xhe1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xhf1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xhg1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xhi1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xhj1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xhk1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xhl1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xhm1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       
                       mxPath(from="xia1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xib1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xic1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xid1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xie1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xif1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xig1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xih1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xij1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xik1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xil1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xim1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       
                       mxPath(from="xja1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xjb1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xjc1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xjd1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xje1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xjf1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xjg1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xjh1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xji1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xjk1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xjl1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xjm1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       
                       mxPath(from="xka1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xkb1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xkc1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xkd1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xke1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xkf1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xkg1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xkh1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xki1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xkj1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xkl1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xkm1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       
                       mxPath(from="xla1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xlb1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xlc1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xld1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xle1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xlf1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xlg1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xlh1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xli1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xlj1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xlk1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xlm1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       
                       mxPath(from="xma1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xmb1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xmc1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xmd1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xme1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xmf1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xmg1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xmh1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xmi1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xmj1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xmk1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       mxPath(from="xml1", arrows=2, free=T,  values=generalRelationshipVar, labels = "rv2"),
                       
                       
                       mxPath(from="xab1", to = "xba1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc2"),
                       mxPath(from="xac1", to = "xca1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc2"),
                       mxPath(from="xad1", to = "xda1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc2"),
                       mxPath(from="xae1", to = "xea1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc2"),
                       mxPath(from="xaf1", to = "xfa1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc2"),
                       mxPath(from="xag1", to = "xga1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc2"),
                       mxPath(from="xah1", to = "xha1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc2"),
                       mxPath(from="xai1", to = "xia1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc2"),
                       mxPath(from="xaj1", to = "xja1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc2"),
                       mxPath(from="xak1", to = "xka1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc2"),
                       mxPath(from="xal1", to = "xla1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc2"),
                       mxPath(from="xam1", to = "xma1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc2"),
                       
                       mxPath(from="xbc1", to = "xcb1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc2"),
                       mxPath(from="xbd1", to = "xdb1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc2"),
                       mxPath(from="xbe1", to = "xeb1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc2"),
                       mxPath(from="xbf1", to = "xfb1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc2"),
                       mxPath(from="xbg1", to = "xgb1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc2"),
                       mxPath(from="xbh1", to = "xhb1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc2"),
                       mxPath(from="xbi1", to = "xib1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc2"),
                       mxPath(from="xbj1", to = "xjb1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc2"),
                       mxPath(from="xbk1", to = "xkb1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc2"),
                       mxPath(from="xbl1", to = "xlb1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc2"),
                       mxPath(from="xbm1", to = "xmb1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc2"),
                       
                       mxPath(from="xcd1", to = "xdc1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc2"),
                       mxPath(from="xce1", to = "xec1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc2"),
                       mxPath(from="xcf1", to = "xfc1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc2"),
                       mxPath(from="xcg1", to = "xgc1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc2"),
                       mxPath(from="xch1", to = "xhc1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc2"),
                       mxPath(from="xci1", to = "xic1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc2"),
                       mxPath(from="xcj1", to = "xjc1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc2"),
                       mxPath(from="xck1", to = "xkc1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc2"),
                       mxPath(from="xcl1", to = "xlc1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc2"),
                       mxPath(from="xcm1", to = "xmc1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc2"),
                       
                       mxPath(from="xde1", to = "xed1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc2"),
                       mxPath(from="xdf1", to = "xfd1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc2"),
                       mxPath(from="xdg1", to = "xgd1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc2"),
                       mxPath(from="xdh1", to = "xhd1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc2"),
                       mxPath(from="xdi1", to = "xid1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc2"),
                       mxPath(from="xdj1", to = "xjd1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc2"),
                       mxPath(from="xdk1", to = "xkd1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc2"),
                       mxPath(from="xdl1", to = "xld1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc2"),
                       mxPath(from="xdm1", to = "xmd1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc2"),
                       
                       mxPath(from="xef1", to = "xfe1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc2"),
                       mxPath(from="xeg1", to = "xge1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc2"),
                       mxPath(from="xeh1", to = "xhe1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc2"),
                       mxPath(from="xei1", to = "xie1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc2"),
                       mxPath(from="xej1", to = "xje1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc2"),
                       mxPath(from="xek1", to = "xke1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc2"),
                       mxPath(from="xel1", to = "xle1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc2"),
                       mxPath(from="xem1", to = "xme1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc2"),
                       
                       mxPath(from="xfg1", to = "xgf1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc2"),
                       mxPath(from="xfh1", to = "xhf1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc2"),
                       mxPath(from="xfi1", to = "xif1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc2"),
                       mxPath(from="xfj1", to = "xjf1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc2"),
                       mxPath(from="xfk1", to = "xkf1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc2"),
                       mxPath(from="xfl1", to = "xlf1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc2"),
                       mxPath(from="xfm1", to = "xmf1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc2"),
                       
                       mxPath(from="xgh1", to = "xhg1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc2"),
                       mxPath(from="xgi1", to = "xig1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc2"),
                       mxPath(from="xgj1", to = "xjg1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc2"),
                       mxPath(from="xgk1", to = "xkg1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc2"),
                       mxPath(from="xgl1", to = "xlg1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc2"),
                       mxPath(from="xgm1", to = "xmg1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc2"),
                       
                       mxPath(from="xhi1", to = "xih1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc2"),
                       mxPath(from="xhj1", to = "xjh1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc2"),
                       mxPath(from="xhk1", to = "xkh1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc2"),
                       mxPath(from="xhl1", to = "xlh1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc2"),
                       mxPath(from="xhm1", to = "xmh1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc2"),
                       
                       mxPath(from="xij1", to = "xji1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc2"),
                       mxPath(from="xik1", to = "xki1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc2"),
                       mxPath(from="xil1", to = "xli1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc2"),
                       mxPath(from="xim1", to = "xmi1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc2"),
                       
                       mxPath(from="xjk1", to = "xkj1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc2"),
                       mxPath(from="xjl1", to = "xlj1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc2"),
                       mxPath(from="xjm1", to = "xmj1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc2"),
                       
                       mxPath(from="xkl1", to = "xlk1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc2"),
                       mxPath(from="xkm1", to = "xmk1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc2"),
                       
                       mxPath(from="xlm1", to = "xml1", arrows=2, free=T,  values=generalRelationshipCov, labels = "rc2")
)



#############################################################################################################################


#other13different

other13different <- mxModel("other13different",
                            type="RAM",
                            manifestVars = c(
                              "xab1",
                              "xac1",
                              "xad1",
                              "xae1",
                              "xaf1",
                              "xag1",
                              "xah1",
                              "xai1",
                              "xaj1",
                              "xak1",
                              "xal1",
                              "xam1",
                              
                              
                              "xba1",
                              "xbc1",
                              "xbd1",
                              "xbe1",
                              "xbf1",
                              "xbg1",
                              "xbh1",
                              "xbi1",
                              "xbj1",
                              "xbk1",
                              "xbl1",
                              "xbm1",
                              
                              
                              "xca1",
                              "xcb1",
                              "xcd1",
                              "xce1",
                              "xcf1",
                              "xcg1",
                              "xch1",
                              "xci1",
                              "xcj1",
                              "xck1",
                              "xcl1",
                              "xcm1",
                              
                              
                              "xda1",
                              "xdb1",
                              "xdc1",
                              "xde1",
                              "xdf1",
                              "xdg1",
                              "xdh1",
                              "xdi1",
                              "xdj1",
                              "xdk1",
                              "xdl1",
                              "xdm1",
                              
                              
                              "xea1",
                              "xeb1",
                              "xec1",
                              "xed1",
                              "xef1",
                              "xeg1",
                              "xeh1",
                              "xei1",
                              "xej1",
                              "xek1",
                              "xel1",
                              "xem1",
                              
                              
                              "xfa1",
                              "xfb1",
                              "xfc1",
                              "xfd1",
                              "xfe1",
                              "xfg1",
                              "xfh1",
                              "xfi1",
                              "xfj1",
                              "xfk1",
                              "xfl1",
                              "xfm1",
                              
                              
                              "xga1",
                              "xgb1",
                              "xgc1",
                              "xgd1",
                              "xge1",
                              "xgf1",
                              "xgh1",
                              "xgi1",
                              "xgj1",
                              "xgk1",
                              "xgl1",
                              "xgm1",
                              
                              
                              "xha1",
                              "xhb1",
                              "xhc1",
                              "xhd1",
                              "xhe1",
                              "xhf1",
                              "xhg1",
                              "xhi1",
                              "xhj1",
                              "xhk1",
                              "xhl1",
                              "xhm1",
                              
                              
                              "xia1",
                              "xib1",
                              "xic1",
                              "xid1",
                              "xie1",
                              "xif1",
                              "xig1",
                              "xih1",
                              "xij1",
                              "xik1",
                              "xil1",
                              "xim1",
                              
                              
                              "xja1",
                              "xjb1",
                              "xjc1",
                              "xjd1",
                              "xje1",
                              "xjf1",
                              "xjg1",
                              "xjh1",
                              "xji1",
                              "xjk1",
                              "xjl1",
                              "xjm1",
                              
                              
                              "xka1",
                              "xkb1",
                              "xkc1",
                              "xkd1",
                              "xke1",
                              "xkf1",
                              "xkg1",
                              "xkh1",
                              "xki1",
                              "xkj1",
                              "xkl1",
                              "xkm1",
                              
                              
                              "xla1",
                              "xlb1",
                              "xlc1",
                              "xld1",
                              "xle1",
                              "xlf1",
                              "xlg1",
                              "xlh1",
                              "xli1",
                              "xlj1",
                              "xlk1",
                              "xlm1",
                              
                              "xma1",
                              "xmb1",
                              "xmc1",
                              "xmd1",
                              "xme1",
                              "xmf1",
                              "xmg1",
                              "xmh1",
                              "xmi1",
                              "xmj1",
                              "xmk1",
                              "xml1"),
                            
                            latentVars = c(
                              "Aon1"     ,      
                              "Bon1"     ,
                              "Con1"     ,
                              "Don1"     ,
                              "Eon1"     ,
                              "Fon1"     ,
                              "Gon1"     ,
                              "Hon1"     ,
                              "Ion1"     ,
                              "Jon1"     ,
                              "Kon1"     ,
                              "Lon1"     ,
                              "Mon1"     ,
                              
                              "onA1"     ,
                              "onB1"     ,
                              "onC1"     ,
                              "onD1"     ,
                              "onE1"     ,
                              "onF1"     ,
                              "onG1"     ,
                              "onH1"     ,
                              "onI1"     ,
                              "onJ1"     ,
                              "onK1"     ,
                              "onL1"     ,
                              "onM1"     ),
                            
                            
                            #person level
                            
                            mxPath(from="Aon1", to= "xab1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Aon1", to= "xac1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Aon1", to= "xad1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Aon1", to= "xae1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Aon1", to= "xaf1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Aon1", to= "xag1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Aon1", to= "xah1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Aon1", to= "xai1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Aon1", to= "xaj1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Aon1", to= "xak1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Aon1", to= "xal1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Aon1", to= "xam1", arrows = 1, free = F, values=1, labels = "la2"),
                            
                            mxPath(from="Bon1", to= "xba1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Bon1", to= "xbc1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Bon1", to= "xbd1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Bon1", to= "xbe1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Bon1", to= "xbf1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Bon1", to= "xbg1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Bon1", to= "xbh1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Bon1", to= "xbi1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Bon1", to= "xbj1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Bon1", to= "xbk1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Bon1", to= "xbl1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Bon1", to= "xbm1", arrows = 1, free = F, values=1, labels = "la2"),
                            
                            mxPath(from="Con1", to= "xca1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Con1", to= "xcb1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Con1", to= "xcd1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Con1", to= "xce1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Con1", to= "xcf1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Con1", to= "xcg1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Con1", to= "xch1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Con1", to= "xci1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Con1", to= "xcj1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Con1", to= "xck1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Con1", to= "xcl1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Con1", to= "xcm1", arrows = 1, free = F, values=1, labels = "la2"),
                            
                            mxPath(from="Don1", to= "xda1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Don1", to= "xdb1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Don1", to= "xdc1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Don1", to= "xde1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Don1", to= "xdf1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Don1", to= "xdg1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Don1", to= "xdh1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Don1", to= "xdi1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Don1", to= "xdj1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Don1", to= "xdk1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Don1", to= "xdl1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Don1", to= "xdm1", arrows = 1, free = F, values=1, labels = "la2"),
                            
                            mxPath(from="Eon1", to= "xea1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Eon1", to= "xeb1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Eon1", to= "xec1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Eon1", to= "xed1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Eon1", to= "xef1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Eon1", to= "xeg1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Eon1", to= "xeh1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Eon1", to= "xei1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Eon1", to= "xej1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Eon1", to= "xek1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Eon1", to= "xel1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Eon1", to= "xem1", arrows = 1, free = F, values=1, labels = "la2"),
                            
                            mxPath(from="Fon1", to= "xfa1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Fon1", to= "xfb1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Fon1", to= "xfc1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Fon1", to= "xfd1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Fon1", to= "xfe1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Fon1", to= "xfg1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Fon1", to= "xfh1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Fon1", to= "xfi1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Fon1", to= "xfj1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Fon1", to= "xfk1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Fon1", to= "xfl1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Fon1", to= "xfm1", arrows = 1, free = F, values=1, labels = "la2"),
                            
                            mxPath(from="Gon1", to= "xga1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Gon1", to= "xgb1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Gon1", to= "xgc1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Gon1", to= "xgd1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Gon1", to= "xge1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Gon1", to= "xgf1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Gon1", to= "xgh1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Gon1", to= "xgi1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Gon1", to= "xgj1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Gon1", to= "xgk1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Gon1", to= "xgl1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Gon1", to= "xgm1", arrows = 1, free = F, values=1, labels = "la2"),
                            
                            mxPath(from="Hon1", to= "xha1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Hon1", to= "xhb1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Hon1", to= "xhc1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Hon1", to= "xhd1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Hon1", to= "xhe1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Hon1", to= "xhf1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Hon1", to= "xhg1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Hon1", to= "xhi1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Hon1", to= "xhj1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Hon1", to= "xhk1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Hon1", to= "xhl1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Hon1", to= "xhm1", arrows = 1, free = F, values=1, labels = "la2"),
                            
                            mxPath(from="Ion1", to= "xia1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Ion1", to= "xib1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Ion1", to= "xic1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Ion1", to= "xid1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Ion1", to= "xie1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Ion1", to= "xif1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Ion1", to= "xig1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Ion1", to= "xih1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Ion1", to= "xij1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Ion1", to= "xik1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Ion1", to= "xil1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Ion1", to= "xim1", arrows = 1, free = F, values=1, labels = "la2"),
                            
                            mxPath(from="Jon1", to= "xja1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Jon1", to= "xjb1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Jon1", to= "xjc1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Jon1", to= "xjd1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Jon1", to= "xje1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Jon1", to= "xjf1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Jon1", to= "xjg1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Jon1", to= "xjh1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Jon1", to= "xji1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Jon1", to= "xjk1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Jon1", to= "xjl1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Jon1", to= "xjm1", arrows = 1, free = F, values=1, labels = "la2"),
                            
                            mxPath(from="Kon1", to= "xka1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Kon1", to= "xkb1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Kon1", to= "xkc1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Kon1", to= "xkd1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Kon1", to= "xke1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Kon1", to= "xkf1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Kon1", to= "xkg1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Kon1", to= "xkh1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Kon1", to= "xki1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Kon1", to= "xkj1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Kon1", to= "xkl1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Kon1", to= "xkm1", arrows = 1, free = F, values=1, labels = "la2"),
                            
                            mxPath(from="Lon1", to= "xla1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Lon1", to= "xlb1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Lon1", to= "xlc1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Lon1", to= "xld1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Lon1", to= "xle1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Lon1", to= "xlf1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Lon1", to= "xlg1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Lon1", to= "xlh1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Lon1", to= "xli1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Lon1", to= "xlj1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Lon1", to= "xlk1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Lon1", to= "xlm1", arrows = 1, free = F, values=1, labels = "la2"),
                            
                            mxPath(from="Mon1", to= "xma1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Mon1", to= "xmb1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Mon1", to= "xmc1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Mon1", to= "xmd1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Mon1", to= "xme1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Mon1", to= "xmf1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Mon1", to= "xmg1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Mon1", to= "xmh1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Mon1", to= "xmi1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Mon1", to= "xmj1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Mon1", to= "xmk1", arrows = 1, free = F, values=1, labels = "la2"),
                            mxPath(from="Mon1", to= "xml1", arrows = 1, free = F, values=1, labels = "la2"),
                            
                            mxPath(from="onA1", to= "xba1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onA1", to= "xca1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onA1", to= "xda1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onA1", to= "xea1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onA1", to= "xfa1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onA1", to= "xga1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onA1", to= "xha1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onA1", to= "xia1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onA1", to= "xja1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onA1", to= "xka1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onA1", to= "xla1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onA1", to= "xma1", arrows = 1, free = F, values=1, labels = "lp2"),
                            
                            mxPath(from="onB1", to= "xab1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onB1", to= "xcb1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onB1", to= "xdb1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onB1", to= "xeb1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onB1", to= "xfb1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onB1", to= "xgb1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onB1", to= "xhb1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onB1", to= "xib1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onB1", to= "xjb1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onB1", to= "xkb1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onB1", to= "xlb1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onB1", to= "xmb1", arrows = 1, free = F, values=1, labels = "lp2"),
                            
                            mxPath(from="onC1", to= "xac1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onC1", to= "xbc1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onC1", to= "xdc1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onC1", to= "xec1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onC1", to= "xfc1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onC1", to= "xgc1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onC1", to= "xhc1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onC1", to= "xic1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onC1", to= "xjc1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onC1", to= "xkc1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onC1", to= "xlc1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onC1", to= "xmc1", arrows = 1, free = F, values=1, labels = "lp2"),
                            
                            mxPath(from="onD1", to= "xad1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onD1", to= "xbd1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onD1", to= "xcd1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onD1", to= "xed1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onD1", to= "xfd1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onD1", to= "xgd1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onD1", to= "xhd1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onD1", to= "xid1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onD1", to= "xjd1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onD1", to= "xkd1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onD1", to= "xld1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onD1", to= "xmd1", arrows = 1, free = F, values=1, labels = "lp2"),
                            
                            mxPath(from="onE1", to= "xae1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onE1", to= "xbe1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onE1", to= "xce1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onE1", to= "xde1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onE1", to= "xfe1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onE1", to= "xge1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onE1", to= "xhe1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onE1", to= "xie1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onE1", to= "xje1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onE1", to= "xke1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onE1", to= "xle1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onE1", to= "xme1", arrows = 1, free = F, values=1, labels = "lp2"),
                            
                            mxPath(from="onF1", to= "xaf1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onF1", to= "xbf1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onF1", to= "xcf1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onF1", to= "xdf1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onF1", to= "xef1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onF1", to= "xgf1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onF1", to= "xhf1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onF1", to= "xif1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onF1", to= "xjf1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onF1", to= "xkf1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onF1", to= "xlf1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onF1", to= "xmf1", arrows = 1, free = F, values=1, labels = "lp2"),
                            
                            mxPath(from="onG1", to= "xag1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onG1", to= "xbg1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onG1", to= "xcg1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onG1", to= "xdg1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onG1", to= "xeg1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onG1", to= "xfg1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onG1", to= "xhg1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onG1", to= "xig1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onG1", to= "xjg1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onG1", to= "xkg1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onG1", to= "xlg1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onG1", to= "xmg1", arrows = 1, free = F, values=1, labels = "lp2"),
                            
                            mxPath(from="onH1", to= "xah1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onH1", to= "xbh1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onH1", to= "xch1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onH1", to= "xdh1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onH1", to= "xeh1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onH1", to= "xfh1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onH1", to= "xgh1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onH1", to= "xih1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onH1", to= "xjh1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onH1", to= "xkh1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onH1", to= "xlh1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onH1", to= "xmh1", arrows = 1, free = F, values=1, labels = "lp2"),
                            
                            mxPath(from="onI1", to= "xai1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onI1", to= "xbi1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onI1", to= "xci1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onI1", to= "xdi1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onI1", to= "xei1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onI1", to= "xfi1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onI1", to= "xgi1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onI1", to= "xhi1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onI1", to= "xji1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onI1", to= "xki1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onI1", to= "xli1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onI1", to= "xmi1", arrows = 1, free = F, values=1, labels = "lp2"),
                            
                            mxPath(from="onJ1", to= "xaj1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onJ1", to= "xbj1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onJ1", to= "xcj1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onJ1", to= "xdj1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onJ1", to= "xej1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onJ1", to= "xfj1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onJ1", to= "xgj1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onJ1", to= "xhj1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onJ1", to= "xij1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onJ1", to= "xkj1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onJ1", to= "xlj1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onJ1", to= "xmj1", arrows = 1, free = F, values=1, labels = "lp2"),
                            
                            mxPath(from="onK1", to= "xak1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onK1", to= "xbk1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onK1", to= "xck1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onK1", to= "xdk1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onK1", to= "xek1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onK1", to= "xfk1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onK1", to= "xgk1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onK1", to= "xhk1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onK1", to= "xik1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onK1", to= "xjk1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onK1", to= "xlk1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onK1", to= "xmk1", arrows = 1, free = F, values=1, labels = "lp2"),
                            
                            mxPath(from="onL1", to= "xal1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onL1", to= "xbl1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onL1", to= "xcl1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onL1", to= "xdl1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onL1", to= "xel1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onL1", to= "xfl1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onL1", to= "xgl1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onL1", to= "xhl1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onL1", to= "xil1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onL1", to= "xjl1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onL1", to= "xkl1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onL1", to= "xml1", arrows = 1, free = F, values=1, labels = "lp2"),
                            
                            mxPath(from="onM1", to= "xam1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onM1", to= "xbm1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onM1", to= "xcm1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onM1", to= "xdm1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onM1", to= "xem1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onM1", to= "xfm1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onM1", to= "xgm1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onM1", to= "xhm1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onM1", to= "xim1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onM1", to= "xjm1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onM1", to= "xkm1", arrows = 1, free = F, values=1, labels = "lp2"),
                            mxPath(from="onM1", to= "xlm1", arrows = 1, free = F, values=1, labels = "lp2"),
                            
                            
                            mxPath(from="Aon1", arrows=2, free=T,  values=differentActorVar, labels = "av2"),
                            mxPath(from="Bon1", arrows=2, free=T,  values=differentActorVar, labels = "av2"),
                            mxPath(from="Con1", arrows=2, free=T,  values=differentActorVar, labels = "av2"),
                            mxPath(from="Don1", arrows=2, free=T,  values=differentActorVar, labels = "av2"),
                            mxPath(from="Eon1", arrows=2, free=T,  values=differentActorVar, labels = "av2"),
                            mxPath(from="Fon1", arrows=2, free=T,  values=differentActorVar, labels = "av2"),
                            mxPath(from="Gon1", arrows=2, free=T,  values=differentActorVar, labels = "av2"),
                            mxPath(from="Hon1", arrows=2, free=T,  values=differentActorVar, labels = "av2"),
                            mxPath(from="Ion1", arrows=2, free=T,  values=differentActorVar, labels = "av2"),
                            mxPath(from="Jon1", arrows=2, free=T,  values=differentActorVar, labels = "av2"),
                            mxPath(from="Kon1", arrows=2, free=T,  values=differentActorVar, labels = "av2"),
                            mxPath(from="Lon1", arrows=2, free=T,  values=differentActorVar, labels = "av2"),
                            mxPath(from="Mon1", arrows=2, free=T,  values=differentActorVar, labels = "av2"),
                            
                            mxPath(from="onA1", arrows=2, free=T,  values=differentPartnerVar, labels = "pv2"),
                            mxPath(from="onB1", arrows=2, free=T,  values=differentPartnerVar, labels = "pv2"),
                            mxPath(from="onC1", arrows=2, free=T,  values=differentPartnerVar, labels = "pv2"),
                            mxPath(from="onD1", arrows=2, free=T,  values=differentPartnerVar, labels = "pv2"),
                            mxPath(from="onE1", arrows=2, free=T,  values=differentPartnerVar, labels = "pv2"),
                            mxPath(from="onF1", arrows=2, free=T,  values=differentPartnerVar, labels = "pv2"),
                            mxPath(from="onG1", arrows=2, free=T,  values=differentPartnerVar, labels = "pv2"),
                            mxPath(from="onH1", arrows=2, free=T,  values=differentPartnerVar, labels = "pv2"),
                            mxPath(from="onI1", arrows=2, free=T,  values=differentPartnerVar, labels = "pv2"),
                            mxPath(from="onJ1", arrows=2, free=T,  values=differentPartnerVar, labels = "pv2"),
                            mxPath(from="onK1", arrows=2, free=T,  values=differentPartnerVar, labels = "pv2"),
                            mxPath(from="onL1", arrows=2, free=T,  values=differentPartnerVar, labels = "pv2"),
                            mxPath(from="onM1", arrows=2, free=T,  values=differentPartnerVar, labels = "pv2"),
                            
                            mxPath(from="Aon1", to = "onA1", arrows=2, free=T,  values=differentActorPartnerCov, labels = "apc2"),
                            mxPath(from="Bon1", to = "onB1", arrows=2, free=T,  values=differentActorPartnerCov, labels = "apc2"),
                            mxPath(from="Con1", to = "onC1", arrows=2, free=T,  values=differentActorPartnerCov, labels = "apc2"),
                            mxPath(from="Don1", to = "onD1", arrows=2, free=T,  values=differentActorPartnerCov, labels = "apc2"),
                            mxPath(from="Eon1", to = "onE1", arrows=2, free=T,  values=differentActorPartnerCov, labels = "apc2"),
                            mxPath(from="Fon1", to = "onF1", arrows=2, free=T,  values=differentActorPartnerCov, labels = "apc2"),
                            mxPath(from="Gon1", to = "onG1", arrows=2, free=T,  values=differentActorPartnerCov, labels = "apc2"),
                            mxPath(from="Hon1", to = "onH1", arrows=2, free=T,  values=differentActorPartnerCov, labels = "apc2"),
                            mxPath(from="Ion1", to = "onI1", arrows=2, free=T,  values=differentActorPartnerCov, labels = "apc2"),
                            mxPath(from="Jon1", to = "onJ1", arrows=2, free=T,  values=differentActorPartnerCov, labels = "apc2"),
                            mxPath(from="Kon1", to = "onK1", arrows=2, free=T,  values=differentActorPartnerCov, labels = "apc2"),
                            mxPath(from="Lon1", to = "onL1", arrows=2, free=T,  values=differentActorPartnerCov, labels = "apc2"),
                            mxPath(from="Mon1", to = "onM1", arrows=2, free=T,  values=differentActorPartnerCov, labels = "apc2"),
                            
                            mxPath(from="one", to = "onA1", arrows=1, free=F,  values=0, labels = "mean2"),
                            mxPath(from="one", to = "onB1", arrows=1, free=F,  values=0, labels = "mean2"),
                            mxPath(from="one", to = "onC1", arrows=1, free=F,  values=0, labels = "mean2"),
                            mxPath(from="one", to = "onD1", arrows=1, free=F,  values=0, labels = "mean2"),
                            mxPath(from="one", to = "onE1", arrows=1, free=F,  values=0, labels = "mean2"),
                            mxPath(from="one", to = "onF1", arrows=1, free=F,  values=0, labels = "mean2"),
                            mxPath(from="one", to = "onG1", arrows=1, free=F,  values=0, labels = "mean2"),
                            mxPath(from="one", to = "onH1", arrows=1, free=F,  values=0, labels = "mean2"),
                            mxPath(from="one", to = "onI1", arrows=1, free=F,  values=0, labels = "mean2"),
                            mxPath(from="one", to = "onJ1", arrows=1, free=F,  values=0, labels = "mean2"),
                            mxPath(from="one", to = "onK1", arrows=1, free=F,  values=0, labels = "mean2"),
                            mxPath(from="one", to = "onL1", arrows=1, free=F,  values=0, labels = "mean2"),
                            mxPath(from="one", to = "onM1", arrows=1, free=F,  values=0, labels = "mean2"),
                            
                            
                            #dyad level
                            
                            mxPath(from="xab1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xac1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xad1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xae1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xaf1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xag1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xah1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xai1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xaj1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xak1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xal1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xam1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            
                            mxPath(from="xba1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xbc1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xbd1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xbe1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xbf1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xbg1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xbh1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xbi1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xbj1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xbk1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xbl1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xbm1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            
                            mxPath(from="xca1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xcb1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xcd1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xce1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xcf1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xcg1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xch1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xci1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xcj1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xck1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xcl1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xcm1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            
                            mxPath(from="xda1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xdb1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xdc1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xde1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xdf1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xdg1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xdh1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xdi1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xdj1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xdk1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xdl1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xdm1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            
                            mxPath(from="xea1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xeb1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xec1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xed1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xef1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xeg1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xeh1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xei1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xej1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xek1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xel1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xem1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            
                            mxPath(from="xfa1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xfb1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xfc1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xfd1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xfe1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xfg1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xfh1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xfi1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xfj1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xfk1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xfl1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xfm1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            
                            mxPath(from="xga1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xgb1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xgc1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xgd1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xge1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xgf1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xgh1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xgi1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xgj1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xgk1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xgl1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xgm1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            
                            mxPath(from="xha1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xhb1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xhc1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xhd1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xhe1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xhf1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xhg1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xhi1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xhj1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xhk1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xhl1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xhm1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            
                            mxPath(from="xia1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xib1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xic1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xid1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xie1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xif1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xig1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xih1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xij1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xik1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xil1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xim1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            
                            mxPath(from="xja1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xjb1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xjc1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xjd1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xje1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xjf1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xjg1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xjh1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xji1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xjk1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xjl1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xjm1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            
                            mxPath(from="xka1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xkb1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xkc1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xkd1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xke1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xkf1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xkg1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xkh1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xki1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xkj1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xkl1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xkm1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            
                            mxPath(from="xla1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xlb1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xlc1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xld1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xle1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xlf1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xlg1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xlh1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xli1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xlj1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xlk1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xlm1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            
                            mxPath(from="xma1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xmb1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xmc1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xmd1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xme1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xmf1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xmg1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xmh1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xmi1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xmj1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xmk1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            mxPath(from="xml1", arrows=2, free=T,  values=differentRelationshipVar, labels = "rv2"),
                            
                            
                            mxPath(from="xab1", to = "xba1", arrows=2, free=T,  values=differentRelationshipCov, labels = "rc2"),
                            mxPath(from="xac1", to = "xca1", arrows=2, free=T,  values=differentRelationshipCov, labels = "rc2"),
                            mxPath(from="xad1", to = "xda1", arrows=2, free=T,  values=differentRelationshipCov, labels = "rc2"),
                            mxPath(from="xae1", to = "xea1", arrows=2, free=T,  values=differentRelationshipCov, labels = "rc2"),
                            mxPath(from="xaf1", to = "xfa1", arrows=2, free=T,  values=differentRelationshipCov, labels = "rc2"),
                            mxPath(from="xag1", to = "xga1", arrows=2, free=T,  values=differentRelationshipCov, labels = "rc2"),
                            mxPath(from="xah1", to = "xha1", arrows=2, free=T,  values=differentRelationshipCov, labels = "rc2"),
                            mxPath(from="xai1", to = "xia1", arrows=2, free=T,  values=differentRelationshipCov, labels = "rc2"),
                            mxPath(from="xaj1", to = "xja1", arrows=2, free=T,  values=differentRelationshipCov, labels = "rc2"),
                            mxPath(from="xak1", to = "xka1", arrows=2, free=T,  values=differentRelationshipCov, labels = "rc2"),
                            mxPath(from="xal1", to = "xla1", arrows=2, free=T,  values=differentRelationshipCov, labels = "rc2"),
                            mxPath(from="xam1", to = "xma1", arrows=2, free=T,  values=differentRelationshipCov, labels = "rc2"),
                            
                            mxPath(from="xbc1", to = "xcb1", arrows=2, free=T,  values=differentRelationshipCov, labels = "rc2"),
                            mxPath(from="xbd1", to = "xdb1", arrows=2, free=T,  values=differentRelationshipCov, labels = "rc2"),
                            mxPath(from="xbe1", to = "xeb1", arrows=2, free=T,  values=differentRelationshipCov, labels = "rc2"),
                            mxPath(from="xbf1", to = "xfb1", arrows=2, free=T,  values=differentRelationshipCov, labels = "rc2"),
                            mxPath(from="xbg1", to = "xgb1", arrows=2, free=T,  values=differentRelationshipCov, labels = "rc2"),
                            mxPath(from="xbh1", to = "xhb1", arrows=2, free=T,  values=differentRelationshipCov, labels = "rc2"),
                            mxPath(from="xbi1", to = "xib1", arrows=2, free=T,  values=differentRelationshipCov, labels = "rc2"),
                            mxPath(from="xbj1", to = "xjb1", arrows=2, free=T,  values=differentRelationshipCov, labels = "rc2"),
                            mxPath(from="xbk1", to = "xkb1", arrows=2, free=T,  values=differentRelationshipCov, labels = "rc2"),
                            mxPath(from="xbl1", to = "xlb1", arrows=2, free=T,  values=differentRelationshipCov, labels = "rc2"),
                            mxPath(from="xbm1", to = "xmb1", arrows=2, free=T,  values=differentRelationshipCov, labels = "rc2"),
                            
                            mxPath(from="xcd1", to = "xdc1", arrows=2, free=T,  values=differentRelationshipCov, labels = "rc2"),
                            mxPath(from="xce1", to = "xec1", arrows=2, free=T,  values=differentRelationshipCov, labels = "rc2"),
                            mxPath(from="xcf1", to = "xfc1", arrows=2, free=T,  values=differentRelationshipCov, labels = "rc2"),
                            mxPath(from="xcg1", to = "xgc1", arrows=2, free=T,  values=differentRelationshipCov, labels = "rc2"),
                            mxPath(from="xch1", to = "xhc1", arrows=2, free=T,  values=differentRelationshipCov, labels = "rc2"),
                            mxPath(from="xci1", to = "xic1", arrows=2, free=T,  values=differentRelationshipCov, labels = "rc2"),
                            mxPath(from="xcj1", to = "xjc1", arrows=2, free=T,  values=differentRelationshipCov, labels = "rc2"),
                            mxPath(from="xck1", to = "xkc1", arrows=2, free=T,  values=differentRelationshipCov, labels = "rc2"),
                            mxPath(from="xcl1", to = "xlc1", arrows=2, free=T,  values=differentRelationshipCov, labels = "rc2"),
                            mxPath(from="xcm1", to = "xmc1", arrows=2, free=T,  values=differentRelationshipCov, labels = "rc2"),
                            
                            mxPath(from="xde1", to = "xed1", arrows=2, free=T,  values=differentRelationshipCov, labels = "rc2"),
                            mxPath(from="xdf1", to = "xfd1", arrows=2, free=T,  values=differentRelationshipCov, labels = "rc2"),
                            mxPath(from="xdg1", to = "xgd1", arrows=2, free=T,  values=differentRelationshipCov, labels = "rc2"),
                            mxPath(from="xdh1", to = "xhd1", arrows=2, free=T,  values=differentRelationshipCov, labels = "rc2"),
                            mxPath(from="xdi1", to = "xid1", arrows=2, free=T,  values=differentRelationshipCov, labels = "rc2"),
                            mxPath(from="xdj1", to = "xjd1", arrows=2, free=T,  values=differentRelationshipCov, labels = "rc2"),
                            mxPath(from="xdk1", to = "xkd1", arrows=2, free=T,  values=differentRelationshipCov, labels = "rc2"),
                            mxPath(from="xdl1", to = "xld1", arrows=2, free=T,  values=differentRelationshipCov, labels = "rc2"),
                            mxPath(from="xdm1", to = "xmd1", arrows=2, free=T,  values=differentRelationshipCov, labels = "rc2"),
                            
                            mxPath(from="xef1", to = "xfe1", arrows=2, free=T,  values=differentRelationshipCov, labels = "rc2"),
                            mxPath(from="xeg1", to = "xge1", arrows=2, free=T,  values=differentRelationshipCov, labels = "rc2"),
                            mxPath(from="xeh1", to = "xhe1", arrows=2, free=T,  values=differentRelationshipCov, labels = "rc2"),
                            mxPath(from="xei1", to = "xie1", arrows=2, free=T,  values=differentRelationshipCov, labels = "rc2"),
                            mxPath(from="xej1", to = "xje1", arrows=2, free=T,  values=differentRelationshipCov, labels = "rc2"),
                            mxPath(from="xek1", to = "xke1", arrows=2, free=T,  values=differentRelationshipCov, labels = "rc2"),
                            mxPath(from="xel1", to = "xle1", arrows=2, free=T,  values=differentRelationshipCov, labels = "rc2"),
                            mxPath(from="xem1", to = "xme1", arrows=2, free=T,  values=differentRelationshipCov, labels = "rc2"),
                            
                            mxPath(from="xfg1", to = "xgf1", arrows=2, free=T,  values=differentRelationshipCov, labels = "rc2"),
                            mxPath(from="xfh1", to = "xhf1", arrows=2, free=T,  values=differentRelationshipCov, labels = "rc2"),
                            mxPath(from="xfi1", to = "xif1", arrows=2, free=T,  values=differentRelationshipCov, labels = "rc2"),
                            mxPath(from="xfj1", to = "xjf1", arrows=2, free=T,  values=differentRelationshipCov, labels = "rc2"),
                            mxPath(from="xfk1", to = "xkf1", arrows=2, free=T,  values=differentRelationshipCov, labels = "rc2"),
                            mxPath(from="xfl1", to = "xlf1", arrows=2, free=T,  values=differentRelationshipCov, labels = "rc2"),
                            mxPath(from="xfm1", to = "xmf1", arrows=2, free=T,  values=differentRelationshipCov, labels = "rc2"),
                            
                            mxPath(from="xgh1", to = "xhg1", arrows=2, free=T,  values=differentRelationshipCov, labels = "rc2"),
                            mxPath(from="xgi1", to = "xig1", arrows=2, free=T,  values=differentRelationshipCov, labels = "rc2"),
                            mxPath(from="xgj1", to = "xjg1", arrows=2, free=T,  values=differentRelationshipCov, labels = "rc2"),
                            mxPath(from="xgk1", to = "xkg1", arrows=2, free=T,  values=differentRelationshipCov, labels = "rc2"),
                            mxPath(from="xgl1", to = "xlg1", arrows=2, free=T,  values=differentRelationshipCov, labels = "rc2"),
                            mxPath(from="xgm1", to = "xmg1", arrows=2, free=T,  values=differentRelationshipCov, labels = "rc2"),
                            
                            mxPath(from="xhi1", to = "xih1", arrows=2, free=T,  values=differentRelationshipCov, labels = "rc2"),
                            mxPath(from="xhj1", to = "xjh1", arrows=2, free=T,  values=differentRelationshipCov, labels = "rc2"),
                            mxPath(from="xhk1", to = "xkh1", arrows=2, free=T,  values=differentRelationshipCov, labels = "rc2"),
                            mxPath(from="xhl1", to = "xlh1", arrows=2, free=T,  values=differentRelationshipCov, labels = "rc2"),
                            mxPath(from="xhm1", to = "xmh1", arrows=2, free=T,  values=differentRelationshipCov, labels = "rc2"),
                            
                            mxPath(from="xij1", to = "xji1", arrows=2, free=T,  values=differentRelationshipCov, labels = "rc2"),
                            mxPath(from="xik1", to = "xki1", arrows=2, free=T,  values=differentRelationshipCov, labels = "rc2"),
                            mxPath(from="xil1", to = "xli1", arrows=2, free=T,  values=differentRelationshipCov, labels = "rc2"),
                            mxPath(from="xim1", to = "xmi1", arrows=2, free=T,  values=differentRelationshipCov, labels = "rc2"),
                            
                            mxPath(from="xjk1", to = "xkj1", arrows=2, free=T,  values=differentRelationshipCov, labels = "rc2"),
                            mxPath(from="xjl1", to = "xlj1", arrows=2, free=T,  values=differentRelationshipCov, labels = "rc2"),
                            mxPath(from="xjm1", to = "xmj1", arrows=2, free=T,  values=differentRelationshipCov, labels = "rc2"),
                            
                            mxPath(from="xkl1", to = "xlk1", arrows=2, free=T,  values=differentRelationshipCov, labels = "rc2"),
                            mxPath(from="xkm1", to = "xmk1", arrows=2, free=T,  values=differentRelationshipCov, labels = "rc2"),
                            
                            mxPath(from="xlm1", to = "xml1", arrows=2, free=T,  values=differentRelationshipCov, labels = "rc2")
)



####################################

#start simulation

############################################################################################################################

set.seed(161103132)



#1vs5 groups, groupsize 5, different, rv

prophvector1vs5gs5different_rv <- rep(NA, 1300)

for (i in 1:1300){
  
  random <- sample(1:100000, size = 6, replace =F)
  
  bothgroups <- mxModel("bothgroups",
                        mxModel(base5, mxData(PopSame[random[1],], type="raw")), mxModel(other5different, mxData(PopDifferent[random[2:6],], type="raw")),
                        mxFitFunctionMultigroup(c("base5", "other5different")))
  
  bothgroupse <- mxModel("bothgroupse",
                         mxModel(base5, mxData(PopSame[random[1],], type="raw")), mxModel(other5different, mxData(PopDifferent[random[2:6],], type="raw")),
                         mxConstraint( (rv1 / (av1 + pv1 + rv1)) == (rv2 / (av2 + pv2 + rv2)) ),
                         mxFitFunctionMultigroup(c("base5", "other5different")))
  
  
  if (has_warning(bothgroupsfit <- mxRunDave(bothgroups)) == T){
    prophvector1vs5gs5different_rv[i] <- -888
  }
  else if (has_warning(bothgroupsefit <- mxRunDave(bothgroupse)) == T){
    if (grepl("optimality conditions to the required accuracy", summary(bothgroupsefit)[39]$npsolMessage) ==T){
      prophvector1vs5gs5different_rv[i] <- summary(bothgroupsefit)[8]$Minus2LogLikelihood - summary(bothgroupsfit)[8]$Minus2LogLikelihood
    }
    else {
      prophvector1vs5gs5different_rv[i] <- -8888}
  }
  else if (length(is(bothgroupsfit)) >  2){
    prophvector1vs5gs5different_rv[i] <- -999
  }
  else if (length(is(bothgroupsefit)) >  2){
    prophvector1vs5gs5different_rv[i] <- -9999
  }
  else {
    prophvector1vs5gs5different_rv[i] <- summary(bothgroupsefit)[8]$Minus2LogLikelihood - summary(bothgroupsfit)[8]$Minus2LogLikelihood
  }
}
save(prophvector1vs5gs5different_rv, file = "prophvector1vs5gs5different_rv.rda")



############################################################################################################################


#1vs5 groups, groupsize 13, different, rv

prophvector1vs5gs13different_rv <- rep(NA, 1300)

for (i in 1:1300){
  
  random <- sample(1:100000, size = 6, replace =F)
  
  bothgroups <- mxModel("bothgroups",
                        mxModel(base13, mxData(PopSame[random[1],], type="raw")), mxModel(other13different, mxData(PopDifferent[random[2:6],], type="raw")),
                        mxFitFunctionMultigroup(c("base13", "other13different")))
  
  bothgroupse <- mxModel("bothgroupse",
                         mxModel(base13, mxData(PopSame[random[1],], type="raw")), mxModel(other13different, mxData(PopDifferent[random[2:6],], type="raw")),
                         mxConstraint( (rv1 / (av1 + pv1 + rv1)) == (rv2 / (av2 + pv2 + rv2)) ),
                         mxFitFunctionMultigroup(c("base13", "other13different")))
  
  
  if (has_warning(bothgroupsfit <- mxRunDave(bothgroups)) == T){
    prophvector1vs5gs13different_rv[i] <- -888
  }
  else if (has_warning(bothgroupsefit <- mxRunDave(bothgroupse)) == T){
    if (grepl("optimality conditions to the required accuracy", summary(bothgroupsefit)[39]$npsolMessage) ==T){
      prophvector1vs5gs13different_rv[i] <- summary(bothgroupsefit)[8]$Minus2LogLikelihood - summary(bothgroupsfit)[8]$Minus2LogLikelihood
    }
    else {
      prophvector1vs5gs13different_rv[i] <- -8888}
  }
  else if (length(is(bothgroupsfit)) >  2){
    prophvector1vs5gs13different_rv[i] <- -999
  }
  else if (length(is(bothgroupsefit)) >  2){
    prophvector1vs5gs13different_rv[i] <- -9999
  }
  else {
    prophvector1vs5gs13different_rv[i] <- summary(bothgroupsefit)[8]$Minus2LogLikelihood - summary(bothgroupsfit)[8]$Minus2LogLikelihood
  }
}
save(prophvector1vs5gs13different_rv, file = "prophvector1vs5gs13different_rv.rda")



############################################################################################################################



#5vs5 groups, groupsize 5, different, rv

prophvector5vs5gs5different_rv <- rep(NA, 1300)

for (i in 1:1300){
  
  random <- sample(1:100000, size = 10, replace =F)
  
  bothgroups <- mxModel("bothgroups",
                        mxModel(base5, mxData(PopSame[random[1:5],], type="raw")), mxModel(other5different, mxData(PopDifferent[random[6:10],], type="raw")),
                        mxFitFunctionMultigroup(c("base5", "other5different")))
  
  bothgroupse <- mxModel("bothgroupse",
                         mxModel(base5, mxData(PopSame[random[1:5],], type="raw")), mxModel(other5different, mxData(PopDifferent[random[6:10],], type="raw")),
                         mxConstraint( (rv1 / (av1 + pv1 + rv1)) == (rv2 / (av2 + pv2 + rv2)) ),
                         mxFitFunctionMultigroup(c("base5", "other5different")))
  
  
  if (has_warning(bothgroupsfit <- mxRunDave(bothgroups)) == T){
    prophvector5vs5gs5different_rv[i] <- -888
  }
  else if (has_warning(bothgroupsefit <- mxRunDave(bothgroupse)) == T){
    if (grepl("optimality conditions to the required accuracy", summary(bothgroupsefit)[39]$npsolMessage) ==T){
      prophvector5vs5gs5different_rv[i] <- summary(bothgroupsefit)[8]$Minus2LogLikelihood - summary(bothgroupsfit)[8]$Minus2LogLikelihood
    }
    else {
      prophvector5vs5gs5different_rv[i] <- -8888}
  }
  else if (length(is(bothgroupsfit)) >  2){
    prophvector5vs5gs5different_rv[i] <- -999
  }
  else if (length(is(bothgroupsefit)) >  2){
    prophvector5vs5gs5different_rv[i] <- -9999
  }
  else {
    prophvector5vs5gs5different_rv[i] <- summary(bothgroupsefit)[8]$Minus2LogLikelihood - summary(bothgroupsfit)[8]$Minus2LogLikelihood
  }
}
save(prophvector5vs5gs5different_rv, file = "prophvector5vs5gs5different_rv.rda")


############################################################################################################################


#5vs5 groups, groupsize 13, different, rv

prophvector5vs5gs13different_rv <- rep(NA, 1300)

for (i in 1:1300){
  
  random <- sample(1:100000, size = 10, replace =F)
  
  bothgroups <- mxModel("bothgroups",
                        mxModel(base13, mxData(PopSame[random[1:5],], type="raw")), mxModel(other13different, mxData(PopDifferent[random[6:10],], type="raw")),
                        mxFitFunctionMultigroup(c("base13", "other13different")))
  
  bothgroupse <- mxModel("bothgroupse",
                         mxModel(base13, mxData(PopSame[random[1:5],], type="raw")), mxModel(other13different, mxData(PopDifferent[random[6:10],], type="raw")),
                         mxConstraint( (rv1 / (av1 + pv1 + rv1)) == (rv2 / (av2 + pv2 + rv2)) ),
                         mxFitFunctionMultigroup(c("base13", "other13different")))
  
  
  if (has_warning(bothgroupsfit <- mxRunDave(bothgroups)) == T){
    prophvector5vs5gs13different_rv[i] <- -888
  }
  else if (has_warning(bothgroupsefit <- mxRunDave(bothgroupse)) == T){
    if (grepl("optimality conditions to the required accuracy", summary(bothgroupsefit)[39]$npsolMessage) ==T){
      prophvector5vs5gs13different_rv[i] <- summary(bothgroupsefit)[8]$Minus2LogLikelihood - summary(bothgroupsfit)[8]$Minus2LogLikelihood
    }
    else {
      prophvector5vs5gs13different_rv[i] <- -8888}
  }
  else if (length(is(bothgroupsfit)) >  2){
    prophvector5vs5gs13different_rv[i] <- -999
  }
  else if (length(is(bothgroupsefit)) >  2){
    prophvector5vs5gs13different_rv[i] <- -9999
  }
  else {
    prophvector5vs5gs13different_rv[i] <- summary(bothgroupsefit)[8]$Minus2LogLikelihood - summary(bothgroupsfit)[8]$Minus2LogLikelihood
  }
}
save(prophvector5vs5gs13different_rv, file = "prophvector5vs5gs13different_rv.rda")


############################################################################################################################



#10vs10 groups, groupsize 5, different, rv

prophvector10vs10gs5different_rv <- rep(NA, 1300)

for (i in 1:1300){
  
  random <- sample(1:100000, size = 20, replace =F)
  
  bothgroups <- mxModel("bothgroups",
                        mxModel(base5, mxData(PopSame[random[1:10],], type="raw")), mxModel(other5different, mxData(PopDifferent[random[11:20],], type="raw")),
                        mxFitFunctionMultigroup(c("base5", "other5different")))
  
  bothgroupse <- mxModel("bothgroupse",
                         mxModel(base5, mxData(PopSame[random[1:10],], type="raw")), mxModel(other5different, mxData(PopDifferent[random[11:20],], type="raw")),
                         mxConstraint( (rv1 / (av1 + pv1 + rv1)) == (rv2 / (av2 + pv2 + rv2)) ),
                         mxFitFunctionMultigroup(c("base5", "other5different")))
  
  
  if (has_warning(bothgroupsfit <- mxRunDave(bothgroups)) == T){
    prophvector10vs10gs5different_rv[i] <- -888
  }
  else if (has_warning(bothgroupsefit <- mxRunDave(bothgroupse)) == T){
    if (grepl("optimality conditions to the required accuracy", summary(bothgroupsefit)[39]$npsolMessage) ==T){
      prophvector10vs10gs5different_rv[i] <- summary(bothgroupsefit)[8]$Minus2LogLikelihood - summary(bothgroupsfit)[8]$Minus2LogLikelihood
    }
    else {
      prophvector10vs10gs5different_rv[i] <- -8888}
  }
  else if (length(is(bothgroupsfit)) >  2){
    prophvector10vs10gs5different_rv[i] <- -999
  }
  else if (length(is(bothgroupsefit)) >  2){
    prophvector10vs10gs5different_rv[i] <- -9999
  }
  else {
    prophvector10vs10gs5different_rv[i] <- summary(bothgroupsefit)[8]$Minus2LogLikelihood - summary(bothgroupsfit)[8]$Minus2LogLikelihood
  }
}
save(prophvector10vs10gs5different_rv, file = "prophvector10vs10gs5different_rv.rda")


############################################################################################################################



#10vs10 groups, groupsize 13, different, rv

prophvector10vs10gs13different_rv <- rep(NA, 1300)

for (i in 1:1300){
  
  random <- sample(1:100000, size = 20, replace =F)
  
  bothgroups <- mxModel("bothgroups",
                        mxModel(base13, mxData(PopSame[random[1:10],], type="raw")), mxModel(other13different, mxData(PopDifferent[random[11:20],], type="raw")),
                        mxFitFunctionMultigroup(c("base13", "other13different")))
  
  bothgroupse <- mxModel("bothgroupse",
                         mxModel(base13, mxData(PopSame[random[1:10],], type="raw")), mxModel(other13different, mxData(PopDifferent[random[11:20],], type="raw")),
                         mxConstraint( (rv1 / (av1 + pv1 + rv1)) == (rv2 / (av2 + pv2 + rv2)) ),
                         mxFitFunctionMultigroup(c("base13", "other13different")))
  
  
  if (has_warning(bothgroupsfit <- mxRunDave(bothgroups)) == T){
    prophvector10vs10gs13different_rv[i] <- -888
  }
  else if (has_warning(bothgroupsefit <- mxRunDave(bothgroupse)) == T){
    if (grepl("optimality conditions to the required accuracy", summary(bothgroupsefit)[39]$npsolMessage) ==T){
      prophvector10vs10gs13different_rv[i] <- summary(bothgroupsefit)[8]$Minus2LogLikelihood - summary(bothgroupsfit)[8]$Minus2LogLikelihood
    }
    else {
      prophvector10vs10gs13different_rv[i] <- -8888}
  }
  else if (length(is(bothgroupsfit)) >  2){
    prophvector10vs10gs13different_rv[i] <- -999
  }
  else if (length(is(bothgroupsefit)) >  2){
    prophvector10vs10gs13different_rv[i] <- -9999
  }
  else {
    prophvector10vs10gs13different_rv[i] <- summary(bothgroupsefit)[8]$Minus2LogLikelihood - summary(bothgroupsfit)[8]$Minus2LogLikelihood
  }
}
save(prophvector10vs10gs13different_rv, file = "prophvector10vs10gs13different_rv.rda")


############################################################################################################################





#20vs20 groups, groupsize 5, different, rv

prophvector20vs20gs5different_rv <- rep(NA, 1300)

for (i in 1:1300){
  
  random <- sample(1:100000, size = 40, replace =F)
  
  bothgroups <- mxModel("bothgroups",
                        mxModel(base5, mxData(PopSame[random[1:20],], type="raw")), mxModel(other5different, mxData(PopDifferent[random[21:40],], type="raw")),
                        mxFitFunctionMultigroup(c("base5", "other5different")))
  
  bothgroupse <- mxModel("bothgroupse",
                         mxModel(base5, mxData(PopSame[random[1:20],], type="raw")), mxModel(other5different, mxData(PopDifferent[random[21:40],], type="raw")),
                         mxConstraint( (rv1 / (av1 + pv1 + rv1)) == (rv2 / (av2 + pv2 + rv2)) ),
                         mxFitFunctionMultigroup(c("base5", "other5different")))
  
  
  if (has_warning(bothgroupsfit <- mxRunDave(bothgroups)) == T){
    prophvector20vs20gs5different_rv[i] <- -888
  }
  else if (has_warning(bothgroupsefit <- mxRunDave(bothgroupse)) == T){
    if (grepl("optimality conditions to the required accuracy", summary(bothgroupsefit)[39]$npsolMessage) ==T){
      prophvector20vs20gs5different_rv[i] <- summary(bothgroupsefit)[8]$Minus2LogLikelihood - summary(bothgroupsfit)[8]$Minus2LogLikelihood
    }
    else {
      prophvector20vs20gs5different_rv[i] <- -8888}
  }
  else if (length(is(bothgroupsfit)) >  2){
    prophvector20vs20gs5different_rv[i] <- -999
  }
  else if (length(is(bothgroupsefit)) >  2){
    prophvector20vs20gs5different_rv[i] <- -9999
  }
  else {
    prophvector20vs20gs5different_rv[i] <- summary(bothgroupsefit)[8]$Minus2LogLikelihood - summary(bothgroupsfit)[8]$Minus2LogLikelihood
  }
}
save(prophvector20vs20gs5different_rv, file = "prophvector20vs20gs5different_rv.rda")


############################################################################################################################



#20vs20 groups, groupsize 13, different, rv

prophvector20vs20gs13different_rv <- rep(NA, 1300)

for (i in 1:1300){
  
  random <- sample(1:100000, size = 40, replace =F)
  
  bothgroups <- mxModel("bothgroups",
                        mxModel(base13, mxData(PopSame[random[1:20],], type="raw")), mxModel(other13different, mxData(PopDifferent[random[21:40],], type="raw")),
                        mxFitFunctionMultigroup(c("base13", "other13different")))
  
  bothgroupse <- mxModel("bothgroupse",
                         mxModel(base13, mxData(PopSame[random[1:20],], type="raw")), mxModel(other13different, mxData(PopDifferent[random[21:40],], type="raw")),
                         mxConstraint( (rv1 / (av1 + pv1 + rv1)) == (rv2 / (av2 + pv2 + rv2)) ),
                         mxFitFunctionMultigroup(c("base13", "other13different")))
  
  
  if (has_warning(bothgroupsfit <- mxRunDave(bothgroups)) == T){
    prophvector20vs20gs13different_rv[i] <- -888
  }
  else if (has_warning(bothgroupsefit <- mxRunDave(bothgroupse)) == T){
    if (grepl("optimality conditions to the required accuracy", summary(bothgroupsefit)[39]$npsolMessage) ==T){
      prophvector20vs20gs13different_rv[i] <- summary(bothgroupsefit)[8]$Minus2LogLikelihood - summary(bothgroupsfit)[8]$Minus2LogLikelihood
    }
    else {
      prophvector20vs20gs13different_rv[i] <- -8888}
  }
  else if (length(is(bothgroupsfit)) >  2){
    prophvector20vs20gs13different_rv[i] <- -999
  }
  else if (length(is(bothgroupsefit)) >  2){
    prophvector20vs20gs13different_rv[i] <- -9999
  }
  else {
    prophvector20vs20gs13different_rv[i] <- summary(bothgroupsefit)[8]$Minus2LogLikelihood - summary(bothgroupsfit)[8]$Minus2LogLikelihood
  }
}
save(prophvector20vs20gs13different_rv, file = "prophvector20vs20gs13different_rv.rda")



############################################################################################################################


