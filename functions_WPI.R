#function to split the data into different groups (criteria),calculate the WPI
#for each group and return both the WPI median and bounds and a new objects that contain
#the raw data split, so it can be further split into other groups
#Arguments: object: a list of 16 elements (one for each site) that contains all the
#species info for that site (occupnacy matrix, guild, redlist status, etc)
#Crit: a character vector with the different levels to split the data into different
# criteria (e.g. Continent, red list status, etc)

SplitdataCalculateWPI<-function(object,crit,nsites,iters,years){
  require(abind)
  require(reshape2)
  nc <- length(years)
  #name of the criteria
  crit.name <- deparse(substitute(crit))
  #sequence<-seq(from=1,to=17000,by=1000)
  wpiMatrix<-array(data=NA,dim=c(nsites*iters,length(years),length(crit)),dimnames=list(NULL,years,as.character(crit)))
  #object to store the new split objects
  newob<-rep(list(object),length(crit))
  names(newob) <- as.character(crit)
  for(j in 1:length(crit)){
    wpi<-numeric()
    for(i in 1:nsites){
      temp<-sapply(object[[i]][[1]], function(x){x[[crit.name]]==crit[j]})
      #this object just contains a list of the species with status "LC"
      if(sum(temp)>0) {
        temp<-subset(object[[i]][[1]],temp==T)
        newob[[j]][[i]][[1]]<-temp
        #SECOND step: Now calculate the WPI for just this group
        #get the matrices out of the list and put them in a 3D array
        occMat<-do.call(abind,list(lapply(temp, "[[", "occMat"),along=3,use.first.dimnames=T))
        #Now calculate WPI
        newwpi<-f.WPI(occMat,object[[i]]$startyear)
        #if(object[[i]]$startyear>1){
        newwpi[,object[[i]]$startyear]<-NA
        #}
        wpi<-rbind(wpi,newwpi)
      }
      else {newwpi<-matrix(NA,nr=iters,nc=length(years))
            wpi<-rbind(wpi,newwpi)
            newob[[j]][[i]][1]<-NULL
            next}#Do this for the next site
    }
    wpiMatrix[,,j]<-wpi
    wpi<-numeric()
  }
  #Extract the wpi statistics
  
  #calculate median
  med<-apply(wpiMatrix,c(2,3),median,na.rm=T)
  #root the wpi values in the first year to 1
  posNA<-apply(med,2,function(x){match(T,!is.na(x))-1})
  coords<-cbind(as.numeric(posNA),1:dim(wpiMatrix)[3])
  med[coords]<-1
  
  #calculate low and hi 80%
  lo80<-apply(wpiMatrix,c(2,3),quantile,0.10,na.rm=T)
  lo80[coords]<-1
  hi80<-apply(wpiMatrix,c(2,3),quantile,0.90,na.rm=T)
  hi80[coords]<-1
  
  #calculate low and hi 50%
  lo50<-apply(wpiMatrix,c(2,3),quantile,0.25,na.rm=T)
  lo50[coords]<-1
  hi50<-apply(wpiMatrix,c(2,3),quantile,0.75,na.rm=T)
  hi50[coords]<-1
  
  fact<-factor(x=c('median','lo50','hi50','lo80','hi80'))
  med<-data.frame(med,stats=levels(fact)[5],year=rownames(med))
  lo80<-data.frame(lo80,stats=levels(fact)[4],year=rownames(med))
  hi80<-data.frame(hi80,stats=levels(fact)[2],year=rownames(med))
  lo50<-data.frame(lo50,stats=levels(fact)[3],year=rownames(med))
  hi50<-data.frame(hi50,stats=levels(fact)[1],year=rownames(med))
  
  res<-rbind(med,lo80,hi80,lo50,hi50)
  res<-melt(res)
  res<-dcast(data=res,formula=year+variable~stats,value.var="value")
  res$year<-as.numeric(as.character(res$year))
  #return both the wpi results and a new list split by crit
  list(wpi=res,splitob=newob)
  #wpiMatrix
}

#function to generate the WPI from the output simulations in JAGS
# psi is a three dimensional matrix with the psi of each species in each year
f.WPI <-function(psi,startyear){
  nsim <- dim(psi)[1]
  nyears <- dim(psi)[2]
  nsp <- dim(psi)[3]
  rel_psi<-numeric()
  wpi<-matrix(NA,nr=nsim,nc=nyears)
  for(i in 1:nsim){
    for(t in startyear:nyears){
      for(s in 1:nsp){
        rel_psi[s] <- log(psi[i,t,s]/psi[i,startyear,s])
      }
      wpi[i,t]<-exp(1/nsp*sum(rel_psi))
    }
  }
  colnames(wpi)<-dimnames(psi)[[2]]
  wpi
  
}

#function to calculate whether a species is increasing or decreasing
# but using lambda
speciesChangesLambda <- function(data){
  require(reshape2)
  mat <- acast(data,iteration~year,value.var="psi")
  lambda <- calc.lambda(mat)
  interv <- dim(lambda)[2]
  median.lam <- median(lambda)
  lo95 <- apply(lambda,2,quantile,0.025)
  hi95 <- apply(lambda,2,quantile,0.975)
  lo80 <- apply(lambda,2,quantile,0.10)
  hi80 <- apply(lambda,2,quantile,0.90)
  
  ind <- numeric(interv)
  ind2 <- numeric(interv)
  for(i in 1:interv){
    if(lo95[i] <= 1 & hi95[i] >=1)
      ind[i] <- 0
    else if(lo95[i] > 1 & hi95[i] > 1)
      ind[i] <- 1
    else
      ind[i] <- -1
  }
  ind95 <- mean(ind)
  for(i in 1:interv){
    if(lo80[i] <= 1 & hi80[i] >=1)
      ind2[i] <- 0
    else if(lo80[i] > 1 & hi80[i] > 1)
      ind2[i] <- 1
    else
      ind2[i] <- -1
  }
  ind80 <- mean(ind2)
  
  sitesp <- data$site.sp[1]
  cont <- data$cont[1]
  rls <- data$rls[1]
  mass <- data$mass[1]
  class <- data$class[1]
  guild <- data$guild[1]
  site <- data$site.abb[1]
  sp <- data$bin[1]
  res <- data.frame(sp,median.lam,ind95,ind80,cont,rls,mass,class,guild,site)
  res
}

speciesChangesLogistic <- function(data){
  #function to calculate a logistic regression of occupancy as a function of time
  require(reshape2)
  require(plyr)
  mat <- acast(data,iteration~year,value.var="psi")
  years <-1:dim(mat)[2]
  coeffs<-apply(mat,1,logRegression,60,years)
  #coeffs
  median.intercept <- median(coeffs[1,])
  median.slope <- median(coeffs[2,])
  ci95 <- quantile(coeffs[2,],c(0.025,0.975))
  ci80 <- quantile(coeffs[2,],c(0.10,0.9))
  
  ind95<-character()
  ind80<-character()
  
  if(ci95[1] <= 0 & ci95[2] >=0)
    ind95 <- "stable"
  else if(ci95[1] > 0 & ci95[2] > 0)
    ind95 <- "increasing"
  else
    ind95 <- "decreasing"
  
  if(ci80[1] <= 0 & ci80[2] >=0)
    ind80 <- "stable"
  else if(ci80[1] > 0 & ci80[2] > 0)
    ind80 <- "increasing"
  else
    ind80 <- "decreasing"
  
  sitesp <- data$site.sp[1]
  cont <- data$cont[1]
  rls <- data$rls[1]
  mass <- data$mass[1]
  class <- data$class[1]
  guild <- data$guild[1]
  site <- data$site.abb[1]
  sp <- data$bin[1]
  nyears <- length(years)
  rare80 <- data$Rare80[1]
  res <- data.frame(sp,site,median.slope,median.intercept,ind95,ind80,cont,rls,mass,class,guild,sitesp,nyears,rare80)
  res

}

logRegression <- function(row,n,years){
  #code to perform a logistic regression on each row of data and return coefficients
  # row is a vector with occupancy values for each year
  # n is the number of camera traps at a site
  y <- round(cbind(row*n,(1-row)*n))
  model<-glm(y~years,family="binomial")
  model$coeff
}

calc.lambda <- function(mat) {
  require(combinat)
  cols.lam <- dim(mat)[2]
  cols.lam <- as.matrix(combn(cols.lam,2))
  nc <- length(cols.lam)/2
  lambda<-matrix(0,nr=dim(mat)[1],nc=nc)
  for(i in 1:nc){
    lambda[,i]<- mat[,cols.lam[2,i]]/mat[,cols.lam[1,i]]
  }
  lambda
}
