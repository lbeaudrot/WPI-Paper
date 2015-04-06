# Load data from WPI system to create figures and analyses
# for the paper

source("functions_WPI.R")
library(reshape2)
library(plyr)
library(dplyr)
library(ggplot2)
library(grid)

# read in the raw data produced by the vertica system
file <- "Corrected Psi For Jorge v2//PsiAll_w_AllCovariates.csv"

#file_bin <- "PsiByModel//psi_species_model_binomial_last1000_v1.csv"
#file_cons <- "PsiByModel//psi_species_model_const-phi-gam-p_last1000_v1.csv"
#file_full <- "PsiByModel//psi_species_simplemodel_last1000_v1.csv"

data <- read.csv(file,h=T)
#d_cons <- read.csv(file_cons,h=T)
#d_full <- read.csv(file_full,h=T)

#d_bin[,"model_type"] <- "binom"
#d_cons[,"model_type"] <- "const"
#d_full[,"model_type"] <- "full"

#data1 <- rbind(d_bin,d_cons,d_full)
#data1 <- read.csv(file,h=T)

# remove iterations 1001 and 1002 from the data (they dont play nicely)
data <- tbl_df(data)
#data$iteration <- as.numeric(as.character(data1$iteration))
data <- filter(data, iteration <= 1000)

#site.key <- read.csv("sitekey.csv",h=T)
#site.key <- site.key[,-1]

# create a new column for species name
#rl.table <- read.csv("taxonomy_red_list_status_20140414.csv",head=T)
#sp.table <- read.csv("taxonomy_scientific_name_wpi_20140414.csv",head=T)
#sp.table[,"bin"] <- paste(sp.table$genus, sp.table$species, sep=" ") 

#get species id
#data1[,"bin"] <- mapvalues(x=data1$species_id,from=sp.table$id,to=sp.table$bin)

#get guild
#data1[,"guild"]<-mapvalues(x=data1$bin,from=sp.table$bin,to=as.character(sp.table$guild))

#get red list species status
#qwe <- mapvalues(x=data1$bin,from=sp.table$bin,to=as.numeric(sp.table$red_list_status_id))
#qwe2 <- mapvalues(x=qwe,from=rl.table$id,to=as.character(rl.table$code))
#data1[,"rls"] <- qwe2

#get body mass
#data1[,"mass"] <- as.numeric(mapvalues(x=data1$bin,from=sp.table$bin,to=sp.table$mass))


#get site name
#site.key[,"site.abb"] <- c('CAX','PSH','CSN','MAN','VB','BIF','BBS','KRP','YAN','COU','UDZ','BCI','NAK','RNF','YAS','NNN')
#data1[,"site_name"] <- mapvalues(x=data1$site_id,from=site.key$site_id,to=as.character(site.key$site_name)) 
#data1[,"site.abb"] <- mapvalues(x=data1$site_id,from=site.key$site_id,to=as.character(site.key$site.abb))

#get continent
#cont <- c('LA','AS','LA','LA','LA','AF','AS','AF','LA','LA','AF','LA','AS','AF','LA','AF')
#site.key[,"cont"] <- cont
#data1[,"cont"] <- mapvalues(x=data1$site_id,from=site.key$site_id,to=as.character(site.key$cont))

#get class
#data1[,"class"] <- mapvalues(x=data1$bin, from=sp.table$bin, to=as.character(sp.table$class))

#get hunting and protection data
#hunt <- read.csv("summary_hunting.csv",h=T)
#names(hunt)[3] <- "orig_poach"
#hunt[,"poaching"] <- c(1,1,1,1,1,2,2,1,2,1,2,1,1,1,1,2)
#hunt[,"poaching_level"]<-ifelse(hunt$poaching==1,"low","high")
#hunt[,"protection_level"] <-ifelse(hunt$protection == 1, "good",ifelse(hunt$protection==2,"fair","poor"))
#hunt[,"site_type"] <- c("settled","extractive","settled","remote","remote","remote","extractive","extractive","settled","remote","settled","settled","extractive","settled","remote","remote")

#data1[,"poaching"] <- mapvalues(x=data1$site.abb, from=hunt$Site, to=hunt$poaching_level)

#data1[,"protection"] <- mapvalues(x=data1$site.abb, from=hunt$Site, to=hunt$protection_level)

#data1[,"sitetype"] <- mapvalues(x=data1$site.abb, from=hunt$Site, to=hunt$site_type)


#data1 <- read.csv("Merged Files_Jan15//Raw_psi_results_Jan_2015_LB_Rare.csv",h=T)

# create a list object to store all the attributes of a species at a site
spOb2<-{list(occMat=matrix(NA,nr=1000,nc=8,dimnames=list(NULL,2007:2014)),
            guild=character(),
            site=character(),
            continent=character(),
            rls=character(),
            class=character(),
            spname=character(),
            protection=character(),
            poaching=character(),
            sitetype=character(),
            popstatus=character(),
            crit="all")}

# figure out how many unique site/species combinations are there
un<-unique(data[,c("sitecode","bin")])
# create this number of lists and put them in a list called allsp
allsp2<-rep(list(spOb2),dim(un)[1])

# Now create a loop that goes through each unique site species combos
# and populate the information - can take 1-2 mins
for(i in 1:dim(un)[1]){
  #i<-1
  un1<-subset(data2,site_name==un[i,1] & bin==un[i,2],drop=T)
  # extract the occupancy matrix
  y<-unique(un1$year)
  a<-dcast(data=un1,formula=iteration~year,value.var='psi')
  a<-a[,-1]
  a<-as.matrix(a)
  allsp2[[i]][[1]][,colnames(a)]<-a
  #site,country,continent,red list, etc info
  allsp2[[i]]$site<-as.character(unique(un1$sitecode))
  #allsp[[i]]$country<-as.character(unique(un1$country_name))
  allsp2[[i]]$continent<-as.character(unique(un1$cont))
  allsp2[[i]]$rls<-as.character(unique(un1$rls))
  allsp2[[i]]$spname<-as.character(unique(un1$bin))
  allsp2[[i]]$guild<-as.character(unique(un1$guild))
  allsp2[[i]]$protection <- as.character(unique(un1$protection))
  allsp2[[i]]$poaching <- as.character(unique(un1$poaching))
  allsp2[[i]]$sitetype <- as.character(unique(un1$sitetype))
  allsp2[[i]]$class <- as.character(unique(un1$class))
  allsp2[[i]]$popstatus <- as.character(unique(un1$Rare80))
  print(paste("Done with i=",i))
  
}

# Create a data frame to populate the site names and startyears--------
list_Sites_StartYears<-data.frame(sites=sort(as.character(unique(un[,1]))),startyear=1)
list_Sites_StartYears[,2]<-c(4,3,3,4,2,5,5,4,3,3,5,4,3,1,5,5)

# Do a loop to create an object that contains 16 objects (each for one site). 
# Inside each of these objects is a list with all the elements of a species
siteob2<-list()
for(i in 1:length(list_Sites_StartYears[,1])){
  temp<-sapply(allsp2, function(x){x$site==as.character(list_Sites_StartYears[i,1])})
  temp<-subset(allsp2,temp==T)
  siteob2[[i]]<-list(temp)
  names(siteob2[[i]])<-as.character(list_Sites_StartYears[i,1])
  siteob2[[i]]$startyear<-list_Sites_StartYears[i,2]
}

# create vectors with all the categories that we would want to split
# these names must match the names inside of the species object list (spOb)
site <- unique(data$sitecode)
spname<-unique(data$bin)
continent <- unique(data$cont)
guild <- unique(data$guild)
rls <- unique(data$rls)
crit<-"all"
protection <- unique(data$protection)
poaching <- unique(data$poaching)
sitetype <- unique(data$sitetype)
class <- unique(data$class)
popstatus <- unique(data$Rare80)
#data1[,"site.sp"] <-paste(data1$bin,"-",data1$site.abb,sep="")

#Do logistic regressions on the modeled occupancy
sp.status2<-ddply(.data=data,.variables="site.sp",.fun=speciesChangesLogistic,.parallel=T)

#Do logistic regressions and get the coefficients for each simulation each population
rawCoeffs<-ddply(.data=data,.variables="site.sp",.fun=speciesChangesLogistic2,.parallel=T)


#get the Land use change data
lucdata<-read.csv("LUC_data.csv",h=T)
#match the site names with the sp.status data
qwe<-match(sp.status2$site,lucdata$Site_Code)
luc<-lucdata[qwe,]
sp.status2 <- cbind(sp.status2,luc[,-c(1:4)])

# load hunting data

qwe <- match(sp.status2$site,hunt$Site)
sp.status2 <- cbind(sp.status2,hunt[qwe,-(1:4)])

write.csv(sp.status2,"Species-site-results_Jan_2015.csv",row.names=F)

##################################################################
#
# Doing some figures for the paper
# Palette for color blind people
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# SUMMARIZE THE DATA ACCORDING TO DIFFERENT CRITERIA
w2 <- SplitdataCalculateWPI(siteob2,crit,nsites=16,iters=1000,years=2007:2014)
pr <- SplitdataCalculateWPI(siteob,protection,nsites=16,iters=1000,years=2007:2014)
h <- SplitdataCalculateWPI(siteob,poaching,nsites=16,iters=1000,years=2007:2014)
st2 <- SplitdataCalculateWPI(siteob2,sitetype,nsites=16,iters=1000,years=2007:2014)
g <- SplitdataCalculateWPI(siteob,guild,nsites=16,iters=1000,years=2007:2014)
g$wpi <- filter(g$wpi, variable!="V5")
cl <- SplitdataCalculateWPI(siteob,class,nsites=16,iters=1000,years=2007:2014)
ct <- SplitdataCalculateWPI(siteob,continent,nsites=16,iters=1000,years=2007:2014)
rl <- SplitdataCalculateWPI(siteob,rls,nsites=16,iters=1000,years=2007:2014)
ps <- SplitdataCalculateWPI(siteob2,popstatus,nsites=16,iters=1000,years=2007:2014)
sit <- SplitdataCalculateWPI(siteob2,site,nsites=16,iters=1000,years=2007:2014)
#graph the results
#Global WPI
library(ggplot2)
w$wpi$year<-as.numeric(as.character(x$wpi$year))
w$wpi[1,-c(1,2)]<-1
p <- ggplot(data=w2$wpi, aes(x=year))
p<-p+geom_line(aes(y=median),size=1.5)
p<-p+geom_ribbon(aes(ymin=lo50,ymax=hi50),alpha=0.2,fill="blue")+geom_ribbon(aes(ymin=lo80,ymax=hi80),alpha=0.1,fill="blue")+xlab("")+ylab("WPI")+geom_hline(yintercept=1,size=0.5,linetype=2)+theme_bw()



WPIFig <- p

#WPI by protection
p <- ggplot(data=pr$wpi, aes(x=year))
p<-p+geom_line(aes(y=median,color=variable),size=2)
p<-p+geom_ribbon(aes(ymin=lo80,ymax=hi80,fill=variable),alpha=0.1)+xlab("Year")+ylab("Wildlife Picture Index")+labs(title="")+geom_hline(yintercept=1,size=0.5,linetype=2) + theme_bw() + ggtitle("Protection Level") + guides(color=guide_legend(NULL),fill=NULL) + scale_fill_discrete(guide = 'none') + theme(legend.position=c(.9,.85), legend.background=element_rect(fill=NA,colour=NA),legend.key.size=unit(0.3,"cm"))

ProtectionFig <- p

#WPI by hunting

p <- ggplot(data=h$wpi, aes(x=year))
p<-p+geom_line(aes(y=median,color=variable),size=2)
p<-p+geom_ribbon(aes(ymin=lo80,ymax=hi80,fill=variable),alpha=0.1)+xlab("Year")+ylab("Wildlife Picture Index")+labs(title="")+geom_hline(yintercept=1,size=0.5,linetype=2) + theme_bw() + ggtitle("Hunting level") + guides(color=guide_legend(NULL),fill=NULL) + scale_fill_discrete(guide = 'none') +  theme(legend.position=c(.9,.85), legend.background=element_rect(fill="white",colour=NA),legend.key.size=unit(0.3,"cm"))

PoachingFig <- p

# by SiteType

p <- ggplot(data=st2$wpi, aes(x=year))
p<-p+geom_line(aes(y=median,color=variable),size=2)
p<-p+geom_ribbon(aes(ymin=lo80,ymax=hi80,fill=variable),alpha=0.1)+xlab("Year")+ylab("Wildlife Picture Index")+labs(title="")+geom_hline(yintercept=1,size=0.5,linetype=2) + theme_bw() + ggtitle("Site type") + guides(color=guide_legend(NULL)) + scale_fill_discrete(guide = 'none') +  theme(legend.position=c(.85,.9), legend.background=element_rect(fill=NA,colour=NA),legend.key.size=unit(0.3,"cm"))

sitetypeFig <- p

# by guild

p <- ggplot(data=g$wpi, aes(x=year))
p<-p+geom_line(aes(y=median,color=variable),size=2)
p<-p+geom_ribbon(aes(ymin=lo80,ymax=hi80,fill=variable),alpha=0.2)+xlab("Year")+ylab("Wildlife Picture Index")+labs(title="Guild")+geom_hline(yintercept=1,size=0.5,linetype=2) + theme_bw() +  guides(color=guide_legend(NULL)) + scale_fill_discrete(guide = 'none') +  theme(legend.position=c(.85,.8), legend.background=element_rect(fill=NA,colour=NA),legend.key.size=unit(0.3,"cm"))

guildFig <- p

#by birds and mammals

p <- ggplot(data=cl$wpi, aes(x=year))
p<-p+geom_line(aes(y=median,colour=variable),size=2)
p<-p+geom_ribbon(aes(ymin=lo80,ymax=hi80,fill=variable),alpha=0.2)+xlab("Year")+ylab("Wildlife Picture Index")+labs(title="Class")+geom_hline(yintercept=1,size=0.5,linetype=2) +  theme_bw() +  guides(color=guide_legend(NULL)) + scale_fill_discrete(guide = 'none') +  theme(legend.position=c(.8,.85), legend.background=element_rect(fill="white",colour=NA),legend.key.size=unit(0.3,"cm"))

classFig <- p

#by continent

p <- ggplot(data=ct$wpi, aes(x=year))
p<-p+geom_line(aes(y=median,colour=variable),size=2)
p<-p+geom_ribbon(aes(ymin=lo80,ymax=hi80,fill=variable),alpha=0.15)+xlab("Year")+ylab("Wildlife Picture Index")+labs(title="Continent")+geom_hline(yintercept=1,size=0.5,linetype=2) + theme_bw() +  guides(color=guide_legend(NULL)) + scale_fill_discrete(guide = 'none') +  theme(legend.position=c(.9,.9), legend.background=element_rect(fill=NA,colour=NA),legend.key.size=unit(0.3,"cm"))

continentFig <- p

#by rls

p <- ggplot(data=rl$wpi, aes(x=year))
p<-p+geom_line(aes(y=median,colour=variable),size=2)
p<-p+geom_ribbon(aes(ymin=lo80,ymax=hi80,fill=variable,colour=variable),alpha=0.1)+xlab("Year")+ylab("Wildlife Picture Index")+labs(title="Red List Status")+geom_hline(yintercept=1,size=0.5,linetype=2) + theme_bw() + guides(color=guide_legend(NULL)) + scale_fill_discrete(guide = 'none')



rlsFig <- p

#by population status
p <- ggplot(data=ps$wpi, aes(x=year))
p<-p+geom_line(aes(y=median,colour=variable),size=2)
p<-p+geom_ribbon(aes(ymin=lo80,ymax=hi80,fill=variable,colour=variable),alpha=0.1)+xlab("Year")+ylab("Wildlife Picture Index")+labs(title="Population Status")+geom_hline(yintercept=1,size=0.5,linetype=2) + theme_bw() + guides(color=guide_legend(NULL)) + scale_fill_discrete(guide = 'none')
psFig <- p


#by site
p <- ggplot(data=sit$wpi, aes(x=year))
p<-p+geom_line(aes(y=median),size=1) + facet_wrap(~variable)
p<-p+xlab("Year")+ylab("Wildlife Picture Index")+labs(title="Site")+geom_hline(yintercept=1,size=0.5,linetype=2) + theme_bw() + guides(color=guide_legend(NULL)) + scale_fill_discrete(guide = 'none')
psFig <- p


# Function to create logistic regression curves for each species

logistic_function <- function(intercept,slope,nyears){
  1/(1+exp(-(intercept+slope*(1:nyears))))
}

# Test the function
logistic_function(-3.29,0.60,3)
results <- matrix(NA, nrow=dim(sp.status2)[1],ncol=max(sp.status2$nyears))
for (i in 1: dim(sp.status2)[1]){
  with(sp.status2,
       results[i,1:nyears[i]]<<- logistic_function(median.intercept[i], median.slope[i], nyears[i])
  )
}

sp.status2 <- data.frame(sp.status2,results)
names(sp.status2)[16:22]<-c('y1','y2','y3','y4','y5','y6','y7')
a<-melt(data = sp.status2,varnames = c(y1,y2,y3,y4,y5,y6,y7))
y1<-filter(a,filter = variable == "y1")
y2<-filter(a,filter = variable == "y2")
y3<-filter(a,filter = variable == "y3")
y4<-filter(a,filter = variable == "y4")
y5<-filter(a,filter = variable == "y5")
y6<-filter(a,filter = variable == "y6")
y7<-filter(a,filter = variable == "y7")
a <- rbind(y1,y2,y3,y4,y5,y6,y7)
a[,"year"] <- numeric()
a[which(a$variable=="y1"),"year"] <- 1
a[which(a$variable=="y2"),"year"] <- 2
a[which(a$variable=="y3"),"year"] <- 3
a[which(a$variable=="y4"),"year"] <- 4
a[which(a$variable=="y5"),"year"] <- 5
a[which(a$variable=="y6"),"year"] <- 6
a[which(a$variable=="y7"),"year"] <- 7

p <- ggplot(a,aes(x=year,y=value))
p <- p + geom_line(aes(color=sitesp,size=ind80),alpha=0.5) + scale_color_discrete(guide='none')
p

#Normalizing the first year
results2 <- results/results[,1]
sp.status2 <- data.frame(sp.status2[,-c(26:32)],results2)
names(sp.status2)[26:32]<-c('y1','y2','y3','y4','y5','y6','y7')
#get rid of sites with two years of data
sp.status3 <- filter(sp.status2,nyears>2)
a<-melt(data = sp.status3,varnames = c(y1,y2,y3,y4,y5,y6,y7))
y1<-filter(a,filter = variable == "y1")
y2<-filter(a,filter = variable == "y2")
y3<-filter(a,filter = variable == "y3")
y4<-filter(a,filter = variable == "y4")
y5<-filter(a,filter = variable == "y5")
y6<-filter(a,filter = variable == "y6")
y7<-filter(a,filter = variable == "y7")
a <- rbind(y1,y2,y3,y4,y5,y6,y7)
a[,"year"] <- numeric()
a[which(a$variable=="y1"),"year"] <- 1
a[which(a$variable=="y2"),"year"] <- 2
a[which(a$variable=="y3"),"year"] <- 3
a[which(a$variable=="y4"),"year"] <- 4
a[which(a$variable=="y5"),"year"] <- 5
a[which(a$variable=="y6"),"year"] <- 6
a[which(a$variable=="y7"),"year"] <- 7

p <- ggplot(a,aes(x=year,y=value))
p <- p + geom_line(aes(group=site.sp,colour=rare80),alpha=0.5,size=1) + scale_color_discrete(guide='none') + guides(color=guide_legend(NULL)) + theme_bw() + xlab("Year") + ylab("Occupancy") + labs(title="Species dynamics") + theme(legend.position=c(.85,.93), legend.background=element_rect(fill=NA,colour=NA),legend.key.size=unit(0.3,"cm"))  #+ scale_y_log10()
p
spDynams <- p

pdf("Figure 2.pdf",width=8,height=11)
pushViewport(viewport(layout=grid.layout(4,2)))
print(WPIFig,vp=viewport(layout.pos.row=1,layout.pos.col=1))
print(spDynams,vp=viewport(layout.pos.row=1,layout.pos.col=2))
print(continentFig,vp=viewport(layout.pos.row=2,layout.pos.col=1))
print(classFig,vp=viewport(layout.pos.row=2,layout.pos.col=2))
print(guildFig,vp=viewport(layout.pos.row=3,layout.pos.col=1))
print(ProtectionFig,vp=viewport(layout.pos.row=3,layout.pos.col=2))
print(PoachingFig,vp=viewport(layout.pos.row=4,layout.pos.col=1))
print(sitetypeFig,vp=viewport(layout.pos.row=4,layout.pos.col=2))

dev.off()
#Figure 2 new version
pdf("Figure 2.1.pdf",width=4,height=11)
pushViewport(viewport(layout=grid.layout(3,1)))
print(WPIFig,vp=viewport(layout.pos.row=1,layout.pos.col=1))
print(spDynams,vp=viewport(layout.pos.row=2,layout.pos.col=1))
#print(continentFig,vp=viewport(layout.pos.row=2,layout.pos.col=1))
#print(classFig,vp=viewport(layout.pos.row=2,layout.pos.col=2))
#print(guildFig,vp=viewport(layout.pos.row=3,layout.pos.col=1))
#print(ProtectionFig,vp=viewport(layout.pos.row=3,layout.pos.col=2))
#print(PoachingFig,vp=viewport(layout.pos.row=4,layout.pos.col=1))
print(sitetypeFig,vp=viewport(layout.pos.row=3,layout.pos.col=1))

dev.off()
