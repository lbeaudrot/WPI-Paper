# Analyze status of TEAM populations using ordinal logistic regression
    # Include TEAM hunted species survey as covariate data
    # Include forest loss for 5 years prior to sampling start as covariate data

library(reshape)
library(ordinal)
library(MuMIn)
library(ggplot2)

################## DATA INPUT #######################
### Load and merge all necessary input data
    # Load main population status data file
    Spdata <- read.csv("Species-site-results_Sep_2014_LB.csv")

    # Remove outdated forest loss data and Francesco's hunting survey results
    Spdata <- data.frame(Spdata[,1:14], site_type=Spdata[,25])
    Spdata <- Spdata[,1:15]

    # If interested in examining Case, bring in object SpdataCase from file Stable_Classification.R and WPI_Analysis.R
    # Spdata <- SpdataCase

    rls2 <- ifelse(Spdata$rls=="DD", "E.DD", ifelse(Spdata$rls=="LC", "D.LC", ifelse(Spdata$rls=="NT", "C.NT", ifelse(Spdata$rls=="VU", "B.VU", ifelse(Spdata$rls=="EN", "A.EN", NA)))))
    Site_cat <- as.factor(ifelse(Spdata$site_type=="remote", "C.Remote", ifelse(Spdata$site_type=="extractive", "B.Extractive", ifelse(Spdata$site_type=="settled", "A.Settled", NA))))
  
    Spdata <- data.frame(Spdata, rls2=rls2, Site_cat=Site_cat)

    # Reclassify YAS site type as extractive
    #Spdata$site_type <- as.factor(ifelse(Spdata$site=="YAS", 1, Spdata$site_type))

    # Read in population level hunting survey data
    #Hunted <- read.csv(file="SpeciesHuntingData.csv")

    # Note alternative hunting data file to include species that are hunted outside of the TEAM core area at COU and YAN
    #Hunted <- read.csv(file="SpeciesHuntingDataB.csv")
    #Hunted <- read.csv(file="SpeciesHuntingDataC.csv") #Classifies all species at BIF, PSH, BCI, UDZ and YAS as hunted b/c of snares & comments; all NNN species as not hunted
    #Hunted <- read.csv(file="SpeciesHuntingDataD.csv") #Classifies all species at NNN as not hunted; all species at BCI & YAS as hunted; all but 1 species at PSH as hunted; adds 2 snared species at BIF
    Hunted <- read.csv(file="SpeciesHuntingDataE.csv") # Classifies all NNN as not hunted; former "no" at BCI to "unknown" to correspond with YAS, UDZ and PSH approaches
    #Hunted <- read.csv(file="SpeciesHuntingDataF.csv") # Classifies all NNN as not hunted; former "no" at BCI and YAS to "yes" based on notes; not PSH b/c of snares
    #Hunted <- read.csv(file="SpeciesHuntingDataG.csv") # Classifies all NNN as not hunted
    Hunted <- Hunted[,3:4]
    Hunted2 <- ifelse(Hunted$Hunted=="No", "No", "Yes")
    Hunted <- data.frame(Hunted, Hunted2=Hunted2)
    # Merge hunting survey data with main population status file
    Spdata <- merge(Spdata, Hunted, by.x="site.sp", by.y="site.sp")
    
    # Manually remove 2 species that site managers indicated were incorrect identifications b/c of IUCN range maps
    # Dasyprocta punctata at YAS (row 116)
    # Muntiacus montanus at PSH (row 280)

    Spdata <- Spdata[-280,]
    Spdata <- Spdata[-116,]

    # Read in site-level forest loss data emailed from Alex on 10/22/2014
    FL <- read.csv("20141021_forest_loss.csv")
    FL <- cast(FL, sitecode ~ aoi, value="loss_pct")

    # Read in site-level survey data
#Survey <- read.csv(file="SiteHuntingData.csv")
    
    # Note alternative site-level survey data that uses BCI-BCNM and VB-La Selva instead of Soberania and BCNP, which have more camera traps
    Survey <- read.csv(file="SiteHuntingDataB.csv")

    # Merge all site-level data
    Sitedata <- merge(FL, Survey, by.x="sitecode", by.y="Site.Code")
    
    # Merge population level and site level data
    WPIdata <- merge(Spdata, Sitedata, by.x="site", by.y="sitecode", all=TRUE)

### Format data for analysis
CV <- function(data){
  sd(data)/mean(data)
}

elev.range <- function(data){
  max(data) - min(data)
}

load("ct_pts_elev.RData")
elevation.data <- ct_pts_elev
names(elevation.data) <- c("Site.Code", "ct_ID", "Elevation")
elevation.mean <- aggregate(elevation.data$Elevation ~ elevation.data$Site.Code, FUN=mean)
elevation.CV <- aggregate(elevation.data$Elevation ~ elevation.data$Site.Code, FUN=CV)[2]
elevation.range <- aggregate(elevation.data$Elevation ~ elevation.data$Site.Code, FUN=elev.range)[2]
elevation <- cbind(elevation.mean, elevation.CV, elevation.range)
colnames(elevation) <- c("Site.Code", "Elev.Mean", "Elev.CV", "Elev.Range")

Area <- read.csv("GFC_Forest_Change_Summary.csv")
names(Area) <- c("Site.Code", "ForestLossSite", "ForestLossZOI", "PA_area", "SA_area", "ZOI_area")
Area <- Area[,-3]
Area <- Area[,-2]

Elev.Area <- merge(elevation, Area, by="Site.Code")

WPIdata <- merge(WPIdata, Elev.Area, by.x="site", by.y="Site.Code")



    # Scale continuous covariates (mass, forest loss)
    WPIdata$mass <- scale(WPIdata$mass)
    WPIdata$CSA <- scale(WPIdata$CSA)
    WPIdata$PA <- scale(WPIdata$PA)
    WPIdata$ZOI <- scale(WPIdata$ZOI)
    WPIdata$ZOIminusPA <- scale(WPIdata$ZOIminusPA)
    WPIdata$Elev.Mean <- scale(WPIdata$Elev.Mean)
    WPIdata$Elev.CV <- scale(WPIdata$Elev.CV)
    WPIdata$Elev.Range <- scale(WPIdata$Elev.Range)
    WPIdata$PA_area <- scale(WPIdata$PA_area)
    WPIdata$SA_area <- scale(WPIdata$SA_area)
    WPIdata$ZOI_area <- scale(WPIdata$ZOI_area)

    # Coerce survey data to factors
    WPIdata$Q1 <- as.factor(WPIdata$Q1)
    WPIdata$Q2 <- as.factor(WPIdata$Q2)
    WPIdata$Q3 <- as.factor(WPIdata$Q3)
    WPIdata$Q4 <- as.factor(WPIdata$Q4)
    WPIdata$Q5 <- as.factor(WPIdata$Q5)
    WPIdata$Q6 <- as.factor(WPIdata$Q6)
    WPIdata$Q7 <- as.factor(WPIdata$Q7)

    # Coerce response variable into ordinal response
    ind80_num <- ifelse(WPIdata$ind80=="decreasing", -1, 
                    ifelse(WPIdata$ind80=="stable", 0, 
                           ifelse(WPIdata$ind80=="increasing", 1, NA)))
    ind80_num <- as.factor(ind80_num)
    ind95_num <- ifelse(WPIdata$ind95=="decreasing", -1, 
                    ifelse(WPIdata$ind95=="stable", 0, 
                           ifelse(WPIdata$ind95=="increasing", 1, NA)))
    ind95_num <- as.factor(ind95_num)
    
    # Explore different combinations of categorical data
    Q1B <- as.factor(ifelse(WPIdata$Q1==2,1,WPIdata$Q1))
    Q1C <- as.factor(ifelse(WPIdata$Q1==3,2,WPIdata$Q1))

    Q2B <- as.factor(ifelse(WPIdata$Q2==5,4,WPIdata$Q2))
    Q2C <- as.factor(ifelse(WPIdata$Q2==5,4,ifelse(WPIdata$Q2==3,2,WPIdata$Q2)))
    Q2D <- as.factor(ifelse(WPIdata$Q2==5,4,ifelse(WPIdata$Q2==2,1,WPIdata$Q2)))
    Q2E <- as.factor(ifelse(WPIdata$Q2==4,3,ifelse(WPIdata$Q2==2,1,WPIdata$Q2)))
    Q2F <- as.factor(ifelse(WPIdata$Q2==4,3,ifelse(WPIdata$Q2==2,1,ifelse(WPIdata$Q2==5,3,WPIdata$Q2))))
    Q2G <- as.factor(ifelse(WPIdata$Q2==3,1,ifelse(WPIdata$Q2==2,1,ifelse(WPIdata$Q2==5,4,WPIdata$Q2))))

    Q3B <- as.factor(ifelse(WPIdata$Q3==4,3,WPIdata$Q3))
    Q3C <- as.factor(ifelse(WPIdata$Q3==3,2,WPIdata$Q3))
    Q3D <- as.factor(ifelse(WPIdata$Q3==3,2,ifelse(WPIdata$Q3==4,2,WPIdata$Q3)))
    Q3E <- as.factor(ifelse(WPIdata$Q3==2,1,ifelse(WPIdata$Q3==4,3,WPIdata$Q3)))
    Q3F <- as.factor(ifelse(WPIdata$Q3==2,1,WPIdata$Q3))
    Q3G <- as.factor(ifelse(WPIdata$Q3==3,1,ifelse(WPIdata$Q3==2,1,WPIdata$Q3)))

    Q6B <- as.factor(ifelse(WPIdata$Q6==4,3,WPIdata$Q6))    
    Q6C <- as.factor(ifelse(WPIdata$Q6==3,2,WPIdata$Q6))    
    Q6D <- as.factor(ifelse(WPIdata$Q6==2,1,WPIdata$Q6))    
    Q6E <- as.factor(ifelse(WPIdata$Q6==2,1,ifelse(WPIdata$Q6==4,3,WPIdata$Q6)))    
    Q6F <- as.factor(ifelse(WPIdata$Q6==3,2,ifelse(WPIdata$Q6==4,2,WPIdata$Q6)))    
    Q6G <- as.factor(ifelse(WPIdata$Q6==3,1,ifelse(WPIdata$Q6==2,1,WPIdata$Q6)))    
    #WPI <- cbind(ind80_num, ind95_num, WPIdata, Q1B, Q1C, Q2B, Q2C, Q2D, Q2E, Q2F, Q2G, Q3B, Q3C, Q3D, Q3E, Q3F, Q3G, Q6B, Q6C, Q6D, Q6E, Q6F, Q6G) 
    
    WPI <- cbind(ind80_num, ind95_num, WPIdata)

################ GOODNESS OF FIT TESTS ######################
table(WPI$Hunted2, WPI$ind80)
g.test(table(WPI$Hunted2, WPI$ind80))
g.test(table(WPI$Hunted2, WPI$ind95))

table(WPI$rls, WPI$ind80)
g.test(table(WPI$rls, WPI$ind80))
g.test(table(WPI$rls, WPI$ind95))

table(WPI$Q1, WPI$ind80)
g.test(table(WPI$Q1, WPI$ind80))
g.test(table(WPI$Q1, WPI$ind95))

table(WPI$Q2, WPI$ind80)
g.test(table(WPI$Q2, WPI$ind80))
g.test(table(WPI$Q2, WPI$ind95))

table(WPI$Q3, WPI$ind80)
g.test(table(WPI$Q3, WPI$ind80))
g.test(table(WPI$Q3, WPI$ind95))

table(WPI$Q6, WPI$ind80)
g.test(table(WPI$Q6, WPI$ind80))
g.test(table(WPI$Q6, WPI$ind95))

table(WPI$site_type, WPI$ind80)
g.test(table(WPI$site_type, WPI$ind80))
g.test(table(WPI$site_type, WPI$ind95))

table(WPI$nyears, WPI$ind80)
g.test(table(WPI$nyears, WPI$ind80))
g.test(table(WPI$nyears, WPI$ind95))

table(WPI$guild, WPI$ind80)
g.test(table(WPI$guild, WPI$ind80))
g.test(table(WPI$guild, WPI$ind95))

table(WPI$Case, WPI$ind80)
g.test(WPI$Case, WPI$ind80)
g.test(WPI$Case, WPI$ind95)
################## ANALYSIS #########################

# Examine a null model and random effects for site, continent and species
m0 <- clmm2(ind80_num ~ 1, data=WPI)
m0.site <- clmm2(ind80_num ~ 1, random=site, data=WPI, Hess=TRUE, nAGQ=10)
m0.cont <- clmm2(ind80_num ~ 1, random=cont, data=WPI, Hess=TRUE, nAGQ=10)
m0.sp <- clmm2(ind80_num ~ 1, random=sp, data=WPI, Hess=TRUE, nAGQ=10)
AIC(m0, m0.site, m0.cont, m0.sp)
# Only the site random effect outperforms the null model. Check top model with and without site effect.

# Examine effects of individual predictors 
# Null model
m0 <- clm(ind80_num ~ 1, data=WPI)

# Species attribute models
m1 <- clm(ind80_num ~ class, data=WPI)
summary(m1)
m2 <- clm(ind80_num ~ mass, data=WPI)
summary(m2)
m3 <- clm(ind80_num ~ guild, data=WPI)
summary(m3)
#m4 <- clm(ind80_num ~ 1, nominal=~rls2, data=WPI)
#summary(m4)
m5 <- clm(ind80_num ~ Hunted, data=WPI)
summary(m5)
m5.2 <- clm(ind80_num ~ Hunted2, data=WPI)
summary(m5.2)

# Site attribute models
m6 <- clm(ind80_num ~ nyears, data=WPI)
summary(m6)
m7 <- clm(ind80_num ~ Q1, data=WPI)
summary(m7)
m8 <- clm(ind80_num ~ Q2, data=WPI)
summary(m8)
m9 <- clm(ind80_num ~ Q3, data=WPI)
summary(m9)
m10 <- clm(ind80_num ~ Q6, data=WPI)
summary(m10)
m11 <- clm(ind80_num ~ PA,  data=WPI) # Fails nominal and scale tests but is continuous, needs scale?
summary(m11)
m12 <- clm(ind80_num ~ ZOIminusPA, data=WPI)
summary(m12)
m13 <- clm(ind80_num ~ Site_cat, data=WPI)
summary(m13)
m14 <- clm(ind80_num ~ cont, data=WPI)
summary(m14)
m15 <- clm(ind80_num ~ site, data=WPI)
summary(m15)
m16 <- clm(ind80_num ~ Elev.CV,  data=WPI)
summary(m16)
m17 <- clm(ind80_num ~ Elev.Range,  data=WPI)
summary(m17)
m18 <- clm(ind80_num ~ PA_area, scale=~PA_area, data=WPI)
summary(m18)
m19 <- clm(ind80_num ~ ZOI_area, data=WPI)
summary(m19)


SinglePredic.Sel <- model.sel(m0, m1, m2, m3, m5, m5.2, m6, m7, m8, m9, m10, m11, m12, m13, m14, m15, m16, m17, m18, m19, rank=AIC)


# Explore combinations of individual predictors that outperformed the null model
m20 <- clm(ind80_num ~ nyears + Site_cat, data=WPI)
summary(m20)
m21 <- clm(ind80_num ~ nyears + site, data=WPI)
summary(m21)
m22 <- clm(ind80_num ~ nyears + Site_cat + Elev.CV, data=WPI)
summary(m22)
m23 <- clm(ind80_num ~ nyears + Elev.CV, data=WPI)
summary(m23)
m24 <- clm(ind80_num ~ nyears, nominal=~Hunted2, data=WPI)
summary(m24)
m25 <- clm(ind80_num ~ nyears + PA_area, data=WPI)
summary(m25)
m26 <- clm(ind80_num ~ nyears, nominal=~cont, data=WPI)
summary(m26)
m27 <- clm(ind80_num ~ nyears, nominal=~Site_cat, data=WPI)
summary(m27)
m28 <- clm(ind80_num ~ nyears + PA, data=WPI)
summary(m28)
m29 <- clm(ind80_num ~ nyears + Q3, data=WPI)
summary(m29)
m30 <- clm(ind80_num ~ nyears + Q1, nominal=~guild, data=WPI)
summary(m30)
m31 <- clm(ind80_num ~ nyears, nominal=~guild + Hunted2, data=WPI)
summary(m31)
#m32 <- clm(ind80_num ~ nyears, nominal=~guild + Hunted, data=WPI)
#summary(m32)

model.sel(m20, m21, m22, m23, m24, m25, m26, m27, m28, m29, m30, m31, m32, m0, rank=AIC)

m33 <- clm(ind80_num ~ nyears + PA_area + Hunted2, data=WPI, link="cloglog")
summary(m33)
m34 <- clm(ind80_num ~ nyears + PA_area + Site_cat, data=WPI)
summary(m34)
m35 <- clm(ind80_num ~ nyears + PA_area + Site_cat + Hunted2, data=WPI)
summary(m35)
m36 <- clm(ind80_num ~ nyears + PA_area + Q1, data=WPI)
summary(m36)

links <- c("logit", "probit", "cloglog", "loglog", "cauchit")
sapply(links, function(link){
  clm(ind80_num ~ nyears + PA_area + Hunted2, data=WPI, link=link)$logLik})

fm4.lgt <- update(m33, link = "logit") ## default
fm4.prt <- update(m33, link = "probit")
fm4.ll <- update(m33, link = "loglog")
fm4.cll <- update(m33, link = "cloglog")
fm4.cct <- update(m33, link = "cauchit")
anova(fm4.lgt, fm4.prt, fm4.ll, fm4.cll, fm4.cct)

# get significance values from list: use summary(model)[[6]]



model.sel(m24, m25, m33, m34, m35)

test <- model.sel(m21, m22, m30)
# Add random effects for site or continent to top model
m33.site <- clmm2(ind80_num ~ nyears + PA_area + Q6, random=site, data=WPI, Hess=TRUE, nAGQ=10)
summary(m33.site)
m33.cont <- clmm2(ind80_num ~ nyears + PA_area + Hunted2, random=cont, data=WPI, Hess=TRUE, nAGQ=10)
summary(m33.cont)

# Test goodness of fit of top model
nominal_test(m33)
scale_test(m33)



x.new1 <- data.frame(site_type="remote", Hunted="Yes", nyears=5)
x.new2 <- data.frame(site_type="remote", Hunted="No", nyears=5)
x.new3 <- data.frame(site_type="remote", Hunted="Unknown", nyears=5)
x.new4 <- data.frame(site_type="extractive", Hunted="Yes", nyears=5)
x.new5 <- data.frame(site_type="extractive", Hunted="No", nyears=5)
x.new6 <- data.frame(site_type="extractive", Hunted="Unknown", nyears=5)
x.new7 <- data.frame(site_type="settled", Hunted="Yes", nyears=5)
x.new8 <- data.frame(site_type="settled", Hunted="No", nyears=5)
x.new9 <- data.frame(site_type="settled", Hunted="Unknown", nyears=5)


predict(m30, x.new1)
?predict
