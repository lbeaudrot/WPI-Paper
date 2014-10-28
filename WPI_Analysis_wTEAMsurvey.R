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

    # Read in population level hunting survey data
    Hunted <- read.csv(file="SpeciesHuntingData.csv")
    Hunted <- Hunted[,3:4]
    # Merge hunting survey data with main population status file
    Spdata <- merge(Spdata, Hunted, by.x="site.sp", by.y="site.sp")
  
    # Read in site-level forest loss data emailed from Alex on 10/22/2014
    FL <- read.csv("20141021_forest_loss.csv")
    FL <- cast(FL, sitecode ~ aoi, value="loss_pct")

    # Read in site-level survey data
    Survey <- read.csv(file="SiteHuntingData.csv")
    
    # Merge all site-level data
    Sitedata <- merge(FL, Survey, by.x="sitecode", by.y="Site.Code")
    
    # Merge population level and site level data
    WPIdata <- merge(Spdata, Sitedata, by.x="site", by.y="sitecode", all=TRUE)

### Format data for analysis

    # Scale continuous covariates (mass, forest loss)
    WPIdata$mass <- scale(WPIdata$mass)
    WPIdata$CSA <- scale(WPIdata$CSA)
    WPIdata$PA <- scale(WPIdata$PA)
    WPIdata$ZOI <- scale(WPIdata$ZOI)
    WPIdata$ZOIminusPA <- scale(WPIdata$ZOIminusPA)

    # Coerce survey data to ordinal predictors
    WPIdata$Q1 <- as.ordered(WPIdata$Q1)
    WPIdata$Q2 <- as.ordered(WPIdata$Q2)
    WPIdata$Q3 <- as.ordered(WPIdata$Q3)
    WPIdata$Q4 <- as.ordered(WPIdata$Q4)
    WPIdata$Q5 <- as.ordered(WPIdata$Q5)
    WPIdata$Q6 <- as.ordered(WPIdata$Q6)
    WPIdata$Q7 <- as.ordered(WPIdata$Q7)

    # Coerce response variable into ordinal response
    ind80_num <- ifelse(WPIdata$ind80=="decreasing", -1, 
                    ifelse(WPIdata$ind80=="stable", 0, 
                           ifelse(WPIdata$ind80=="increasing", 1, NA)))
    ind80_num <- as.factor(ind80_num)
    ind95_num <- ifelse(WPIdata$ind95=="decreasing", -1, 
                    ifelse(WPIdata$ind95=="stable", 0, 
                           ifelse(WPIdata$ind95=="increasing", 1, NA)))
    ind95_num <- as.factor(ind95_num)
    WPI <- cbind(ind80_num, ind95_num, WPIdata) 


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
m4 <- clm(ind80_num ~ rls, data=WPI)
summary(m4)
m5 <- clm(ind80_num ~ Hunted, data=WPI)
summary(m5)

# Site attribute models
m6 <- clm(ind80_num ~ nyears, data=WPI)
summary(m6)
m7 <- clm(ind80_num ~ Q2, data=WPI)
summary(m7)
m8 <- clm(ind80_num ~ Q3, data=WPI)
summary(m8)
m9 <- clm(ind80_num ~ Q6, data=WPI)
summary(m9)
m10 <- clm(ind80_num ~ PA, data=WPI)
summary(m10)
m11 <- clm(ind80_num ~ ZOIminusPA, data=WPI)
summary(m1)
m12 <- clm(ind80_num ~ site_type, data=WPI)
summary(m12)
m13 <- clm(ind80_num ~ cont, data=WPI)
summary(m13)
m14 <- clm(ind80_num ~ site, data=WPI)
summary(m14)

SinglePredic.Sel <- model.sel(m0, m1, m2, m3, m4, m5, m6, m7, m8, m9, m10, m11, m12, m13, m14, rank=AIC)

m15.site <- clmm2(ind80_num ~ nyears + Q2, random=site, data=WPI, Hess=TRUE, nAGQ=10)
summary(m15.site)

