# Conduct regional defaunation analysis suggested by Tim O'Brien
library(fields)

#Bring in data from multiple sources
# First, all data extracted from 50 km radius of TEAM site based on Map of Life IUCN range maps for all mammals and birds
MoL <- read.csv(file="Map_of_Life_Data_ALL.csv")

# Second, a list of all species monitored by the WPI for each site that also contains Family data
TEAM <- read.csv(file="TEAM_Species_Site_Summaries.csv")

# Third, select species from the Map of Life data for each site that are in the same families as the WPI
TEAM_fams <- subset(MoL,MoL$Family %in% TEAM$Family)

# Calculate % of regional species pool monitored by TEAM at each site based on monitored species
Region_Pool <- as.data.frame(table(TEAM_fams$Site))
names(Region_Pool) <- c("Site.Code", "SpPool")
Region_Percent_obs <- as.data.frame((table(TEAM$SITE)/table(TEAM_fams$Site))*100)
names(Region_Percent_obs) <- c("Site.Code", "Region_Per_obs")
#write.csv(Region_Percent, file="Region_Percent_Monitored.csv", row.names=FALSE)

TEAM_covs <- read.csv(file="TEAM_site_covariates.csv")
TEAM_covs <- merge(TEAM_covs, Region_Percent_obs, by.x="sitecode", by.y="Site.Code")
 

# Bring in estimated species richness and calculate regional species pool based on estimates
CT_averages <- read.csv(file="CTaverages_overall.csv")
Region_Per_est <- as.data.frame((CT_averages$CT.median/Region_Pool$SpPool)*100)
names(Region_Per_est) <- c("Region_Per_est")
TEAM_covs <- merge(TEAM_covs, CT_averages, by.x="sitecode", by.y="Site.Code")
TEAM_covs <- merge(TEAM_covs, Region_Pool, by.x="sitecode", by.y="Site.Code")
TEAM_covs <- data.frame(TEAM_covs, Region_Per_est=Region_Per_est)

#### Examine relationship between coverage of regional species pool and protection level
# Use observed species richness
boxplot(TEAM_covs$Region_Per_obs ~ TEAM_covs$protection_level, ylab="Percent Regional Pool (observed)", xlab="Protection Level")
summary(lm(TEAM_covs$Region_Per_obs ~ TEAM_covs$protection_level))

# Use estimated species richness
boxplot(TEAM_covs$Region_Per_est ~ TEAM_covs$protection_level, ylab="Percent Regional Pool (estimated)", xlab="Protection Level")
summary(lm(TEAM_covs$Region_Per_est ~ TEAM_covs$protection_level))



#### Examine relationship between coverage of regional species pool and site type
# Use observed species richness
boxplot(TEAM_covs$Region_Per_obs ~ TEAM_covs$site_type, ylab="Percent of Regional Pool (observed richness)", xlab="Site Type")
summary(lm(TEAM_covs$Region_Per_obs ~ TEAM_covs$site_type))

# Use estimated species richness
boxplot(TEAM_covs$Region_Per_est ~ TEAM_covs$site_type, ylab="Percent of Regional Pool (estimated richness)", xlab="Site Type")
summary(lm(TEAM_covs$Region_Per_est ~ TEAM_covs$site_type))


# Check PA areas for species richness effects
Areas <- read.csv("GFC_Forest_Change_Summary.csv")
TEAM_covs <- merge(TEAM_covs, Areas, by.x="sitecode", by.y="Site_Code")
boxplot(PA_ha~protection_level, data=TEAM_covs, xlab="Protection Level", ylab="Protected Area Size (ha)")
boxplot(PA_ha~site_type, data=TEAM_covs, xlab="Site Type", ylab="Protected Area Size (ha)")


summary(lm(TEAM_covs$Region_Per_est ~ TEAM_covs$site_type + log(TEAM_covs$PA_ha)))
area_rm <- resid(lm(TEAM_covs$Region_Per_est ~ TEAM_covs$PA_ha))

set.panel(2,4)
boxplot(area_rm ~ TEAM_covs$protection_level, xlab="Protection Level", ylab="Residuals")
summary(lm(area_rm ~ TEAM_covs$protection_level))
boxplot(area_rm ~ TEAM_covs$site_type, xlab="Site Type", ylab="Residuals")
summary(lm(area_rm ~ TEAM_covs$site_type))

