# Conduct regional defaunation analysis request by Tim O'Brien

MoL <- read.csv(file="Map_of_Life_Data_ALL.csv")

# Need a list of all species monitored by the WPI for each site that also contains Family data

TEAM <- read.csv(file="TEAM_Species_Site_Summaries.csv")

# Then select species from the Map of Life data for each site that are in the same families as the WPI

subset(splist,splist$Unique_Name %in% sitelist & splist$Include==1)

TEAM_fams <- subset(MoL,MoL$Family %in% TEAM$Family)

# Calculate % of regional species pool monitored by TEAM at each site
Region_Pool <- as.data.frame(table(TEAM_fams$Site))
names(Region_Pool) <- c("Site.Code", "SpPool")
Region_Percent <- as.data.frame((table(TEAM$SITE)/table(TEAM_fams$Site))*100)
names(Region_Percent) <- c("Site.Code", "Region_Per")
#write.csv(Region_Percent, file="Region_Percent_Monitored.csv", row.names=FALSE)

TEAM_covs <- read.csv(file="TEAM_site_covariates.csv")
TEAM_covs <- merge(TEAM_covs, Region_Percent, by.x="sitecode", by.y="Site.Code")


# Use observed species richness
boxplot(TEAM_covs$Region_Per ~ TEAM_covs$protection_level, ylab="Percent Regional Pool (observed)", xlab="Protection Level")
summary(lm(TEAM_covs$Region_Per ~ TEAM_covs$protection_level))

# Use estimated species richness 
CT_averages <- read.csv(file="CTaverages_overall.csv")
TEAM_covs <- merge(TEAM_covs, CT_averages, by.x="sitecode", by.y="Site.Code")
TEAM_covs <- merge(TEAM_covs, Region_Pool, by.x="sitecode", by.y="Site.Code")
boxplot((TEAM_covs$CT.median/TEAM_covs$SpPool)*100 ~ TEAM_covs$protection_level, ylab="Percent Regional Pool (estimated)", xlab="Protection Level")
summary(lm((TEAM_covs$CT.median/TEAM_covs$SpPool)*100 ~ TEAM_covs$protection_level))

