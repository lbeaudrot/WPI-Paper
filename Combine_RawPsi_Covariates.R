# WPI Analysis with binomial case treated as "rare" categorization
rm(list=ls())
library(reshape)
source("g.test.R")

# Load Summarized WPI output (produced from code below)
#newdata <- read.csv("Species-site-results_Jan_2015_LB.csv")

############ COMBINE ALL VARIABLES TO RAW FILE FOR WPI FIGURE & CREATE A SUMMARIZED FILE FOR WPI ANALYSIS ############

# Read in the last 1000 iterations of the WPI and add column for model_type
WPIbinomial <- read.csv("psi_species_model_binomial_last1000_v1.csv")
WPIbinomial <- cbind(WPIbinomial, model_type=rep("binom", dim(WPIbinomial)[1]))
WPIsimple <- read.csv("psi_species_simplemodel_last1000_v1.csv")
WPIsimple <- cbind(WPIsimple, model_type=rep("simple", dim(WPIsimple)[1]))
WPIconst <- read.csv("psi_species_model_const-phi-gam-p_last1000_v1.csv")
WPIconst <- cbind(WPIconst, model_type=rep("simple", dim(WPIconst)[1]))

WPIall <- rbind(WPIbinomial, WPIsimple, WPIconst)

# Read in the taxonomy data that corresponds with the ids used in the WPI analysis
taxonomy <- read.csv("taxonomy_scientific_name_wpi_20140414_LB2.csv")

# Fill in missing guild and mass values
# Arborophila chloropus - omnivore
# Crypturellus bartletti - omnivore
# Nothocrax urumutum - herbivore (100% fruit diet)
# Trichys fasciculata - omnivore

taxonomy$guild <- ifelse(as.character(taxonomy$bin=="Arborophila chloropus")==TRUE, "Omnivore", 
                  ifelse(as.character(taxonomy$bin=="Crypturellus bartletti")==TRUE, "Omnivore",
                  ifelse(as.character(taxonomy$bin=="Nothocrax urumutum")==TRUE, "Herbivore", 
                  ifelse(as.character(taxonomy$bin=="Trichys fasciculata")==TRUE, "Omnivore", as.character(taxonomy$guild)))))

taxonomy$mass <- ifelse(as.character(taxonomy$bin=="Lophotibis cristata")==TRUE, 850, taxonomy$mass)

# Match taxonomy to species id values
WPIall <- merge(WPIall, taxonomy, by.x="species_id", by.y="id")

# Match site code to site id values
sitekey <- read.csv("sitecode_key.csv")
key <- data.frame(sitecode=as.character(sitekey$sitecode), siteid=sitekey$database_code, cont_long=as.character(sitekey$continent_long))
WPIall <- merge(WPIall, key, by.x="site_id", by.y="siteid")

# Create site.sp variable
site.sp <- paste(WPIall$bin, WPIall$sitecode, sep="-")
WPIall <- data.frame(WPIall, site.sp)


# EXTRACT UNIQUE POPULATIONS for summarized file
WPIunique <- WPIall[duplicated(WPIall$site.sp)==FALSE,]

# Merge ind80 and ind95 data
status <- read.csv(file="Species-site-results_Jan_2015.csv")
status <- data.frame(status[,1:8], nyears=status[,13], site_type=status[,25])
WPIunique <- merge(WPIunique, status, by="site.sp")

# Add new 95% and 80% columns that also have "rare" category
Rare95 <- ifelse(WPIunique$model_type=="binom", "rare", as.character(WPIunique$ind95))
Rare80 <- ifelse(WPIunique$model_type=="binom", "rare", as.character(WPIunique$ind80))
WPIunique <- data.frame(WPIunique, Rare95, Rare80)
NewRare80 <-ifelse(WPIunique$ind80=="decreasing" & WPIunique$Rare80=="rare", "decreasing", ifelse(WPIunique$ind80=="increasing" & WPIunique$Rare80=="rare", "increasing", as.character(WPIunique$Rare80)))
NewRare95 <-ifelse(WPIunique$ind95=="decreasing" & WPIunique$Rare95=="rare", "decreasing", ifelse(WPIunique$ind95=="increasing" & WPIunique$Rare95=="rare", "increasing", as.character(WPIunique$Rare95)))
WPIunique <- cbind(WPIunique, NewRare80, NewRare95)

# Remove excess columns
WPIunique <- WPIunique[,-27] # Remove sp (duplicate of bin)
WPIunique <- WPIunique[,-24] # Remove wpi_include (all are true)
WPIunique <- WPIunique[,-21] # Remove red_list_status_id (rls provides acronyms)
# Remove iteration, year, task_id, job_id, mammal_id, spcrecid
WPIunique <- data.frame(WPIunique[, 1:3], psi=WPIunique$psi, model_type=WPIunique$model_type, WPIunique[,12:36])
WPIunique <- WPIunique[,-19] # Remove site because it duplicates sitecode and has MAN for MAS

# Remove unused site levels
WPIunique$sitecode <- factor(WPIunique$sitecode)



# Add species level hunting data
Hunted <- read.csv(file="SpeciesHuntingDataE.csv") # Classifies all NNN as not hunted; former "no" at BCI to "unknown" to correspond with YAS, UDZ and PSH approaches
Hunted <- Hunted[,3:4]
Hunted2 <- ifelse(Hunted$Hunted=="No", "No", "Yes")
Hunted <- data.frame(Hunted, Hunted2=Hunted2)

WPIdata <- merge(WPIunique, Hunted, by.x="site.sp", by.y="site.sp", all=TRUE)

# Inspect NA values cbind(WPIdata[,1], WPIdata[,29:31])
# There are several mismatches between species hunting survey name list and WPIunique name list
# See file "Known_Species_Counts_10.31.2014Database.csv"
# Populations with WPI output available but missing hunting data

#load("ct_data2014-10-31.gzip")
#load("ct_data2015-03-05")
head(cam_trap_data)
alldata <- cam_trap_data
Site.Code <- substr(alldata$Sampling.Unit.Name,4,6)
alldata <- cbind(alldata, Site.Code)
site.sp <- paste(alldata$Genus, alldata$Species, alldata$Site.Code, sep=" ")
alldata <- cbind(alldata, site.sp)


names <- as.character(unique(alldata$site.sp))
count <- table(alldata$site.sp)
length(names(count))
length(count)
#write.csv(as.data.frame(count), file="Known_Species_Counts_03.05.2015.csv")


# MANUALLY FIX PROBLEMATIC SPECIES
# REMOVE THE FOLLOWING POPULATIONS FOR MISIDENTIFICATIONS
WPIdata <- WPIdata[-523,] # Dasyprocta punctata at YAS misidentification from Hunting Survey; Out of IUCN range
WPIdata <- WPIdata[-520,] # Xenoperdix udzungwensis at UDZ based on misidentification comment from Francesco (row XX); Out of IUCN range
WPIdata <- WPIdata[-287,] # Muntiacus montanus at PSH misidentification; Out of IUCN range
WPIdata <- WPIdata[-132,] # Name Change Needed Dendrohyrax arboreus-UDZ; Rovero says genus should be validus; Exclude because of concern over misidentification 

# UPDATE HUNTING STATUS FOR THE FOLLOWING POPULATIONS
WPIdata$Hunted <- ifelse(WPIdata$site.sp=="Chalcophaps indica-NAK", "Unknown", WPIdata$Hunted) # Change to Unknown bc excluded from survey
WPIdata$Hunted2 <- ifelse(WPIdata$site.sp=="Chalcophaps indica-NAK", "Yes", WPIdata$Hunted2) 

WPIdata$Hunted <- ifelse(WPIdata$site.sp=="Paradoxurus hermaphroditus-BBS", "Unknown", WPIdata$Hunted) # Change to Unknown bc excluded from survey
WPIdata$Hunted2 <- ifelse(WPIdata$site.sp=="Paradoxurus hermaphroditus-BBS", "Yes", WPIdata$Hunted2) 

WPIdata$Hunted <- ifelse(WPIdata$site.sp=="Tragulus javanicus-NAK", "Unknown", WPIdata$Hunted) # Change to Unknown bc excluded from survey
WPIdata$Hunted2 <- ifelse(WPIdata$site.sp=="Tragulus javanicus-NAK", "Yes", WPIdata$Hunted2) 

WPIdata$Hunted <- ifelse(WPIdata$site.sp=="Tragulus napu-PSH", "Unknown", WPIdata$Hunted) # Change to Unknown bc excluded from survey
WPIdata$Hunted2 <- ifelse(WPIdata$site.sp=="Tragulus napu-PSH", "Yes", WPIdata$Hunted2) 

WPIdata$Hunted <- ifelse(WPIdata$site.sp=="Tragulus kanchil-BBS", "Unknown", WPIdata$Hunted) # Change to Unknown bc excluded from survey
WPIdata$Hunted2 <- ifelse(WPIdata$site.sp=="Tragulus kanchil-BBS", "Yes", WPIdata$Hunted2) 

WPIdata$Hunted <- ifelse(WPIdata$site.sp=="Proechimys semispinosus-COU", "No", WPIdata$Hunted) # Change to No because no hunting at Cashu
WPIdata$Hunted2 <- ifelse(WPIdata$site.sp=="Proechimys semispinosus-COU", "No", WPIdata$Hunted2) 

WPIdata$Hunted <- ifelse(WPIdata$site.sp=="Sciurus igniventris-YAS", "Yes", WPIdata$Hunted)  # Change to Yes because marked as hunted on survey
WPIdata$Hunted2 <- ifelse(WPIdata$site.sp=="Sciurus igniventris-YAS", "Yes", WPIdata$Hunted2) 

# REMOVE POPULATIONS LACKING WPI

WPIdata <- WPIdata[is.na(WPIdata$NewRare80)==FALSE,]

# Add site level survey data
# Get ride of extra site categories in survey data
Survey <- read.csv(file="SiteHuntingDataB.csv")
Survey <- Survey[-14,] # Remove VB-BCNP
Survey <- Survey[-3,] # Remove BCI-Soberania

WPIdata <- merge(WPIdata, Survey, by.x="sitecode", by.y="Site.Code", all=TRUE)
WPIdata$Q1 <- as.factor(WPIdata$Q1)
WPIdata$Q2 <- as.factor(WPIdata$Q2)
WPIdata$Q3 <- as.factor(WPIdata$Q3)
WPIdata$Q4 <- as.factor(WPIdata$Q4)
WPIdata$Q5 <- as.factor(WPIdata$Q5)
WPIdata$Q6 <- as.factor(WPIdata$Q6)
WPIdata$Q7 <- as.factor(WPIdata$Q7)

WPI <- WPIdata

#Provide body mass categories
mass_cat <- ifelse(WPI$mass<1001, "1", ifelse(WPI$mass>1000 & WPI$mass<5001, "2", ifelse(WPI$mass>5000, "3", NA)))
mass_cat <- factor(mass_cat)
WPI <- cbind(WPI, mass_cat=mass_cat)

# New FL data
FL <- read.csv(file="20141004_forest_loss.csv")
FL <- melt(FL)
FL <- FL[,-3]
NewFL <- cast(FL, sitecode ~ aoi)

ZOI_cat <- ifelse(NewFL$ZOI<1, "1", ifelse(NewFL$ZOI>=1 & NewFL$ZOI<5, "2", ifelse(NewFL$ZOI>=5, "3", NA)))
ZOIminusPA_cat <- ifelse(NewFL$ZOIminusPA<1, "1", ifelse(NewFL$ZOIminusPA>=1 & NewFL$ZOIminusPA<5, "2", ifelse(NewFL$ZOIminusPA>=5, "3", NA)))
PA_cat <- ifelse(NewFL$PA<0.1, "1", ifelse(NewFL$PA>=0.1 & NewFL$PA<1, "2", ifelse(NewFL$PA>=1, "3", NA)))

NewFL <- cbind(NewFL, ZOIminusPA_cat, PA_cat, ZOI_cat)


# Create Forest Loss Categories

# Merge NewFL with summarized table
WPI <- merge(WPI, NewFL, by="sitecode")

# Write summarized table
write.csv(WPI, file="Species-site-results_March_2015_LB_SurveyData.csv")


# Combine summarized table with raw psi values for Jorge's use in WPI figure (WPIall)

# Reduce WPIall columns for merge to avoid duplicates
head(WPIall)
names(WPIall) # remove covariates already present in WPI
All.reduced <- data.frame(WPIall[,1:7], site.sp=WPIall$site.sp)

# Reduce WPI columns for merge to avoid duplicates
names(WPI) # remove site_id, species_id, psi 
WPI.subset <- data.frame(WPI[,1:2], WPI[,6:47])

# Merge
PsiAll_w_AllCovariates <- merge(All.reduced, WPI.subset, by="site.sp")

write.csv(PsiAll_w_AllCovariates, file="PsiAll_w_AllCovariates.csv")

