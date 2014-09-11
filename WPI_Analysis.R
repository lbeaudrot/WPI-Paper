# Descriptive Results for WPI Paper

# Read in the last 1000 iterations of the WPI for the 3 types of modeling cases
WPIbinomial <- read.csv("psi_species_model_binomial_last1000.csv")
WPIsimple <- read.csv("psi_species_simplemodel_last1000.csv")
WPIconst <- read.csv("psi_species_model_const-phi-gam-p_last1000.csv")

# Identify unique species identifications
binomial <- unique(WPIbinomial$species_id)
const <- unique(WPIconst$species_id)
simple <- unique(WPIsimple$species_id)

# Read in the taxonomy data that corresponds with the ids used in the WPI analysis
taxonomy <- read.csv("taxonomy_scientific_name_wpi_20140414_LB2.csv")

# Create objects to identify species lists for each of 3 cases of modeling for WPI analysis (i.e. naieve (binomial), constant (<8% detection) or (simple) covariate models)
binomial_sp <- taxonomy[match(binomial, taxonomy$id),]
binomial_sp <- cbind(binomial_sp, binomial=rep(1, dim(binomial_sp)[1]))
const_sp <- taxonomy[match(const, taxonomy$id),]
const_sp <- cbind(const_sp, const=rep(1, dim(const_sp)[1]))
simple_sp <- taxonomy[match(simple, taxonomy$id),]
simple_sp <- cbind(simple_sp, simple=rep(1, dim(simple_sp)[1]))


# Add column to the 3 objects above pre-merge so that case can be identified in merged object
Merge1 <- merge(binomial_sp, const_sp, by.x="bin", by.y="bin", all=TRUE)
Merge2 <- merge(Merge1, simple_sp, by.x="bin", by.y="bin", all=TRUE)
ModelType <- data.frame(bin=as.character(Merge2$bin), binomial=Merge2$binomial, const=Merge2$const, simple=Merge2$simple)
ModelType <- ModelType[1:248,]
rownames(ModelType) <- ModelType$bin
ModelType <- ModelType[,2:4]

# Identify number for each of the cases of modeling for WPI analysis (i.e. naieve (binomial), constant (<8% detection) or (simple) covariate models)
colSums(ModelType, na.rm=TRUE)

#Identify the number of species that used 1, 2, or 3 types of models (i.e. for different populations)
table(apply(ModelType, MARGIN=1, FUN=sum, na.rm=TRUE))

# Read in population trend results from output of WPI (removing duplicate columns pre-merge) 
wpidata <- read.csv("Species-site-results_reduced_9Sept2014.csv")
wpidata <- wpidata[,-29]
wpidata <- wpidata[,-28]
wpidata <- wpidata[,-9]

# Merge WPI output data with taxonomic and trait information from taxonomy file (removing duplicate columns pre-merge)
# Use taxonomy file rather than masterlist file because taxonomy file is more up to date for spelling and guild values
taxonomy <- taxonomy[,1:12]
alldata <- merge(wpidata, taxonomy, by.x="sp", by.y="bin", all=FALSE)

# Merge WPI output data with site characteristics from Francesco's survey results
newdata <- read.csv("protection_hunting_TEAM sites.csv")
alldata <- merge(alldata, newdata, by.x="site", by.y="Site", all=TRUE)

# Examine WPI output
table(alldata$site, alldata$class)
length(unique(alldata$sp))
hist(table(alldata$sp))
table(table(alldata$sp))
