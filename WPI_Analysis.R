# Descriptive Results for WPI Paper
# Make test change here for commit/push/pull
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

# Examine distribution of species modeled with covariates compared to overall species
# load function from file g.test.R
# Note that function likelihood.test in package "Deducer" and function g.test produce identical results; can use either; g.test will work on server
library(Deducer)
table(simple_sp$guild)
table(taxonomy$guild)
Sprop.table <- cbind(table(simple_sp$guild), table(taxonomy$guild))
g.test(Sprop.table)

table(taxonomy$guild)/sum(table(taxonomy$guild)) # overall proportion of species in each guild
table(simple_sp$guild)/table(taxonomy$guild) # proportion of species in each guild that were modeled with covariates

Cprop.table <- cbind(table(const_sp$guild), table(taxonomy$guild))
g.test(Cprop.table)
table(taxonomy$guild)/sum(table(taxonomy$guild))
table(const_sp$guild)/table(taxonomy$guild)

Bprop.table <- cbind(table(binomial_sp$guild), table(taxonomy$guild))
g.test(Bprop.table)
table(taxonomy$guild)/sum(table(taxonomy$guild))
table(binomial_sp$guild)/table(taxonomy$guild)

# i.e. fewer carnivores (10% less), more herbivores (7% more), ~ same insectivores and fewer omnivores (14% less) were modeled with covariates

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















#############  Analysis with ordinal logistic regression from 24 June 2014 #############

# Analysis of predictors of WPI output
library(ordinal)
library(MuMIn)
library(ggplot2)

# Read in WPI output table and format data for ordinal logistic modeling
#WPI <- read.csv("Species-site-results_reduced.csv")
WPI <- alldata
ind80_num <- ifelse(WPI$ind80=="decreasing", -1, 
                    ifelse(WPI$ind80=="stable", 0, 
                           ifelse(WPI$ind80=="increasing", 1, NA)))
ind80_num <- as.factor(ind80_num)
ind95_num <- ifelse(WPI$ind95=="decreasing", -1, 
                    ifelse(WPI$ind95=="stable", 0, 
                           ifelse(WPI$ind95=="increasing", 1, NA)))
ind95_num <- as.factor(ind95_num)
WPI <- cbind(ind80_num, ind95_num, WPI)
#WPI <- cbind(WPI[,1:17], WPI[,35:39])
names(WPI)
WPI$protection <- as.ordered(WPI$protection)
WPI$Category <- as.ordered(WPI$Category)
WPI$poaching <- as.ordered(WPI$poaching)

# Exclude Pasoh to examine effects of outlier on forest loss effect
# WPI <- subset(WPI, site!="PSH")
# Exclude Volcan Barva to examine effects of potential outlier on number of years
# WPI <- subset(WPI, site!="VB")

# Ordered logistic models

# Examine a null model and random effects for site, continent and species
m0 <- clmm2(ind80_num ~ 1, data=WPI)
m1.site <- clmm2(ind80_num ~ 1, random=site, data=WPI, Hess=TRUE, nAGQ=10)
m1.cont <- clmm2(ind80_num ~ 1, random=cont, data=WPI, Hess=TRUE, nAGQ=10)
m1.sp <- clmm2(ind80_num ~ 1, random=sp, data=WPI, Hess=TRUE, nAGQ=10)
AIC(m0, m1.site, m1.cont, m1.sp)
# Only the site random effect outperforms the null model. Check top model with and without site effect.

# Examine effects of individual predictors 
m0 <- clm(ind80_num ~ 1, data=WPI)
m2 <- clm(ind80_num ~ nyears, data=WPI)
summary(m2)
m3 <- clm(ind80_num ~ protection_level, data=WPI)
summary(m3)
m4 <- clm(ind80_num ~ Category, data=WPI)
summary(m4)
m5 <- clm(ind80_num ~ log(T75_Loss_SampleArea_Pct), data=WPI)
summary(m5)
m6 <- clm(ind80_num ~ log(mass), data=WPI)
summary(m6)
m7 <- clm(ind80_num ~ guild, data=WPI)
summary(m7)
m8 <- clm(ind80_num ~ class, data=WPI)
summary(m8)
m9 <- clm(ind80_num ~ rls, data=WPI[-343,])
summary(m9)
m10 <- clm(ind80_num ~ poaching_level, data=WPI)
summary(m10)
m11 <- clm(ind80_num ~ land_use, data=WPI)
summary(m10)

# Use model selection to compare ranking of individual predictors
SinglePredic.Sel <- model.sel(m0, m2, m3, m4, m5, m6, m7, m8, m9, m10, m11, rank=AIC)
#write.table(SinglePredic.Sel, file="SinglePredic.Sel.csv", sep=",")

# Combine predictors that outperformed the null model into a single model
# Compare model with and without a random effect for site
m11 <- clmm2(ind80_num ~ nyears + protection + T75_Loss_SampleArea_Pct, random=site, data=WPI, Hess=TRUE, nAGQ=10)
summary(m11)
m11.cont <- clmm2(ind80_num ~ nyears + protection + T75_Loss_SampleArea_Pct, random=cont, data=WPI, Hess=TRUE, nAGQ=10)
summary(m11.cont)
m11.sp <- clmm2(ind80_num ~ nyears + protection + T75_Loss_SampleArea_Pct, random=sp, data=WPI, Hess=TRUE, nAGQ=10)
summary(m11.sp)

m12 <- clmm2(ind80_num ~ nyears + protection + T75_Loss_SampleArea_Pct, data=WPI)
summary(m12)
exp(cbind(odds=coef(m12)[3:6], confint(m12)))

m12.95 <- clm(ind95_num ~ nyears + protection + T75_Loss_SampleArea_Pct, data=WPI)
summary(m12.95)
exp(cbind(odds=coef(m12.95)[3:6], confint(m12.95)))

m12clm <- clm(ind80_num ~ nyears + protection + T75_Loss_SampleArea_Pct, data=WPI)
TopPredic.Sel <- model.sel(m0, m2, m3, m5, m12clm, rank=AIC)
#write.table(TopPredic.Sel, file="TopPredic.Sel.csv", sep=",")


# Top model (m12) includes the number of years, protection status and a continuous measure of forest loss but no site random effect
# Results are consisent for status defintions based on 80% (m12) or 95% (m12.95) credible intervals

# Create output table 
m12.table <- round(exp(cbind(odds=coef(m12)[3:6], confint(m12))), 2)
rownames(m12.table) <- c("N years", "Protection (medium)", "Protection (high)", "Forest loss (%)")

m12.95table <- round(exp(cbind(odds=coef(m12.95)[3:6], confint(m12.95))), 2)
rownames(m12.95table) <- c("N years", "Protection (medium)", "Protection (high)", "Forest loss (%)")

#write.table(rbind(m12.table, m12.95table), "m12output.csv", sep=",")

#write.table(summary(m12)[6], "m12estimates.csv", sep=",")
#write.table(summary(m12.95)[5], "m12.95estimates.csv", sep=",")

# Visualize significant predictors

wpi <- read.csv("Species-site-results_reduced.csv")

nyears.raw <- ftable(wpi$ind80, wpi$nyears)
nyears.total <- cbind(
  rep(colSums(nyears.raw)[1], 3),
  rep(colSums(nyears.raw)[2], 3),
  rep(colSums(nyears.raw)[3], 3),
  rep(colSums(nyears.raw)[4], 3),
  rep(colSums(nyears.raw)[5], 3))
nyears.per <- (nyears.raw/nyears.total)*100
barplot(nyears.per, horiz=FALSE)
hold <- rbind(nyears.per[2:3,], nyears.per[1,])
colnames(hold) <- names(table(wpi$nyears))
rownames(hold) <- c("Increasing", "Stable", "Decreasing")
barplot(hold, horiz=TRUE, col=c("green4", "gray", "red4"), ylab="Monitoring Duration (years)")


protec.raw <- ftable(wpi$ind80, wpi$protection)
protec.total <- cbind(
  rep(colSums(protec.raw)[1], 3),
  rep(colSums(protec.raw)[2], 3),
  rep(colSums(protec.raw)[3], 3))
protec.per <- (protec.raw/protec.total)*100
barplot(protec.per, horiz=FALSE)
holda <- rbind(protec.per[2:3,], protec.per[1,])
colnames(holda) <- levels(WPI$protection)
hold1 <- cbind(holda[,3], holda[,2], holda[,1])
colnames(hold1) <- c("Low", "Medium", "High")
rownames(hold1) <- c("Increasing", "Stable", "Decreasing")
barplot(hold1, horiz=TRUE, col=c("green4", "gray", "red4"), ylab="Protection")


categ.raw <- ftable(WPI$ind80, WPI$Category)
categ.total <- cbind(
  rep(colSums(categ.raw)[1], 3),
  rep(colSums(categ.raw)[2], 3),
  rep(colSums(categ.raw)[3], 3),
  rep(colSums(categ.raw)[4], 3))
categ.per <- (categ.raw/categ.total)*100
hold2 <- rbind(categ.per[2:3,], categ.per[1,])
colnames(hold2) <- levels(WPI$Category)
rownames(hold2) <- c("Increasing", "Stable", "Decreasing")
barplot(categ.per, horiz=FALSE)
barplot(hold2, horiz=TRUE, col=c("green4", "gray", "red4"), ylab="Forest Loss Severity")

library(fields)
pdf(file="Barplot_WPIeffects.pdf")
set.panel(3,1)
par(mar=c(1,7,1,1), oma=c(7,3,1,1))
barplot(hold, horiz=TRUE, las=1, mgp=c(5, 1, 0), cex.lab=1.4, col=c("green4", "gray", "red4"), ylab="Monitoring (years)")
barplot(hold2, horiz=TRUE, las=1, mgp=c(5, 1, 0), cex.lab=1.4, col=c("green4", "gray", "red4"), ylab="Forest Loss Severity")
barplot(hold1, horiz=TRUE, las=1, mgp=c(5, 1, 0), cex.lab=1.4, cex.axis=1, col=c("green4", "gray", "red4"), ylab="Protection Level")
mtext("                Populations (%)", side=1, line=3, outer=TRUE, cex=1.5)
mtext("                     Increasing", adj=0, side=1, line=2, outer=TRUE, col="green4")
mtext("Decreasing   ", adj=1, side=1, line=2, outer=TRUE, col="red4")
mtext("Site Characteristics", side=2, line=1, outer=TRUE, cex=1.5)
dev.off()
