PLGS Analysis of TEAM data WPI output

library(ape) 
library(caper) 
library(phytools)
library(dummies)

#load unresolved phylo (FROM JP LESSARD)
phylo<-read.nexus("mammalST_MSW05_all.tre")
phylo<-phylo[[1]]


# Change names to match phylogeny
# "Pardofelis_temminckii” to "Catopuma_temminckii"
# "Phataginus_tricuspis” to "Manis_tricuspis"

Spdata$sp <- gsub("Pardofelis temminckii", "Catopuma temminckii", Spdata$sp)
Spdata$sp <- gsub("Phataginus tricuspis", "Manis tricuspis", Spdata$sp)
Spdata$sp <- gsub("Callosciurus erythraeus", "Callosciurus inornatus", Spdata$sp)
Spdata$sp <- gsub("Caracal aurata", "Caracal caracal", Spdata$sp)
Spdata$sp <- gsub("Eulemur rufifrons", "Eulemur rufus", Spdata$sp)
Spdata$sp <- gsub("Nesotragus moschatus", "Neotragus moschatus", Spdata$sp)
Spdata$sp <- gsub("Smutsia gigantea", "Manis gigantea", Spdata$sp)


# Add dummy variables for categorical predictors

G <- dummy(Spdata$guild)
R <- dummy(Spdata$rls)
S <- dummy(Spdata$site_type)
H <- dummy(Spdata$Hunted)
H2 <- dummy(Spdata$Hunted2)

Spdata <- cbind(Spdata, G, R, S, H, H2)





# Separate mammal and bird populations. Will need to do load another phylogeny and run analysis separately for birds
Mtraits <- Spdata[Spdata$class=="MAMMALIA",]
Btraits <- Spdata[Spdata$class=="AVES",]


# Add column to trait data with species name formatting with "_" instead of space 
spname <- sub(" ", "_", Mtraits$sp)
spname2 <- paste(spname, Mtraits$site, sep="_")
Mtraits <- data.frame(spname=spname, spname2=spname2, Mtraits)


# 1. "Prune" the large mammal tree to include only those in trait data
tree<-multi2di(drop.tip(phylo, phylo$tip.label[!phylo$tip.label %in% Mtraits$spname]))

# Now add in duplicate species to include all populations with corresponding site codes
# Add a new tip for each row in spname2

hold <- list()
tree_pops <- tree

for(i in 1:length(spname2)){
  hold <- spname2[i]
  tree_pops <- add.species.to.genus(tree_pops, hold, where="root")
}

# Then remove original species without site codes
# Drop all tips from spname

tree2 <- drop.tip(tree_pops, tip=spname)


# 2. Look at the phylogeny
plot.phylo(tree2,type="phylogram",cex=0.07, adj=1, label.offset=3)

# Explore traits across phylogeny
plotBranchbyTrait(z$phy, log(z$data$mass), mode="edges", cex=0.07, adj=1, label.offset=3)
plot.phylo(z$phy,type="phylogram",cex=0.07, adj=1, label.offset=3, tip.color=unclass(as.factor(z$data$guild)))
plot.phylo(z$phy,type="phylogram",cex=0.07, adj=1, label.offset=3, tip.color=c("green", "black", "red")[unclass(as.factor(z$data$Hunted))])
plot.phylo(z$phy,type="phylogram",cex=0.07, adj=1, label.offset=3, tip.color=c("green", "red")[unclass(as.factor(z$data$Hunted2))])
plot.phylo(z$phy,type="phylogram",cex=0.07, adj=1, label.offset=3, tip.color=c("red", "green", "black")[unclass(as.factor(z$data$ind80))])


# 3. Set up a "Comparative Data Object" ####

#We use the function "comparative.data' in the caper package to build an object that combines our phylogeny with our dataset.
#Before this stage it is important to make sure that the names of species in your phylogeny match those in the dataset.
#You don't need to worry about this now, as we have done the matching already, but its something to keep in mind for future projects.

z<-comparative.data(tree2,data=Mtraits, names="spname2",na.omit=F)

# This object containts the original data (accesed through "z$data"), the phylogeny (accessed through "z$phy"),
# and some other important information such as any species that were dropped when merging the phylogeny with the trait data (accessed thruogh "z$dropped")
# These "dropped tips" are likely the result of missing trait data for species in the tree, or mismatches between the species names in the tree and the dataset.

# Check to see if any species were dropped
z$dropped


# Data Exploration


plot(median.coeff ~ log(mass), data=z$data, pch=20,xlab= "log Body Mass (g)", ylab="Median Coefficient")
plot(median.coeff ~ guild, data=z$data, pch=20,xlab= "log Body Mass (g)", ylab="Median Coefficient")
plot(median.coeff ~ Hunted2, data=z$data, pch=20,xlab= "log Body Mass (g)", ylab="Median Coefficient")
plot(median.coeff ~ nyears, data=z$data, pch=20,xlab= "log Body Mass (g)", ylab="Median Coefficient")


## PGLS with the Maximum Likelihood estimate of Lambda
pgls1<-pgls(median.coeff ~ log(mass), data=z, lambda="ML")
summary(pgls1)
abline(pgls1, col=3, lwd=2, lty=1)

pgls2<-pgls(median.coeff ~ guild, data=z, lambda="ML")
summary(pgls2)
abline(pgls2, col=3, lwd=2, lty=1)

pgls3<-pgls(median.coeff ~ Hunted2, data=z, lambda="ML")
summary(pgls3)
abline(pgls3, col=3, lwd=2, lty=1)

pgls4<-pgls(median.coeff ~ nyears + Hunted2No + guildCarnivore + guildHerbivore + guildInsectivore, data=z, lambda="ML")
summary(pgls4)
abline(pgls4, col=3, lwd=2, lty=1)


# Examine phylogenetic signal in residuals of regression model following Revell 2010
test <- phyl.resid(z$phy, x=as.vector(z$data$nyears), Y=as.vector(z$data$median.coeff), method="lambda")
phylosig(z$phy, test$resid, method="lambda")

IND80 <- z$data$ind80_num
IND80 <- as.numeric(as.character(IND80))
names(IND80) <- z$phy$tip.label


MC <- as.vector(z$data$median.coeff)
names(MC) <- z$phy$tip.label

NYEARS <- as.vector(z$data$nyears)
names(NYEARS) <- z$phy$tip.label

GUILDS <- as.matrix(cbind(z$data$guildCarnivore, z$data$guildHerbivore, z$data$guildInsectivore))
names(GUILDS) <- z$phy$tip.label

Hunting2 <- as.vector(z$data$Hunted2No)
names(Hunting) <- z$phy$tip.label

test1 <- phyl.resid(z$phy, x=NYEARS, Y=MC, method="lambda")
phylosig(z$phy, test1$resid, method="K")

test2 <- phyl.resid(z$phy, x=GUILDS, Y=MC, method="lambda")
phylosig(z$phy, test2$resid, method="K")

test3 <- phyl.resid(z$phy, x=Hunting2, Y=MC, method="lambda")
phylosig(z$phy, test3$resid, method="K")

test4 <- phyl.resid(z$phy, x=log(z$data$mass), Y=MC, method="lambda")
phylosig(z$phy, test4$resid, method="K")

test5 <- phyl.resid(z$phy, x=cbind(NYEARS, GUILDS, Hunting2), Y=MC, method="lambda")
phylosig(z$phy, test5$resid, method="K")










######## From Max Farrell


## Importing data 

#To conduct comparative methods, one needs a matrix of trait values per species
#and a phylogenetic tree with branch lengths - in our case the branch lengths represent geologic time

# Phylogeny for all mammals complied by Fritz et al.2009 
# Available in the Supplementary Material - Ecol Lett. 2009 (6):538-49. (doi: 10.1111/j.1461-0248.2009.01307.x.)
mammal_tree<-read.tree("mammals.tre")

# Brain and body sizes for 94 mammal species from Navette et al. 2011. Nautre 480, 91–93 (doi:10.1038/nature10629)
# Available from http://www.nature.com/nature/journal/v480/n7375/abs/nature10629.html
# Life history traits from the panTHERIA database
# Jones et al. 2009. Ecology 90:9. (doi: 10.1890/08-1494.1)
# Available from http://www.esapubs.org/Archive/ecol/E090/184/default.htm
traits<-read.csv("traits.csv", header=T)

# Have a look at the dataset to see the format:
# each species is given its own row and traits are in the columns.

head(traits) #or View(traits)

# Aside from the brain mass and body mass measurements, other life history traits have been included for you to explore
# These include Max Longevity in months (lifespan), Neonate Body Mass in grams, Geographical Range Area in kilometres squared, 
# and the mean annual precipitation for the species range in mm.

#########################################################################

#### B. Pruning and Plotting Trees ####

# The data we are working with is only for a few mammal species.
# To see the phylogeny for these species, we have to "prune" the full mammal tree to include only those in the dataset.

# 1. "Prune" the large mammal tree to include only those in "traits$Species"
tree<-multi2di(drop.tip(mammal_tree,mammal_tree$tip.label[!mammal_tree$tip.label %in%traits$Species]))

# 2. Have a look at the phylogeny
plot.phylo(tree,type="phylogram",cex=0.3, adj=1, label.offset=2)

# NOTE: There are lots of fun ways to plot trees 
# try type=  "phylogram" (the default), "fan", "cladogram", "unrooted", or "radial"


##################################################################


#### C. Setting up a "Comparative Data Object" ####

#We use the function "comparative.data' in the caper package to build an object that combines our phylogeny with our dataset.
#Before this stage it is important to make sure that the names of species in your phylogeny match those in the dataset.
#You don't need to worry about this now, as we have done the matching already, but its something to keep in mind for future projects.

z<-comparative.data(tree,data=traits, names="Species",na.omit=F)

# This object containts the original data (accesed through "z$data"), the phylogeny (accessed through "z$phy"),
# and some other important information such as any species that were dropped when merging the phylogeny with the trait data (accessed thruogh "z$dropped")
# These "dropped tips" are likely the result of missing trait data for species in the tree, or mismatches between the species names in the tree and the dataset.

# Check to see if any species were dropped
z$dropped


#################################################################


#### D. Data Exploration ####

# Plot the data
plot(log(Brain_mass_g) ~ log(BodyMass_g), data=z$data, pch=20,xlab= "log Body Mass (g)", ylab="log Brain Mass (g)")

# Add a basic trend line - this regression assumes statistical independence of the data points - i.e. no phylogenetic information
lm_1<-lm(log(Brain_mass_g) ~ log(BodyMass_g), data=z$data)
abline(lm_1, lwd=1.5,col=1)


# Now look at the data with the points coloured by the taxonomix order of each species
plot(log(Brain_mass_g) ~ log(BodyMass_g), data=z$data, col=1:10,pch=20,xlab= "log Body Mass (g)", ylab="log Brain Mass (g)")
legend("bottomright",legend=unique(z$data$Order),col=1:10,cex=0.8,pch=20)

# To give an idea of why non-independence is important we can plot trend lines 
# for a few taxonomic orders separately:

# RODENTIA
lm_rodentia<-lm(log(Brain_mass_g) ~ log(BodyMass_g), data=z$data[z$data$Order=="Rodentia",])
abline(lm_rodentia, lwd=2,col=1)

# LAGOMORPHA
lm_lagomorpha<-lm(log(Brain_mass_g) ~ log(BodyMass_g), data=z$data[z$data$Order=="Lagomorpha",])
abline(lm_lagomorpha, lwd=2,col=2)

# CARNIVORA
lm_carnivora<-lm(log(Brain_mass_g) ~ log(BodyMass_g), data=z$data[z$data$Order=="Carnivora",])
abline(lm_carnivora, lwd=2,col=6)

# PRIMATES
lm_primates<-lm(log(Brain_mass_g) ~ log(BodyMass_g), data=z$data[z$data$Order=="Primates",])
abline(lm_primates, lwd=2,col=3)

# The slopes of these regressions show that the relationship between brain size and body size differs based 
# taxonomic order. This is interesting in its own right, but we are interested in getting an estimate of this
# relationship across all mammals. To do this we need to take into account these differences due to shared
# evolutionary history.


####################################################################


####  E. Phylogenetically Independent Contrasts  ####

# The APE package has a function to automatically compute the contrasts
# Since we are dealing with continuous data, we use the "crunch" function. 
# If you use PIC with categorical data, you should use the "brunch" function.

pic1<-crunch(log(Brain_mass_g)~log(BodyMass_g), data=z)

# You should get an error message that reads 
# "The phylogeny contains either negative or zero branch lengths and crunch has not been set to use equal branch lengths"

# This is because PIC takes the difference in trait values between pairs of species that are the most closely related.
# The tree we are using has a few nodes where more than two species split from the same common ancestor.
# As a result, these species are not paired and the contrasts cannot be computed.
# These are the nodes that have the "zero branch lengths" that the function doesn't like.

# Have a look at the phylogeny and see if you can spot these problem nodes - usually reffered to as "polytomies"
plot.phylo(tree,type="phylogram",cex=0.3, adj=1, label.offset=2)


# There are a few ways to get around this problem, but the easiest is to add a small length to every branch in the tree.
# This effectivley replaces the zero branch lengths with a very small number.
z2<-z
z2$phy$edge.length<-z2$phy$edge.length + 0.00001

pic1<-crunch(log(Brain_mass_g)~log(BodyMass_g), data=z2)
summary(pic1)

#To display the contrast regression graphically:
contrasts1<-caic.table(pic1)
plot(contrasts1[,1] ~ contrasts1[,2], xlab = "Contrasts in log Body Mass (g)", ylab = "Contrasts in log Brain Mass (g)",pch=20)
abline(pic1)

#This indicates there is a strong correlation between changes in body mass and changes in brain mass among species.
#Notice that the x-axis contrasts are positive - this is because contrasts can be in either direction, to so help with
#interpretation of the contrasts, the x-axis is usually forced to be positive by subtracting the smaller x-axis trait value 
#from the larger for each species pair.


######################################################################

#### F. Linear Regression Vs. Phylogenetic Generalized Least Squares ####

## Plot Data
## The trend lines for different models will be drawn on this plot

plot(log(Brain_mass_g) ~ log(BodyMass_g), data=z$data, pch=20,xlab= "log Body Mass (g)", ylab="log Brain Mass (g)")
legend("topleft", c("lm","PGLS Lambda = 0", "PGLS Lambda = 1", "PGLS Lambda = ML"),lty=c(1,3,1,1),lwd=c(1,3,1,1),col=c(1,2,4,3), cex=0.8)


## Linear Model

# Again, this is the typical way comparative analyses used to be conducted - treating all species as independent.
lm_1<-lm(log(Brain_mass_g) ~ log(BodyMass_g), data=z$data)
summary(lm_1)
abline(lm_1, lwd=1.5,col=1)


## Phylogenetic Generalized Least Squares (PGLS)

#NOTE: PGLS is much less demanding that PIC - it can handle unresolved trees (ones with polytomies)
#      With PGLS we don't have to modify the branch lengths like we did before running the PIC.


## PGLS with Lambda set to 0
pgls1<-pgls(log(Brain_mass_g) ~ log(BodyMass_g), data=z, lambda=0.000001) #PGLS doesn't like lambda set exactly to zero
summary(pgls1)
abline(pgls1, col=2, lwd=3,lty=3)

# Q: How does this compare to the basic linear regression (lm_1)?



## PGLS with Lambda set to 1
pgls2<-pgls(log(Brain_mass_g) ~ log(BodyMass_g), data=z, lambda=1)
summary(pgls2)
abline(pgls2, col=4, lwd=2,lty=1)

# Q:How does this compare to the slope found using PICs? Hint: check "summary(pic1)"


## PGLS with the Maximum Likelihood estimate of Lambda
pgls3<-pgls(log(Brain_mass_g) ~ log(BodyMass_g), data=z, lambda="ML")
summary(pgls3)
abline(pgls3, col=3, lwd=2, lty=1)

# Q: What is the maximum likelihood estimate of lambda? 
# Q: What does this tell you about the phylogenetic non-independence of these data?


## Lambda Likelihood Profile

#Generate the likelihood profile for lambda
lambda_profile=pgls.profile(pgls3, which="lambda")

#Plot the likelihood profile
#PGLS automatically determines the maximum likelihood estimate and plots it
#along with the 95% confidence interval.
plot(lambda_profile)


## Model Comparison
#The PGLS function also automatically calulates AIC values to allow for model comparison - accessed through model_name$aic

pgls1$aic
pgls2$aic
pgls3$aic

#Which model is best according to its AIC value?



# If you are interested, other life history traits have been included in "traits.csv"
# Try running the same script but replacing the variables.

