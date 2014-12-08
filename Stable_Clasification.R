# Determine whether populations classified as "stable" in early WPI analyses can be considered stable or need to be reclassified as unknown
# Visualize posterior distributions for all TEAM populations

library(denstrip)

d <- read.csv(file="raw_logistic_coeffs.csv")
# Remove populations from Manaus
d <- d[d$site!="MAN",]
# Manually remove 2 species that site managers indicated were incorrect identifications b/c of IUCN range maps
    # Dasyprocta punctata at YAS (row 116) Dasyprocta punctata-YAS
    # Muntiacus montanus at PSH (row 280) Muntiacus montanus-PSH
d <- d[d$site.sp!="Dasyprocta punctata-YAS",]
d <- d[d$site.sp!="Muntiacus montanus-PSH",]
# NB there are 507 populations represented in d, but only 506 in other WPI data. Need to determine additional population. 

# Connect attribute data
# Load main population status data file
    Spdata <- read.csv("Species-site-results_Sep_2014_LB.csv")

    # Remove outdated forest loss data and Francesco's hunting survey results
    Spdata <- data.frame(Spdata[,1:14], site_type=Spdata[,25])

    # Add data on case used to model each population (see object "Cases" from file "WPI_Analysis.R")
    SpdataCase <- merge(Spdata, Cases, by.x="site.sp", by.y="site.sp")

dCase <- merge(d, SpdataCase, by.x="site.sp", by.y="site.sp")
d_simple <- dCase[dCase$Case=="simple",]
d_constant <- dCase[dCase$Case=="constant",]
d_binomial <- dCase[dCase$Case=="binomial",]

d_80inc <- dCase[dCase$ind80=="increasing",]
d_80dec <- dCase[dCase$ind80=="decreasing",]
d_80sta <- dCase[dCase$ind80=="stable",]

# Create a list where each element in the list is the posterior distribution of a single population
# Then loop through list to create denstrip figure
# Sample code from Miguel:
    #You can plot an empty plot and then use
    #denstrip(posterior1,at=.1,width=0.013)
    #denstrip(posterior2,at=.12,width=0.013)
    #you have to play with "at" and "width" to get the structure that you want; "at" will depend on the scale of the y axis.
 
# Use to split the input data into a list based on site.sp
d.list <- dlply(d, "site.sp")
names(d.list) <- unique(d$site.sp)

# Separate out by case
d.listS <- dlply(d_simple, "site.sp")
names(d.listS) <- unique(d_simple$site.sp)

d.listC <- dlply(d_constant, "site.sp")
names(d.listC) <- unique(d_constant$site.sp)

d.listB <- dlply(d_binomial, "site.sp")
names(d.listB) <- unique(d_binomial$site.sp)

# Separate out by status
d.list80Inc <- dlply(d_80inc, "site.sp")
names(d.list80Inc) <- unique(d_80inc$site.sp)

d.list80Dec <- dlply(d_80dec, "site.sp")
names(d.list80Dec) <- unique(d_80dec$site.sp)

d.list80Sta <- dlply(d_80sta, "site.sp")
names(d.list80Sta) <- unique(d_80sta$site.sp)


# Plot using denstrip
# Plot a single density strip
x <- d.list[[1]][,7]
plot(x, xlim=c(-3,3), ylim=c(-3,3), xlab="x", ylab="y", type="n")
denstrip(x, at=0)

# Loop over all elements in list
use.list <- d.listB

plot(x, xlim=c(-3,3), ylim=c(0,length(use.list)), xlab="x", ylab="y", type="n")
for(i in 1:length(use.list)){
  denstrip(use.list[[i]][,7], at=i, width=1)
}

# See example code 
## Lattice example data: heights of singing voice types

bwplot(voice.part ~ height, data=singer, xlab="Height (inches)",
       panel=panel.violin, xlim=c(50,80))
bwplot(voice.part ~ height, data=singer, xlab="Height (inches)",
       panel = function(x, y) {
           xlist <- split(x, factor(y))
           for (i in seq(along=xlist))
               panel.denstrip(x=xlist[[i]], at=i)
       },
       xlim=c(50,80)
       )

# Mimic example code
bwplot(site.sp ~ coeffs, data=d, xlab="Logistic Coefficients",
       panel=panel.violin, xlim=c(-4,4))
bwplot(site.sp ~ coeffs, data=d, xlab="Logistic Coefficients",
       panel = function(x, y) {
           xlist <- split(x, factor(y))
           for (i in seq(along=xlist))
               panel.denstrip(x=xlist[[i]], at=i)
       },
       xlim=c(-4,4)
       )


# Examine posterior distributions based on model case (i.e. 1, 2, or 3) used with expectation that rare species will have widest distributions
# Examine based on current classification of decreasing/stable/increasing
# Examine after ordering them based on median of distribution; add red line to represent median
