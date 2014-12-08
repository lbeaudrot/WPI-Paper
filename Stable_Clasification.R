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

hist(d.list[[1]][,7])
plot()
denstrip(d.list[[1]][,7], at=1, width=0.1)


# Examine posterior distributions based on model case (i.e. 1, 2, or 3) used with expectation that rare species will have widest distributions

