# Examine population level summary results for all WPI populations
rm(list=ls())
source("g.test.R")


WPI <- read.csv("Species-site-results_March_06_2015_LB_SurveyData.csv")
Threat <- as.factor(ifelse(WPI$rls=="CR" | WPI$rls=="EN" | WPI$rls=="VU", "TH", as.character(WPI$rls)))
WPI <- cbind(WPI, Threat)

WPI.hold <- WPI
WPI.all <- WPI
# Explore significantly changing rare populations (N=77; 47 decreasing; 30 increasing)

RC <- WPI[WPI$NewRare80!="rare" & WPI$Rare80=="rare",]
RC$NewRare80 <- factor(RC$NewRare80)
table(RC$NewRare80)
table(RC$NewRare80, RC$site_type)
table(RC$NewRare80, RC$Hunted)


Rdec <- RC[RC$NewRare80=="decreasing",]
Rinc <- RC[RC$NewRare80=="increasing",]

round((table(Rdec$order_team)/table(WPI$order_team))*100,0)




#WPI <- WPI[WPI$Rare80!="rare",]

WPI <- WPI[WPI$NewRare80!="rare",]
WPI$model_type <- factor(WPI$model_type)
WPI$Rare80 <- factor(WPI$Rare80)
WPI$Rare95 <- factor(WPI$Rare95)
WPI$NewRare80 <- factor(WPI$NewRare80)






g.test(table(WPI$model_type, WPI$NewRare80))
g.test(table(WPI$model_type, WPI$Rare80))
g.test(table(WPI$model_type, WPI$Rare95))

g.test(table(WPI$site_type, WPI$NewRare80))
g.test(table(WPI$site_type, WPI$Rare80))
g.test(table(WPI$site_type, WPI$Rare95))

g.test(table(WPI$nyears, WPI$NewRare80))
g.test(table(WPI$nyears, WPI$Rare80))
g.test(table(WPI$nyears, WPI$Rare95))

g.test(table(WPI$guild, WPI$NewRare80))
g.test(table(WPI$guild, WPI$Rare80))
g.test(table(WPI$guild, WPI$Rare95))

g.test(table(WPI$Q1, WPI$NewRare80))
g.test(table(WPI$Q1, WPI$Rare80))
g.test(table(WPI$Q1, WPI$Rare95))

g.test(table(WPI$Q2, WPI$NewRare80))
g.test(table(WPI$Q2, WPI$Rare80))
g.test(table(WPI$Q2, WPI$Rare95))

g.test(table(WPI$Q3, WPI$NewRare80))
g.test(table(WPI$Q3, WPI$Rare80))
g.test(table(WPI$Q3, WPI$Rare95))

g.test(table(WPI$Q6, WPI$NewRare80))
g.test(table(WPI$Q6, WPI$Rare80))
g.test(table(WPI$Q6, WPI$Rare95))

g.test(table(WPI$cont, WPI$NewRare80))
g.test(table(WPI$cont, WPI$Rare80))
g.test(table(WPI$cont, WPI$Rare95))

g.test(table(WPI$Hunted2, WPI$NewRare80))
g.test(table(WPI$Hunted2, WPI$Rare80))
g.test(table(WPI$Hunted2, WPI$Rare95))

g.test(table(WPI$Hunted, WPI$NewRare80))
g.test(table(WPI$Hunted, WPI$Rare80))
g.test(table(WPI$Hunted, WPI$Rare95))

g.test(table(WPI$rls, WPI$NewRare80))
g.test(table(WPI$rls, WPI$Rare80))
g.test(table(WPI$rls, WPI$Rare95))

g.test(table(WPI$mass_cat, WPI$NewRare80))
g.test(table(WPI$ZOIminusPA_cat, WPI$NewRare80))
g.test(table(WPI$ZOI_cat, WPI$NewRare80))
g.test(table(WPI$PA_cat, WPI$NewRare80))




Rare80num <- as.factor(ifelse(WPI$Rare80=="decreasing", "-1", ifelse(WPI$Rare80=="increasing", "1", ifelse(WPI$Rare80=="stable",0,NA))))
ST_num <- as.factor(ifelse(WPI$site_type=="remote", "1", ifelse(WPI$site_type=="settled", "3", 2)))
WPI <- data.frame(WPI, Rare80num=Rare80num, ST_num=ST_num)
barplot(t(table(WPI$ST_num, WPI$Rare80num)/rowSums(table(WPI$ST_num, WPI$Rare80num))*100), 
        col=c("red3", "gray", "green4"),
        xlab="Site Type", las=1, cex.lab=1,
        ylab="Percent of Populations", cex.axis=1,
        names=c("Remote", "Extractive", "Settled"))
mtext("   Decreasing", side=2, line=0.25, adj=0, col="red3", cex=0.7)
mtext("Stable", side=2, line=0.25, adj=0.5, col="gray55", cex=0.7)
mtext("Increasing   ", side=2, line=0.5, adj=1, col="green4", cex=0.7)


barplot(t(table(WPI$STnum, WPI$Rare80num)), col=c("red", "gray", "green4"), beside=TRUE)

threat <- ifelse(WPI$rls=="CR"|WPI$rls=="EN"|WPI$rls=="VU", "TH", as.character(WPI$rls))
table(threat)
WPI <- data.frame(WPI, threat=threat)

WPI.all <- 


#Per.Rare <- (table(WPI$sitecode, WPI$Rare80)[,3]/table(WPI$sitecode))*100
#Per.Dec <-  (table(WPI$sitecode, WPI$Rare80)[,1]/table(WPI$sitecode[WPI$model_type!="binom"]))*100

Per.Rare <- (table(WPI.all$sitecode, WPI.all$NewRare80)[,3]/table(WPI.all$sitecode))*100
Per.Dec <-  (table(WPI$sitecode, WPI$NewRare80)[,1]/table(WPI$sitecode[WPI$NewRare80!="rare"]))*100
Per.No <- (table(WPI$Hunted, WPI$sitecode)[1,]/table(WPI$sitecode))*100
Per.Unknown <- (table(WPI$Hunted, WPI$sitecode)[2,]/table(WPI$sitecode))*100
Per.Yes <- (table(WPI$Hunted, WPI$sitecode)[3,]/table(WPI$sitecode))*100

Per.PA1 <- (table(WPI$PA_cat, WPI$NewRare80)[1,]/table(WPI$sitecode))*100
Per.PA2
Per.PA3
Per.ZOIminusPA1
Per.ZOIminusPA2
Per.ZOIminusPA3

Type <- c("S", "E", "S", "R", "R", "R", "E", NA, "S", "R", "S", "S", "E", "S", "R", "R", NA, NA)
Type <- Type[is.na(Type)!=TRUE]
NUM <- ifelse(Type=="S", 3, ifelse(Type=="E", 2, ifelse(Type=="R", 1, NA)))
NUM <- NUM[is.na(NUM)!=TRUE]
TYPE <- data.frame(Per.Rare=Per.Rare, Per.Dec=Per.Dec, Per.No, Per.Unknown, Per.Yes, Type=Type, Num=NUM)
#TYPE <- TYPE[1:16,]
#TYPE <- TYPE[-8,]

library(fields)
set.panel(2,1)
par(mar=c(1,4,0,0), oma=c(4,1,1,1))
boxplot(TYPE$Per.Dec.Freq ~ TYPE$Num, las=1, 
        names=c("", "", ""), 
        xlab="", ylab="Declining Populations (%)")
boxplot(TYPE$Per.Rare.Freq ~ TYPE$Num, las=1, 
        names=c("Remote", "Extractive", "Settled"), 
        xlab="Site Type", ylab="Rare Populations (%)")
mtext("Site Type", side=1, line=2, outer=TRUE, adj=0.6, cex=1.5)

library(vioplot)
S.Per.Rare <- TYPE$Per.Rare.Freq[TYPE$Type=="S"]
E.Per.Rare <- TYPE$Per.Rare.Freq[TYPE$Type=="E"]
R.Per.Rare <- TYPE$Per.Rare.Freq[TYPE$Type=="R"]

S.Dec.Rare <- TYPE$Per.Dec.Freq[TYPE$Type=="S"]
E.Dec.Rare <- TYPE$Per.Dec.Freq[TYPE$Type=="E"]
R.Dec.Rare <- TYPE$Per.Dec.Freq[TYPE$Type=="R"]

S.No <- TYPE$Freq[TYPE$Type=="S"]
E.No <- TYPE$Freq[TYPE$Type=="E"]
R.No <- TYPE$Freq[TYPE$Type=="R"]

S.Unknown <- TYPE$Freq.1[TYPE$Type=="S"]
E.Unknown <- TYPE$Freq.1[TYPE$Type=="E"]
R.Unknown <- TYPE$Freq.1[TYPE$Type=="R"]

S.Yes <- TYPE$Freq.2[TYPE$Type=="S"]
E.Yes <- TYPE$Freq.2[TYPE$Type=="E"]
R.Yes <- TYPE$Freq.2[TYPE$Type=="R"]


vioplot(R.Dec.Rare, E.Dec.Rare, S.Dec.Rare, col=c("chocolate2"), names=c("", "", ""))
mtext("Declining Populations (%)", side=2, line=3)
vioplot(R.Per.Rare, E.Per.Rare, S.Per.Rare, col=c("gray55"), names=c("Remote", "Extractive", "Settled"))
mtext("Rare Populations (%)", side=2, line=3)
mtext("Site Type", side=1, line=2, outer=TRUE, adj=0.6, cex=1.2)


set.panel(3,1)
vioplot(R.No, E.No, S.No, col=c("chocolate2"), names=c("", "", ""))
mtext("Not Hunted Populations (%)", side=2, line=3, cex=0.6)
vioplot(R.Unknown, E.Unknown, S.Unknown, col=c("chocolate2"), names=c("", "", ""))
mtext("Unknown Populations (%)", side=2, line=3, cex=0.6)
vioplot(R.Yes, E.Yes, S.Yes, col=c("chocolate2"), names=c("Remote", "Extractive", "Settled"))
mtext("Hunted Populations (%)", side=2, line=3, cex=0.6)
mtext("Site Type", side=1, line=2, outer=TRUE, adj=0.6, cex=1.2)




# Examine the consistency in population status for species found at multiple sites

multiples <- table(WPI$sp,WPI$ind80)
multiples <- cbind(multiples, tot=rowSums(multiples))
multiples <- as.data.frame(multiples)
multiples <- multiples[multiples$tot>1,]

multiples <- match(WPI$sp, rownames(multiples))
multiples <- cbind(WPI, multiples)
multiples <- na.omit(multiples)
multiples <- data.frame(site.sp=as.character(multiples$site.sp), 
                        sp=as.character(multiples$sp), 
                        site=as.character(multiples$site),
                        Rare80=as.character(multiples$Rare80), 
                        Rare95=as.character(multiples$Rare95), 
                        site_type=as.character(multiples$site_type))
#write.csv(multiples, "MultiplePopulations2.csv") # Manually assign whether a species has multiple site types (Y/N) in column "types" 
#multiples <- read.csv("MultiplePopulationComparison2.csv")

#compare <- multiples[multiples$types=="Yes",]
#write.csv(compare, "Compare.csv") # Manually assign whether a species fares worse at a settled site
#CP <- read.csv("ComparePopulations.csv")

CP <- read.csv("MultiplePopulationComparison2.csv")

settled <- CP[CP$site_type=="settled" & CP$types=="Yes",]
table(settled$Compare80)
settled <- settled[settled$Settled80!="Rare",]
