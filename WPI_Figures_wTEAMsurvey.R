# Create histogram displaying coefficients by site type (Remote/Extractive/Settled) and number of years (3/4/5+)
library(ggplot2)
library(ggthemes)

# Add column to data for 5+ years of sampling
Year_cat <- ifelse(WPI$nyears==5, "5+", ifelse(WPI$nyears==6,"5+", ifelse(WPI$nyears==7, "5+", WPI$nyears)))
Site_cat <- as.factor(ifelse(WPI$site_type=="remote", "I - Remote", ifelse(WPI$site_type=="extractive", "II - Extractive", ifelse(WPI$site_type=="settled", "III - Settled", NA))))
WPI <- data.frame(WPI, Year_cat, Site_cat)

ggplot(WPI, aes(median.coeff)) +
  geom_histogram() +
  facet_grid(nyears~.) +
  xlim(-3, 3) +
  xlab("Number of Years") +
  ylab("Frequency")


ggplot(WPI, aes(median.coeff)) +
  geom_histogram() +
  facet_grid(Site_cat~Year_cat) +
  xlim(-3, 3) +
  xlab("Median Coefficient") +
  ylab("Number of Populations")

WPIplot1 <- ggplot(WPI, aes(median.coeff, fill=factor(ind80))) +
  geom_histogram() +
  facet_grid(Year_cat~Hunted2) +
  xlim(-2, 2) +
  xlab("Median coefficient") +
  ylab("Terrestrial Vertebrate Populations") 
  

#test + theme_bw() + theme(legend.position = c(0.1, 0.9)) + guides(fill=guide_legend(title=NULL))
WPIplot1 + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(legend.position = c(0.1, 0.87)) + scale_fill_manual(values=c( "red", "blue","light gray"),
                                                                             name="Population Status",
                                                                             breaks=c("decreasing", "stable", "increasing"),
                                                                             labels=c("Decreasing", "Stable", "Increasing"))



WPIplot2 <- ggplot(WPI, aes(x=median.coeff, fill=Hunted2)) + geom_histogram(binwidth=.5, alpha=.5, position="identity") +
  facet_grid(nyears~.) +
  xlab("Median coefficient") +
  ylab("Terrestrial Vertebrate Populations") 

WPIplot2 + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(legend.position = c(0.1, 0.87)) + scale_fill_manual(values=c( "dark blue", "brown"),
                                                                                                                                                       name="Hunted", 
                                                                                                                                                       breaks=c("No", "Yes"),
                                                                                                                                                       labels=c("No", "Yes"))



ggplot(WPI, aes(x=median.coeff, fill=Hunted2)) + geom_density(alpha=.3) + 
  facet_grid(nyears~.) +
  geom_vline(aes(xintercept=mean(WPI$median.coeff[WPI$Hunted2=="No"])), colour="black") +
  geom_vline(aes(xintercept=mean(WPI$median.coeff[WPI$Hunted2=="Yes"])), colour="red", linewidth=3)

  geom_hline(aes(yintercept=0)) +
  geom_line(stat="vline", xintercept="mean", y=10)
  
  geom_line(stat="vline", xintercept="mean")

