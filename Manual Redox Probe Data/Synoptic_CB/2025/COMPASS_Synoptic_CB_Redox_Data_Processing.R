#######################################################################################
####### COMPASS Synoptic 
####### Chesapeake Bay Sites
####### Data Analysis Code: Manual Redox Profiles 
##############################################################################################


############################# Information ##################################
#Author: Stephanie J. Wilson
#Edited: 20220727

#Measurements taken with SWAP Redox Probes & Calomel Ref probe
#Protocol: https://docs.google.com/document/d/1dL2ZPhLFBAZJmVBCfKxup_KOYY96qA5IWfJmQen50pY/edit
#Units = mV 
############################################################################

############################ Read in data ##################################
#set working directory
setwd("S:/Biogeochemistry/People/Wilson (Steph)/Data/Redox")

#read in data 
dat <- read.csv("COMPASS_Synoptic_CB_Redox_Data_Raw.csv")

#Quick look at dataframe
head(dat)

#packages:
library(ggplot2)
library(matrixStats)
library(gridExtra)
library(ggpubr)
library(grid)
library(dplyr)
library(data.table)

std.error <- function (x, na.rm=FALSE){
  sqrt(var(x, na.rm = na.rm) / sum(!is.na(x)))}
############################################################################

############################ Adjust Raw Readings before avg. all ##################################

#need to adjust for pH: 
# redox = redoxi + [(ph - 5)59 + 224] This is the equation from Megonigal, Patrick, & Faulkner 1993
# They used a calomel electrode and used 224 to adjust to H2, we use a KCl ref probe, so 
# We need to add 202 instead 

soildat <- read.csv("CB_Soil_pH.csv")
head(soildat)

head(dat)

prelong <- dat[ , -c(4) ]
head(prelong)

#switch to long format 
datlong <-  melt(setDT(prelong), id.vars = c("Date","Month", "Year", "Site","Zone", "Depth"), variable.name = "RP_mV")
head(datlong)

#subset by site and transect location then adjust 
#MSM
MSM <- as.data.frame(subset(datlong, Site == "MSM", select = Date:value))
head(MSM)

MSMup <- as.data.frame(subset(MSM, Zone == "Upland", select = Date:value))
MSMup$value_corr <- (MSMup$value + ((4.38 - 5)*59) + 202) #4.38 is the pH of the soil in that zone
head(MSMup)

MSMtr <- as.data.frame(subset(MSM, Zone == "Transition", select = Date:value))
MSMtr$value_corr <- (MSMtr$value + ((5.31 - 5)*59) + 202) #5.31 is the pH of the soil in that zone
head(MSMtr)

MSMwc <- as.data.frame(subset(MSM, Zone == "Wetland", select = Date:value))
MSMwc$value_corr <- (MSMwc$value + ((4.72 - 5)*59) + 202) #5.31 is the pH of the soil in that zone
head(MSMwc)

#GWI
GWI <- as.data.frame(subset(datlong, Site == "GWI ", select = Date:value))
head(GWI)

GWIup <- as.data.frame(subset(GWI, Zone == "Upland", select = Date:value))
GWIup$value_corr <- (GWIup$value + ((4.38 - 5)*59) + 202) #4.38 is the pH of the soil in that zone
head(GWIup)

GWItr <- as.data.frame(subset(GWI, Zone == "Transition", select = Date:value))
GWItr$value_corr <- (GWItr$value + ((5.31 - 5)*59) + 202) #5.31 is the pH of the soil in that zone
head(GWItr)

GWIwc <- as.data.frame(subset(GWI, Zone == "Wetland", select = Date:value))
GWIwc$value_corr <- (GWIwc$value + ((4.72 - 5)*59) + 202) #5.31 is the pH of the soil in that zone
head(GWIwc)

#Gcrew
Gcrew <- as.data.frame(subset(datlong, Site == "Gcrew", select = Date:value))
head(Gcrew)

Gcrewup <- as.data.frame(subset(Gcrew, Zone == "Upland", select = Date:value))
Gcrewup$value_corr <- (Gcrewup$value + ((4.38 - 5)*59) + 202) #4.38 is the pH of the soil in that zone
head(Gcrewup)

Gcrewtr <- as.data.frame(subset(Gcrew, Zone == "Transition", select = Date:value))
Gcrewtr$value_corr <- (Gcrewtr$value + ((5.31 - 5)*59) + 202) #5.31 is the pH of the soil in that zone
head(Gcrewtr)

Gcrewwc <- as.data.frame(subset(Gcrew, Zone == "Wetland", select = Date:value))
Gcrewwc$value_corr <- (Gcrewwc$value + ((4.72 - 5)*59) + 202) #5.31 is the pH of the soil in that zone
head(Gcrewwc)


#SWH
SWH <- as.data.frame(subset(datlong, Site == "SWH", select = Date:value))
head(SWH)

SWHupcon <- as.data.frame(subset(SWH, Zone == "Upland", select = Date:value))
SWHupcon$value_corr <- (SWHupcon$value + ((4.9 - 5)*59) + 202) #4.9 is the pH of the soil in that zone - using same pH as upland here
head(SWHupcon)

SWHup <- as.data.frame(subset(SWH, Zone == "Swamp", select = Date:value))
SWHup$value_corr <- (SWHup$value + ((4.9 - 5)*59) + 202) #4.9 is the pH of the soil in that zone
head(SWHup)

SWHtrA <- as.data.frame(subset(SWH, Zone == "Transition-Hummock", select = Date:value))
SWHtrA$value_corr <- (SWHtrA$value + ((6.2 - 5)*59) + 202) #6.2 is the pH of the soil in that zone
head(SWHtrA)

SWHtrB <- as.data.frame(subset(SWH, Zone == "Transition-Hollow", select = Date:value))
SWHtrB$value_corr <- (SWHtrB$value + ((6.2 - 5)*59) + 202) #6.2 is the pH of the soil in that zone
head(SWHtrB)

SWHtrC <- as.data.frame(subset(SWH, Zone == "Transition", select = Date:value))
SWHtrC$value_corr <- (SWHtrC$value + ((6.2 - 5)*59) + 202) #6.2 is the pH of the soil in that zone
head(SWHtrC)

SWHwc <- as.data.frame(subset(SWH, Zone == "Wetland", select = Date:value))
SWHwc$value_corr <- (SWHwc$value + ((5.9 - 5)*59) + 202) #5.9 is the pH of the soil in that zone
head(SWHwc)


#Bring all the data back together 
all_corr <- rbind(MSMup, MSMtr, MSMwc, GWIup, GWItr, GWIwc, Gcrewup, Gcrewtr, Gcrewwc, SWHupcon, SWHup, SWHtrA, SWHtrB, SWHwc)
head(all_corr)

red <- all_corr %>%
  group_by(Site, Zone, Depth) %>%
  summarise(Date = mean(Date), Avg_mV = mean(value_corr, na.rm=TRUE), Std_dv = sd(value_corr, na.rm=TRUE), Std_err = std.error(value_corr, na.rm=TRUE))
head(red)

red1 <- all_corr %>%
  group_by(Site, Zone, Month, Year, Depth) %>%
  summarise(Date = mean(Date), Avg_mV = mean(value_corr, na.rm=TRUE), Std_dv = sd(value_corr, na.rm=TRUE), Std_err = std.error(value_corr, na.rm=TRUE))
head(red1)

red2 <- all_corr %>%
  group_by(Site, Zone, Year, Depth) %>%
  summarise(Date = mean(Date), Avg_mV = mean(value_corr, na.rm=TRUE), Std_dv = sd(value_corr, na.rm=TRUE), Std_err = std.error(value_corr, na.rm=TRUE))
head(red2)


#pull out data from during N cycle experiment
#TenCm <-  as.data.frame(subset(all_corr, Depth == "10", select = Month:value_corr))
#head(TenCm)

#Aug <- as.data.frame(subset(TenCm, Month == "August", select = Month:value_corr))
#September <- as.data.frame(subset(TenCm, Month == "September", select = Month:value_corr))
#out <- rbind(Aug, September)

#write.csv(out, "COMPASS_Synoptic_Aug_Sept_Redox_10cm.csv")

########################################################################################

############################ Write Out Adjusted & Avg'd Data ##################################

#Look at the data 
head(red1)


#write a .csv of the new data frame 
write.csv(red1, "COMPASS_Synoptic_CB_Redox_Data_Processed.csv")

############################################################################

########################### Plot all data together - Correct Colors here #####################################
#rename sweet hall zones 
red <- red %>% 
  mutate(Type = case_when(
    str_detect(Site, "SWH" & Zone, "Upland") ~ "Swamp",
    str_detect(Site, "SWH" & Zone, "Upland Control") ~ "Upland"))


head(red)

rall <- ggplot(red, aes(x=Depth, y= Avg_mV, col=Zone, fill=Zone))+
  geom_point(size=3) +  geom_line(size=1.5) +
  scale_x_reverse( lim=c(45,5)) + ylim(c(-150, 700)) +
  coord_flip() + theme_classic() +
  facet_grid(~factor(Site, levels=c("MSM", "GWI ", "Gcrew", "SWH"))) +
  labs(x="Depth (cm)", y="Redox Potential (mV)", title ="May 2022") +
  scale_colour_manual(values = c("darkblue","turquoise4", "darkgreen", "black", "purple", "orange")) +
  theme(panel.spacing = unit(.5, "lines"),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        strip.background = element_rect(color = "black", size = 1),
        plot.margin = margin(t=0.1, r=0.1, l=0.1, b=0.1, "cm"),
        plot.title = element_text(hjust = 0.0),
        axis.text=element_text(size=12),axis.title=element_text(size=12),
        legend.text=element_text(size=12), legend.title=element_blank(),
        legend.position = "bottom")+
  geom_errorbar(aes(ymin=Avg_mV-Std_err,
                    ymax=Avg_mV+Std_err),width=0.3,position=position_dodge(.1))
rall

#pretty steady redox across all the months so maybe lets try averaging across the months? 
#head(dat)
#dat1 <- dat[ , -c(1,2,3,4)]
#head(dat1)

#long <- melt(setDT(dat1), id.vars = c("Site","Zone", "Depth"), variable.name = "RP_mV")
#head(long)
#
#allavg <- long %>%
 # group_by(Site, Zone, Depth) %>%
#  summarise(
    #Avg = mean(value, na.rm=TRUE),
   # sd = sd(value, na.rm=TRUE),
  #  n = n(),
 #   se = sd / sqrt(n))
#allavg

#allavg <- as.data.frame(allavg)
#head(allavg)

head(all_corr)

allavg <- all_corr %>%
  group_by(Site, Zone, Depth) %>%
  summarise(
    Avg = mean(value, na.rm=TRUE),
    sd = sd(value, na.rm=TRUE),
    n = n(),
    se = sd / sqrt(n))
allavg

allavg <- as.data.frame(allavg)
head(allavg)

allavg$Site <- factor(allavg$Site, levels=c("Gcrew", "GWI ","MSM", "SWH"))
levels(allavg$Site)

#allavg$Site <- recode_factor(allavg$Site, Gcrew = "GCrew", GWI  = "Goodwin Island", MSM = "MoneyStump")
#levels(allavg$Site) <- list(Gcrew = "GCrew", GWI  = "Goodwin Island", MSM = "MoneyStump")
#allavg$Site <- recode_factor(allavg$Site, GCrew = "Gcrew", Goodwin_Island  = "GWI ", MoneyStump = "MSM")

levels(allavg$Site)[levels(allavg$Site) == "Gcrew"]  <- "GCReW"
levels(allavg$Site)[levels(allavg$Site) == "GWI "] <- "Goodwin Islands"
levels(allavg$Site)[levels(allavg$Site) == "MSM"]  <- "Moneystump"
levels(allavg$Site)[levels(allavg$Site) == "SWH"]  <- "Sweet Hall"

compass_coolors <- c("#3B5D32", "#E8DDB5", "#777469")
compass_coolors1 <- c("#3B5D32", "#CFB863", "#777469", "black")
compass_coolors2 <- c("#726E75", "#F2AF29",   "darkorange","#F2AF29", "black", "#5A9D43")
compass_coolors3 <- c("black", "#F2AF29",  "#5A9D43")


# CFB863

avall <- ggplot(allavg, aes(x=Depth, y= Avg, fill=Zone, group=Zone))+
  geom_line(size=2, aes(col=Zone)) + 
  geom_point(size=4.5, shape = 21, col= "black" ) +
  scale_x_reverse( lim=c(38,8)) + #ylim(c(-350, 500)) +
  coord_flip() + theme_classic() +
  facet_grid(~factor(Site))+ #, levels=c("GCrew", "Goodwin Island","MoneyStump"))) +
  labs(x="Depth (cm)", y="Redox Potential (mV)", title =" ") +
  scale_colour_manual(values = compass_coolors2) +
  scale_fill_manual(values = compass_coolors2) +
   #scale_fill_manual(values = c("#777469","#E8DDB5", "#3B5D32")) +
  theme(panel.spacing = unit(.5, "lines"),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        strip.background = element_rect(color = "black", size = 1),
        strip.text.x = element_text(size = 16),
        plot.margin = margin(t=0.1, r=0.1, l=0.1, b=0.1, "cm"),
        plot.title = element_text(hjust = 0.0),
        axis.text=element_text(size=16),axis.title=element_text(size=16),
        legend.text=element_text(size=16), legend.title=element_blank(),
        legend.position = "bottom")+
  geom_errorbar(aes(ymin=(Avg-(2*se)), ymax=(Avg+(2*se))),width=0,position=position_dodge(.1))  #2 times the se here for scale 
avall

ggsave(plot = avall, "Redox_Corr_Average_All_20221104.png", h = 6, w = 12, type = "cairo-png")


ravall <- ggplot(allavg, aes(x=Depth, y= Avg, shape=Zone))+
  geom_rect(aes(xmin=5, xmax=8, ymin=300, ymax=600), fill="green4", alpha=0.05)+
  geom_rect(aes(xmin=5, xmax=8, ymin=10, ymax=350), fill="mediumpurple3", alpha=0.05)+
  geom_rect(aes(xmin=5, xmax=8, ymin=-150, ymax=220), fill="orange1", alpha=0.05)+
  geom_rect(aes(xmin=5, xmax=8, ymin=-150, ymax=-100), fill="tan4", alpha=0.05)+
  geom_point(size=3) +  geom_line(size=1.5) +
  scale_x_reverse( lim=c(45,5)) + ylim(c(-300, 700)) +
  coord_flip() + theme_classic() +
  facet_grid(~factor(Site, levels=c("MSM", "GWI ", "Gcrew"))) +
  labs(x="Depth (cm)", y="Redox Potential (mV)", title =" ") +
  scale_shape_manual(values = c(15, 16, 17)) +
  theme(panel.spacing = unit(.5, "lines"),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        strip.background = element_rect(color = "black", size = 1),
        plot.margin = margin(t=0.1, r=0.1, l=0.1, b=0.1, "cm"),
        plot.title = element_text(hjust = 0.0),
        axis.text=element_text(size=12),axis.title=element_text(size=12),
        legend.text=element_text(size=12), legend.title=element_blank(),
        legend.position = "none")+
  geom_errorbar(aes(ymin=Avg-se, ymax=Avg+se),width=0.3,position=position_dodge(.1))
ravall 

ggsave(plot = ravall, "Redox_Average_Corr_All_20220915.png", h = 6, w = 12, type = "cairo-png")

#########################################################################################

######################## Box plot of redox in zones ########################################

boxall <- ggplot(all_corr, aes(x=Zone, y= value_corr, fill=Zone, group=Zone))+
  geom_boxplot()+ 
  theme_classic() +
  facet_grid(~factor(Site))+ #, levels=c("GCrew", "Goodwin Island","MoneyStump"))) +
  labs(x="Depth (cm)", y="Redox Potential (mV)", title =" ") +
  scale_colour_manual(values = compass_coolors2) +
  scale_fill_manual(values = compass_coolors2) +
  theme(panel.spacing = unit(.5, "lines"),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        strip.background = element_rect(color = "black", size = 1),
        strip.text.x = element_text(size = 12),
        plot.margin = margin(t=0.1, r=0.1, l=0.1, b=0.1, "cm"),
        plot.title = element_text(hjust = 0.0),
        axis.text=element_text(size=12),axis.title=element_text(size=12),
        legend.text=element_text(size=12), legend.title=element_blank(),
        legend.position = "bottom")
boxall

#pull out SWH and plot separately for GCW, MSM, GWI and SWH

swh_all_corr <- all_corr %>% 
  mutate(Site = "SWH")

other_all_corr <- all_corr %>% 
  filter(!Site == "SWH")

swh_all_corr$Zone <- as.character(swh_all_corr$Zone)
swh_all_corr$Zone <- factor(swh_all_corr$Zone, levels=c("Upland", "Swamp", "Transition","Transition-Hummock", "Transition-Hollow", "Wetland"))
levels(swh_all_corr$Zone)

swhall <- ggplot(swh_all_corr, aes(x=Zone, y= value_corr, fill=Zone, group=Zone))+
  geom_boxplot()+ 
  theme_classic() +
  #facet_grid(~factor(Site))+ #, levels=c("GCrew", "Goodwin Island","MoneyStump"))) +
  labs(x="Depth (cm)", y="Redox Potential (mV)", title =" ") +
  scale_colour_manual(values = compass_coolors2) +
  scale_fill_manual(values = compass_coolors2) +
  theme(panel.spacing = unit(.5, "lines"),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        strip.background = element_rect(color = "black", size = 1),
        strip.text.x = element_text(size = 12),
        plot.margin = margin(t=0.1, r=0.1, l=0.1, b=0.1, "cm"),
        plot.title = element_text(hjust = 0.0),
        axis.text=element_text(size=12),axis.title=element_text(size=12),
        legend.text=element_text(size=12), legend.title=element_blank(),
        legend.position = "bottom")
swhall

other_all_corr$Zone <- as.character(other_all_corr$Zone)
other_all_corr$Zone <- factor(other_all_corr$Zone, levels=c("Upland", "Transition","Wetland"))
levels(other_all_corr$Zone)

otherall <- ggplot(data=other_all_corr, aes(x=Zone, y= value_corr, fill=Zone, group=Zone))+
  geom_boxplot()+ 
  theme_bw() +
  facet_grid(~factor(Site))+ #, levels=c("GCrew", "Goodwin Island","MoneyStump"))) +
  labs(x=" ", y="Redox Potential (mV)", title =" ") +
  scale_colour_manual(values = compass_coolors3) +
  scale_fill_manual(values = compass_coolors3) +
  ylim(-500,1000)+
  theme(panel.spacing = unit(.5, "lines"),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        strip.background = element_rect(color = "black", size = 1),
        strip.text.x = element_text(size = 12),
        plot.margin = margin(t=0.1, r=0.1, l=0.1, b=0.1, "cm"),
        plot.title = element_text(hjust = 0.0),
        axis.text=element_text(size=12),axis.title=element_text(size=12),
        legend.text=element_text(size=12), legend.title=element_blank(),
        legend.position = "bottom")
otherall


##########################################################################################################

head(all_corr)

######################### For ESS PI Meeting ############################################
head(red)
#just pull out MSM data 
msm_redox <- as.data.frame(subset(red, Site == 'MSM', select = Site:Std_err))
head(msm_redox)

#grab Fausto Data 
setwd("S:/Biogeochemistry/People/Wilson (Steph)/Conferences/ESS PI 2024")

le_redox <- read.csv("df_PTR_REDOX_PI_meeting_2024_RBP_b.csv")
head(le_redox)

le_red <- le_redox %>%
  group_by(Site, Zone, Depth) %>%
  summarise(Date = mean(Date), Avg_mV = mean(Ehc_mV, na.rm=TRUE), Std_dv = sd(Ehc_mV, na.rm=TRUE), Std_err = std.error(Ehc_mV, na.rm=TRUE))
head(le_red)

setwd("S:/Biogeochemistry/People/Wilson (Steph)/Data/Redox")

fig_red <- rbind(msm_redox, le_red)
head(fig_red)
fig_red$Zone <- factor(fig_red$Zone, levels=c('Upland', 'Transition', 'Wetland'))

redfig <- ggplot(fig_red, aes(x=Depth, y= Avg_mV, col=Zone, fill=Zone))+
  geom_point(size=4) +  geom_line(size=1) +
  scale_x_reverse( lim=c(40,0)) + ylim(c(-150, 700)) +
  coord_flip() + theme_classic() +
  facet_grid(~factor(Site, levels=c("MSM", "PTR"))) +
  labs(x="Depth (cm)", y="Redox Potential (mV)", title =" ") +
  scale_colour_manual(values = c( "#726E75","#F2AF29",  "#5A9D43")) +
  theme(panel.spacing = unit(.5, "lines"),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        strip.background = element_rect(color = "black", size = 1),
        strip.text.x = element_text(size = 16),
        plot.margin = margin(t=0.1, r=0.1, l=0.1, b=0.1, "cm"),
        plot.title = element_text(hjust = 0.0),
        axis.text=element_text(size=16),axis.title=element_text(size=16),
        legend.text=element_text(size=16), legend.title=element_blank(),
        legend.position = "bottom")+
  geom_errorbar(aes(ymin=Avg_mV-Std_err,
                    ymax=Avg_mV+Std_err),width=0.3,position=position_dodge(.1))
redfig

##########################################################################################



#other ways to plot and look at the data: 

############################ OLD Plot Data with Shapes ##################################

# Create a text
grob <- grobTree(textGrob("Denitrification", x=30,  y=200, 
                          gp=gpar(col="green4", fontsize=12)))
# Plot
sp2 + annotation_custom(grob)

rmay1 <- ggplot(RedMay, aes(x=Depth, y= Avg, shape=Zone))+
  geom_rect(aes(xmin=5, xmax=8, ymin=300, ymax=600), fill="green4", alpha=0.05)+
  geom_rect(aes(xmin=5, xmax=8, ymin=10, ymax=350), fill="mediumpurple3", alpha=0.05)+
  geom_rect(aes(xmin=5, xmax=8, ymin=-150, ymax=220), fill="orange1", alpha=0.05)+
  geom_rect(aes(xmin=5, xmax=8, ymin=-150, ymax=-100), fill="tan4", alpha=0.05)+
  geom_point(size=3) +  geom_line(size=1.5) +
  scale_x_reverse( lim=c(45,5)) + ylim(c(-150, 700)) +
  coord_flip() + theme_classic() +
  facet_grid(~factor(Site, levels=c("MSM", "GWI ", "Gcrew"))) +
  #facet_grid(cols = vars(RedMay$Site)) +
  labs(x="Depth (cm)", y="Redox Potential (mV)", title ="May 2022") +
  scale_shape_manual(values = c(15, 16, 17)) +
  theme(panel.spacing = unit(.5, "lines"),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        strip.background = element_rect(color = "black", size = 1),
        plot.margin = margin(t=0.1, r=0.1, l=0.1, b=0.1, "cm"),
        plot.title = element_text(hjust = 0.0),
        axis.text=element_text(size=12),axis.title=element_text(size=12),
        legend.text=element_text(size=12), legend.title=element_blank(),
        legend.position = "none")+
  geom_errorbar(aes(ymin=Avg-Std_err,
                    ymax=Avg+Std_err),width=0.3,position=position_dodge(.1))
rmay1 


rjune1 <- ggplot(RedJune,  aes(x=Depth, y= Avg, shape=Zone))+
  geom_rect(aes(xmin=5, xmax=8, ymin=300, ymax=600), fill="green4", alpha=0.05)+
  geom_rect(aes(xmin=5, xmax=8, ymin=10, ymax=350), fill="mediumpurple3", alpha=0.05)+
  geom_rect(aes(xmin=5, xmax=8, ymin=-150, ymax=220), fill="orange1", alpha=0.05)+
  geom_rect(aes(xmin=5, xmax=8, ymin=-150, ymax=-100), fill="tan4", alpha=0.05)+
  geom_point(size=3) +  geom_line(size=1.5) +
  scale_x_reverse( lim=c(45,5)) + ylim(c(-150, 700)) +
  coord_flip() + theme_classic() +
  #facet_grid(cols = vars(RedJune$Site)) +
  facet_grid(~factor(Site, levels=c("MSM", "GWI ", "Gcrew"))) +
  labs(x="Depth (cm)", y="Redox Potential (mV)", title="June 2022") +
  #  title =" ") +
  scale_shape_manual(values = c(15, 16, 17)) +
  theme(panel.spacing = unit(.5, "lines"),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        strip.background = element_rect(color = "black", size = 1),
        plot.margin = margin(t=0.1, r=0.1, l=0.1, b=0.1, "cm"),
        plot.title = element_text(hjust = 0.0),
        axis.text=element_text(size=12),axis.title=element_text(size=12),
        legend.text=element_text(size=12), legend.title=element_blank(),
        legend.position = "none")+
  geom_errorbar(aes(ymin=Avg-Std_err,
                    ymax=Avg+Std_err),width=0.3,position=position_dodge(.1))
rjune1

rjuly1 <- ggplot(RedJuly,  aes(x=Depth, y= Avg, shape=Zone))+
  geom_rect(aes(xmin=5, xmax=8, ymin=300, ymax=600), fill="green4", alpha=0.05)+
  geom_rect(aes(xmin=5, xmax=8, ymin=10, ymax=350), fill="mediumpurple3", alpha=0.05)+
  geom_rect(aes(xmin=5, xmax=8, ymin=-150, ymax=220), fill="orange1", alpha=0.05)+
  geom_rect(aes(xmin=5, xmax=8, ymin=-150, ymax=-100), fill="tan4", alpha=0.05)+
  geom_point(size=3) +  geom_line(size=1.5) +
  scale_x_reverse( lim=c(45,5)) + ylim(c(-150, 700)) +
  coord_flip() + theme_classic() +
  #facet_grid(cols = vars(RedJune$Site)) +
  facet_grid(~factor(Site, levels=c("MSM", "GWI ", "Gcrew"))) +
  labs(x="Depth (cm)", y="Redox Potential (mV)", title="July 2022") +
  #  title =" ") +
  scale_shape_manual(values = c(15, 16, 17)) +
  theme(panel.spacing = unit(.5, "lines"),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        strip.background = element_rect(color = "black", size = 1),
        plot.margin = margin(t=0.1, r=0.1, l=0.1, b=0.1, "cm"),
        plot.title = element_text(hjust = 0.0),
        axis.text=element_text(size=12),axis.title=element_text(size=12),
        legend.text=element_text(size=12), legend.title=element_blank(),
        legend.position = "bottom")+
  geom_errorbar(aes(ymin=Avg-Std_err,
                    ymax=Avg+Std_err),width=0.3,position=position_dodge(.1))
rjuly1

#plot these together 
ggarrange(rmay1, rjune1, rjuly1, ncol=1, nrow=3, common.legend = TRUE, legend="bottom")
############################################################################################

############################ Subset Data ##################################
head(red1)

#subset by Month 
RedMay <- as.data.frame(subset(red1, Month == "May", select = Date:Std_err))
head(RedMay)

#RedMay$Site <- factor(RedMay$Site, levels=c("MSM", "GWI"))
#levels(RedMay$Site) <- c("MSM", "GWI")

RedJune <- as.data.frame(subset(red, Month == "June", select = Date:Std_err))
head(RedJune)

#levels(RedJune$Site) <- c("MSM", "GWI", "Gcrew" )
#RedJune$Site = factor(RedJune$Site, levels=c("MSM", "GWI", "GCrew"), labels=c("MSM", "GWI", "GCrew"))

RedJuly <- as.data.frame(subset(red, Month == "July", select = Date:Std_err))
head(RedJuly)

#######################################################################################

############################ Plot Data ##################################
rmay <- ggplot(RedMay, aes(x=Depth, y= Avg, col=Zone, fill=Zone))+
  geom_point(size=3) +  geom_line(size=1.5) +
  scale_x_reverse( lim=c(45,5)) + ylim(c(-150, 700)) +
  coord_flip() + theme_classic() +
  facet_grid(~factor(Site, levels=c("MSM", "GWI ", "Gcrew"))) +
  #facet_grid(cols = vars(RedMay$Site)) +
  labs(x="Depth (cm)", y="Redox Potential (mV)", title ="May 2022") +
  scale_colour_manual(values = c("darkblue","turquoise4", "darkgreen")) +
  theme(panel.spacing = unit(.5, "lines"),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        strip.background = element_rect(color = "black", size = 1),
        plot.margin = margin(t=0.1, r=0.1, l=0.1, b=0.1, "cm"),
        plot.title = element_text(hjust = 0.0),
        axis.text=element_text(size=12),axis.title=element_text(size=12),
        legend.text=element_text(size=12), legend.title=element_blank(),
        legend.position = "none")+
  geom_errorbar(aes(ymin=Avg-Std_err,
                    ymax=Avg+Std_err),width=0.3,position=position_dodge(.1))
rmay 


rjune <- ggplot(RedJune, aes(x=Depth, y= Avg, col=Zone, fill=Zone))+
  geom_point(size=3) +  geom_line(size=1.5) +
  scale_x_reverse( lim=c(45,5)) + ylim(c(-150, 700)) +
  coord_flip() + theme_classic() +
  #facet_grid(cols = vars(RedJune$Site)) +
  facet_grid(~factor(Site, levels=c("MSM", "GWI ", "Gcrew"))) +
  labs(x="Depth (cm)", y="Redox Potential (mV)", title="June 2022") +
  #  title =" ") +
  scale_colour_manual(values = c("darkblue","turquoise4", "darkgreen")) +
  theme(panel.spacing = unit(.5, "lines"),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        strip.background = element_rect(color = "black", size = 1),
        plot.margin = margin(t=0.1, r=0.1, l=0.1, b=0.1, "cm"),
        plot.title = element_text(hjust = 0.0),
        axis.text=element_text(size=12),axis.title=element_text(size=12),
        legend.text=element_text(size=12), legend.title=element_blank(),
        legend.position = "none")+
  geom_errorbar(aes(ymin=Avg-Std_err,
                    ymax=Avg+Std_err),width=0.3,position=position_dodge(.1))
rjune

rjuly <- ggplot(RedJuly, aes(x=Depth, y= Avg, col=Zone, fill=Zone))+
  geom_point(size=3) +  geom_line(size=1.5) +
  scale_x_reverse( lim=c(45,5)) + ylim(c(-150, 700)) +
  coord_flip() + theme_classic() +
  #facet_grid(cols = vars(RedJune$Site)) +
  facet_grid(~factor(Site, levels=c("MSM", "GWI ", "Gcrew"))) +
  labs(x="Depth (cm)", y="Redox Potential (mV)", title="July 2022") +
  #  title =" ") +
  scale_colour_manual(values = c("darkblue","turquoise4", "darkgreen")) +
  theme(panel.spacing = unit(.5, "lines"),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        strip.background = element_rect(color = "black", size = 1),
        plot.margin = margin(t=0.1, r=0.1, l=0.1, b=0.1, "cm"),
        plot.title = element_text(hjust = 0.0),
        axis.text=element_text(size=12),axis.title=element_text(size=12),
        legend.text=element_text(size=12), legend.title=element_blank(),
        legend.position = "bottom")+
  geom_errorbar(aes(ymin=Avg-Std_err,
                    ymax=Avg+Std_err),width=0.3,position=position_dodge(.1))
rjuly

#plot these together 
ggarrange(rmay, rjune, rjuly, ncol=1, nrow=3, common.legend = TRUE, legend="bottom")

################################

#### END
