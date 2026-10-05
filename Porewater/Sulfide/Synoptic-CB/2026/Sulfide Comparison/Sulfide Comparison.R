library(dplyr)
library(ggplot2)


####Comparison of values from each month####
July <- read.csv("COMPASS_SynopticCB_H2S_202607_processed.csv")
July <- July %>% 
  rename(H2S_mean_uM_July = H2S_mean_uM) %>% 
  dplyr::select(Site, Zone, Replicate, Depth_cm, H2S_mean_uM_July)


Sept <- read.csv("COMPASS_SynopticCB_H2S_202609_MD_processed.csv")
Sept <- Sept %>% 
  rename(H2S_mean_uM_Sept = H2S_mean_uM) %>% 
  dplyr::select(Site, Zone, Replicate, Depth_cm, H2S_mean_uM_Sept)

all <- merge(July, Sept, by = c("Site", "Zone", "Replicate", "Depth_cm"))


ggplot(all)+
  geom_point(aes(x=H2S_mean_uM_July,y=H2S_mean_uM_Sept),size=3,alpha=0.8)+
  theme_light() + 
  labs(x="July Sulfide (uM)",y="September Sulfide (uM)")+
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed", linewidth = 1)




####Comparison ofsamples from July reran in September####
Orig <- read.csv("COMPASS_SynopticCB_H2S_202607_processed.csv")
Orig <- Orig %>% 
  rename(H2S_mean_uM_Orig = H2S_mean_uM) %>% 
  dplyr::select(Site, Zone, Replicate, Depth_cm, H2S_mean_uM_Orig)


Reruns <- read.csv("July_samples.csv")
Reruns <- Reruns %>% 
  rename(H2S_mean_uM_Reruns = H2S_mean_uM) %>% 
  mutate(Sample_ID = gsub("_JULY", "", Sample_ID)) %>% 
  separate(
    col = Sample_ID,
    sep = "_",
    into = c("Site", "Zone", "Replicate", "Depth_cm"),
    remove = FALSE)  %>% 
  mutate(Depth_cm = ifelse(Zone == "SW", "0", Depth_cm)) %>% 
  dplyr::select(Site, Zone, Replicate, Depth_cm, H2S_mean_uM_Reruns)

all <- merge(Orig, Reruns, by = c("Site", "Zone", "Replicate", "Depth_cm"))


ggplot(all)+
  geom_point(aes(x=H2S_mean_uM_Orig,y=H2S_mean_uM_Reruns),size=3,alpha=0.8)+
  theme_light() + 
  labs(x="Original Sulfide (7/20/26)",y="Rerun Sulfide (10/01/26)")+
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed", linewidth = 1)+ 
  coord_cartesian(xlim = c(0, 3500), ylim = c(0, 3500))
