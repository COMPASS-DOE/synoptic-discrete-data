library(dplyr)
library(ggplot2)

July <- read.csv("COMPASS_SynopticCB_H2S_202607_processed.csv")
July <- July %>% 
  rename(H2S_mean_uM_July = H2S_mean_uM) %>% 
  select(Site, Zone, Replicate, Depth_cm, H2S_mean_uM_July)


Sept <- read.csv("COMPASS_SynopticCB_H2S_202609_MD_processed.csv")
Sept <- Sept %>% 
  rename(H2S_mean_uM_Sept = H2S_mean_uM) %>% 
  select(Site, Zone, Replicate, Depth_cm, H2S_mean_uM_Sept)

all <- merge(July, Sept, by = c("Site", "Zone", "Replicate", "Depth_cm"))


ggplot(all)+
  geom_point(aes(x=H2S_mean_uM_July,y=H2S_mean_uM_Sept),size=3,alpha=0.8)+
  theme_light() + 
  labs(y="July Sulfide",x="Sept Sulfide")+
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed", linewidth = 1)

