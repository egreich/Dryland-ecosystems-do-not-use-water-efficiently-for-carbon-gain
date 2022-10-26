### This script file will read in the SAM model output across sites
### and make summary graphs

# Load packages
library(tidyverse)
library(ggforce) # for facet grids
library(cowplot) # for ggsave
library(ggpubr) # for 1:1 graphs stat_cor

# Create necessary folders if they do not already exist
if(!file.exists("plots")) { dir.create("plots")}

path_out = "./plots/" # set save path



### Look at stuff  ###
#load("./output_coda/zc_seg.RData")


sum_IN_seg = read.csv("./output_dfs/df_sum_seg.csv")
sum_IN_ses = read.csv("./output_dfs/df_sum_ses.csv")
sum_IN_wjs = read.csv("./output_dfs/df_sum_wjs.csv")
sum_IN_mpj = read.csv("./output_dfs/df_sum_mpj.csv")
sum_IN_vcp1 = read.csv("./output_dfs/df_sum_vcp1.csv")
sum_IN_vcp2 = read.csv("./output_dfs/df_sum_vcp2.csv")
# sum_IN_vcp = data.frame(X = sum_IN_vcp1$X, var = sum_IN_vcp1$var, 
#                         ID1 = sum_IN_vcp1$ID1, ID2 = sum_IN_vcp1$ID2,
#                         mean = (sum_IN_vcp1$mean + sum_IN_vcp2$mean)/2,
#                         median = (sum_IN_vcp1$median + sum_IN_vcp2)/2,
#                         sd = (sum_IN_vcp1$sd + sum_IN_vcp2$sd)/2,
#                         pc2.5 = (sum_IN_vcp1$pc2.5 + sum_IN_vcp2$pc2.5)/2,
#                         pc97.5 = (sum_IN_vcp1$pc97.5 + sum_IN_vcp2$pc97.5)/2,
#                         site = "vcp")
sum_IN_vcm1 = read.csv("./output_dfs/df_sum_vcm1.csv")
sum_IN_vcm2 = read.csv("./output_dfs/df_sum_vcm2.csv")
sum_IN_vcs = read.csv("./output_dfs/df_sum_vcs.csv")

sum_IN = rbind(sum_IN_seg, sum_IN_ses, sum_IN_wjs, sum_IN_mpj, 
               sum_IN_vcp1, sum_IN_vcp2, sum_IN_vcm1, sum_IN_vcm2, sum_IN_vcs)
sum_IN$ID1 = as.numeric(sum_IN$ID1)
sum_IN$site <- factor(sum_IN$site, levels = c("seg", "ses", "wjs", "mpj", "vcp1", "vcp2", "vcm1", "vcm2", "vcs"))

# Graph
p <- sum_IN %>%
  filter(var %in% c("wP","wSd","wSs","wT","wV",
                    "wP.monthly","wP.weekly")) %>%
  ggplot(aes(x=ID1, y=mean, color=site)) + 
  geom_errorbar(aes(ymax = pc97.5, ymin = pc2.5), position = "dodge") +
  geom_point(pch=21, size = .5) +
  scale_color_manual(values=c("red", "tomato2", "sienna2", "orange","light blue","light blue","blue","blue", "dark blue")) +
  facet_grid(site ~ var) +
  labs(title = NULL, y = "weight", x = NULL)+
  theme_bw() +
  theme(legend.position = "none",
        legend.text=element_text(size=14),
        text = element_text(size=14),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5))
p

ggsave2("weights_grid.png", plot = p, path = path_out)


df_p <- sum_IN %>%
  filter(var %in% c("beta1")) #%>%
#filter(site %in% c("seg","ses", "wjs", "mpj")) #%>%
#filter(site %in% c("vcp1","vcp2", "vcm1", "vcm2", "vcs"))
df_p$ID1 <-factor(df_p$ID1 , levels = c("1","2","3","4", "5", "6"), labels = c("VPD", "T", "P", "PAR", "Sshall", "Sdeep"))

p <- ggplot(data = df_p, aes(x=ID1, y=mean, color=site)) + 
  geom_errorbar(aes(ymax = pc97.5, ymin = pc2.5), position = "dodge") +
  geom_point(pch=21, size = .5) +
  scale_color_manual(values=c("red", "tomato2", "sienna2", "orange","light blue","light blue","blue","blue", "dark blue")) +
  facet_grid(var ~ site) +
  labs(title = NULL, y = "main effects", x = NULL)+
  theme_bw() +
  theme(legend.position = "none",
        legend.text=element_text(size=14),
        text = element_text(size=14),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5))
p

ggsave2("main_effects_grid.png", plot = p, path = path_out)

#interactions
#X2  = cbind(X1[,1]*X1[,2], X1[,1]*X1[,3], X1[,1]*X1[,5], X1[,1]*X1[,6], X1[,2]*X1[,3], X1[,2]*X1[,5], X1[,2]*X1[,6], 
            #X1[,3]*X1[,5], X1[,3]*X1[,6], X1[,5]*X1[,6])
#X[1,i] <- VPDant[i]       ## Also included as squared term
#X[2,i] <- TAant[i]        ## Also included as squared term
#X[3,i] <- PPTant[i]
#X[4,i] <- PAR[Yday[i]]    ## not included in interactions
#X[5,i] <- Sshall_ant[i]
#X[6,i] <- Sdeep_ant[i]
df_p <- sum_IN %>%
  filter(var %in% c("beta2"))
df_p$ID1 <-factor(df_p$ID1 , levels = c("1","2","3","4", "5", "6", "7","8","9","10"), 
                  labels = c("VPD*T", "VPD*P", "VPD*Sshall", "VPD*Sdeep", "T*P", "T*Sshall", "T*Sdeep", "P*Sshall", "P*Sdeep", "Sshall*Sdeep"))

p <- ggplot(data = df_p, aes(x=ID1, y=mean, color=site)) + 
  geom_errorbar(aes(ymax = pc97.5, ymin = pc2.5), position = "dodge") +
  geom_point(pch=21, size = .5) +
  scale_color_manual(values=c("red", "tomato2", "sienna2", "orange","light blue","light blue","blue","blue", "dark blue")) +
  facet_grid(var ~ site) +
  labs(title = NULL, y = "interactive effects", x = NULL)+
  theme_bw() +
  theme(legend.position = "none",
        legend.text=element_text(size=14),
        text = element_text(size=14),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5))
p

ggsave2("int_effects_grid.png", plot = p, path = path_out)


p1 <- sum_IN %>%
  filter(var %in% c("Y", "Y.rep")) %>%
  filter(site %in% c("seg","ses", "wjs", "mpj")) %>%
  #filter(site %in% c("vcp1","vcp2", "vcm1", "vcm2", "vcs")) %>%
  pivot_wider(id_cols = c(site,ID1), names_from = var, values_from = c(mean,median,sd,pc2.5,pc97.5)) %>%
  ggplot(aes(x = mean_Y, y= mean_Y.rep)) +
  #geom_point() +
  geom_pointrange(aes(ymin=pc2.5_Y, ymax=pc97.5_Y), alpha=0.5)+
  geom_smooth(method="lm", se = F, color = "red") +
  geom_abline(slope=1, intercept=0, lty=2, col="blue", size=1.25)+
  stat_cor(aes(label = ..rr.label..), color = "red", label.x = 0.5, size = 3) +
  facet_col("site", strip.position = "right") +
  labs(title = NULL, x="Y", y="Y.rep") +
  #theme_classic(base_size = 12)+
  theme(legend.position = "right",
        legend.text=element_text(size=14),
        text = element_text(size=14),
        legend.title = element_blank(),
        panel.background = element_rect(fill="white"),
        axis.line = element_line(color = "black"),
        axis.text.x = element_text(colour="black"),
        plot.title = element_text(hjust = 0.5))

p1

ggsave2("fit_low.png", plot = p1, path = path_out)

p2 <- sum_IN %>%
  filter(var %in% c("Y", "Y.rep")) %>%
  #filter(site %in% c("seg","ses", "wjs", "mpj")) %>%
  filter(site %in% c("vcp1","vcp2", "vcm1", "vcm2", "vcs")) %>%
  pivot_wider(id_cols = c(site,ID1), names_from = var, values_from = c(mean,median,sd,pc2.5,pc97.5)) %>%
  ggplot(aes(x = mean_Y, y= mean_Y.rep)) +
  #geom_point() +
  geom_pointrange(aes(ymin=pc2.5_Y, ymax=pc97.5_Y), alpha=0.5)+
  geom_smooth(method="lm", se = F, color = "red") +
  geom_abline(slope=1, intercept=0, lty=2, col="blue", size=1.25)+
  stat_cor(aes(label = ..rr.label..), color = "red", label.x = 0.5, size = 3) +
  facet_col("site", strip.position = "right") +
  labs(title = NULL, x="Y", y="Y.rep") +
  #theme_classic(base_size = 12)+
  theme(legend.position = "right",
        legend.text=element_text(size=14),
        text = element_text(size=14),
        legend.title = element_blank(),
        panel.background = element_rect(fill="white"),
        axis.line = element_line(color = "black"),
        axis.text.x = element_text(colour="black"),
        plot.title = element_text(hjust = 0.5))

p2

ggsave2("fit_high.png", plot = p2, path = path_out)

#dYdX[i,1] <- dYdVPD[i]
#dYdX[i,2] <- dYdT[i]
#dYdX[i,3] <- dYdP[i]
#dYdX[i,4] <- dYdSs[i]
#dYdX[i,5] <- dYdSd[i]

############## Net sensitivities

sum_seg = sum_IN_seg %>%
  filter(var %in% c("dYdX")) %>%
  pivot_wider(id_cols = c(site,ID1,ID2), names_from = var, values_from = c(mean,median,sd,pc2.5,pc97.5)) %>%
  filter(ID2 == 1)
dataIN_T = read.csv("./clean_data/d_B_seg.csv")
dataIN_WUE = read.csv("./clean_data/d_B_wue_seg.csv")
dataIN <- left_join(dataIN_T, dataIN_WUE)
sum_seg = sum_seg %>%
  mutate(date = dataIN[184:nrow(dataIN),1])


sum_IN_ses = read.csv("./output_dfs/df_sum_ses.csv")
sum_IN_wjs = read.csv("./output_dfs/df_sum_wjs.csv")
sum_IN_mpj = read.csv("./output_dfs/df_sum_mpj.csv")
sum_IN_vcp1 = read.csv("./output_dfs/df_sum_vcp1.csv")
sum_IN_vcp2 = read.csv("./output_dfs/df_sum_vcp2.csv")
sum_IN_vcm1 = read.csv("./output_dfs/df_sum_vcm1.csv")
sum_IN_vcm2 = read.csv("./output_dfs/df_sum_vcm2.csv")
sum_IN_vcs = read.csv("./output_dfs/df_sum_vcs.csv")

df_p <- sum_IN %>%
  filter(var %in% c("dYdX")) %>%
  #filter(site %in% c("vcp1", "vcp2", "vcm1", "vcm2", "vcs")) %>%
  pivot_wider(id_cols = c(site,ID1,ID2), names_from = var, values_from = c(mean,median,sd,pc2.5,pc97.5))

df_p$ID2 <-factor(df_p$ID2 , levels = c("1","2","3","4", "5"), labels = c("VPD", "T", "P", "Sshall", "Sdeep"))

p <- ggplot(data = df_p, aes(x=ID1, y=mean_dYdX, color=site)) + 
  geom_errorbar(aes(ymax = pc97.5_dYdX, ymin = pc2.5_dYdX), position = "dodge") +
  geom_point(pch=21, size = .5) +
  scale_color_manual(values=c("red", "tomato2", "sienna2", "orange","light blue","light blue","blue","blue", "dark blue")) +
  facet_grid(ID2 ~ site) +
  labs(title = NULL, y = "dWUE/dX", x = "Timeseries index")+
  theme_bw() +
  theme(legend.position = "none",
        legend.text=element_text(size=14),
        text = element_text(size=14),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5))
p

ggsave2("net_sens_grid.png", plot = p, path = path_out)




