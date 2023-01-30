### This script file will read in the SAM model output across sites
### and make summary graphs

# Load packages
library(tidyverse)
library(gridExtra)
library(ggforce) # for facet grids
library(cowplot) # for ggsave
library(ggpubr) # for 1:1 graphs stat_cor

# Load self-made functions
source("./scripts/coda_functions.R")

# Create necessary folders if they do not already exist
if(!file.exists("plots")) { dir.create("plots")}

path_out = "./plots/" # set save path

sum_OUT <- list()
for(i in 1:5){
  
  sum_IN_seg = read.csv(paste("./output_dfs/df_sum_seg_v", i,".csv", sep=""))
  sum_IN_ses = read.csv(paste("./output_dfs/df_sum_ses_v", i,".csv", sep=""))
  sum_IN_wjs = read.csv(paste("./output_dfs/df_sum_wjs_v", i,".csv", sep=""))
  sum_IN_mpj = read.csv(paste("./output_dfs/df_sum_mpj_v", i,".csv", sep=""))
  sum_IN_vcp = read.csv(paste("./output_dfs/df_sum_vcp_v", i,".csv", sep=""))
  sum_IN_vcm1 = read.csv(paste("./output_dfs/df_sum_vcm1_v", i,".csv", sep=""))
  sum_IN_vcm2 = read.csv(paste("./output_dfs/df_sum_vcm2_v", i,".csv", sep=""))
  sum_IN_vcs = read.csv(paste("./output_dfs/df_sum_vcs_v", i,".csv", sep=""))
  
  sum_IN = rbind(sum_IN_seg, sum_IN_ses, sum_IN_wjs, sum_IN_mpj, 
                 sum_IN_vcp, sum_IN_vcm1, sum_IN_vcm2, sum_IN_vcs)
  sum_IN$ID1 = as.numeric(sum_IN$ID1)
  sum_IN$ID2 = as.numeric(sum_IN$ID2)
  sum_IN$site <- factor(sum_IN$site, levels = c("seg", "ses", "wjs", "mpj", "vcp", "vcm1", "vcm2", "vcs"),
                        labels = c("US-Seg", "US-Ses", "US-Wjs", "US-Mpj", "US-Vcp", "US-Vcm1", "US-Vcm2", "US-Vcs"))
  sum_IN$modelv <- i
  sum_OUT[[i]] <- sum_IN
}
sum_IN <- bind_rows(sum_OUT)

##### model fit metrics/ model selection ######
# model 5 seems to be the best

# read in DIC and pD output, calculate WAIC, get Dinf, and make table
sum_list <- list()
k=1 # counter for out list
for(i in 1:5){ # model level
  for(j in 1:8){ # site level
    
    if(j == 1){
      key="seg"
      sitename = "US-Seg"
    } else if(j == 2){
      key="ses"
      sitename = "US-Ses"
    } else if(j == 3){
      key="wjs"
      sitename = "US-Wjs"
    } else if(j == 4){
      key="mpj"
      sitename = "US-Mpj"
    } else if(j == 5){
      key="vcp"
      sitename = "US-Vcp"
    } else if(j == 6){
      key="vcm1"
      sitename = "US-Vcm1"
    } else if(j == 7){
      key="vcm2"
      sitename = "US-Vcm2"
    } else if(j == 8){
      key="vcs"
      sitename = "US-Vcs"
    }
    
    # DIC
    dic_IN = read.csv(paste("./output_dfs/table.dic_", key,"_v", i,".csv", sep=""))
    
    # WAIC
    sum_dx <- sum_IN %>%
      filter(site == sitename) %>%
      filter(modelv == i) %>%
      filter(var == "dx")
    sum_ldx <- sum_IN %>%
      filter(site == sitename) %>%
      filter(modelv == i) %>%
      filter(var == "ldx")
    lppd = sum(log(sum_dx$mean)) # Only get rows associated with “pd”
    pwaic = sum(sum_ldx$sd^2) # Only get rows associated with “lpd” variance
    waic = lppd + pwaic
    
    #Dinf
    sum_dsum <- sum_IN %>%
      filter(site == sitename) %>%
      filter(modelv == i) %>%
      filter(var == "Dsum")
    Dinf = sum_dsum$mean
    
    # R2
    sum_R2 <- sum_IN %>%
      filter(site == sitename) %>%
      filter(modelv == i) %>%
      filter(var == "R2")
    R2 = sum_R2$mean
    
    sum_temp <- dic_IN %>%
      mutate(WAIC = waic, Dinf = Dinf, R2 = R2)
    
    sum_list[[k]] <- sum_temp
    k = k + 1 # update counter
  }
}
sum_stats <- bind_rows(sum_list)
sum_stats$site <- factor(sum_stats$site, levels = c("seg", "ses", "wjs", "mpj", "vcp", "vcm1", "vcm2", "vcs"),
                         labels = c("US-Seg", "US-Ses", "US-Wjs", "US-Mpj", "US-Vcp", "US-Vcm1", "US-Vcm2", "US-Vcs"))
sum_stats <- sum_stats %>% 
  select(-1) %>%
  arrange(site)

write.csv(sum_stats, file = "output_dfs/Model_comparison.csv")  # save for model comparisons

### observed vs predicted
p1 <- sum_IN %>%
  filter(modelv==5) %>%
  filter(var %in% c("ET", "ET.pred")) %>%
  filter(site %in% c("US-Seg","US-Ses", "US-Wjs", "US-Mpj")) %>%
  #filter(site %in% c("US-Vcp", "US-Vcm1", "US-Vcm2", "US-Vcs")) %>%
  pivot_wider(id_cols = c(site,ID1), names_from = var, values_from = c(mean,median,sd,pc2.5,pc97.5)) %>%
  ggplot(aes(x = mean_ET, y= mean_ET.pred)) +
  #geom_point() +
  geom_pointrange(aes(ymin=pc2.5_ET.pred, ymax=pc97.5_ET.pred), alpha=0.5)+
  geom_smooth(method="lm", se = F, color = "red") +
  geom_abline(slope=1, intercept=0, lty=2, col="blue", size=1.25)+
  stat_cor(aes(label = ..rr.label..), color = "red", geom = "label") +
  #stat_cor(aes(label = ..rr.label..), color = "red", label.x = 0.5, size = 3) +
  ylim(0,6) + xlim(0,6) +
  facet_row("site", strip.position = "top") +
  labs(title = NULL, x="observed ET", y="predicted ET") +
  #theme_classic(base_size = 12)+
  theme(legend.position = "right",
        legend.text=element_text(size=14),
        text = element_text(size=14),
        legend.title = element_blank(),
        panel.background = element_rect(fill="white"),
        axis.line = element_line(color = "black"),
        axis.text.x = element_text(colour="black"),
        aspect.ratio=1,
        plot.title = element_text(hjust = 0.5))

p1

ggsave2("fit_low.png", plot = p1, path = path_out)

p2 <- sum_IN %>%
  filter(modelv==5) %>%
  filter(var %in% c("ET", "ET.pred")) %>%
  #filter(site %in% c("seg","ses", "wjs", "mpj")) %>%
  filter(site %in% c("US-Vcp", "US-Vcm1", "US-Vcm2", "US-Vcs")) %>%
  pivot_wider(id_cols = c(site,ID1), names_from = var, values_from = c(mean,median,sd,pc2.5,pc97.5)) %>%
  ggplot(aes(x = mean_ET, y= mean_ET.pred)) +
  #geom_point() +
  geom_pointrange(aes(ymin=pc2.5_ET.pred, ymax=pc97.5_ET.pred), alpha=0.5)+
  geom_smooth(method="lm", se = F, color = "red") +
  geom_abline(slope=1, intercept=0, lty=2, col="blue", size=1.25)+
  stat_cor(aes(label = ..rr.label..), color = "red", geom = "label") +
  #stat_cor(aes(label = ..rr.label..), color = "red", label.x = 0.5, size = 3) +
  ylim(0,6) + xlim(0,6) +
  facet_row("site", strip.position = "top") +
  labs(title = NULL, x="observed ET", y="predicted ET") +
  #theme_classic(base_size = 12)+
  theme(legend.position = "right",
        legend.text=element_text(size=14),
        text = element_text(size=14),
        legend.title = element_blank(),
        panel.background = element_rect(fill="white"),
        axis.line = element_line(color = "black"),
        axis.text.x = element_text(colour="black"),
        aspect.ratio=1,
        plot.title = element_text(hjust = 0.5))
p2
ggsave2("fit_high.png", plot = p2, path = path_out)

p <- grid.arrange(p1,p2, nrow=2)
p

ggsave2("fit.png", plot = p, path = path_out)


################### Determine pairwise differences across seasons ####################

# weights: ID1 = season, ID2 = lag time
pair_list <- list()
k = 1 # counter
for(i in c(4:5)){ # model level
  #print(paste("i:", i, sep=""))
  for( j in 1:8){ # site level
   # print(paste("j:", j, sep=""))
    
    if(j == 1){
      key="seg"
      sitename = "US-Seg"
    } else if(j == 2){
      key="ses"
      sitename = "US-Ses"
    } else if(j == 3){
      key="wjs"
      sitename = "US-Wjs"
    } else if(j == 4){
      key="mpj"
      sitename = "US-Mpj"
    } else if(j == 5){
      key="vcp"
      sitename = "US-Vcp"
    } else if(j == 6){
      key="vcm1"
      sitename = "US-Vcm1"
    } else if(j == 7){
      key="vcm2"
      sitename = "US-Vcm2"
    } else if(j == 8){
      key="vcs"
      sitename = "US-Vcs"
    }
    
    temp_df <- sum_IN %>%
      filter(modelv == i) %>%
      filter(site == sitename) %>%
      filter(var %in% c("wP", "wV", "wSs", "wT"))

    w_pair <- bd.check_similarity(data = temp_df,
                    params = "var",
                    which.params = c("wP","wV", "wSs", "wT"),
                    m="mean",
                    l="pc2.5",
                    u="pc97.5",
                    anchor = "ID2", # anchor by weight lag time
                    ids = "ID1") # id is season (comparing across seasons)

    w_pair <- w_pair %>%
      mutate(modelv = i, site = sitename)


    pair_list[[k]] <- w_pair
    k = k + 1 # update counter
  }
}

pair_list <- bind_rows(pair_list)

sum_w_pair <- sum_IN %>%
  filter(var %in% c("wP", "wV", "wSs", "wT")) %>%
  left_join(pair_list, by = c("ID1", "ID2", "var", "site", "modelv")) %>%
  rename(w_pair = letters)



# main effects: ID1 = parameter type, ID2 = season
pair_list <- list()
k = 1 # counter
for(i in c(2,5)){ # model level
  for( j in 1:8){ # site level
  
  if(j == 1){
    key="seg"
    sitename = "US-Seg"
  } else if(j == 2){
    key="ses"
    sitename = "US-Ses"
  } else if(j == 3){
    key="wjs"
    sitename = "US-Wjs"
  } else if(j == 4){
    key="mpj"
    sitename = "US-Mpj"
  } else if(j == 5){
    key="vcp"
    sitename = "US-Vcp"
  } else if(j == 6){
    key="vcm1"
    sitename = "US-Vcm1"
  } else if(j == 7){
    key="vcm2"
    sitename = "US-Vcm2"
  } else if(j == 8){
    key="vcs"
    sitename = "US-Vcs"
  }
    temp_df <- sum_IN %>%
      filter(modelv == i) %>%
      filter(site == sitename) %>%
      filter(var %in% c("beta1"))
    
    m_pair <- bd.check_similarity(data = temp_df,
                                  params = "var",
                                  which.params = c("beta1"),
                                  m="mean",
                                  l="pc2.5",
                                  u="pc97.5",
                                  anchor = "ID1", # anchor by parameter type
                                  ids = "ID2") # id is season (comparing across seasons)
    m_pair <- m_pair %>%
      mutate(modelv = i, site = sitename)
    
    
    pair_list[[k]] <- m_pair
    k = k + 1 # update counter
  }
}

pair_list <- bind_rows(pair_list)

sum_m_pair <- sum_IN %>%
  filter(var %in% c("beta1")) %>%
  left_join(pair_list, by = c("ID1", "ID2", "var", "site", "modelv")) %>%
  rename(m_pair = letters)

  




# Graph Outputs

# weights
# 1- Spring 2-Summer 3-Fall 4-Winter
df_p <- sum_w_pair %>%
  filter(modelv==5)
df_p$ID1 <-factor(df_p$ID1 , levels = c("1","2","3","4"), labels = c("Spring", "Summer", "Fall", "Winter"))
#old
# C1 = c(0, 1, 2, 4, 6, 13, 20, 27, 55, 83, 111), #stop times for covariates
# C2 = c(0, 1, 2, 3, 5, 7, 14, 21, 28, 56, 84), #start times for covariates
# P1 = c(6, 13, 20, 27, 55, 83, 111, 139, 167), #stop times for precip
# P2 = c(0, 7, 14, 21, 28, 56, 84, 112, 140)) #start times for precip
#new
# C1 = c(0, 1, 2, 4, 6, 13, 20, 27) #stop times for covariates
# C2 = c(0, 1, 2, 3, 5, 7, 14, 21) #start times for covariates
# P1 = c(6, 13, 20, 27, 55, 83, 111, 139, 167) #stop times for precip
# P2 = c(0, 7, 14, 21, 28, 56, 84, 112, 140) #start times for precip
df_p1 <- df_p %>%
  filter(var %in% c("wSd","wSs","wT","wV"))
df_p1$ID2 <-factor(df_p1$ID2 , levels = c("1","2","3","4","5","6","7","8"), 
                  labels = c("day of","previous day","2 days ago","days 3-4 ago","days 5-6 ago","previous week","2 weeks ago","3 weeks ago"))
df_p2 <- df_p %>%
  filter(var %in% c("wP"))
df_p2$ID2 <-factor(df_p2$ID2 , levels = c("1","2","3","4","5","6","7","8","9"), 
                  labels = c("day of","previous week","2 weeks ago","3 weeks ago","1 month ago","2 months ago","3 months ago", "4 months ago", "5 months ago"))
df_p <- rbind(df_p1, df_p2)

p <- df_p %>%
  filter(var %in% c("wP","wSd","wSs","wT","wV")) %>%
  ggplot() + 
  geom_pointrange(aes(x=ID2, y=mean, color=ID1, shape=ID1, ymax = pc97.5, ymin = pc2.5), position = position_dodge(width=0.7), fatten = .5) +
  #geom_point(pch=21, size = .5) +
  scale_color_brewer(palette = "Dark2") +
  facet_grid(site ~ var, scales="free_x") +
  geom_text(aes(x=ID2, y=mean, group=ID1), label=df_p$w_pair, position = position_dodge(width = .7), vjust = -1) +
  labs(title = NULL, y = "weight", x = NULL)+
  theme_bw() +
  theme(legend.position = "top",
        legend.text=element_text(size=12),
        text = element_text(size=12),
        axis.text.x = element_text(size = 12, colour="black", angle = 90, vjust = 0.6, hjust = 1),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5))
p

ggsave2("weights_grid_pairwise.png", plot = p, path = path_out, width = 14, height = 10)

#install.packages("ggh4x")
library(ggh4x)
p <- df_p %>%
  filter(var %in% c("wP","wSd","wSs","wT","wV")) %>%
  ggplot() + 
  geom_pointrange(aes(x=ID2, y=mean, color=ID1, ymax = pc97.5, ymin = pc2.5), position = position_dodge(width = 1), fatten = .5) +
  #geom_point(pch=21, size = .5) +
  scale_color_brewer(palette = "Dark2") +
  #facet_grid(site ~ var) +
  facet_nested(site ~ var + ID1, scales="free_x") + 
  geom_label(aes(x=ID2, y=mean), label=df_p$w_pair, nudge_x = 0.25, nudge_y = 0.25, check_overlap = T) +
  labs(title = NULL, y = "weight", x = NULL)+
  theme_bw() +
  theme(legend.position = "top",
        legend.text=element_text(size=12),
        text = element_text(size=12),
        axis.text.x = element_text(size = 11, angle = 90, vjust = 0.6, hjust = 1),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5))
p

ggsave2("weights_grid2_pairwise.png", plot = p, path = path_out, width = 25, height = 10)

df_p <- sum_m_pair %>%
  filter(modelv==5) %>%
  filter(site!="US-Vcm1") %>%
  filter(var %in% c("beta1")) #%>%
#filter(site %in% c("seg","ses", "wjs", "mpj")) #%>%
#filter(site %in% c("vcp1","vcp2", "vcm1", "vcm2", "vcs"))
df_p$ID1 <-factor(df_p$ID1 , levels = c("1","2","3","4", "5", "6"), labels = c("VPD", "T", "P", "PAR", "Sshall", "Sdeep"))
df_p$ID2 <-factor(df_p$ID2 , levels = c("1","2","3","4"), labels = c("Spring", "Summer", "Fall", "Winter"))

p <- ggplot(data = df_p) + 
  geom_pointrange(aes(x=mean, y=ID1, color=site, xmax = pc97.5, xmin = pc2.5), position = position_dodge(width = 1), fatten = .5) +
  geom_vline(xintercept = 0, linetype = 2) + 
  scale_color_manual(values=c("red", "tomato2", "sienna2", "orange","light blue","blue", "dark blue")) +
  #facet_grid(var ~ site, scales = "free_x") +
  facet_nested(site ~ var + ID2, scales="free_x") + 
  geom_text(aes(x=mean, y=ID1, group=site), label=df_p$m_pair, position = position_dodge(width = .7)) +
  labs(title = NULL, y = "main effects", x = NULL)+
  theme_bw() +
  theme(legend.position = "none",
        legend.text=element_text(size=14),
        text = element_text(size=14),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5))
p

ggsave2("main_effects_grid.png", plot = p, path = path_out, width = 25, height = 10)

p <- ggplot(data = df_p) + 
  geom_pointrange(aes(x=mean, y=ID1, color=ID2, xmax = pc97.5, xmin = pc2.5), position = position_dodge(width = 1), fatten = .5) +
  geom_vline(xintercept = 0, linetype = 2) +
  scale_color_brewer(palette = "Dark2") +
  #scale_color_manual(values=c("red", "tomato2", "sienna2", "orange","light blue","blue", "dark blue")) +
  #facet_grid(var ~ site, scales = "free_x") +
  facet_nested(site ~ var + ID2, scales="free_x") +
  geom_text(aes(x=mean, y=ID1, group=ID2), label=df_p$m_pair, position = position_dodge(width = .7), hjust=-1.2) +
  labs(title = NULL, y = "main effects", x = NULL)+
  theme_bw() +
  theme(legend.position = "none",
        legend.text=element_text(size=14),
        text = element_text(size=14),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5))
p

ggsave2("main_effects_grid_season_pairwise.png", plot = p, path = path_out, width = 25, height = 10)

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
  filter(modelv==5) %>%
  filter(site!="US-Vcm1") %>%
  filter(var %in% c("beta2"))
df_p$ID1 <-factor(df_p$ID1 , levels = c("1","2","3","4", "5", "6", "7","8","9","10"), 
                  labels = c("VPD*T", "VPD*P", "VPD*Sshall", "VPD*Sdeep", "T*P", "T*Sshall", "T*Sdeep", "P*Sshall", "P*Sdeep", "Sshall*Sdeep"))

p <- ggplot(data = df_p) + 
  geom_pointrange(aes(x=mean, y=ID1, color=site, xmax = pc97.5, xmin = pc2.5), position = position_dodge(width = 1), fatten = .5) +
  geom_vline(xintercept = 0, linetype = 2) +
  scale_color_manual(values=c("red", "tomato2", "sienna2", "orange","light blue","blue", "dark blue")) +
  facet_grid(var ~ site, scales = "free_x") +
  #xlim(-0.3,0.3) +
  labs(title = NULL, y = "interactive effects", x = NULL)+
  theme_bw() +
  theme(legend.position = "none",
        legend.text=element_text(size=14),
        text = element_text(size=14),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5))
p

ggsave2("int_effects_grid.png", plot = p, path = path_out, width = 25, height = 10)


#dYdX[i,1] <- dYdVPD[i]
#dYdX[i,2] <- dYdT[i]
#dYdX[i,3] <- dYdP[i]
#dYdX[i,4] <- dYdSs[i]
#dYdX[i,5] <- dYdSd[i]

############## Net sensitivities

df_p <- sum_IN %>%
  filter(var %in% c("dYdX")) %>%
  pivot_wider(id_cols = c(site,ID1,ID2), names_from = var, values_from = c(mean,median,sd,pc2.5,pc97.5))

df_p$ID2 <-factor(df_p$ID2 , levels = c("1","2","3","4", "5"), labels = c("VPD", "T", "P", "Sshall", "Sdeep"))

p <- ggplot(data = df_p, aes(x=ID1, y=mean_dYdX, color=site)) + 
  geom_errorbar(aes(ymax = pc97.5_dYdX, ymin = pc2.5_dYdX), position = "dodge") +
  geom_point(pch=21, size = .5) +
  scale_color_manual(values=c("red", "tomato2", "sienna2", "orange","light blue","blue","blue", "dark blue")) +
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



############## WUE uncertainty

p <- sum_IN %>%
  filter(modelv==5) %>%
  filter(site %in% c("US-Seg", "US-Ses", "US-Wjs", "US-Mpj")) %>%
  filter(var == "WUE.pred") %>%
  ggplot(aes(x = ID1, y= mean)) +
  geom_pointrange(aes(ymin=pc2.5, ymax=pc97.5), alpha=0.5, fatten =.2)+
  facet_col("site", strip.position = "right") +
  labs(title = NULL, x="day", y="predicted WUE") +
  ylim(0,20) +
  xlim(250,850) +
  #theme_classic(base_size = 12)+
  theme(legend.position = "right",
        legend.text=element_text(size=14),
        text = element_text(size=14),
        legend.title = element_blank(),
        panel.background = element_rect(fill="white"),
        axis.line = element_line(color = "black"),
        axis.text.x = element_text(colour="black"),
        plot.title = element_text(hjust = 0.5))

p

ggsave2("WUE_uncertainty_low.png", plot = p, path = path_out)


p <- sum_IN %>%
  filter(modelv==5) %>%
  filter(site %in% c("US-Vcp", "US-Vcm2", "US-Vcs")) %>%
  filter(var == "WUE.pred") %>%
  ggplot(aes(x = ID1, y= mean)) +
  geom_pointrange(aes(ymin=pc2.5, ymax=pc97.5), alpha=0.5, fatten =.2)+
  facet_col("site", strip.position = "right") +
  labs(title = NULL, x="day", y="predicted WUE") +
  ylim(0,20) +
  xlim(250,850) +
  #theme_classic(base_size = 12)+
  theme(legend.position = "right",
        legend.text=element_text(size=14),
        text = element_text(size=14),
        legend.title = element_blank(),
        panel.background = element_rect(fill="white"),
        axis.line = element_line(color = "black"),
        axis.text.x = element_text(colour="black"),
        plot.title = element_text(hjust = 0.5))

p

ggsave2("WUE_uncertainty_high.png", plot = p, path = path_out)


p <- sum_IN %>%
  filter(modelv==5) %>%
  filter(var == "WUE.pred") %>%
  ggplot(aes(site, mean)) + 
  ggdist::stat_halfeye(#aes(fill=Season),
                         alpha = 0.45,
                         ## custom bandwidth
                         adjust = .5, 
                         ## adjust height
                         width = 1.2, 
                         ## move geom to the right
                         justification = -.1, 
                         ## remove slab interval
                         .width = 0, 
                         point_colour = NA) +
    geom_boxplot(
      width = .15) +
    scale_fill_brewer(palette = "Dark2")+
  ylim(0,30) +
    labs(title = NULL, y = expression(paste("WUE (g C / mm ", H[2], "O)")), x = NULL, fill = NULL) +
    theme(legend.position = c(.17,0.8),
          legend.text=element_text(size=16),
          text = element_text(size=16),
          panel.background = element_rect(fill="white"),
          axis.line = element_line(color = "black"),
          axis.text.x = element_text(size = 16, colour="black", angle = 12, vjust = 0.7, hjust = 0.6),
          plot.title = element_text(hjust = 0.5))
p

ggsave2("p_WUE_raincloud.png", plot = p, path = path_out)







