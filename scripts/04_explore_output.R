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


### Create date dataframes ###

date_list <- list()
for(j in 1:8){ # number of sites
  
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

# Load data for the correct site/key
load(paste("./clean_data/dataIN_",key,".RData",sep=""))

# define df names based on key
date_IN <- get(paste("dataIN_",key,sep="")) # daily time series

date_IN <- date_IN %>%
  mutate(Tperiod = ifelse(Season == "Spring", 1, NA)) %>%
  mutate(Tperiod = ifelse(Season == "Summer", 2, Tperiod)) %>%
  mutate(Tperiod = ifelse(Season == "Fall", 3, Tperiod)) %>%
  mutate(Tperiod = ifelse(Season == "Winter", 4, Tperiod)) %>%
  mutate(site = sitename) %>%
  rowid_to_column("dayind") %>%
  select(date, year, month, day, month_name, water_year, Season, dayind, P, ET, GPP, VPD,Tair,S,Smid,Sdeep,P_acc,PPFD_IN, site)

# Growing season starts and stops for each year
Gstart = date_IN %>%
  filter(month == 4 & day == 1) %>%
  select(dayind)

# If running the burned site, start later to accommodate January data start
if(key == "vcm2"){
  Gstart = date_IN %>%
    filter(month == 7 & day == 1) %>%
    select(dayind)}

# Get rid of pre-April starts (Oct and Nov)
YIN = date_IN[Gstart$dayind[1]:nrow(date_IN),]
# Filter out winter
YIN <- YIN %>% filter(Season != "Winter")
# Make dayind match output
YIN$dayind = c(Gstart$dayind[1]:(nrow(YIN) + Gstart$dayind[1] - 1))


date_list[[j]] <- YIN
}
date_OUT <- bind_rows(date_list)

### Read in model output ###

version_num <- c(1,2,2.5,3,4,5,5.5,6,6.5,7)
version_length <- length(version_num)
sum_OUT <- list()
for(i in 1:version_length){ # number of model versions
  
  sum_IN_seg = read.csv(paste("./output_dfs/df_sum_seg_v", version_num[i],".csv", sep=""))
  sum_IN_ses = read.csv(paste("./output_dfs/df_sum_ses_v", version_num[i],".csv", sep=""))
  sum_IN_wjs = read.csv(paste("./output_dfs/df_sum_wjs_v", version_num[i],".csv", sep=""))
  sum_IN_mpj = read.csv(paste("./output_dfs/df_sum_mpj_v", version_num[i],".csv", sep=""))
  sum_IN_vcp = read.csv(paste("./output_dfs/df_sum_vcp_v", version_num[i],".csv", sep=""))
  sum_IN_vcm1 = read.csv(paste("./output_dfs/df_sum_vcm1_v", version_num[i],".csv", sep=""))
  sum_IN_vcm2 = read.csv(paste("./output_dfs/df_sum_vcm2_v", version_num[i],".csv", sep=""))
  sum_IN_vcs = read.csv(paste("./output_dfs/df_sum_vcs_v", version_num[i],".csv", sep=""))
  
  sum_IN = rbind(sum_IN_seg, sum_IN_ses, sum_IN_wjs, sum_IN_mpj, 
                 sum_IN_vcp, sum_IN_vcm1, sum_IN_vcm2, sum_IN_vcs)
  sum_IN$ID1 = as.numeric(sum_IN$ID1)
  sum_IN$ID2 = as.numeric(sum_IN$ID2)
  sum_IN$site <- factor(sum_IN$site, levels = c("seg", "ses", "wjs", "mpj", "vcp", "vcm1", "vcm2", "vcs"),
                        labels = c("US-Seg", "US-Ses", "US-Wjs", "US-Mpj", "US-Vcp", "US-Vcm1", "US-Vcm2", "US-Vcs"))
  sum_IN$modelv <- version_num[i]
  sum_OUT[[i]] <- sum_IN
}
sum_IN <- bind_rows(sum_OUT)

##### model fit metrics/ model selection ######

# read in DIC and pD output, calculate WAIC, get Dinf, and make table
sum_list <- list()
k=1 # counter for out list
for(i in 1:version_length){ # model version level
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
    dic_IN = read.csv(paste("./output_dfs/table.dic_", key,"_v", version_num[i],".csv", sep=""))
    
    # WAIC
    sum_dx <- sum_IN %>%
      filter(site == sitename) %>%
      filter(modelv == version_num[i]) %>%
      filter(var == "dx")
    sum_ldx <- sum_IN %>%
      filter(site == sitename) %>%
      filter(modelv == version_num[i]) %>%
      filter(var == "ldx")
    lppd = sum(log(sum_dx$mean)) # Only get rows associated with “pd”
    pwaic = sum(sum_ldx$sd^2) # Only get rows associated with “lpd” variance
    waic = lppd + pwaic
    
    #Dinf
    sum_dsum <- sum_IN %>%
      filter(site == sitename) %>%
      filter(modelv == version_num[i]) %>%
      filter(var == "Dsum")
    Dinf = sum_dsum$mean
    
    # R2
    sum_R2 <- sum_IN %>%
      filter(site == sitename) %>%
      filter(modelv == version_num[i]) %>%
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
model_version=3
p1 <- sum_IN %>%
  filter(modelv==model_version) %>%
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

#ggsave2("fit_low.png", plot = p1, path = path_out)

p2 <- sum_IN %>%
  filter(modelv==model_version) %>%
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
#ggsave2("fit_high.png", plot = p2, path = path_out)

p <- grid.arrange(p1,p2, nrow=2)
p

ggsave2(paste("fit","_v", model_version,".png",sep=""), plot = p, path = path_out)


################### Determine pairwise differences across seasons ####################

# weights: ID1 = season, ID2 = lag time
version_num <- c(4,5,5.5)
version_length <- length(version_num)
pair_list <- list()
k = 1 # counter
for(i in version_num){ # model level
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

pair_list1 <- bind_rows(pair_list)

version_num <- c(4,5,5.5,6,6.5,7)
version_length <- length(version_num)
pair_list <- list()
k = 1 # counter
for(i in version_num){ # model level
  print(paste("i:", i, sep=""))
  for( j in 1:8){ # site level
    print(paste("j:", j, sep=""))
    
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
    
    if(i %in% c(4,5,5.5)){
      w_pair <- bd.check_similarity(data = temp_df,
                                    params = "var",
                                    which.params = c("wP","wV", "wSs", "wT"),
                                    m="mean",
                                    l="pc2.5",
                                    u="pc97.5",
                                    anchor = "ID1", # anchor by season
                                    ids = "ID2") # id weight lag time (comparing across lag times)
    }else{
      w_pair <- bd.check_similarity(data = temp_df,
                                    params = "var",
                                    which.params = c("wP","wV", "wSs", "wT"),
                                    m="mean",
                                    l="pc2.5",
                                    u="pc97.5",
                                    anchor = "site", # anchor by season
                                    ids = "ID1") # id weight lag time (comparing across lag times)
    }
  
    
    w_pair <- w_pair %>%
      mutate(modelv = i, site = sitename)
    
    
    pair_list[[k]] <- w_pair
    k = k + 1 # update counter
  }
}

pair_list2 <- bind_rows(pair_list)

sum_w_pair <- sum_IN %>%
  filter(var %in% c("wP", "wV", "wSs", "wT")) %>%
  left_join(pair_list1, by = c("ID1", "ID2", "var", "site", "modelv")) %>%
  rename(w_season_letters = letters) %>%
  left_join(pair_list2, by = c("ID1", "ID2", "var", "site", "modelv")) %>%
  rename(w_letters = letters)

sum_w_pair$site <- factor(sum_w_pair$site, levels = c("US-Seg", "US-Ses", "US-Wjs", "US-Mpj", "US-Vcp", "US-Vcm1", "US-Vcm2", "US-Vcs"))



# main effects: ID1 = parameter type, ID2 = season
version_num <- c(2,2.5,5,5.5,6,6.5,7)
version_length <- length(version_num)
pair_list <- list()
k = 1 # counter
for(i in version_num){ # model level
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
    
    if(i %in% c(4,5,5.5)){
    m_pair <- bd.check_similarity(data = temp_df,
                                  params = "var",
                                  which.params = c("beta1"),
                                  m="mean",
                                  l="pc2.5",
                                  u="pc97.5",
                                  anchor = "ID2", # anchor by season
                                  ids = "ID1") # id is parameter 
    } else{
      m_pair <- bd.check_similarity(data = temp_df,
                                    params = "var",
                                    which.params = c("beta1"),
                                    m="mean",
                                    l="pc2.5",
                                    u="pc97.5",
                                    anchor = "ID1", # anchor by parameter type
                                    ids = "ID2") # id is season (comparing across seasons)
    }
    m_pair <- m_pair %>%
      mutate(modelv = i, site = sitename)
    
    
    pair_list[[k]] <- m_pair
    k = k + 1 # update counter
  }
}

pair_list1 <- bind_rows(pair_list)


sum_m_pair <- sum_IN %>%
  filter(var %in% c("beta1")) %>%
  left_join(pair_list1, by = c("ID1", "ID2", "var", "site", "modelv")) %>%
  rename(m_season_letters = letters)

sum_m_pair$site <- factor(sum_m_pair$site, levels = c("US-Seg", "US-Ses", "US-Wjs", "US-Mpj", "US-Vcp", "US-Vcm1", "US-Vcm2", "US-Vcs"))





# Graph Outputs
model_version=3

# weights
# 1- Spring 2-Summer 3-Fall 4-Winter
df_p <- sum_w_pair %>%
  filter(modelv==model_version) %>%
  filter(site!="US-Vcm1") %>%
  filter(var %in% c("wP","wSd","wSs","wT","wV"))
  
if(model_version %in% c(4,5,5.5)){
  df_p$ID1 <-factor(df_p$ID1 , levels = c("1","2","3"), labels = c("Spring", "Summer", "Fall"))
}

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
if(model_version %in% c(4,5,5.5)){
df_p1$ID2 <-factor(df_p1$ID2 , levels = c("1","2","3","4","5","6","7","8"), 
                  labels = c("day of","previous day","2 days ago","3-4 days ago","5-6 days ago","previous week","2 weeks ago","3 weeks ago"))
}else{
  df_p1$ID1 <-factor(df_p1$ID1 , levels = c("1","2","3","4","5","6","7","8"), 
                     labels = c("day of","previous day","2 days ago","3-4 days ago","5-6 days ago","previous week","2 weeks ago","3 weeks ago"))
}
df_p2 <- df_p %>%
  filter(var %in% c("wP"))
if(model_version %in% c(4,5,5.5)){
df_p2$ID2 <-factor(df_p2$ID2 , levels = c("1","2","3","4","5","6","7","8","9"), 
                  labels = c("day of","previous week","2 weeks ago","3 weeks ago","1 month ago","2 months ago","3 months ago", "4 months ago", "5 months ago"))
}else{
  df_p2$ID1 <-factor(df_p2$ID1 , levels = c("1","2","3","4","5","6","7","8","9"), 
                     labels = c("day of","previous week","2 weeks ago","3 weeks ago","1 month ago","2 months ago","3 months ago", "4 months ago", "5 months ago"))
}
df_p <- rbind(df_p1, df_p2)

df_p$var <- factor(df_p$var, levels = c("wP","wSd","wSs","wT","wV"))

if(model_version %in% c(4,5,5.5)){
p <- df_p %>%
  filter(var %in% c("wP","wSs","wT","wV")) %>%
  ggplot(aes(x=ID2, y=mean, group=ID1)) + 
  geom_pointrange(aes(color=ID1, ymax = pc97.5, ymin = pc2.5), position = position_dodge(width=0.7), fatten = .5) +
  scale_color_brewer(palette = "Dark2") +
  facet_grid(site ~ var, scales="free_x") +
  geom_text(data = df_p, aes(label=w_season_letters, color=ID1), position = position_dodge(width = .7), vjust = -3, size=3) +
  geom_text(data = df_p, aes(label=w_letters), position = position_dodge(width = .7), vjust = -4, size=3) +
  labs(title = NULL, y = "weight", x = NULL)+
  ylim(0,1.4) +
  theme_bw() +
  theme(legend.position = "top",
        legend.text=element_text(size=12),
        text = element_text(size=12),
        axis.text.x = element_text(size = 12, colour="black", angle = 90, vjust = 0.6, hjust = 1),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5))
}else{
  #df_p <- df_p %>%
    #filter(var %in% c("wT","wV"))
  p <- df_p %>%
    ggplot(aes(x=ID1, y=mean)) + 
    geom_pointrange(aes(color=site, ymax = pc97.5, ymin = pc2.5), fatten = .5) +
    scale_color_manual(values=c("red", "tomato2", "sienna2", "orange","light blue","blue", "dark blue")) +
    #scale_color_brewer(palette = "Dark2") +
    facet_grid(site ~ var, scales="free_x") +
    geom_text(data = df_p, aes(label=w_letters),vjust = -0.2, hjust = -0.6, size=3) +
    labs(title = NULL, y = "weight", x = NULL)+
    ylim(0,1) +
    theme_bw() +
    theme(legend.position = "none",
          legend.text=element_text(size=12),
          text = element_text(size=12),
          axis.text.x = element_text(size = 12, colour="black", angle = 90, vjust = 0.6, hjust = 1),
          legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5))
}

p

ggsave2(paste("weights_grid_pairwise","_v", model_version,".png",sep=""), plot = p, path = path_out, width = 10, height = 7)

#install.packages("ggh4x")
library(ggh4x)
if(model_version %in% c(4,5,5.5)){
p <- df_p %>%
  filter(var %in% c("wP","wSd","wSs","wT","wV")) %>%
  ggplot(aes(x=ID2, y=mean, group=ID1)) + 
  geom_pointrange(aes(color=ID1, ymax = pc97.5, ymin = pc2.5), position = position_dodge(width = 1), fatten = .5) +
  #geom_point(pch=21, size = .5) +
  scale_color_brewer(palette = "Dark2") +
  #facet_grid(site ~ var) +
  facet_nested(site ~ var + ID1, scales="free_x") + 
  geom_text(data = df_p, aes(label=w_letters), position = position_dodge(width = .7), vjust = -3, size=3) +
  #geom_label(data = df_p, aes(label=w_letters), position = position_dodge(width = .7), size=3) +
  labs(title = NULL, y = "weight", x = NULL)+
  theme_bw() +
  theme(legend.position = "top",
        legend.text=element_text(size=12),
        text = element_text(size=12),
        axis.text.x = element_text(size = 11, angle = 90, vjust = 0.6, hjust = 1),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5))
p

ggsave2(paste("weights_grid2_pairwise","_v", model_version,".png",sep=""), plot = p, path = path_out, width = 25, height = 10)
}

df_p <- sum_m_pair %>%
  filter(modelv==model_version) %>%
  filter(site!="US-Vcm1") %>%
  filter(var %in% c("beta1")) #%>%
#filter(site %in% c("seg","ses", "wjs", "mpj")) #%>%
#filter(site %in% c("vcp1","vcp2", "vcm1", "vcm2", "vcs"))
df_p$ID1 <-factor(df_p$ID1 , levels = c("1","2","3","4", "5", "6"), labels = c("VPD", "Tair", "P", "PAR", "Sshall", "Sdeep"))
df_p$ID2 <-factor(df_p$ID2 , levels = c("1","2","3"), labels = c("Spring", "Summer", "Fall"))


p <- ggplot(data = df_p, aes(x=mean, y=ID1, group=site)) + 
  geom_pointrange(aes(color=site, xmax = pc97.5, xmin = pc2.5, alpha = ifelse(pc2.5 <= 0 & pc97.5 >= 0, .9, 1)), position = position_dodge(width = .7), fatten = .5) +
  geom_vline(xintercept = 0, linetype = 2) + 
  scale_color_manual(values=c("red", "tomato2", "sienna2", "orange","light blue","blue", "dark blue")) +
  #facet_grid(var ~ site, scales = "free_x") +
  facet_nested(site ~ var + ID2, scales="free_x") + 
  #geom_text(aes(x=mean, y=ID1, group=site), label=df_p$m_pair, position = position_dodge(width = .7)) +
  geom_text(data = df_p, aes(label=m_season_letters), position = position_dodge(width = .7), vjust = -.5, size=3) +
  #xlim(-10,10) +
  labs(title = NULL, y = "main effects", x = NULL)+
  theme_bw() +
  theme(legend.position = "none",
        legend.text=element_text(size=10),
        text = element_text(size=10),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5))
p

ggsave2(paste("main_effects_grid","_v", model_version,".png",sep=""), plot = p, path = path_out, width = 10, height = 8)

#df_p <- df_p %>%
#  filter(ID1 != "PAR") 

p <- df_p %>%
  ggplot(aes(x=mean,y=var,group=ID2)) + 
  geom_pointrange(aes(color=ID2, xmax = pc97.5, xmin = pc2.5, alpha = ifelse(pc2.5 <= 0 & pc97.5 >= 0, "nonsignificant", "significant")), position = position_dodge(width = 0.5), fatten = .5) +
  scale_alpha_discrete(range = c(.3,1)) +
  geom_vline(xintercept = 0, linetype = 2) +
  scale_color_brewer(palette = "Dark2") +
  #scale_color_manual(values=c("red", "tomato2", "sienna2", "orange","light blue","blue", "dark blue")) +
  facet_grid(ID1 ~ site, scales = "free") +
  #facet_nested(site ~ var + ID2, scales="free_x") +
  geom_text(data = df_p, aes(label=m_season_letters, color=ID2, alpha = ifelse(pc2.5 <= 0 & pc97.5 >= 0, "nonsignificant", "significant")), position = position_dodge(width = .7), hjust=-1.2) +
  #xlim(-2,2) +
  labs(title = NULL, y = "main effects", x = NULL)+
  theme_bw() +
  theme(legend.position = "top",
        legend.text=element_text(size=14),
        text = element_text(size=14),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.text.y=element_blank(),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5))
p

ggsave2(paste("main_effects_grid_season_pairwise","_v", model_version,".png",sep=""), plot = p, path = path_out, width = 10, height = 6)

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
  filter(modelv==model_version) %>%
  filter(site!="US-Vcm1") %>%
  filter(var %in% c("beta2"))
df_p$ID1 <-factor(df_p$ID1 , levels = c("1","2","3","4", "5", "6", "7","8","9","10"), 
                  labels = c("VPD*T", "VPD*P", "VPD*Sshall", "VPD*Sdeep", "T*P", "T*Sshall", "T*Sdeep", "P*Sshall", "P*Sdeep", "Sshall*Sdeep"))

p <- ggplot(data = df_p) + 
  geom_pointrange(aes(x=mean, y=ID1, xmax = pc97.5, xmin = pc2.5, color = ifelse(pc2.5 <= 0 & pc97.5 >= 0, "nonsignificant", "significant")), position = position_dodge(width = 1), fatten = .5) +
  geom_vline(xintercept = 0, linetype = 2) +
  scale_color_manual(values=c("gray", "black")) +
  facet_row("site", scales = "free_x") +
  #xlim(-0.3,0.3) +
  labs(title = NULL, y = "interactive effects", x = NULL)+
  theme_bw() +
  theme(legend.position = "top",
        legend.text=element_text(size=14),
        text = element_text(size=14),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5))
p

ggsave2(paste("int_effects_grid","_v", model_version,".png",sep=""), plot = p, path = path_out, width = 8, height = 6)


#dYdX[i,1] <- dYdVPD[i]
#dYdX[i,2] <- dYdT[i]
#dYdX[i,3] <- dYdP[i]
#dYdX[i,4] <- dYdSs[i]
#dYdX[i,5] <- dYdSd[i]

############## Net sensitivities

date_temp <- as.data.frame(date_OUT[which(date_OUT$site !="US-Vcm1"),])
date_temp <- date_temp %>%
  rename(ID1 = dayind) %>%
  filter(Season != "Winter")

df_p <- sum_IN %>%
  filter(modelv==model_version) %>%
  filter(site!="US-Vcm1") %>%
  filter(var %in% c("dYdX")) %>%
  pivot_wider(id_cols = c(site,ID1,ID2), names_from = var, values_from = c(mean,median,sd,pc2.5,pc97.5))


temp_WUE <- sum_IN %>%
  filter(modelv==model_version) %>%
  filter(site!="US-Vcm1") %>%
  filter(var %in% c("WUE.pred")) %>%
  pivot_wider(id_cols = c(site,ID1), names_from = var, values_from = c(mean,median,sd,pc2.5,pc97.5))

df_p <- left_join(df_p, temp_WUE, by = c("ID1", "site"))
#df_p$mean_dYdX <- df_p$mean_dYdX * df_p$mean_WUE.pred
#df_p$pc2.5_dYdX <- df_p$pc2.5_dYdX * df_p$pc2.5_WUE.pred
#df_p$pc97.5_dYdX <- df_p$pc97.5_dYdX * df_p$pc97.5_WUE.pred
#

df_p$ID2 <-factor(df_p$ID2 , levels = c("1","2","3","4", "5", "6"), labels = c("VPD", "Tair", "P", "Sshall", "Sdeep", "PAR"))

df_p <- left_join(df_p, date_temp, by = c("ID1", "site"))

#df_p <- df_p %>%
#  filter(Season != "Winter") 

df_p$date <- as.Date(df_p$date)
df_p$site <- factor(df_p$site, levels = c("US-Seg", "US-Ses", "US-Wjs", "US-Mpj", "US-Vcp", "US-Vcm2", "US-Vcs"))
df_p$Season <- factor(df_p$Season, levels = c("Spring", "Summer", "Fall"))


k=0.5
df_p$mean_dYdX_opt <- (df_p$GPP/df_p$ET)*k*df_p$VPD^(k-1)
k=0.9
df_p$mean_dYdX_subopt <- (df_p$GPP/df_p$ET)*k*df_p$VPD^(k-1)


p <- df_p %>%
  filter(year>2017) %>%
  filter(year<2020) %>%
  ggplot(aes(x=date, y=mean_dYdX)) +
  geom_pointrange(aes(ymax = pc97.5_dYdX, ymin = pc2.5_dYdX, color = ifelse(pc2.5_dYdX <= 0 & pc97.5_dYdX >= 0, "nonsignificant", "significant")), position = position_dodge(width = 1), fatten = .5) +
  scale_color_manual(values = c("gray", "black")) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_rect(aes(xmin = date, xmax = dplyr::lead(date), ymin = -Inf, ymax = Inf, fill = Season), alpha = 0.1) +
  scale_fill_brewer(palette = "Dark2") +
  facet_grid(ID2 ~ site, scales = "free") +
  labs(title = NULL, y = "dWUE/dX", x = "Year")+
  theme_bw() +
  theme(legend.position = "top",
        legend.text=element_text(size=14),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        text = element_text(size=14),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5))
p

ggsave2(paste("net_sens_grid_supp","_v", model_version,".png",sep=""), plot = p, path = path_out, width = 25, height = 15)


p <- df_p %>%
  filter(ID2 %in% c("VPD", "Tair")) %>%
  ggplot(aes(x=date, y=mean_dYdX)) +
  geom_pointrange(aes(ymax = pc97.5_dYdX, ymin = pc2.5_dYdX, color = ifelse(pc2.5_dYdX <= 0 & pc97.5_dYdX >= 0, "nonsignificant", "significant")), position = position_dodge(width = 1), fatten = .5) +
  scale_color_manual(values = c("gray", "black")) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_rect(aes(xmin = date, xmax = dplyr::lead(date), ymin = -Inf, ymax = Inf, fill = Season), alpha = 0.1) +
  scale_fill_brewer(palette = "Dark2") +
  facet_grid(site ~ ID2) +
  labs(title = NULL, y = "dWUE/dX", x = "Year")+
  #ylim(-100,100) +
  theme_bw() +
  theme(legend.position = "top",
        legend.text=element_text(size=14),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        text = element_text(size=14),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5))
p

ggsave2(paste("net_sens_grid","_v", model_version,".png",sep=""), plot = p, path = path_out, width = 25, height = 15)

# net sens graph for presentations
library(scales)
library(bdscale)
library(patchwork)
library(grid)
site_list <- c("US-Ses", "US-Wjs", "US-Vcp")
p1 <- list()
for(i in c(1:3)){
df_p1 <- df_p %>%
  filter(ID2 == "VPD") %>%
  filter(year > 2016) %>%
  filter(year < 2018) %>%
  filter(site %in% c(site_list[i]))

p1[[i]] <- df_p1 %>%
  ggplot(aes(x=date)) +
  geom_pointrange(aes(y = mean_dYdX, ymax = pc97.5_dYdX, ymin = pc2.5_dYdX, color = ifelse(pc2.5_dYdX <= 0 & pc97.5_dYdX >= 0, "nonsignificant", "significant")), position = position_dodge(width = 1), fatten = .5, alpha=1) +
  #geom_line(aes(y = pc97.5_dYdX, color = ifelse(pc2.5_dYdX <= 0 & pc97.5_dYdX >= 0, "nonsignificant", "significant")), alpha=0.2) +
  #geom_point(aes(y = mean_dYdX, color = ifelse(pc2.5_dYdX <= 0 & pc97.5_dYdX >= 0, "nonsignificant", "significant")), alpha=1) +
  #geom_line(aes(y = pc2.5_dYdX, color = ifelse(pc2.5_dYdX <= 0 & pc97.5_dYdX >= 0, "nonsignificant", "significant")), alpha=0.2) +
  geom_point(aes(y=mean_dYdX_opt, color = "optimal net sensitivity"), alpha=1, size=.5) +
  geom_point(aes(y=mean_dYdX_subopt, color = "suboptimal net sensitivity"), alpha=1, size=.5) +
  scale_color_manual(values = c("black", "gray","red", "blue"), breaks = c("significant", "nonsignificant", "optimal net sensitivity", "suboptimal net sensitivity")) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_rect(aes(xmin = date, xmax = dplyr::lead(date), ymin = -Inf, ymax = Inf, fill = Season), alpha = 0.1) +
  scale_fill_brewer(palette = "Dark2") +
  facet_col("site", scales = "free_y", strip.position = "right") +
  scale_x_bd(business.dates=df_p1$date, max.major.breaks=3, labels = date_format("%Y")) + # labels=date_format("%Y") # labels = function(x) { x <- as.character(x); x[!1:0] <- " "; x}
  labs(title = NULL, y = NULL, x = NULL)+
  ylim(-1,10) +
  theme_bw() +
  theme(legend.position = "right",
        legend.text=element_text(size=10),
        #axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.title.x = element_blank(),
        text = element_text(size=14),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5))

}


layout <- '
A
B
C
'
p <- wrap_plots(A = p1[[1]], B = p1[[2]], C = p1[[3]], design = layout)
p <- p + plot_layout(guides = "collect") & theme(legend.position = 'right')
p <- grid.arrange(patchworkGrob(p), left = textGrob("dWUE/dVPD", rot = 90, gp=gpar(fontsize=14)), bottom = textGrob("Year", gp=gpar(fontsize=14)))


ggsave2(paste("net_sens_VPD1","_v", model_version,".png",sep=""), plot = p, path = path_out, width = 7, height = 5)



site_list <- c("US-Ses", "US-Wjs", "US-Vcp")
p1 <- list()
for(i in c(1:3)){
  df_p1 <- df_p %>%
    filter(ID2 == "VPD") %>%
    filter(year > 2016) %>%
    filter(site %in% c(site_list[i]))
  
  p1[[i]] <- df_p1 %>%
    ggplot(aes(x=date)) +
    geom_pointrange(aes(y = mean_dYdX, ymax = pc97.5_dYdX, ymin = pc2.5_dYdX, color = ifelse(pc2.5_dYdX <= 0 & pc97.5_dYdX >= 0, "nonsignificant", "significant")), position = position_dodge(width = 1), fatten = .5, alpha=0.2) +
    #geom_line(aes(y = pc97.5_dYdX, color = ifelse(pc2.5_dYdX <= 0 & pc97.5_dYdX >= 0, "nonsignificant", "significant")), alpha=0.2) +
    #geom_point(aes(y = mean_dYdX, color = ifelse(pc2.5_dYdX <= 0 & pc97.5_dYdX >= 0, "nonsignificant", "significant")), alpha=1) +
    #geom_line(aes(y = pc2.5_dYdX, color = ifelse(pc2.5_dYdX <= 0 & pc97.5_dYdX >= 0, "nonsignificant", "significant")), alpha=0.2) +
    geom_point(aes(y=mean_dYdX_opt, color = "optimal net sensitivity"), alpha=0.2, size=.5) +
    geom_point(aes(y=mean_dYdX_subopt, color = "suboptimal net sensitivity"), alpha=0.2, size=.5) +
    scale_color_manual(values = c("black", "gray","red", "blue"), breaks = c("significant", "nonsignificant", "optimal net sensitivity", "suboptimal net sensitivity")) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_rect(aes(xmin = date, xmax = dplyr::lead(date), ymin = -Inf, ymax = Inf, fill = Season), alpha = 0.1) +
    scale_fill_brewer(palette = "Dark2") +
    facet_col("site", scales = "free_y", strip.position = "right") +
    scale_x_bd(business.dates=df_p1$date, max.major.breaks=7, labels = date_format("%Y")) + # labels=date_format("%Y") # labels = function(x) { x <- as.character(x); x[!1:0] <- " "; x}
    labs(title = NULL, y = NULL, x = NULL)+
    ylim(-1,10) +
    theme_bw() +
    theme(legend.position = "top",
          legend.text=element_text(size=10),
          #axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          axis.title.x = element_blank(),
          text = element_text(size=14),
          legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5))
  
}


layout <- '
A
B
C
'
p <- wrap_plots(A = p1[[1]], B = p1[[2]], C = p1[[3]], design = layout)
p <- p + plot_layout(guides = "collect") & theme(legend.position = 'top')
p <- grid.arrange(patchworkGrob(p), left = textGrob("dWUE/dVPD", rot = 90, gp=gpar(fontsize=14)), bottom = textGrob("Year", gp=gpar(fontsize=14)))


ggsave2(paste("net_sens_VPD2","_v", model_version,".png",sep=""), plot = p, path = path_out, width = 10, height = 5)


site_list <- c("US-Seg", "US-Ses", "US-Wjs", "US-Mpj", "US-Vcp", "US-Vcm2", "US-Vcs")
p1 <- list()
for(i in c(1:length(site_list))){
  df_p1 <- df_p %>%
    filter(ID2 == "VPD") %>%
    filter(year > 2016) %>%
    filter(year < 2018) %>%
    filter(site %in% c(site_list[i]))
  
  maxy <- ifelse(site_list[i] == "US-Ses", 6, 2)
  
  p1[[i]] <- df_p1 %>%
    ggplot(aes(x=date)) +
    geom_pointrange(aes(y = mean_dYdX, ymax = pc97.5_dYdX, ymin = pc2.5_dYdX, color = ifelse(pc2.5_dYdX <= 0 & pc97.5_dYdX >= 0, "nonsignificant", "significant")), position = position_dodge(width = 1), fatten = .5, alpha=1) +
    #geom_point(aes(y=mean_dYdX_opt, color = "optimal net sensitivity"), alpha=1, size=.5) +
    #geom_point(aes(y=mean_dYdX_subopt, color = "suboptimal net sensitivity"), alpha=1, size=.5) +
    geom_hline(aes(yintercept = .5, color = "optimal net sensitivity"), linetype = 2) +
    geom_hline(aes(yintercept = 1, color = "suboptimal net sensitivity"), linetype = 2) +
    scale_color_manual(values = c("black", "gray","red", "blue"), breaks = c("significant", "nonsignificant", "optimal net sensitivity", "suboptimal net sensitivity")) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_rect(aes(xmin = date, xmax = dplyr::lead(date), ymin = -Inf, ymax = Inf, fill = Season), alpha = 0.1) +
    scale_fill_brewer(palette = "Dark2") +
    facet_col("site", scales = "free_y", strip.position = "right") +
    scale_x_bd(business.dates=df_p1$date, max.major.breaks=3, labels = date_format("%b-%Y")) + # labels=date_format("%Y") # labels = function(x) { x <- as.character(x); x[!1:0] <- " "; x}
    labs(title = NULL, y = NULL, x = NULL)+
    ylim(-1,maxy) +
    theme_bw() +
    theme(legend.position = "right",
          legend.text=element_text(size=10),
          #axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          axis.title.x = element_blank(),
          text = element_text(size=14),
          legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5))
  
}


layout <- '
AB
CD
EG
'
p <- wrap_plots(A = p1[[1]], B = p1[[2]], C = p1[[3]],
                D = p1[[4]], E = p1[[5]], G = p1[[7]],
                design = layout)
p <- p + plot_layout(guides = "collect") & theme(legend.position = 'right')
p <- grid.arrange(patchworkGrob(p), left = textGrob("dWUE/dVPD", rot = 90, gp=gpar(fontsize=14)), bottom = textGrob("Year", gp=gpar(fontsize=14)))


ggsave2(paste("net_sens_log_VPD1","_v", model_version,".png",sep=""), plot = p, path = path_out, width = 7, height = 5)



site_list <- c("US-Ses", "US-Wjs", "US-Vcp")
p1 <- list()
for(i in c(1:3)){
  df_p1 <- df_p %>%
    filter(ID2 == "VPD") %>%
    filter(year > 2016) %>%
    filter(site %in% c(site_list[i]))
  
  p1[[i]] <- df_p1 %>%
    ggplot(aes(x=date)) +
    geom_pointrange(aes(y = mean_dYdX, ymax = pc97.5_dYdX, ymin = pc2.5_dYdX, color = ifelse(pc2.5_dYdX <= 0 & pc97.5_dYdX >= 0, "nonsignificant", "significant")), position = position_dodge(width = 1), fatten = .5, alpha=1) +
    #geom_point(aes(y=mean_dYdX_opt, color = "optimal net sensitivity"), alpha=1, size=.5) +
    #geom_point(aes(y=mean_dYdX_subopt, color = "suboptimal net sensitivity"), alpha=1, size=.5) +
    geom_hline(aes(yintercept = .5, color = "optimal net sensitivity"), linetype = 2) +
    geom_hline(aes(yintercept = 1, color = "suboptimal net sensitivity"), linetype = 2) +
    scale_color_manual(values = c("black", "gray","red", "blue"), breaks = c("significant", "nonsignificant", "optimal net sensitivity", "suboptimal net sensitivity")) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_rect(aes(xmin = date, xmax = dplyr::lead(date), ymin = -Inf, ymax = Inf, fill = Season), alpha = 0.1) +
    scale_fill_brewer(palette = "Dark2") +
    facet_col("site", scales = "free_y", strip.position = "right") +
    scale_x_bd(business.dates=df_p1$date, max.major.breaks=7, labels = date_format("%b-%Y")) + # labels=date_format("%Y") # labels = function(x) { x <- as.character(x); x[!1:0] <- " "; x}
    labs(title = NULL, y = NULL, x = NULL)+
    ylim(-1,10) +
    theme_bw() +
    theme(legend.position = "top",
          legend.text=element_text(size=10),
          #axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          axis.title.x = element_blank(),
          text = element_text(size=14),
          legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5))
  
}


layout <- '
A
B
C
'
p <- wrap_plots(A = p1[[1]], B = p1[[2]], C = p1[[3]], design = layout)
p <- p + plot_layout(guides = "collect") & theme(legend.position = 'top')
p <- grid.arrange(patchworkGrob(p), left = textGrob("dWUE/dVPD", rot = 90, gp=gpar(fontsize=14)), bottom = textGrob("Year", gp=gpar(fontsize=14)))


ggsave2(paste("net_sens_log_VPD2","_v", model_version,".png",sep=""), plot = p, path = path_out, width = 10, height = 5)


df_p


site_list <- c("US-Ses", "US-Wjs", "US-Vcp")
p1 <- list()
for(i in c(1:3)){
  df_p1 <- df_p %>%
    filter(ID2 == "VPD") %>%
    filter(site %in% site_list[i])
  
  df_p1$dist <- abs(df_p1$mean_dYdX-.5)
  library(splines)
  library(ggdist)
  
  p1[[i]] <- df_p1 %>%
    ggplot(aes(x=Season, y = dist)) +
    PupillometryR::geom_flat_violin(position = position_nudge(x = .25, y = 0),adjust =2, trim = TRUE)+
    geom_point(aes(color = VPD), position = position_jitter(width = .15), size = .25)+
    geom_boxplot(outlier.shape = NA, alpha = 0.3, width = .1) +
    #geom_point(aes(color = ifelse(pc2.5_dYdX <= 0 & pc97.5_dYdX >= 0, "nonsignificant", "significant"))) +
    #geom_smooth(aes(y = dist)) +
    #stat_smooth(method = 'lm', formula = y ~ ns(x, 10)) +
    #geom_hline(aes(yintercept = .5, color = "optimal net sensitivity"), linetype = 2) +
    #geom_hline(aes(yintercept = 1, color = "suboptimal net sensitivity"), linetype = 2) +
    #scale_color_manual(values = c("black", "gray"), breaks = c("significant", "nonsignificant")) +
    scale_colour_viridis_c(option = "plasma") +
    geom_hline(yintercept = 0, linetype = 2) +
    facet_col("site", scales = "free_y", strip.position = "right") +
    labs(title = NULL, y = NULL, x = NULL)+
    coord_flip() +
    theme_bw() +
    theme(legend.position = "top",
          legend.text=element_text(size=10),
          #axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          axis.title.x = element_blank(),
          text = element_text(size=14),
          legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5))
  
}


layout <- '
A
B
C
'
p <- wrap_plots(A = p1[[1]], B = p1[[2]], C = p1[[3]], design = layout)
p <- p + plot_layout(guides = "collect") & theme(legend.position = 'top')
p <- grid.arrange(patchworkGrob(p), left = textGrob("dWUE/dVPD", rot = 90, gp=gpar(fontsize=14)), bottom = textGrob("Year", gp=gpar(fontsize=14)))


df_p1 %>%
    ggplot(aes(x=site, y = dist)) +
    PupillometryR::geom_flat_violin(position = position_nudge(x = .25, y = 0),adjust =2, trim = TRUE)+
    geom_point(aes(color = S), position = position_jitter(width = .15), size = .25)+
    geom_boxplot(outlier.shape = NA, alpha = 0.3, width = .1) +
    #geom_point(aes(color = ifelse(pc2.5_dYdX <= 0 & pc97.5_dYdX >= 0, "nonsignificant", "significant"))) +
    #geom_smooth(aes(y = dist)) +
    #stat_smooth(method = 'lm', formula = y ~ ns(x, 10)) +
    #geom_hline(aes(yintercept = .5, color = "optimal net sensitivity"), linetype = 2) +
    #geom_hline(aes(yintercept = 1, color = "suboptimal net sensitivity"), linetype = 2) +
    #scale_color_manual(values = c("black", "gray"), breaks = c("significant", "nonsignificant")) +
    scale_colour_viridis_c(option = "plasma") +
    #scale_colour_brewer(palette = "Dark2") +
    geom_hline(yintercept = 0, linetype = 2) +
    #facet_col("site", scales = "free_y", strip.position = "right") +
    labs(title = NULL, y = NULL, x = NULL)+
    coord_flip() +
    theme_bw() +
    theme(legend.position = "top",
          legend.text=element_text(size=10),
          #axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          axis.title.x = element_blank(),
          text = element_text(size=14),
          legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5))

df_p1 %>%
  #filter(site == "US-Seg") %>%
  ggplot(aes(x=Tair, y = dist)) +
  geom_point(aes(color = S), position = position_jitter(width = .15), size = .25)+
  geom_smooth() +
  scale_colour_viridis_c(option = "plasma") + 
  geom_hline(yintercept = 0, linetype = 2) +
  facet_grid(site~Season, scales = "free") +
  #labs(title = NULL, y = NULL, x = NULL)+
  #coord_flip() +
  theme_bw() +
  theme(legend.position = "top",
        text = element_text(size=14),
        plot.title = element_text(hjust = 0.5))


df_p1 <- df_p %>%
  filter(ID2 == "VPD")

dist <- abs(df_p1$mean_dYdX-.5)

df_p1 <- df_p %>%
  filter(ID2 == "Sshall")
df_p1$dist <- dist

p <- df_p1 %>%
  ggplot(aes(x=mean_dYdX, y = dist)) +
  geom_point(aes(color = S), size = .25)+
  scale_colour_viridis_c(option = "plasma") + 
  #scale_colour_viridis_d(option = "plasma") +
  geom_hline(yintercept = 0, linetype = 2) +
  labs(color = "shallow SWC") +
  xlab("Net sensitivity of WUE to shallow SWC") + ylab("Distance from optimal sensitivity to VPD") +
  theme_bw() +
  theme(legend.position = "top",
        text = element_text(size=14),
        plot.title = element_text(hjust = 0.5))
p
ggsave2(paste("optimaldist_vs_netsensSWC","_v", model_version,".png",sep=""), plot = p, path = path_out, width = 5, height = 5)


p <- df_p1 %>%
  ggplot(aes(x=mean_dYdX, y = dist)) +
  geom_point(aes(color = S), size = .25)+
  scale_colour_viridis_c(option = "plasma") + 
  #scale_colour_viridis_d(option = "plasma") +
  geom_hline(yintercept = 0, linetype = 2) +
  facet_col("site", scales = "free", strip.position = "right") +
  labs(color = "shallow SWC") +
  xlab("Net sensitivity of WUE to shallow SWC") + ylab("Distance from optimal sensitivity to VPD") +
  theme_bw() +
  theme(legend.position = "top",
        text = element_text(size=14),
        plot.title = element_text(hjust = 0.5))
p
ggsave2(paste("optimaldist_vs_netsensSWC_facet","_v", model_version,".png",sep=""), plot = p, path = path_out, width = 5, height = 9)

  

df_p1 <- df_p %>%
  filter(ID2 == "VPD")
dist <- abs(df_p1$mean_dYdX-.5)
df_p1$dist <- dist
df_p1 <- df_p1 %>%
  select(site,ID1,dist)
df_p1 <- left_join(df_p, df_p1)

library(ggh4x)
p <- df_p1 %>%
  ggplot(aes(x=mean_dYdX, y = dist)) +
  geom_point(aes(color = Season), size = .25, alpha = .5)+
  #scale_colour_viridis_c(option = "plasma") + 
  scale_colour_viridis_d(option = "plasma") +
  geom_hline(yintercept = 0, linetype = 2) +
  facet_nested(site~ID2+Season, scales = "free") +
  xlab("Net sensitivity of log WUE") + ylab("Distance from optimal sensitivity to VPD") +
  theme_bw() +
  theme(legend.position = "top",
        text = element_text(size=14),
        plot.title = element_text(hjust = 0.5))
p
ggsave2(paste("optimaldist_vs_netsens2","_v", model_version,".png",sep=""), plot = p, path = path_out, width = 15, height = 8)



######## end log scale graphs

site_list <- c("US-Seg","US-Ses", "US-Wjs","US-Mpj","US-Vcp","US-Vcm2","US-Vcs")
p1 <- list()
for(i in c(1:length(site_list))){
df_p1 <- df_p %>%
  filter(ID2 == "VPD") %>%
  filter(year > 2013) %>%
  filter(pc2.5_dYdX > 0) %>% # filter out nonsignificant values
  filter(site %in% c(site_list[i]))

p1[[i]] <- df_p1 %>%
  ggplot(aes(x=mean_dYdX_opt, y = mean_dYdX, group = Season, color = Season)) +
  geom_pointrange(aes( ymax = pc97.5_dYdX, ymin = pc2.5_dYdX, color = Season), fatten = .5, alpha=0.2) +
  scale_color_brewer(palette = "Dark2") +
  geom_smooth(method = lm, size = 1, color = "black", se = FALSE, alpha=0.5) + 
  geom_smooth(method="lm", se = F, size = 0.5, alpha=0.5) +
  geom_abline(slope=1, intercept=0, lty=2, col="black", size=1.25)+
  facet_col("site", scales = "free_y", strip.position = "top") +
  labs(title = NULL, y = NULL, x = NULL)+
  ylim(0,10) + xlim(0,10) +
  theme_bw() +
  theme(legend.position = c(.8,.2),
        legend.text=element_text(size=10),
        text = element_text(size=14),
        aspect.ratio=1,
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5))
}

layout <- '
ABCD
EFG#
'
p <- wrap_plots(A = p1[[1]], B = p1[[2]], C = p1[[3]],D = p1[[4]],E = p1[[5]],F = p1[[6]],G = p1[[7]], design = layout)
p <- p + plot_layout(guides = "collect") #& theme(legend.position = c(0.8,0.2))
p <- grid.arrange(patchworkGrob(p), left = textGrob("dWUE/dVPD", rot = 90, gp=gpar(fontsize=14)), bottom = textGrob("optimal response", gp=gpar(fontsize=14)))


ggsave2(paste("net_sens_VPD_opt_comp","_v", model_version,".png",sep=""), plot = p, path = path_out, width = 10, height = 5)


# just looking
df_p %>%
  filter(ID2 == "VPD") %>%
  filter(site == "US-Wjs") %>%
  filter(year>2016) %>%
  ggplot(aes(x=date, y=mean_dYdX)) +
  geom_point(aes(y=mean_dYdX), alpha=0.5) +
  #geom_pointrange(aes(ymax = pc97.5_dYdX, ymin = pc2.5_dYdX, color = ifelse(pc2.5_dYdX <= 0 & pc97.5_dYdX >= 0, "nonsignificant", "significant")), position = position_dodge(width = 1), fatten = .5, alpha=0.9) +
  geom_point(aes(y=mean_dYdX_opt, color = "optimal net sensitivity"), alpha=0.5, size=.3) +
  #geom_point(aes(y=mean_dYdX_subopt, color = "suboptimal net sensitivity"), alpha=0.5, size=.5) +
  scale_color_manual(values = c("black", "gray","red", "blue"), breaks = c("significant", "nonsignificant", "optimal net sensitivity", "suboptimal net sensitivity")) +
  geom_hline(yintercept = 0, linetype = 2) +
  ylim(-1,10) +
  facet_col("site") +
  theme_bw() +
  theme(legend.position = "none",
        legend.text=element_text(size=14),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        text = element_text(size=14),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5))



############## WUE uncertainty

date_temp <- as.data.frame(date_OUT[which(date_OUT$site !="US-Vcm1"),])
date_temp <- date_temp %>%
  rename(ID1 = dayind)

df_p <- sum_IN %>%
  filter(modelv==model_version) %>%
  filter(var == "WUE.pred") %>%
  filter(site!="US-Vcm1") %>%
  pivot_wider(id_cols = c(site,ID1,ID2), names_from = var, values_from = c(mean,median,sd,pc2.5,pc97.5))

df_p <- left_join(df_p, date_temp, by = c("ID1", "site"))
df_p$date <- as.Date(df_p$date)
df_p$site <- factor(df_p$site, levels = c("US-Seg", "US-Ses", "US-Wjs", "US-Mpj", "US-Vcp", "US-Vcm2", "US-Vcs"))

df_p <- df_p %>%
  filter(year>2015)

p <- df_p %>%
  ggplot(aes(x = date, y= mean_WUE.pred)) +
  geom_pointrange(aes(ymin=pc2.5_WUE.pred, ymax=pc97.5_WUE.pred), alpha=0.5, fatten =.2)+
  geom_rect(aes(xmin = date, xmax = dplyr::lead(date), ymin = -Inf, ymax = Inf, fill = Season), alpha = 0.1) +
  scale_fill_brewer(palette = "Dark2") +
  facet_col("site", strip.position = "right", scales = "free_y") +
  labs(title = NULL, x="date", y="predicted WUE") +
  theme_classic(base_size = 12)+
  theme(legend.position = "top",
        legend.text=element_text(size=14),
        text = element_text(size=14),
        legend.title = element_blank(),
        panel.background = element_rect(fill="white"),
        axis.line = element_line(color = "black"),
        axis.text.x = element_text(colour="black"),
        plot.title = element_text(hjust = 0.5))

p

ggsave2(paste("WUE_uncertainty","_v", model_version,".png",sep=""), plot = p, path = path_out)




p <- sum_IN %>%
  filter(modelv==5) %>%
  filter(site %in% c("US-Seg", "US-Ses", "US-Wjs", "US-Mpj", "US-Vcp", "US-Vcm2", "US-Vcs")) %>%
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



