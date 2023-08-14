#### IMPORT PACKAGES #### 
library(ggpubr)
library(rstatix)
library(readr)
library(dplyr)
library(lme4)
library(lmerTest)
library(lmtest)
library(reghelper)
library(sjPlot)


#### READ RAW DATA ####
df <- readRDS(file = "...raw_data.rds") # TO ADAPT


#Drop all the rows where „firstResponse” has a value of 0.
#drop all the rows where „tripletType” has a value of „X”, „T”, or „R”
#drop all the rows where „isPractice”  has a value of 1
df<-subset(df, firstResponse!="0" )
df <- subset (df, tripletType!="X" & tripletType!= "T" & tripletType!= "R")
df <- subset (df, isPractice!= "")

df <- subset(df, trial_type== "serial-reaction-time")


# Calculate median reaction times for every subject (35 subjects), every session (3 sessions), every epoch (4 epochs) and for both both tripletTypes (high and low frequency triplets)
# This results in: 35 x 3 x 4 x 2 = 840 data rows:
rtmedian <- aggregate(x = df$rt,                # Specify data column
                      by = list(df$subject,df$session,df$epoch,df$tripletType,df$group,df$MoCA_BL,df$age),              # Specify group indicator
                      FUN = median)                           # Specify function (i.e.mean)
colnames(rtmedian) <- c("subject","session","epoch","tripletType","group","MoCA","age","rt")


#### GENERATE RESULTS DESCRIBED IN SECTION 3.1 ####
rtmedian$epoch = rtmedian$epoch -1 # start epoch numbering from 0 to 3 instead of 1 to 4 (important for correct lmm interpretation as epoch will be inserted as interaction term)
# generate LMM 
lmm_res_3.1 <- lmer(rt ~ group*session+tripletType+epoch*session+(1+epoch|subject),data = rtmedian)
summary(lmm_res_3.1)

#### GENERATE RESULTS DESCRIBED IN SECTION 3.2 ####
rtmedian$MoCA = rtmedian$MoCA - mean(rtmedian$MoCA) # mean-centered MoCA scores
rtmedian$age = rtmedian$age - mean(rtmedian$age)    # mean-centered age

# generate LMM 
lmm_res_3.2 <- lmer(rt ~ group*session*MoCA+tripletType+session*age*group+epoch*session*MoCA+(1+epoch|subject),data = rtmedian)   # actually better
summary(lmm_res_3.2)

#### GENERATE FIGURES 3 & 4 ####
# For visualization purposes we need to re-calculate median reaction times excluding epochs
rtmedian <- aggregate(x = df$rt,                # Specify data column
                      by = list(df$subject,df$session,df$group,df$MoCA_BL,df$age,df$condition),              # Specify group indicator
                      FUN = median)                           # Specify function (i.e.mean)
colnames(rtmedian) <- c("subject","session","group","MoCA_BL","age","Cond","rt")

# BASELINE CORRECTION 
rtmedian <- rtmedian %>%
  arrange(subject, session) %>%
  group_by(subject) %>%
  mutate(rt_Imp = rt - rt[1L]) %>%
  ungroup()

# Adding Age column
rtmedian <-  rtmedian %>% mutate(Age =
                                   case_when(age <=68 ~ "≤68",
                                             age > 68 ~ ">68"))

rtmedian$Age <- as.factor(rtmedian$Age)

levels(rtmedian$Cond) <- c("≥26", "<26")
levels(rtmedian$session) <- c("Baseline", "Post-treatment", "3 months later")
rtmedian$MoCA <- rtmedian$Cond

# FIGURE 3
# ggline
vml_line <- ggline(rtmedian, x = "session", y = "rt", add = c("mean_se"),
                   shape = "group",
                   color = "group", palette = c("orangered1", "steelblue"))+
  theme_classic()+
  theme(axis.title = element_text(size=15), axis.text = element_text(size=15), strip.text.x = element_text(size=15))+
  theme(plot.title = element_blank(), legend.text = element_text(size=15),legend.title = element_blank(),
        legend.position="bottom",axis.title.x = element_blank())+
  ylab("Reaction time in ms") 

# ggbar
rtmedian_BC_plot <- subset(rtmedian,session!="Baseline")
rtmedian_BC_plot$rt_Imp <- rtmedian_BC_plot$rt_Imp*(-1)

levels(rtmedian_BC_plot$session) <- c("Baseline", "Direct effect", "Long-term effect")

vml_offline <- ggbarplot(rtmedian_BC_plot, x = "session", y = "rt_Imp", add = c("mean_se", "jitter"),
                         shape = "group",
                         color = "group", palette = c("orangered1", "steelblue"), 
                         position = position_dodge(0.8))+
  theme_classic()+
  theme(axis.title = element_text(size=20), axis.text = element_text(size=15), strip.text.x = element_text(size=15))+
  theme(plot.title = element_text(hjust=0.5,size=20), legend.text = element_text(size=15),legend.title = element_blank(),
        legend.position="bottom",axis.title.x = element_blank())+
  ylab("Offline visuomotor learning in ms")

ggarrange(vml_line, vml_offline,
          labels = c("a", "b"),
          font.label = list(size = 20),
          ncol = 2, 
          nrow = 1,
          widths = c(2, 1))

# FIGURE 4
# SESSION X GROUP X AGE 
rtmedian$rt_Imp <- rtmedian$rt_Imp*(-1)
vml_age_line <- ggline(rtmedian, x = "session", y = "rt_Imp", add = "mean_se",shape = "Age",
                       color = "Age", palette = c("blueviolet", "cadetblue"),facet.by = "group")+
  theme_classic()+
  theme(axis.title = element_text(size=15), axis.text = element_text(size=15), strip.text.x = element_text(size=15))+
  theme(plot.title = element_text(size=20,hjust = 0.5),legend.text = element_text(size=15),legend.title = element_text(size=15),
        legend.position="bottom",axis.title.x = element_blank())+
  xlab("Session") + ylab("Offline visuomotor learning in ms") 

# SESSION X GROUP X MoCA 
vml_MoCA_line <- ggline(rtmedian, x = "session", y = "rt_Imp", add = "mean_se",shape = "MoCA",
                        color = "MoCA", palette = c("darkgreen", "orange"),facet.by = "group")+
  theme_classic()+
  theme(axis.title = element_text(size=15), axis.text = element_text(size=15), strip.text.x = element_text(size=15))+
  theme(plot.title = element_text(size=20,hjust = 0.5),legend.text = element_text(size=15),legend.title = element_text(size=15),
        legend.position="bottom",axis.title.x = element_blank())+
  xlab("Session") + ylab("Offline visuomotor learning in ms") 

ggarrange(vml_age_line,vml_MoCA_line,
          labels = c("a", "b"),
          font.label = list(size = 20),
          ncol = 1, 
          nrow = 2
)

#### GENERATE RESULTS DESCRIBED IN SUPPLEMENTARY MATERIALS SECTION S3.1 ####

# Subset baseline data only
df1 <- subset(df, session== 1)


# Calculate median reaction times for every subject (35 subjects), every epoch (4 epochs) and for both both tripletTypes (high and low frequency triplets)
# This results in: 35 x 4 x 2 = 280 data rows:
rtmedian <- aggregate(x = df1$rt,                # Specify data column
                      by = list(df1$subject,df1$MoCA_BL,df1$age,df1$epoch,df1$tripletType),              # Specify group indicator
                      FUN = median)                           # Specify function (i.e.mean)
colnames(rtmedian) <- c("subject","MoCA","age","epoch","tripletType","rt")

rtmedian$epoch = rtmedian$epoch -1 # start epoch numbering from 0 to 3 instead of 1 to 4 (important for correct lmm interpretation as epoch will be inserted as interaction term)
rtmedian$MoCA = rtmedian$MoCA - mean(rtmedian$MoCA) # mean-centered MoCA scores
rtmedian$age = rtmedian$age - mean(rtmedian$age)    # mean-centered age

# generate LMM 
lmm_res_S3.1 <- lmer(rt ~ epoch*MoCA+tripletType*age+(1+epoch|subject), data=rtmedian) 
summary(lmm_res_S3.1)

# generate figure S.3.1
# For visualization purposes we need to re-generate the lmm using non mean-centered MoCA and age values
rtmedian <- aggregate(x = df1$rt,                # Specify data column
                      by = list(df1$subject,df1$MoCA_BL,df1$age,df1$epoch,df1$tripletType),              # Specify group indicator
                      FUN = median)                           # Specify function (i.e.mean)
colnames(rtmedian) <- c("subject","MoCA","age","epoch","tripletType","rt")

# generate LMM 
lmm_res_S3.1_plot <- lmer(rt ~ epoch*MoCA+tripletType*age+(1+epoch|subject), data=rtmedian) 

figS3.1_a <- plot_model(lmm_res_S3.1_plot,type="pred",terms = c("age","tripletType"))+
  theme_classic()+
  theme(axis.title = element_text(size=20), axis.text = element_text(size=15), strip.text.x = element_text(size=15))+
  theme(plot.title = element_blank(), legend.text = element_text(size=15),legend.title = element_text(size=15),
        legend.position="bottom")+
  xlab("Age") + ylab("Reaction time in ms") + labs(color = "trial type")

figS3.1_b <- plot_model(lmm_res_S3.1_plot,type="pred",terms = c("epoch","tripletType","MoCA [22,28]"))+
  theme_classic()+
  theme(axis.title = element_text(size=20), axis.text = element_text(size=15), strip.text.x = element_text(size=15))+
  theme(plot.title = element_text(size=20,hjust = 0.5), legend.text = element_text(size=15),legend.title = element_text(size=15),
        legend.position="bottom")+
  xlab("Epoch") + ylab("Reaction time in ms") + labs(color = "trial type")

ggarrange(figS3.1_a+ rremove("legend"), figS3.1_b, 
          labels = c("a", "b"),
          font.label = list(size = 20),
          ncol = 2, nrow = 1)

