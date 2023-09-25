library(tidyverse)
library(emmeans)
library(multcomp)
library(nlme)
library(MuMIn)
library(car)
library(corrplot)

# Set working directory to location where script is saved
setwd("C:/Users/dburruss/Documents/GitHub/xscale/")
#dir()

# Import the data
cross.import <- read_csv("Dylans_figures/data/OH_data_combined_block.csv") %>%
  mutate(treatment=factor(treatment,levels = c("control","plant-scale","patch-scale","plant- & patch-scale"))) 

# Graph of all the responses
ggplot(cross.import, aes(x = year, y = cover_pct, col = treatment)) +
  geom_point() +
  geom_line() +
  facet_grid(block ~ type) +
  theme(legend.position = "top") +
  scale_color_manual(values = c("#F8766D", "#7CAE00", "#00BFC4", "#C77CFF")) +
  scale_x_continuous(breaks = 2013:2021)

#Define responsone variable
respvar <- "Perennial Grasses - Alive" #Annual Grasses , Perennial Grasses - Alive , Forb_live

# Graph of Perennial Grass cover only
# The variability increases a lot beginning in 2017
cross.import %>% 
  dplyr::filter(type==respvar) %>% 
  mutate(blk_trt = paste0(block,"-",treatment)) %>%
  ggplot(aes(x = year, y = cover_pct, col = treatment, group = blk_trt)) +
  geom_point() +
  geom_line() +
  #facet_wrap(~ block, ncol = 3) +
  theme(legend.position = "top") +
  scale_color_manual(values = c("#F8766D", "#7CAE00", "#00BFC4", "#C77CFF")) +
  scale_x_continuous(breaks = 2013:2021)

# Adjust data set for modeling
cross.import.analysis <- cross.import %>%
  mutate(year=factor(year),
         block=factor(block),
         type=factor(type)) %>%
  mutate(blk_trt = paste0(block,"-",treatment),
         time = as.numeric(year)) %>%
  rename(Ann_grass_t0 = `Annual Grasses_t0`,
         Ann_grass_tm1 = `Annual Grasses_t-1`,
         P_grass_t0 = `Perennial Grasses - Alive_t0`,
         P_grass_tm1 = `Perennial Grasses - Alive_t-1`,
         Forb_t0 = `Forb_live_t0`,
         Forb_tm1 = `Forb_live_t-1`,
         Litter_tm1 = `Litter_t-1`)

#ANOVA for Treatment

pg.gls <- gls(cover_pct ~ treatment, 
                            correlation = corAR1(form = ~ time | blk_trt),
                            data=cross.import.analysis %>% dplyr::filter(type==respvar))
anova(pg.gls) # Block is not significant
qqnorm(pg.gls, abline = c(0,1)) 


anova(lme(cover_pct ~ block, 
          random=~1 | treatment, 
          method="ML", 
          data=cross.import.analysis %>% dplyr::filter(type==respvar)))

# # Perennial Grass Analysis with Fixed blocks
# pg.gls <- gls(cover_pct ~ block + treatment + year + treatment:year, 
#               correlation = corAR1(form = ~ time | blk_trt),
#               data=cross.import.analysis %>% dplyr::filter(type==respvar))
# 
# anova(pg.gls) # Block is not significant
# qqnorm(pg.gls, abline = c(0,1)) # Residuals are not great





# First look at year*treatment LSMeans
trt.year.lsmeans <- emmeans(pg.gls, ~ treatment+year, lmer.df = "kenward-roger") %>%
  multcomp::cld( Letters = LETTERS, decreasing = FALSE) %>%
  as_tibble() %>%
  arrange(desc(emmean)) %>%
  mutate(.group = trimws(.group))

pd <- position_dodge(0.5)
ggplot(trt.year.lsmeans, aes(x = year, y = emmean, group = treatment, col = treatment)) +  
  theme_bw() +
  geom_point(position=pd) +
  geom_line(position=pd) +
  geom_errorbar(aes(ymin = emmean - SE , ymax = emmean + SE), width=0.2, position=pd) +
  geom_text(aes(label = .group, y = emmean + SE), 
            position = pd, vjust = -0.5, col = "black", size = 2) +
  theme(legend.position = "top") +
  ggtitle(paste0(respvar," Cover: Treatment * Year Least Squares Means"),
          subtitle = "Means with the same letter are not different at alpha = 0.05 (Tukey method)") +
  ylab(paste0("Model - Based ", respvar, " cover S.E."))
#ggsave(paste0(respvar, " Trt x Year LSMeans.jpg"), height = 6, width = 10, units = "in", dpi = 300)

# Treatment main effect
trt.lsmeans <- emmeans(pg.gls, ~ treatment, lmer.df = "kenward-roger") %>%
  multcomp::cld(Letters = LETTERS, decreasing = FALSE) %>%
  as_tibble() %>%
  arrange(desc(emmean)) %>%
  mutate(.group = trimws(.group))

ggplot(trt.lsmeans, aes(x = treatment, y = emmean, fill = treatment)) + 
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = emmean - SE , ymax = emmean + SE), width=0.2) +
  geom_text(aes(label = .group, y = emmean + SE), 
            vjust = -0.5, col = "black", size = 4) 


################################################################################
# # See if model that adjusts for control fits better (y changes to cover_pct_d_yearly)
# pg.gls.adj <- gls(cover_pct_d_yearly ~ block + treatment + year + treatment:year, 
#               correlation = corAR1(form = ~ time | blk_trt),
#               data=cross.import.analysis %>% dplyr::filter(type==respvar & treatment != "control"))
# 
# anova(pg.gls.adj) # Block is not significant
# qqnorm(pg.gls.adj, abline = c(0,1)) # SLightly better bust still not great
# 
# # Because of singular fit, lmer.df = "kenward-roger" does not work; Use approximate Satterhwaite instead
# adj.trt.year.lsmeans <- emmeans(pg.gls.adj, ~ treatment*year, mode = "appx-satterthwaite") %>%
#   multcomp::cld(Letters = LETTERS, decreasing = FALSE) %>%
#   as_tibble() %>%
#   arrange(desc(emmean)) %>%
#   mutate(.group = trimws(.group))
# 
# pd <- position_dodge(0.5)
# ggplot(adj.trt.year.lsmeans, aes(x = year, y = emmean, group = treatment, col = treatment)) +  
#   theme_bw() +
#   geom_point(position=pd) +
#   geom_line(position=pd) +
#   geom_errorbar(aes(ymin = emmean - SE , ymax = emmean + SE), width=0.2, position=pd) +
#   geom_text(aes(label = .group, y = emmean + SE), 
#             position = pd, vjust = -0.5, col = "black", size = 4) +
#   theme(legend.position = "top") +
#   ggtitle("Control-Adjusted Perennial Grass Cover: Treatment * Year Least Squares Means",
#           subtitle = "Means with the same letter are not different at alpha = 0.05 (Tukey method)") +
#   ylab("Model - Based Perennial Grass cover ? S.E.")
# 
# # Same basic results as model with Control

# use gls() to run a repeated measures analysis

################################################################################
#test for correlation among predictors
res <- cor(cross.import.analysis[, c(
  "ppt_cm_gs",
  "ppt_cm_wy",
  "Bare_2011",
  "Shrub_2011",
  "accum_bsne",
  "rel_bsne",
  "Ann_grass_t0",
  #"Ann_grass_tm1",
  "Forb_t0",
  #"Forb_tm1",
  "P_grass_t0",
  #"Perennial Grasses - Alive_t-1",
  "Litter_t0"#,"Litter_tm1"
)])

corrplot(res,
         is.corr = FALSE, method = "square")

#no apparent correlation issues

#test for colinarity among predictor variables

colmod <- lm(cover_pct ~ 
               treatment +
               ppt_6mo_prior +
               ppt_cm_gs +
               ppt_cm_wy +
               Bare_2011 +
               Shrub_2011 +
               #accum_bsne +
               rel_bsne +
               #Ann_grass_t0 +
               #Ann_grass_tm1 +
               #Forb_t0 +
               #Forb_tm1 +
               #P_grass_t0 +
               #P_grass_tm1 +
               Litter_t0 +
               Litter_tm1,
             data=cross.import.analysis)

vif(colmod) #values > 5 suggest there is significant collinearity between vars.
#collinearity between vars does not appear to be a problem

#create 2016
dat2016 <- cross.import.analysis %>% 
  filter(year=="2016",
         type=="Perennial Grasses - Alive") %>%
  mutate(forb2016 = Forb_t0) %>%
  select(block, treatment, forb2016)

#estimate cover for significant groups
combined_trt <- cross.import.analysis %>%
  left_join(dat2016) %>%
  filter(treatment =="control", #plant- & patch-scale
         type == "Perennial Grasses - Alive",
         time > 3)   #year greater than 2015 [year 3]


Mod1 <- gls(cover_pct ~ 
             #treatment +  
             Shrub_2011:ppt_6mo_prior +
             #Litter_t0:ppt_cm_wy +
             #ppt_cm_gs +
             #ppt_cm_wy +
             ppt_6mo_prior +
             # Bare_2011 +
             Shrub_2011 +
             # accum_bsne +
             # rel_bsne +
             # Ann_grass_t0 +
             # Forb_t0 +
             # Litter_t0 ,
             Ann_grass_tm1 +
             Forb_tm1 ,
             # Litter_tm1 ,
             #P_grass_tm1,
             
           correlation = corAR1(form = ~ time | block), #blk_trt
           method = "ML",
           data=combined_trt)

#model summary
summary(Mod1)



#create a scatterplot the observed versus fitted line for all values
plot(Mod1, which=1)

#Create a scatterplot of the observed versus fitted line by block
plot(Mod1, fitted(.) ~ cover_pct | as.factor(round(Shrub_2011,digits=2)), 
     abline = c(0,1), 
     xlim=c(0,60), ylim=c(0,60),
     text(labels=))

d0 <- dredge(Mod1,
             subset = !(rel_bsne && accum_bsne) && !(ppt_6mo_prior && ppt_cm_gs) && !(Bare_2011 && Shrub2011),
             extra = c("R^2","adjR^2"),
             trace = 2)

# use model averaging to combine multiple models
# summary(model.avg(d0, subset=delta < 4))

# clean up dredge results
dredge_out <- d0 %>%
  as_tibble()


# Treatment main effect using only significant variables in Mod1
model.lsmeans <- emmeans(Mod1, specs = c("ppt_cm_gs",
                                           "ppt_cm_wy",
                                           "Bare_2011",
                                           "Shrub_2011",
                                           "accum_bsne",
                                           "rel_bsne",
                                           "Ann_grass_t0",
                                           "Forb_t0",
                                           "Litter_t0")) %>%
  as_tibble() %>%
  pivot_longer(cols = -c(emmean, SE, df, lower.CL, upper.CL),
               names_to = "vars",
               values_to = "values")

ggplot(model.lsmeans, aes(x = vars, y = values, fill = vars)) + 
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = values - SE , ymax = values + SE), width=0.2) 



