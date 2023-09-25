library(tidyverse)
library(nlme)
library(MuMIn)


setwd("C:/Users/dburruss/Documents/GitHub/xscale/")

#read in the data
oh <- read.csv("Dylans_figures/data/OH_data_combined_block.csv") %>%
  as_tibble() %>%
  mutate(treatment=factor(treatment,levels = c("control","plant-scale","patch-scale","plant- & patch-scale"))) %>%
  filter(type=="Perennial Grasses - Alive") %>% #Perennial Grasses - Alive  Annual Grasses   Forb_live
  filter(year>2015) %>% #to accommodate year-1 variables
  
  mutate(year=factor(year),
         block=factor(block),
         type=factor(type)) %>%
  #dplyr::select(year, block, treatment, cover_pct_d_yearly, Litter_t0, Annual.Grasses_t0, Forb_live_t0,Litter_t.1, Annual.Grasses_t.1, Forb_live_t.1, ppt_cm_wy) %>%
  mutate(blk_trt = paste0(block,"-",treatment),
         time = as.numeric(year))

####Dredging to calculate model weight

###plant- & patch-scale

oh_comb <- oh %>%
  filter(treatment=="plant- & patch-scale")  #"control","plant-scale","patch-scale","plant- & patch-scale"

gls_comb <- gls(cover_pct ~ 
             Annual.Grasses_t.1 +
             pg_2017 +
             ppt_cm_wy + 
             Shrub_2011 +
             vert_acc_2017,

           correlation = corAR1(form = ~ time | blk_trt),
           method = "ML",
           data=oh_comb)

summary(gls_comb)

# calc model weight
d_comb <- dredge(gls_comb,
             extra = c("R^2","adjR^2"),
             trace = 2)

#write out
write.xlsx(d_comb, paste0("C:/Users/dburruss/Documents/GitHub/xscale/scripts/Dylan_working/OH_final_model_dredge_RIV.xls"),
           append = TRUE,
           sheetName = "Combined")

###patch-scale

oh_patch <- oh %>%
  filter(treatment=="patch-scale")  #"control","plant-scale","patch-scale","plant- & patch-scale"

gls_patch <- gls(cover_pct ~ 
                   Litter_t.1 +
                   Litter_t0 +
                   pg_2017 + 
                   Shrub_2011 ,
                
                correlation = corAR1(form = ~ time | blk_trt),
                method = "ML",
                data=oh_patch)

summary(gls_patch)

# calc model weight
d_patch <- dredge(gls_patch,
                 extra = c("R^2","adjR^2"),
                 trace = 2)

#write out
write.xlsx(d_patch, paste0("C:/Users/dburruss/Documents/GitHub/xscale/scripts/Dylan_working/OH_final_model_dredge_RIV.xls"),
           append = TRUE,
           sheetName = "Patch")

###plant-scale

oh_plant <- oh %>%
  filter(treatment=="plant-scale")  #"control","plant-scale","patch-scale","plant- & patch-scale"

gls_plant <- gls(cover_pct ~ 
                   Forb_live_t.1 +
                   Litter_t.1 +
                   pg_2017 + 
                   ppt_cm_wy ,
                 
                 correlation = corAR1(form = ~ time | blk_trt),
                 method = "ML",
                 data=oh_plant)

summary(gls_plant)

# calc model weight
d_plant <- dredge(gls_plant,
                  extra = c("R^2","adjR^2"),
                  trace = 2)

#write out
write.xlsx(d_plant, paste0("C:/Users/dburruss/Documents/GitHub/xscale/scripts/Dylan_working/OH_final_model_dredge_RIV.xls"),
           append = TRUE,
           sheetName = "Plant")

###control

oh_control <- oh %>%
  filter(treatment=="control")  #"control","plant-scale","patch-scale","plant- & patch-scale"

gls_control <- gls(cover_pct ~ 
                     ppt_cm_wy +
                     pg_2017 ,
                 
                 correlation = corAR1(form = ~ time | blk_trt),
                 method = "ML",
                 data=oh_control)

summary(gls_control)

# calc model weight
d_control <- dredge(gls_control,
                  extra = c("R^2","adjR^2"),
                  trace = 2)

#write out
write.xlsx(d_control, paste0("C:/Users/dburruss/Documents/GitHub/xscale/scripts/Dylan_working/OH_final_model_dredge_RIV.xls"),
           append = TRUE,
           sheetName = "Control")
