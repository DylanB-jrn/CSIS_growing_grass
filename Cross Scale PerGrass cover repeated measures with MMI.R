
library(tidyverse)
library(emmeans)
library(multcomp)
library(nlme)

# Set working directory to location where script is saved
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
dir()

# Import the data
cross.import <- read_csv("OH_data_combined_block.csv") %>%
  mutate(treatment=factor(treatment,levels = c("control","plant-scale","patch-scale","plant- & patch-scale"))) 

# Graph of all the responses
ggplot(cross.import, aes(x = year, y = cover_pct, col = treatment)) +
  geom_point() +
  geom_line() +
  facet_grid(block ~ type) +
  theme(legend.position = "top") +
  scale_color_manual(values = c("#F8766D", "#7CAE00", "#00BFC4", "#C77CFF")) +
  scale_x_continuous(breaks = 2013:2021)

# Graph of Perennial Grass cover only
# The variability increases a lot beginning in 2017
cross.import %>% 
  dplyr::filter(type=="Perennial Grasses - Alive") %>%
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
         time = as.numeric(year))

# Perennial Grass Analysis with Fixed blocks
pg.gls <- gls(cover_pct ~ block + treatment + year + treatment:year, 
              correlation = corAR1(form = ~ time | blk_trt),
              data=cross.import.analysis %>% dplyr::filter(type=="Perennial Grasses - Alive"))

anova(pg.gls) # Block is not significant
qqnorm(pg.gls, abline = c(0,1)) 

# First look at year*treatment LSMeans
trt.year.lsmeans <- emmeans(pg.gls, ~ treatment*year, lmer.df = "kenward-roger") %>%
  multcomp::cld(Letters = LETTERS, decreasing = FALSE) %>%
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
  ggtitle("Perennial Grass Cover: Treatment * Year Least Squares Means",
          subtitle = "Means with the same letter are not different at alpha = 0.05 (Tukey method)") +
  ylab("Model - Based Perennial Grass cover ? S.E.")
ggsave("Pgrass Trt x Year LSMeans.jpg", height = 6, width = 10, units = "in", dpi = 300)


# Model different treatments separately

library(AICcmodavg) # For calculating AICc
library(drc) # For fitting standard 3-parameter logistic function
library(segmented) # For fitting segmented regression


# Use "plant- & patch-scale" as first example
# Cover needs to be between 0 and 1 for logistic model to work
# Year needs to be numeric
grass.mmi.df <- trt.year.lsmeans %>%
  dplyr::filter(treatment=="plant- & patch-scale") %>%
  mutate(cover = emmean/100,
         year = as.character(year) %>% as.numeric())

ggplot(grass.mmi.df, aes(x = year, y = cover)) +
  geom_point() +
  scale_x_continuous(breaks = unique(grass.mmi.df$year)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  ggtitle("plant- & patch-scale model-based estimates")

# Null model: slope of 0
null.lm <- lm(cover ~ 1, data = grass.mmi.df)

# Single line
single.line.lm <- lm(cover ~ year, data = grass.mmi.df)

# logistic
# Use drm() with fct = L.3() argument
logistic.lm <- drm(cover ~ year, data = grass.mmi.df, fct = L.3())
AIC(logistic.lm)
AICc(logistic.lm) # Does not work
# One measure of R-squared: squared correlation of fitted values to actual values
cor(fitted(logistic.lm), grass.mmi.df$cover)^2

grass.mmi.df %>%
  mutate(logistic = fitted(logistic.lm)) %>%
  ggplot() +
  geom_point(aes(x = year, y = cover)) +
  geom_line(aes(x = year, y = logistic))

# Look for breakpoint
segmented.lm <- segmented(single.line.lm, npsi = 1) 
summary(segmented.lm) # breakpoint is 2015.226


null.stats <- data.frame(variable = "cover_pct", model = "null",
                            AIC = AIC(null.lm),
                            AICc = AICc(null.lm),
                            rSquared = summary(null.lm)$r.squared,
                            Adj_rSquared = summary(null.lm)$adj.r.squared,
                            stringsAsFactors = FALSE)

single.line.stats <- data.frame(variable = "cover_pct", model = "single line",
                                 AIC = AIC(single.line.lm),
                                 AICc = AICc(single.line.lm),
                                 rSquared = summary(single.line.lm)$r.squared,
                                 Adj_rSquared = summary(single.line.lm)$adj.r.squared,
                                 stringsAsFactors = FALSE)

logistic.stats <- data.frame(variable = "cover_pct", model = "logistic",
                              AIC = AIC(logistic.lm),
                              AICc = as.numeric(NA),
                              rSquared = as.numeric(NA),
                              Adj_rSquared = cor(fitted(logistic.lm), grass.mmi.df$cover)^2,
                              stringsAsFactors = FALSE)

segmented.stats <- data.frame(variable = "cover_pct", model = "split line",
                              breakpoint = summary(segmented.lm)$psi[2],
                              AIC = AIC(segmented.lm),
                              AICc = AICc(segmented.lm),
                              rSquared = summary(segmented.lm)$r.squared,
                              Adj_rSquared = summary(segmented.lm)$adj.r.squared,
                              stringsAsFactors = FALSE)

model.results <- null.stats %>%
  bind_rows(single.line.stats) %>%
  bind_rows(logistic.stats) %>%
  bind_rows(segmented.stats) %>%
  mutate(aic_min = min(AIC),
         delta_AIC =  AIC - min(AIC)) %>%
  arrange(delta_AIC)

model.results

write_csv(model.results, "perennial grass MMI on repeated measures LSMeans.csv")
           
# Graph of models
model.mmi.fig <- grass.mmi.df %>%
  dplyr::select(treatment, year, cover) %>%
  mutate(null = fitted(null.lm),
         single_line = fitted(single.line.lm),
         logistic = fitted(logistic.lm),
         segmented = fitted(segmented.lm)) %>%
  gather(model, prediction, null, single_line, logistic, segmented) %>%
  ggplot(aes(x= year, y = prediction, col = model)) +
  geom_point(aes(x = year, y = cover)) +
  geom_line() +
  scale_x_continuous(breaks = unique(grass.mmi.df$year)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  ggtitle("plant- & patch-scale model-based estimates", subtitle = "using repeated measures model LSMeans") +
  ylab("Model - Based Perennial Grass cover")
model.mmi.fig


# Repeat with raw data, ignoring the blocks----

data.model <- cross.import %>%
  dplyr::filter(type=="Perennial Grasses - Alive" & treatment=="plant- & patch-scale" ) %>%
  mutate(cover = cover_pct / 100)

# null model
null.lm.dat <- lm(cover ~ 1, data = data.model)

# Single line
single.line.lm.dat <- lm(cover ~ year, data = data.model)

# logistic
# Use drm() with fct = L.3() argument
logistic.lm.dat <- drm(cover ~ year, data = data.model, fct = L.3())

# Look for breakpoint
segmented.lm.dat <- segmented(single.line.lm.dat, npsi = 1) 
summary(segmented.lm.dat)$psi[2] # breakpoint is 2015.226


null.stats.dat <- data.frame(variable = "cover_pct", model = "null",
                         AIC = AIC(null.lm.dat),
                         AICc = AICc(null.lm.dat),
                         rSquared = summary(null.lm.dat)$r.squared,
                         Adj_rSquared = summary(null.lm.dat)$adj.r.squared,
                         stringsAsFactors = FALSE)

single.line.stats.dat <- data.frame(variable = "cover_pct", model = "single line",
                                AIC = AIC(single.line.lm.dat),
                                AICc = AICc(single.line.lm.dat),
                                rSquared = summary(single.line.lm.dat)$r.squared,
                                Adj_rSquared = summary(single.line.lm.dat)$adj.r.squared,
                                stringsAsFactors = FALSE)

logistic.stats.dat <- data.frame(variable = "cover_pct", model = "logistic",
                             AIC = AIC(logistic.lm.dat),
                             AICc = as.numeric(NA),
                             rSquared = as.numeric(NA),
                             Adj_rSquared = cor(fitted(logistic.lm.dat), data.model$cover)^2,
                             stringsAsFactors = FALSE)

segmented.stats.dat <- data.frame(variable = "cover_pct", model = "split line",
                              breakpoint = summary(segmented.lm.dat)$psi[2],
                              AIC = AIC(segmented.lm.dat),
                              AICc = AICc(segmented.lm.dat),
                              rSquared = summary(segmented.lm.dat)$r.squared,
                              Adj_rSquared = summary(segmented.lm.dat)$adj.r.squared,
                              stringsAsFactors = FALSE)

dat.model.results <- null.stats.dat %>%
  bind_rows(single.line.stats.dat) %>%
  bind_rows(logistic.stats.dat) %>%
  bind_rows(segmented.stats.dat) %>%
  mutate(aic_min = min(AIC),
         delta_AIC =  AIC - min(AIC)) %>%
  arrange(delta_AIC)

dat.model.results
write_csv(dat.model.results, "perennial grass MMI on raw data.csv")

# Graph of models
data.mmi.fig <- data.model %>%
  dplyr::select(treatment, year, cover) %>%
  mutate(null = fitted(null.lm.dat),
         single_line = fitted(single.line.lm.dat),
         logistic = fitted(logistic.lm.dat),
         segmented = fitted(segmented.lm.dat)) %>%
  gather(model, prediction, null, single_line, logistic, segmented) %>%
  ggplot(aes(x= year, y = prediction, col = model)) +
  geom_point(aes(x = year, y = cover)) +
  geom_line() +
  scale_x_continuous(breaks = unique(data.model$year)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  ggtitle("plant- & patch-scale model-based estimates", subtitle = "using raw data without accounting for blocks") +
  ylab("Perennial Grass cover (raw data)")

data.mmi.fig

library(patchwork)
model.mmi.fig / data.mmi.fig

ggsave("Regression estimates.jpg", height = 10, width = 10, units = "in", dpi = 300)

  


