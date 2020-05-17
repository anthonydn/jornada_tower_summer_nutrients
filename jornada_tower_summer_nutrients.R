# R version 4.0.0 (2020-04-24)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 10 x64 (build 16299)
#
# Relevant packages:
# treateffect_0.3.1 -- devtools::install_github("anthonydn/treateffect@075602b")
# lubridate_1.7.8
# emmeans_1.4.6
# nlme_3.1-147
# dplyr_0.8.99.9003 (this is the dev version that should basically be 1.0)
# readr_1.3.1
# tidyr_1.0.3
# tibble_3.0.1
# ggplot2_3.3.0
# tidyverse_1.3.0

library(treateffect)
library(nlme)
library(emmeans)
library(lubridate)
theme_te()

jorn15 <- read_csv("jornada_tower_summer_nutrients.csv") %>%
  mutate(Type = ordered(Type), doy = `Day of Year`, plot_num = factor(plot_num),
    soil_moisture = (1 - percent_dry_mass) * 100)

jorn15_scale <- function() scale_x_continuous(
  breaks = c(155, 168, 175, 182, 189, 196),
  limits = c(150, 198),
  labels = c("Jun 4", "Jun 17", "Jun 24", "Jul 1", "Jul 8", "Jul 15"))

#average response variables over summer
jorn15s.te <-
  treateffect(jorn15, percentC + percentN + MBC + MBN + nh4 + no3 + po4 +
  BGLUC + NAG + CELLOBIO + XYLO + AGLUC + PHOS + Perox + PhenolOx +
  soil_moisture + respiration ~ Type,
  block = "plot_num", average_subsamples = TRUE,
  comp_function = zhouRR)

#Response ratios among cover types based on Zhou intervals
jorn15s.te$comparisons <- jorn15s.te$comparisons %>%
  mutate(response = recode_factor(response, "    CELLOBIO" = "CELLOBIO",
    "    respiration" = "respiration") %>%
    fct_relevel(c("soil_moisture", "percentC", "percentN", "nh4", "no3", "po4",
      "MBC", "MBN", "AGLUC", "BGLUC", "CELLOBIO", "XYLO", "NAG", "PHOS", 
      "PhenolOx", "Perox", "respiration"))) %>%
  arrange(response) %>%
  mutate(response = fct_recode(response, "Soil moisture" = "soil_moisture",
    "Bulk soil C" = "percentC", "Bulk soil N" = "percentN",
    "NH4+" = "nh4", "NO3-" = "no3", "PO43-" = "po4",
    "Microbial biomass C" = "MBC", "Microbial biomass N" = "MBN",
    "a-glucosidase" = "AGLUC", "B-glucosidase" = "BGLUC",
    "Cellobiohydrolase" = "CELLOBIO", "Xylosidase" = "XYLO",
    "N-acetyl-glucosaminidase" = "NAG", "Phosphatase" = "PHOS",
    "Phenol oxidase" = "PhenolOx", "Peroxidase" = "Perox",
    "Soil respiration" = "respiration"),
    comparison = fct_recode(comparison, "Creosote / Bare" = "Creosote - Bare",
    "Mesquite / Bare" = "Mesquite - Bare",
    "Mesquite / Creosote" = "Mesquite - Creosote"))

plotRR(jorn15s.te)

#Response ratios wet vs. dry averaged across cover types using Zhou intervals
jorn15dw0 <-
  jorn15 %>%
  mutate(drywet = case_when(doy %in% c(175,182) ~ "dry", doy %in% c(189,196) ~ "wet"),
    cover = case_when(Type == "Bare" ~ 0.88, Type == "Creosote" ~ 0.11, 
      Type == "Mesquite" ~ 0.01)) %>%
  filter(!is.na(drywet)) %>%
  select(doy, cover, plot_num, drywet, MBC, MBN, nh4, no3, po4, BGLUC, NAG, 
    CELLOBIO, XYLO, AGLUC, PHOS, soil_moisture)
jorn15dw0[5:16] <- jorn15dw0[5:16] * t(jorn15dw0$cover)

jorn15dw.te <- jorn15dw0 %>%
  group_by(doy, plot_num, drywet) %>%
  summarize_at(c("MBC", "MBN", "nh4", "no3", "po4", "BGLUC", "NAG", "CELLOBIO",
    "XYLO", "AGLUC", "PHOS", "soil_moisture"), .funs = list(sum = sum)) %>%
treateffect(MBC_sum + MBN_sum + nh4_sum + no3_sum + po4_sum +
    BGLUC_sum + NAG_sum + CELLOBIO_sum + XYLO_sum + AGLUC_sum + PHOS_sum + 
      soil_moisture_sum ~ drywet,
    block = "plot_num", average_subsamples = TRUE,
    comp_function = zhouRR)

jorn15dw.te$comparisons <- jorn15dw.te$comparisons %>%
  mutate(response = recode_factor(response, "    NAG_sum" = "NAG_sum",
      "    soil_moisture_sum" = "soil_moisture_sum") %>%
    fct_relevel(c("soil_moisture_sum", "MBC_sum", "MBN_sum", "nh4_sum", "no3_sum", 
    "po4_sum", "AGLUC_sum", "BGLUC_sum", "CELLOBIO_sum", "XYLO_sum", "NAG_sum"))) %>%
  arrange(response) %>%
  mutate(response = fct_recode(response, "Soil moisture" = "soil_moisture_sum",
    "NH4+" = "nh4_sum", "NO3-" = "no3_sum", "PO43-" = "po4_sum",
    "Microbial biomass C" = "MBC_sum", "Microbial biomass N" = "MBN_sum",
    "a-glucosidase" = "AGLUC_sum", "B-glucosidase" = "BGLUC_sum",
    "Cellobiohydrolase" = "CELLOBIO_sum", "Xylosidase" = "XYLO_sum",
    "N-acetyl-glucosaminidase" = "NAG_sum", "Phosphatase" = "PHOS_sum"),
    comparison = fct_recode(comparison, "Wet / Dry" = "wet - dry"))

plotRR(jorn15dw.te)

#gravimetric soil moisture
treateffect(jorn15, soil_moisture ~ Type, times = "doy",
    block = "plot_num", pool_variance = c("doy", "Type")) %>%
plot(dodge = 2) +
  scale_x_continuous(
    breaks = c(121, 155, 168, 175, 182, 189, 196),
    limits = c(120, 198),
    labels = c("May 1", "Jun 4", "Jun 17", "Jun 24", "Jul 1", "Jul 8", "Jul 15"))

# tower weather data
precip_2015 <-
    read_csv("JER_Biomet_daily_MayJuneJuly_2015_20190906.csv") %>%
  mutate(doy = yday(date))

#rainfall
precip_2015 %>%
  ggplot(aes(doy, P_rain_1_1_1)) +
  geom_bar(stat="identity", fill = "navyblue") +
  scale_x_continuous(
    breaks = c(121, 155, 168, 175, 182, 189, 196),
    labels = c("May 1", "Jun 4", "Jun 17", "Jun 24", "Jul 1", "Jul 8", "Jul 15"))

#TDR volumetric soil moisture
precip_2015 %>%
  ggplot(aes(doy, SWC_1_1_1)) +
  geom_line() +
  scale_x_continuous(
    breaks = c(121, 155, 168, 175, 182, 189, 196),
    labels = c("May 1", "Jun 4", "Jun 17", "Jun 24", "Jul 1", "Jul 8", "Jul 15"))

#NEE
nee <-
  read_csv("JER_NEE_Reco_GPP_MJJ_2015_20190830.csv")

nee %>%
  mutate(GPP_daily = -GPP_daily) %>%
  pivot_longer(cols = c("NEE_daily", "GPP_daily", "Reco_daily")) %>%
  ggplot(aes(DoY, value, col = name)) +
  geom_point() +
  geom_line() +
  geom_hline(yintercept = 0) +
  scale_color_manual(values = c("forestgreen", "black", "brown")) +
  scale_x_continuous(
    breaks = c(121, 155, 168, 175, 182, 189, 196),
    limits = c(120, 198),
    labels = c("May 1", "Jun 4", "Jun 17", "Jun 24", "Jul 1", "Jul 8", "Jul 15"))

#bulk C and N
jorn15 %>% treateffect(percentC + percentN + CN_ratio ~ Type,
  block = "plot_num", average_subsamples = TRUE) %>%
plot(treatcol = c("black", "red", "orange"))

#respiration
jorn15 %>%
  filter(doy == 189 | doy == 196) %>%
  treateffect(respiration ~ Type | doy, block = "plot_num") %>%
  plot(treatcol = c("black", "red", "orange"))

#MBC, NH4+, NO3-
fig1.te <- treateffect(jorn15, MBC + nh4 + no3 ~ Type, times = "doy")
fig1.te$data$y[fig1.te$data$response == "MBC" & fig1.te$data$y > 400] <- NA
fig1.te$data$y[fig1.te$data$response == "no3" & fig1.te$data$y > 10] <- NA
plot(fig1.te, dodge = 2) + jorn15_scale()

#Exoenzyme activity (Fig. 4)
enz.te <- treateffect(jorn15, AGLUC + BGLUC + CELLOBIO + XYLO + NAG + 
  PHOS + PhenolOx + Perox ~ Type, times = "doy")
plot(enz.te, dodge = 2) + jorn15_scale()

#MBN - SUPP
mbn.te <- treateffect(jorn15, MBN ~ Type, times = "doy")
mbn.te$data$y[mbn.te$data$response == "MBN" & mbn.te$data$y > 40] <- NA
plot(mbn.te, dodge = 2) + jorn15_scale()

#Phosphate - SUPP
phos.te <- treateffect(jorn15, po4 ~ Type, times = "doy")
plot(phos.te, dodge = 2) + jorn15_scale()

#splom - SUPP
library(lattice)
pdf("splom.pdf", width = 26, height = 26, useDingbats = FALSE)
  splom(~b[c(7:12, 14:15, 17:25, 27)], groups = Type, data = b,
  col = treatcol_default(3))
dev.off()
