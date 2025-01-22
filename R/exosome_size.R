##<<<<<<<<<<<<<<<<<<<<<<<<<HEAD

## load packages
library(tidyverse)
library(ggpubr)
library(ggprism)
library(ggh4x)

## load data
dat <- read_csv("Data/exosome_size.csv")

## plot data
p <- dat %>% 
  ggplot(aes(x= dia_nm, y = conc_part_ml, color = mouse)) + 
  geom_line() +
  scale_x_continuous(guide = "axis_minor") + 
  scale_color_manual(values = c("#2166ac", "#b2182b", "black")) +
  theme_bw(base_size = 14) +
  theme(legend.position = "none") +
  xlab("Vesicle size (nm)") +
  ylab("Concentration (particles/ml)")

## save plot 
svg("Results/Figure1/exosome_size.svg", width = 4, height = 4)
print(p)
dev.off()

