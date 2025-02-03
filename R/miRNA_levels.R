##<<<<<<<<<<<<<<<<<<<<<<<<<HEAD
## This script contains the data analysis of miRNA serum levels in a humanized tumor mouse (HTM) model. This proof of principle investigation 
## was conducted to establish the measurement of exosomal miRNAs in HTM after immunotherapy with or without supplementing radiotherapy (abscopal effect)

## load packages
library(tidyverse)
library(ggpubr)
library(ggprism)

## load data
dat <- read_csv("Data/miRNA_serum_levels_190224.csv")

## convert data to tidy format
dat <- dat %>% gather("miRNA", "expression", -c(mouse, treatment)) 
dat$plotid <- paste(dat$mouse, dat$treatment)
dat$plotid <- factor(dat$plotid, levels = c("MDA-MB231 HTM 6Gy_d55", "MDA-MB231 HTM 6Gy_d76",  "MDA-MB231 HTM ctrl",         
                                            "JIMT-1 HTM 6Gy_d34", "JIMT-1 HTM 6Gy_d34_anti_PDL1", "JIMT-1 HTM ctrl",           
                                            "NSG background ctrl"))
dat$mouse <- factor(dat$mouse, levels = c("MDA-MB231 HTM", "JIMT-1 HTM", "NSG background"))

## function to plot data and save as png
save_plots <- function(data, mir.input){
  
  dat <- data[data$miRNA == mir.input,]
  background <- dat[dat$mouse == "NSG background",]$expression
  dat <- dat[dat$mouse != "NSG background",]
  ## plot by cell line and treatment
  p <- dat %>% 
    ggplot(aes(plotid, expression, fill = mouse)) +
    geom_bar(stat = "identity", width = 0.75, color = "black", size = 1)+
    geom_hline(yintercept = background, color = "red", size = 1, lty = 3) +
    theme_prism() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1,face = "bold"),
          legend.position = "none",
          plot.margin = margin(10, 10, 10, 60)) +
    scale_x_discrete(guide = "prism_offset") +
    scale_y_continuous(guide = "prism_offset_minor", expand = c(0.03,0)) +
    xlab("") +
    ylab("MiRNA expression (a.u.)") +
    labs(title = str_replace_all(mir.input, "hsa-mir", "miR")) +
    scale_fill_manual(values = c("white", "dimgrey")) 
  
  svg(paste("Results/", mir.input, ".svg", sep =""), width = 4, height = 6)
  print(p)
  dev.off()
}

## apply function over all miRNAs
lapply(unique(dat$miRNA) , save_plots, data = dat)





