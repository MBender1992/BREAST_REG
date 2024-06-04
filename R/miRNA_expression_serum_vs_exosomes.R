##<<<<<<<<<<<<<<<<<<<<<<<<<HEAD

## load packages
library(tidyverse)
library(ggpubr)
library(ggprism)

## load data
dat <- read_csv("Data/Fig-5-miR-Serum-vs-exosome.csv")

## convert data to tidy format and wrangle data
dat <- dat %>% 
  gather("miRNA", "expression", -c(cell_line, treatment, type)) %>%
  mutate(expression = as.numeric(expression)) %>%
  mutate(unrel = ifelse(is.na(expression), 1, NA)) %>%
  mutate(mismatch = ifelse(stringr::str_detect(miRNA, "^\\*m"), "One",NA)) %>%
  mutate(mismatch = ifelse(stringr::str_detect(miRNA, "\\*\\*"), "Two",mismatch)) %>%
  mutate(mismatch = ifelse(is.na(mismatch), "Zero",mismatch)) %>%
  mutate(treatment = factor(treatment, levels = c("ctrl", "6Gy_34d", "6Gy_34d_anti_PDL1"), 
                            labels = c("Untreated", "6 Gy after 34 days", "6 Gy after 34 days + anti PD-L1"))) %>%
  mutate(miRNA = str_replace_all(.$miRNA, "mir", "miR"))

## filter data for small plots
dat_small <- dat %>% 
  filter(miRNA %in% c("miR-29b-3p","miR-34c-5p", "miR-382-5p","miR-203a-3p", "**miR-378g")) %>%
  mutate(miRNA = factor(miRNA, levels = c("miR-29b-3p","miR-34c-5p", "miR-382-5p","miR-203a-3p", "**miR-378g")))


##************************************
##define functions

## set theme for small plots
theme_small_plot <- function(...){
  theme_prism(...) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
          axis.line.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(), 
          legend.position = "none")
}

## function to generate small plots (miRNAs which are overexpressed in comparison to serum levels)
small_plots <- function(data, trt){
  data %>% 
    filter(treatment == trt & miRNA %in% c("miR-29b-3p","miR-34c-5p", "miR-382-5p","miR-203a-3p", "**miR-378g")) %>%
    mutate(miRNA = factor(miRNA, levels = c("miR-29b-3p","miR-34c-5p", "miR-382-5p","miR-203a-3p", "**miR-378g"))) %>%
    ggplot(aes(miRNA, expression, fill = type)) +
    geom_bar(stat = "identity", position = position_dodge(), color = "black") +
    geom_point(aes(y = unrel*50), position = position_dodge(0.9), show.legend = FALSE) +
    theme_prism() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), 
          legend.position = "none") +
    scale_x_discrete(guide = "prism_offset", expand = c(0.15,0.021)) +
    scale_y_continuous(guide = "prism_offset_minor", expand = c(0.03,0)) +
    xlab("") +
    ylab("miRNA expression (a.u.)") +
    scale_fill_manual(values = c("white", "dimgrey")) 
}

save_small_plots <- function(data, trt){
  lapply(levels(data$miRNA), function(i){
    dat_plot <- data[data$miRNA == i & data$treatment == trt,]
    p <- dat_plot %>% 
      ggplot(aes(miRNA, expression, fill = type)) +
      geom_bar(stat = "identity", position = position_dodge(), color = "black", width =0.4) +
      geom_point(aes(y = unrel*5), position = position_dodge(0.9), show.legend = FALSE) +
      theme_prism() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), 
            legend.position = "none") +
      scale_x_discrete(guide = "prism_offset", expand = c(0.15,0.021)) +
      scale_y_continuous(guide = "prism_offset_minor", expand = c(0.03,0)) +
      xlab("") +
      ylab("miRNA expression (a.u.)") +
      scale_fill_manual(values = c("white", "dimgrey")) 
    
    png(paste("Results/", unique(data$cell_line), trt, str_remove(i, "\\*+"), "high_exo_ratio", ".png", sep ="_"), units = "in", width = 3, height = 4, res = 600)
    print(p)
    dev.off()
  })
}

## function to save plots based on treatment condition
save_plots <- function(data, trt){
  dat_plot <- data %>%
    filter(treatment == trt & mismatch == "Zero")
  p1 <- dat_plot %>%
    ggplot(aes(factor(miRNA), expression, fill = type)) +
    geom_bar(stat = "identity", position = position_dodge(), color = "black") +
    geom_point(aes(y = unrel*50), position = position_dodge(0.9), show.legend = FALSE) +
    theme_prism() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), 
          legend.position = c(0.05, 0.95), legend.justification = c(0.05, 0.95)) +
    scale_x_discrete(guide = "prism_offset") +
    scale_y_continuous(guide = "prism_offset_minor", expand = c(0.01,0), limits = c(0,20000)) +
    xlab("") +
    ylab("miRNA expression (a.u.)") +
    scale_fill_manual(values = c("white", "dimgrey")) 
  
  dat_plot <- data %>%
    filter(treatment == trt & mismatch == "One") 
  p2 <- dat_plot %>%
    ggplot(aes(factor(miRNA), expression, fill = type)) +
    geom_bar(stat = "identity", position = position_dodge(), color = "black") +
    geom_point(aes(y = unrel*50), position = position_dodge(0.9), show.legend = FALSE) +
    theme_small_plot() + 
    scale_x_discrete(guide = "prism_offset") +
    scale_y_continuous(guide = "prism_offset_minor", expand = c(0.01,0), limits = c(0,20000)) +
    xlab("") +
    ylab("") +
    scale_fill_manual(values = c("white", "dimgrey")) 
   
  
  dat_plot <- data %>%
    filter(treatment == trt & mismatch == "Two")
  p3 <- dat_plot %>%
    ggplot(aes(factor(miRNA), expression, fill = type)) +
    geom_bar(stat = "identity", position = position_dodge(), color = "black") +
    geom_point(aes(y = unrel*50), position = position_dodge(0.9), show.legend = FALSE) +
    theme_small_plot() + 
    scale_x_discrete(guide = "prism_offset") +
    scale_y_continuous(guide = "prism_offset_minor", expand = c(0.01,0), limits = c(0,20000)) +
    xlab("") +
    ylab("") +
    scale_fill_manual(values = c("white", "dimgrey"))
    
  
  p4 <- ggpubr::ggarrange(p1,p2, p3, ncol = 3, nrow = 1,  widths = c(0.76, 0.13, 0.11))
  p4 <- annotate_figure(p4, top = text_grob(paste(dat_plot$cell_line, dat_plot$treatment), 
                                      color = "black", face = "bold", size = 14))
  
  png(paste("Results/", dat_plot$cell_line, dat_plot$treatment, "serum_exosome_miRNA", ".png", sep ="_"), units = "in", width = 18, height = 6, res = 600)
  print(p4)
  dev.off()
}

##************************************
## analysis and saving of plots
lapply(levels(dat$treatment) , save_plots, data = dat)

## saving small plots
lapply(levels(dat$treatment) , save_small_plots, data = dat_small)









