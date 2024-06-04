##<<<<<<<<<<<<<<<<<<<<<<<<<HEAD

## load packages
library(tidyverse)
library(ggpubr)
library(ggprism)

## load data
dat <- read_csv("Data/Mouse-all-190224.csv")

## convert data to tidy format
dat <- dat %>% gather("miRNA", "expression", -c(cell_line, treatment)) 
dat$plotid <- paste(dat$cell_line, dat$treatment)
dat$plotid <- factor(dat$plotid, levels = c("MDA-MB231 HTM 6Gy_d55", "MDA-MB231 HTM 6Gy_d76",  "MDA-MB231 HTM ctrl",         
                                            "JIMT-1 HTM 6Gy_d34", "JIMT-1 HTM 6Gy_d34_anti_PDL1", "JIMT-1 HTM ctrl",           
                                            "NSG background ctrl"))
dat$cell_line <- factor(dat$cell_line, levels = c("MDA-MB231 HTM", "JIMT-1 HTM", "NSG background"))

## function to plot data and save as png
save_plots <- function(data, mir.input){
  
  dat <- data[data$miRNA == mir.input,]
  background <- dat[dat$cell_line == "NSG background",]$expression
  dat <- dat[dat$cell_line != "NSG background",]
  ## plot by cell line and treatment
  p1 <- dat %>% 
    ggplot(aes(plotid, expression, fill = cell_line)) +
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
  
  # ## aggregated plot
  # p2 <- dat %>% 
  #   ggbarplot(x = "cell_line", y = "expression", fill = "cell_line", width = 0.6, size = 1,
  #             add = "mean_sd", add.params = list(width = 0.2, size = 0.8)) +
  #   geom_point(shape = 1, size = 3)+
  #   geom_hline(yintercept = background, color = "red", size = 1, lty = 3) +
  #   theme_bw() +
  #   theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, face = "bold"),
  #         axis.line = element_line(color = 'black', linewidth = 1),
  #         panel.grid = element_blank(), 
  #         panel.border = element_blank(),
  #         plot.margin = margin(10, 10, 10, 15)) +
  #   scale_y_continuous(expand = c(0.1,0)) +
  #   xlab("") +
  #   ylab("MiRNA expression (a.u.)") +
  #   labs(title = paste(str_replace_all(mir.input, "hsa-mir", "miR"), "- aggregated")) +
  #   scale_fill_manual(values = c("white", "dimgrey")) +
  #   guides(fill=guide_legend(title="Cell line"))
  # 
  # p3 <- ggarrange(p1,p2,nrow = 1, align = "h")
  svg(paste("Results/", mir.input, ".svg", sep =""), width = 4, height = 6)
  print(p1)
  dev.off()
}

## apply function over all miRNAs
lapply(unique(dat$miRNA) , save_plots, data = dat)





