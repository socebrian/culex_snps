#sonia cebrian camison
#02/02/2025
#scebrian27@gmail.com sonia.cebrian@ebd.csic.es
#script to plot LD heatmap
#not figure worthy, just as checking

library(ggplot2)
library(tidyr)
library(dplyr)
library(tibble)
#input the data

ld_data <- read.table("2202pipiens123_1_maf05.ld", header=T)
ld_data <- arrange(ld_data, R2)
maint <- "chromosome 1 MAF>5%"
f = "2202pipiens123_1_maf05_LD50130.png"

fig <- ggplot(data=ld_data, aes(x=BP_A/1e6, y=BP_B/1e6, color=R2)) +
  geom_point(alpha=0.5) +
  labs(title= maint) +
  xlab("Position Mbp") +
  ylab("Position Mbp") +
  scale_color_viridis_c()+
  theme_bw()+
  xlim(50,130)+
  ylim(50,130)
  
fig
ggsave(f, device = "png", units = "in", plot=fig, width=5, height=4)


