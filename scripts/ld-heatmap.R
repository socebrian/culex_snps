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

ld_data <- read.table("pipiens1234_3_maf03.ld", header=T)
ld_data <- arrange(ld_data, R2)
maint <- "chromosome 3 MAF>3%"
f = "pipiens1234_3_maf03_hor.png"

fig <- ggplot(data=ld_data, aes(x=BP_A/1e6, y=BP_B/1e6, color=R2)) +
  geom_point(alpha=0.5) +
  labs(title= maint) +
  xlab("Position Mbp") +
  ylab("Position Mbp") +
  scale_color_viridis_c()+
  theme_bw()
 # xlim(50,130)+
#  ylim(50,130)
  
fig
ggsave(f, device = "png", units = "in", plot=fig, width=5, height=4)



#plot with ohorizontal genome position

ld_data <- ld_data %>%
  mutate(
    mid = (BP_A + BP_B) / 2 / 1e6,           # midpoint in Mbp for x-axis
    dist = abs(BP_A - BP_B) / 1e6            # distance between SNPs for y-axis
  )

fig2 <- ggplot(ld_data, aes(x = mid, y = dist, color = R2)) +
  geom_point(alpha = 0.5, size = 0.6) +
  scale_color_viridis_c(name = "R²") +
  labs(
    title = "Horizontally projected LD triangle",
    x = "Genomic position (Mb)",
    y = "Distance between SNPs (Mb)"
  ) +
  theme_bw()
ggsave(f, device = "png", units = "in", plot=fig2, width=11, height=4)
