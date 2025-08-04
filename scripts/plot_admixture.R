## Sonia cebrian camison
## proyecto culex snps
# 25 junio 2025
# Script to plot admixture analysis

## libraires
library(ggplot2)
library(tidyr)
library(dplyr)

# fuiles
cv_table <- read.table("cv_table_wholegen.txt", col.names = c("K", "CVerr"))
samples <- read.table("pipiens1234_all_PCA.fam", sep = " ",
                      col.names = c("FID", "IID", "loc", "mother", "sex", "pheno"))[, c("IID", "loc")]
qvals <- read.table(paste0("pipiens1234_all_PCA.2(1).Q"),
                    col.names = c("Q1", "Q2"))
metadata <- read.table("metadata_pipiens1234.txt", sep="\t", header=T)


# plot cv values ----
(cvplot <- ggplot() +
  geom_line(data = cv_table,
            aes(x = K, y = CVerr)) +
  geom_point(data = cv_table,
             aes(x = K, y = CVerr),
             fill = "grey", shape = 21, size = 2.5) +
  xlab(paste0("Number of Populations")) +
  ylab(paste0("Cross Validation Error")) +
  scale_x_continuous(n.breaks = 8) +
  theme_bw())
ggsave(filename = "wholegen_0miss.cvplot.png", plot = cvplot)


#admixture plots
# loop through k 3 to 6 (ORDER BY BELONGING TO CLUSTERS) ----
my_colors = c("Q2"="#00BFC4", "Q1"="#F8766D", "Q3"= "#ffbc42", "Q4" = "#85AE09", "Q5"="#335C67", "Q6"="grey")
my_colors = c("Q2"="#44B7C2", "Q1"="#D38217")
maint= " " 
file = "SEG25pipiens1234.plot.2K.png"
# loop through k 3 to 6 (ORDER BY SAMPLE ORDER IN K2 ----
# Step 1: Get sample order from K = 2
# STEP 1: Determine sample order based on K = 2
k_order <- 2
qvals_k2 <- read.table(paste0("pipiens1234_0miss_5maf_PCA.", k_order, ".Q"),
                       col.names = paste0("Q", seq(1:k_order)))
qvals_k2$IID <- samples$IID

admix_k2 <- left_join(samples, qvals_k2, by = "IID") %>%
  gather(Q, value, starts_with("Q"))

# Order by Q1 membership (descending)
ordering_df <- admix_k2 %>%
  filter(Q == "Q1") %>%
  arrange(desc(value))

# Save this fixed sample order
IID_order <- ordering_df$IID

# STEP 2: Loop through all desired K values using this fixed order
for (k in 2:6) {
  qvals <- read.table(paste0("pipiens1234_0miss_5maf_PCA.", k, ".Q"),
                      col.names = paste0("Q", seq(1:k)))
  qvals$IID <- samples$IID
  
  admix <- left_join(samples, qvals, by = "IID") %>%
    gather(Q, value, starts_with("Q"))
  
  # Apply fixed sample order
  admix$IID <- factor(admix$IID, levels = IID_order)
  
  admix_plot <- ggplot(admix, aes(x = IID, y = value, fill = factor(Q))) +
    geom_bar(stat = "identity", position = "stack") +
    xlab("Sample") + ylab("Ancestry") + labs(title = paste("Whole genome Cx. pipiens (MISS=0, MAF=>5%), K =", k)) +
    scale_fill_manual(values = my_colors) +
    theme_bw() + theme(axis.text.x = element_text(angle = 60, hjust = 1))
  
  ggsave(filename = paste0("wholegen0miss_p1234.plot.K", k, ".order.png"),
         plot = admix_plot, width = 35, height = 10, units = "cm")
}

#update metadata

metadata_2p <- merge(qvals_k2, metadata_pipiens1234, by.x = "IID", by.y = "ID", all.x = T)
str(metadata_2p)
metadata_2p <- metadata_2p[,1:10]

write.table(metadata_2p, "metadata_pipiens1234_0miss5maf.txt", row.names = F )


