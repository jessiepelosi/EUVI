################################
## PopulationStructure.R      ##
## Jessie Pelosi              ##
## Last updated July, 2025    ##
##                            ##
################################

library(vcfR)
library(adegenet)
library(dplyr)
library(ggplot2)
# Using data that includes LAFI samples: 

thinned_1SDremoved <- read.vcfR("populations.snps.thinned.1SD-removed.noLAFI.maf0.01.NoExHet.recode.vcf")

popmap3 <- read.delim("samples_to_keep.txt", sep = " ")

genind <- vcfR2genind(thinned_1SDremoved)
pop(genind) <- popmap3$V4
genind

# make a normal PCA of the genomic data 

# replace missing data with mean genotype value 
genind.nomissing <- tab(genind, freq=TRUE, NA.method="mean")

pca.test <- dudi.pca(genind.nomissing, center = T, scale = F, nf = 50, scannf=F)

cols <- c("skyblue2", "dodgerblue3", "orange", "gold")

s.class(pca.test$li, fac=pop(genind), clabel = 0.5, col=cols, axesell = F)
add.scatter.eig(pca.test$eig[1:50], 3, 1, 2, ratio = 0.3)



pca.df <- pca.test$li
pca.df$sample <- row.names(pca.test$li)
pca.df$pop <- popmap3$V2
pca.df$region <- popmap3$V3
pca.df$species <- popmap3$V4

ggplot(data = pca.df, mapping = aes(x = Axis1, y = Axis2, color = species)) + geom_point() + theme_bw() +
  scale_color_manual(values = c("orange", "skyblue2", "gold", "dodgerblue2")) + 
  xlab(paste0("PC1 (", signif(pca.test$eig[1], 3), "% of variation)")) + 
  ylab(paste0("PC1 (", signif(pca.test$eig[2], 3), "% of variation)")) +
  stat_ellipse() +
  theme(axis.text = element_text(size = 12), legend.text = element_text(size = 12), 
        axis.title = element_text(size = 14), legend.title = element_text(size = 14)) +
  labs(color = "Region")

ggsave("EUVI_PCA.pdf", height = 5, width = 8)
  
# use DAPC -
# Note that for DAPC, we don't need to replace missing data, there is an internal function that deals with it
grp <- find.clusters(genind, max.n.clust = 15, n.pca = 150)
dapc1 <- dapc(genind, grp$grp)
scatter(dapc1)
compoplot(dapc1)

library(ggh4x)
postprobs <- as.data.frame(round(dapc1$posterior, 4))
clusters <- tibble::rownames_to_column(postprobs, var = "ind") %>% 
  mutate(site = popmap3$V4)
clusters_long <- reshape2::melt(clusters)
ggplot(data = clusters_long, mapping = aes(x = ind, y = value, fill = factor(variable))) + 
  geom_col(color = "gray", size=0.01) + 
  facet_nested(~site, switch = "x", nest_line = element_line(size = 1, lineend = "round"),
               scales = "free", space = "free") +
  theme_minimal() +
  labs(x = "Samples", y = "Membership Probability") +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_discrete(expand = expansion(add = 0.5)) +
  theme(panel.spacing.x = unit(0.18, "lines"), 
        panel.grid = element_blank(), axis.text.x = element_text(angle = 90))

ggsave('dapc-k=2.png', width = 40, height = 5, limitsize = F, dpi = 300)

##

grp <- find.clusters(genind, max.n.clust = 15, n.pca = 300)
dapc1 <- dapc(genind, grp$grp)
scatter(dapc1)
compoplot(dapc1)

library(ggh4x)
postprobs <- as.data.frame(round(dapc1$posterior, 4))
clusters <- tibble::rownames_to_column(postprobs, var = "ind") %>% 
  mutate(site = popmap2$STRATA)
clusters_long <- reshape2::melt(clusters)
ggplot(data = clusters_long, mapping = aes(x = ind, y = value, fill = factor(variable))) + 
  geom_col(color = "gray", size=0.01) + 
  facet_nested(~site, switch = "x", nest_line = element_line(size = 1, lineend = "round"),
               scales = "free", space = "free") +
  theme_minimal() +
  labs(x = "Samples", y = "Membership Probability") +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_discrete(expand = expansion(add = 0.5)) +
  theme(panel.spacing.x = unit(0.18, "lines"), 
        panel.grid = element_blank(), axis.text.x = element_text(angle = 90))

ggsave('dapc-k=3.png', width = 40, height = 5, limitsize = F, dpi = 300)

### Plot ADMIXTURE RESULTS

library(StructuRly)
library(popkin)

crossVal <- read.delim("CVs.txt", header = F, sep = "\t")
ggplot(data = crossVal, mapping = aes(x=V1, y = V2)) + geom_point() + geom_line() +
  theme_bw() + xlab("K") + ylab("Cross Validation Error")

ggsave("CV_error.pdf", height = 5, width = 7)


k2 <- read.table("population_structure/populations.snps.thinned.1SD-removed.noLAFI.maf0.01.NoExHet.admix.2.Q")
k3 <- read.table("population_structure/populations.snps.thinned.1SD-removed.noLAFI.maf0.01.NoExHet.admix.3.Q")

#popmap <- read.delim("population_map.txt", header = F)
LAFI_samples <- read.delim("LAFI_samples.txt", sep = "\t", header = F)
popmap2 <- read.delim("population_map.1SD-removed-edited.txt") 
#popmap2 <- subset(popmap, !(V1 %in% STACKS_remove$sample))
popmap3 <- subset(popmap2, !(INDIVIDUALS %in% LAFI_samples$V1))

k2_pops <- cbind(popmap3, k2)

neworder <- c("California", "PNW", "Greece", "Turkey")
ordered <- data.frame(matrix(nrow=0,ncol=6))
for (i in neworder) {
  ordered <- rbind(ordered,k2_pops[k2_pops$REGION==i,])
}
Q <- ordered[,5:6]
samples <- ordered[,1:4]

plot_admix(Q, col=c("skyblue2","orange"), labs = samples$REGION, 
           labs_cex = 1, labs_sep=TRUE,labs_lwd = 2, labs_col="black", leg_omit = TRUE)




k3_pops <- cbind(popmap3, k3)

neworder <- c("California", "PNW", "Greece", "Turkey")
ordered <- data.frame(matrix(nrow=0,ncol=7))
for (i in neworder) {
  ordered <- rbind(ordered,k3_pops[k3_pops$REGION==i,])
}
Q <- ordered[,5:7]
samples <- ordered[,1:4]

plot_admix(Q, col=c("orange", "skyblue2","darkorange3"), labs =samples$REGION, 
           labs_cex = 1, labs_sep = TRUE, labs_lwd=2, labs_col = "black", xlab_line = 2, leg_omit = TRUE)

