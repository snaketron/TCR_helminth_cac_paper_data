.libPaths(new = "/mnt/nfs/simo/rpack_4_4/")
require(ggplot2)
require(reshape2)
require(patchwork)
require(IgGeneUsage)
require(ggrepel)
require(ggforce)
require(scales)

#### Prepare usage data ####

# consider only genes with pmax>=0.99 and |delta|>=0.01
pmax <- 0.99
delta <- 0.01

# load data
mixcr_clones <- get(load("data/mixcr_clones.RData"))
mixcr_clones$sample_id <- paste0(mixcr_clones$experiment, "|", mixcr_clones$sample, "|", mixcr_clones$locus)

# get the V and J gene segments
mixcr_clones$v <- unlist(lapply(X = mixcr_clones$allVHitsWithScore, FUN = function(x) {
  v <- unlist(strsplit(x = x, split = ","))[1]
  v <- unlist(strsplit(x = v, split = "\\*"))[1]
  return(v)
}))
mixcr_clones$j <- unlist(lapply(X = mixcr_clones$allJHitsWithScore, FUN = function(x) {
  v <- unlist(strsplit(x = x, split = ","))[1]
  v <- unlist(strsplit(x = v, split = "\\*"))[1]
  return(v)
}))

# count usage
v <- aggregate(readCount~sample_id+locus+group+v, data = mixcr_clones, FUN = sum)
j <- aggregate(readCount~sample_id+locus+group+j, data = mixcr_clones, FUN = sum)

# split by gene segment
va <- v[v$locus=="TRA",c("sample_id", "group", "v", "readCount")]
vb <- v[v$locus=="TRB",c("sample_id", "group", "v", "readCount")]
ja <- j[j$locus=="TRA",c("sample_id", "group", "j", "readCount")]
jb <- j[j$locus=="TRB",c("sample_id", "group", "j", "readCount")]

colnames(va) <- c("individual_id", "condition", "gene_name", "gene_usage_count")
colnames(vb) <- c("individual_id", "condition", "gene_name", "gene_usage_count")
colnames(ja) <- c("individual_id", "condition", "gene_name", "gene_usage_count")
colnames(jb) <- c("individual_id", "condition", "gene_name", "gene_usage_count")

#### DGU with IgGeneUsage ####

# 4 x DGU 
dgu_ja <- IgGeneUsage::DGU(ud = ja,
                           mcmc_warmup = 1500,
                           mcmc_steps = 3000,
                           mcmc_chains = 4,
                           mcmc_cores = 4,
                           hdi_lvl = 0.95,
                           adapt_delta = 0.99,
                           max_treedepth = 11,
                           paired = FALSE)
dgu_jb <- IgGeneUsage::DGU(ud = jb,
                           mcmc_warmup = 1500,
                           mcmc_steps = 3000,
                           mcmc_chains = 4,
                           mcmc_cores = 4,
                           hdi_lvl = 0.95,
                           adapt_delta = 0.99,
                           max_treedepth = 11,
                           paired = FALSE)
dgu_va <- IgGeneUsage::DGU(ud = va,
                           mcmc_warmup = 1500,
                           mcmc_steps = 3000,
                           mcmc_chains = 4,
                           mcmc_cores = 4,
                           hdi_lvl = 0.95,
                           adapt_delta = 0.99,
                           max_treedepth = 11,
                           paired = FALSE)
dgu_vb <- IgGeneUsage::DGU(ud = vb,
                           mcmc_warmup = 1500,
                           mcmc_steps = 3000,
                           mcmc_chains = 4,
                           mcmc_cores = 4,
                           hdi_lvl = 0.95,
                           adapt_delta = 0.99,
                           max_treedepth = 11,
                           paired = FALSE)

dgu <- list(dgu_va = dgu_va, dgu_vb = dgu_vb, 
            dgu_ja = dgu_ja, dgu_jb = dgu_jb,
            data_va = va, data_vb = vb, 
            data_ja = ja, data_jb = jb)
save(dgu, file = "output/dgu_data.RData", compress = TRUE)




#### Figures ####
dgu <- get(load("output/dgu_data.RData"))


# figures
x <- merge(x = merge(x = dgu$dgu_ja$gu[dgu$dgu_ja$gu$condition == "CAC",], 
                     y = dgu$dgu_ja$gu[dgu$dgu_ja$gu$condition == "CAC+Hp",], by = "gene_name"), 
           y = dgu$dgu_ja$dgu_prob, by = "gene_name")

fig_traj <- ggplot(data = x)+
  geom_errorbar(aes(x = prob_mean.x, y = prob_mean.y, ymin = prob_L.y, ymax = prob_H.y))+
  geom_errorbarh(aes(x = prob_mean.x, y = prob_mean.y, xmin = prob_L.x, xmax = prob_H.x))+
  geom_point(aes(x = prob_mean.x, y = prob_mean.y, col = pmax))+
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", col = "darkgray")+
  scale_color_distiller(name = expression(pi), palette = "Spectral", limits = c(0, 1))+
  theme_bw(base_size = 10)+
  scale_x_continuous(name = "CAC", limits = c(0, 0.065))+
  scale_y_continuous(name = "CAC+Hpb", limits = c(0, 0.065))+
  ggtitle(label = "TRAJ")+
  theme(legend.position = "top")


x <- merge(x = merge(x = dgu$dgu_jb$gu[dgu$dgu_jb$gu$condition == "CAC",], 
                     y = dgu$dgu_jb$gu[dgu$dgu_jb$gu$condition == "CAC+Hp",], by = "gene_name"), 
           y = dgu$dgu_jb$dgu_prob, by = "gene_name")

fig_trbj <- ggplot(data = x)+
  geom_errorbar(aes(x = prob_mean.x, y = prob_mean.y, ymin = prob_L.y, ymax = prob_H.y))+
  geom_errorbarh(aes(x = prob_mean.x, y = prob_mean.y, xmin = prob_L.x, xmax = prob_H.x))+
  geom_point(aes(x = prob_mean.x, y = prob_mean.y, col = pmax))+
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", col = "darkgray")+
  scale_color_distiller(name = expression(pi), palette = "Spectral", limits = c(0, 1))+
  theme_bw(base_size = 10)+
  scale_x_continuous(name = "CAC", limits = c(0, 0.22))+
  scale_y_continuous(name = "CAC+Hpb", limits = c(0, 0.22))+
  ggtitle(label = "TRBJ")+
  theme(legend.position = "top")



x <- merge(x = merge(x = dgu$dgu_va$gu[dgu$dgu_va$gu$condition == "CAC",], 
                     y = dgu$dgu_va$gu[dgu$dgu_va$gu$condition == "CAC+Hp",], 
                     by = "gene_name"), 
           y = dgu$dgu_va$dgu_prob, by = "gene_name")

y <- x[abs(x$es_prob_mean)>delta & x$pmax>=pmax,]
y$gene_name_short <- gsub(pattern = "TRA|TRB", replacement = '', x = y$gene_name)

fig_trav <- ggplot()+
  geom_errorbar(data = x, aes(x = prob_mean.x, y = prob_mean.y, ymin = prob_L.y, ymax = prob_H.y))+
  geom_errorbarh(data = x, aes(x = prob_mean.x, y = prob_mean.y, xmin = prob_L.x, xmax = prob_H.x))+
  geom_point(data = x, aes(x = prob_mean.x, y = prob_mean.y, col = pmax))+
  geom_text_repel(data = y, aes(x = prob_mean.x, y = prob_mean.y, label = gene_name_short), 
                  size = 2.5, min.segment.length = 0, max.overlaps = 100)+
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", col = "darkgray")+
  scale_color_distiller(name = expression(pi), palette = "Spectral", limits = c(0, 1))+
  theme_bw(base_size = 10)+
  scale_x_continuous(name = "CAC", limits = c(0, 0.17))+
  scale_y_continuous(name = "CAC+Hpb", limits = c(0, 0.17))+
  ggtitle(label = "TRAV")+
  theme(legend.position = "top")




x <- merge(x = merge(x = dgu$dgu_vb$gu[dgu$dgu_vb$gu$condition == "CAC",], 
                     y = dgu$dgu_vb$gu[dgu$dgu_vb$gu$condition == "CAC+Hp",], by = "gene_name"), 
           y = dgu$dgu_vb$dgu_prob, by = "gene_name")
y <- x[abs(x$es_prob_mean)>delta & x$pmax>=pmax,]
y$gene_name_short <- gsub(pattern = "TRA|TRB", replacement = '', x = y$gene_name)

fig_trbv <- ggplot(data = x)+
  geom_errorbar(aes(x = prob_mean.x, y = prob_mean.y, ymin = prob_L.y, ymax = prob_H.y))+
  geom_errorbarh(aes(x = prob_mean.x, y = prob_mean.y, xmin = prob_L.x, xmax = prob_H.x))+
  geom_point(aes(x = prob_mean.x, y = prob_mean.y, col = pmax))+
  geom_text_repel(data = y, aes(x = prob_mean.x, y = prob_mean.y, label = gene_name_short), 
                  size = 2.5, min.segment.length = 0, max.overlaps = 100)+
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", col = "darkgray")+
  scale_color_distiller(name = expression(pi), palette = "Spectral", limits = c(0, 1))+
  theme_bw(base_size = 10)+
  scale_x_continuous(name = "CAC", limits = c(0, 0.24))+
  scale_y_continuous(name = "CAC+Hpb", limits = c(0, 0.24))+
  ggtitle(label = "TRBV")+
  theme(legend.position = "top")


rm(x)
# patch
fig_usage <- (fig_trav|fig_trbv)
fig_usage_supplement <- (fig_traj|fig_trbj)





x <- rbind(merge(x = merge(x = dgu$dgu_ja$gu[dgu$dgu_ja$gu$condition == "CAC",], 
                           y = dgu$dgu_ja$gu[dgu$dgu_ja$gu$condition == "CAC+Hp",], by = "gene_name"), 
                 y = dgu$dgu_ja$dgu_prob, by = "gene_name"),
           merge(x = merge(x = dgu$dgu_jb$gu[dgu$dgu_jb$gu$condition == "CAC",], 
                           y = dgu$dgu_jb$gu[dgu$dgu_jb$gu$condition == "CAC+Hp",], by = "gene_name"), 
                 y = dgu$dgu_jb$dgu_prob, by = "gene_name"),
           merge(x = merge(x = dgu$dgu_va$gu[dgu$dgu_va$gu$condition == "CAC",], 
                           y = dgu$dgu_va$gu[dgu$dgu_va$gu$condition == "CAC+Hp",], by = "gene_name"), 
                 y = dgu$dgu_va$dgu_prob, by = "gene_name"),
           merge(x = merge(x = dgu$dgu_vb$gu[dgu$dgu_vb$gu$condition == "CAC",], 
                           y = dgu$dgu_vb$gu[dgu$dgu_vb$gu$condition == "CAC+Hp",], by = "gene_name"), 
                 y = dgu$dgu_vb$dgu_prob, by = "gene_name"))


igdata <- rbind(dgu$dgu_va$ppc$ppc_rep, dgu$dgu_vb$ppc$ppc_rep, dgu$dgu_ja$ppc$ppc_rep, dgu$dgu_jb$ppc$ppc_rep)
igdata <- igdata[igdata$gene_name %in% x$gene_name[abs(x$es_prob_mean)>delta & x$pmax>=pmax], ]
igdata$group <- substr(igdata$gene_name, start = 1, stop = 4)
igdata <- igdata[igdata$group %in% c("TRAV", "TRBV", "TRAJ", "TRBJ"),]

fig_dags <- ggplot(data = igdata)+
  facet_grid(~group, scales = "free_x", space = "free_x")+
  geom_sina(aes(x = gene_name, y = observed_prop, col = condition))+
  theme_bw(base_size = 10)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "top")+
  scale_color_manual(name = "Group", values = c("CAC"="black", "CAC+Hp"="darkgray"))+
  ylab(label = "Relative gene usage")+
  xlab(label = "")


fig_manuscript <- (fig_usage/fig_dags)+patchwork::plot_layout(heights = 1, 1)

ggsave(filename = "output/Fig_usage.pdf",
       plot = fig_manuscript,
       device = "pdf",
       width = 6,
       height = 7.5)

ggsave(filename = "output/SFig_usage_J.pdf",
       plot = fig_usage_supplement,
       device = "pdf",
       width = 6.5,
       height = 4)




#### PPC ####

dgu <- get(load("output/dgu_data.RData"))

sfig_ppc <- 
  ((ggplot(data = dgu$dgu_va$ppc$ppc_rep)+
      facet_wrap(facets = ~sample_id)+
      geom_errorbar(aes(x = observed_count, y = ppc_mean_count, ymin = ppc_L_count, ymax = ppc_H_count), col = "gray")+
      geom_point(aes(x = observed_count, y = ppc_mean_count), size = 0.5)+
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", col = "darkgray")+
      scale_x_log10("Observed usage [count]",
                    breaks = 10^(0:6),
                    labels = trans_format("log10", math_format(10^.x)))+
      scale_y_log10("Predicted usage [count]",
                    breaks = 10^(0:6),
                    labels = trans_format("log10", math_format(10^.x)))+
      annotation_logticks(base = 10, sides = "bl")+
      theme_bw(base_size = 10)+
      theme(strip.text = element_text(size = 9, margin = margin()))+
      ggtitle(label = "TRAV"))|
     (ggplot(data = dgu$dgu_ja$ppc$ppc_rep)+
        facet_wrap(facets = ~sample_id)+
        geom_errorbar(aes(x = observed_count, y = ppc_mean_count, ymin = ppc_L_count, ymax = ppc_H_count), col = "gray")+
        geom_point(aes(x = observed_count, y = ppc_mean_count), size = 0.5)+
        geom_abline(slope = 1, intercept = 0, linetype = "dashed", col = "darkgray")+
        scale_x_log10("Observed usage [count]",
                      breaks = 10^(0:6),
                      labels = trans_format("log10", math_format(10^.x)))+
        scale_y_log10("Predicted usage [count]",
                      breaks = 10^(0:6),
                      labels = trans_format("log10", math_format(10^.x)))+
        annotation_logticks(base = 10, sides = "bl")+
        theme_bw(base_size = 10)+
        theme(strip.text = element_text(size = 9, margin = margin()))+
        ggtitle(label = "TRAJ")))/
  ((ggplot(data = dgu$dgu_vb$ppc$ppc_rep)+
      facet_wrap(facets = ~sample_id)+
      geom_errorbar(aes(x = observed_count, y = ppc_mean_count, ymin = ppc_L_count, ymax = ppc_H_count), col = "gray")+
      geom_point(aes(x = observed_count, y = ppc_mean_count), size = 0.5)+
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", col = "darkgray")+
      scale_x_log10("Observed usage [count]",
                    breaks = 10^(0:6),
                    labels = trans_format("log10", math_format(10^.x)))+
      scale_y_log10("Predicted usage [count]",
                    breaks = 10^(0:6),
                    labels = trans_format("log10", math_format(10^.x)))+
      annotation_logticks(base = 10, sides = "bl")+
      theme_bw(base_size = 10)+
      theme(strip.text = element_text(size = 9, margin = margin()))+
      ggtitle(label = "TRBV"))|
     (ggplot(data = dgu$dgu_jb$ppc$ppc_rep)+
        facet_wrap(facets = ~sample_id)+
        geom_errorbar(aes(x = observed_count, y = ppc_mean_count, ymin = ppc_L_count, ymax = ppc_H_count), col = "gray")+
        geom_point(aes(x = observed_count, y = ppc_mean_count), size = 0.5)+
        geom_abline(slope = 1, intercept = 0, linetype = "dashed", col = "darkgray")+
        scale_x_log10("Observed usage [count]",
                      breaks = 10^(0:6),
                      labels = trans_format("log10", math_format(10^.x)))+
        scale_y_log10("Predicted usage [count]",
                      breaks = 10^(0:6),
                      labels = trans_format("log10", math_format(10^.x)))+
        annotation_logticks(base = 10, sides = "bl")+
        theme_bw(base_size = 10)+
        theme(strip.text = element_text(size = 9, margin = margin()))+
        ggtitle(label = "TRBJ")))
sfig_ppc  


ggsave(filename = "output/SFig_PPC_IgGeneUsage.pdf",
       plot = sfig_ppc,
       device = "pdf",
       width = 8.6,
       height = 8)
