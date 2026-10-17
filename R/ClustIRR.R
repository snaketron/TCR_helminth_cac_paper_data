# ClustIRR analysis

#### Clustering ####


# make sure blast is installed via conda (see ClustIRR documentation)
Sys.setenv(PATH = paste(Sys.getenv("PATH"), "/miniconda3/bin/", sep=.Platform$path.sep))
require(ggplot2)
require(reshape2)
require(patchwork)
require(ggforce)
require(ClustIRR)

# load data
mixcr_clones <- get(load("data/mixcr_clones.RData"))
mixcr_clones$sample_id <- paste0(mixcr_clones$experiment, "|", mixcr_clones$sample, "|", mixcr_clones$locus)

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
mixcr_clones$clonotype_id <- paste0(mixcr_clones$v, "+", mixcr_clones$j, "+", mixcr_clones$aaSeqCDR3)
mixcr_clones <- aggregate(readCount~clonotype_id+sample_id+group+aaSeqCDR3+v+j+locus, 
                          data = mixcr_clones, FUN = sum)

# split data by chain
tra <- mixcr_clones[mixcr_clones$locus == "TRA",]
trb <- mixcr_clones[mixcr_clones$locus == "TRB",]

# do clustirr analysis separately
# 1. TRB
cb <- clustirr(s = data.frame(CDR3b = trb$aaSeqCDR3, sample = trb$sample_id, clone_size = trb$readCount),
               meta = trb, control = list(blast_gmi = 0.9, blast_cores = 40))
save(cb, file = "manuscript_CD8/clustirr/cb.RData", compress = TRUE)

gb <- detect_communities(graph = cb$graph, algorithm = "leiden", metric = "average", 
                         resolution = 1, iterations = 1000, chains = "CDR3b")
save(gb, file = "manuscript_CD8/clustirr/gb.RData", compress = TRUE)


# 2. TRA
ca <- clustirr(s = data.frame(CDR3a = tra$aaSeqCDR3, sample = tra$sample_id, clone_size = tra$readCount),
               meta = tra, control = list(blast_gmi = 0.9, blast_cores = 40))
save(ca, file = "manuscript_CD8/clustirr/ca.RData", compress = TRUE)

ga <- detect_communities(graph = ca$graph, algorithm = "leiden", metric = "average", 
                         resolution = 1, iterations = 1000, chains = "CDR3a")
save(ga, file = "manuscript_CD8/clustirr/ga.RData", compress = TRUE)



#### Modeling ####
require(ggplot2)
require(reshape2)
require(patchwork)
require(ggforce)
require(ClustIRR)


# DCO of CDR3a
ga <- get(load("output/ga.RData"))
meta <- tra[,c("sample_id", "group")]
meta <- meta[duplicated(meta)==FALSE,]
meta$group_num <- ifelse(test = meta$group=="CAC", yes = 1, no = 2)
com <- ga$community_occupancy_matrix
com <- com[,meta$sample_id]
colnames(com)==meta$sample_id
rm(ga, mixcr_clones);gc();gc();gc()


# slimify com: if 0 in 8 out of 9 samples -> reduce
i <- which(apply(com, MARGIN = 1, FUN = function(x) {
  sum(x==0)
})>=8)
z <- com[i,]
z <- colSums(z)
nz <- com[-i,]
nz <- rbind(nz, z)


# DCO
dcoa_h <- dco(community_occupancy_matrix = nz,
              groups = meta$group_num,
              mcmc_control = list(mcmc_warmup = 500,
                                  mcmc_iter = 1750,
                                  mcmc_chains = 3,
                                  mcmc_cores = 3,
                                  adapt_delta = 0.9,
                                  max_treedepth = 10))

save(dcoa_h, file = "output/dcoa_h.RData", compress = TRUE)
ps <- dcoa_h$posterior_summary
save(ps, file = "output/psa_h.RData", compress = TRUE)






# DCO of CDR3b
gb <- get(load("output/gb.RData"))

meta <- trb[,c("sample_id", "group")]
meta <- meta[duplicated(meta)==FALSE,]
meta$group_num <- ifelse(test = meta$group=="CAC", yes = 1, no = 2)
com <- gb$community_occupancy_matrix
com <- com[,meta$sample_id]
colnames(com)==meta$sample_id
rm(gb, mixcr_clones);gc();gc();gc()


# slimify com
i <- which(apply(com, MARGIN = 1, FUN = function(x) {
  sum(x==0)
})>=8)
z <- com[i,]
z <- colSums(z)
nz <- com[-i,]
nz <- rbind(nz, z)

# DCO
dcob_h <- dco(community_occupancy_matrix = nz,
              groups = meta$group_num,
              mcmc_control = list(mcmc_warmup = 500,
                                  mcmc_iter = 1750,
                                  mcmc_chains = 3,
                                  mcmc_cores = 3,
                                  adapt_delta = 0.9,
                                  max_treedepth = 10))
save(dcob_h, file = "output/dcob_h.RData", compress = TRUE)
ps <- dcob_h$posterior_summary
save(ps, file = "output/psb_h.RData", compress = TRUE)






#### Figures ####
require(ggplot2)
require(ggforce)
require(ClustIRR)
require(patchwork)

# minimum probability supporting differentially 
# abundant communities on the joint graph
pmax <- 0.99



# TRB analysis
ps <- get(load("output/psb_h.RData"))
ep <- ps$epsilon
ep <- ep[ep$contrast == "1-2",]
which(ep$pmax>=0.99)
# none found at this stage


# volcano plot
g_volcano <- ggplot(data = ep)+
  geom_hline(yintercept = 0.99, linetype = "dashed")+
  geom_point(aes(x = mean*100, y = pmax), shape = 21, 
             col = "white", fill = "black", stroke = 0.1)+
  theme_bw(base_size = 10)+
  xlab(label = expression(epsilon*" [%]"))+
  ylab(label = expression(pi))


# save volcano plot
ggsave(filename = "output/CDR3b_ClustIRR_volcano.pdf",
       plot = g_volcano,
       device = "pdf",
       width = 4,
       height = 4)



# TRA analysis
ps <- get(load("output/psa_h.RData"))
ep <- ps$epsilon
ep <- ep[ep$contrast == "1-2",]
which(ep$pmax>=0.99)
# hits!

# volcano plot
g_volcano <- ggplot(data = ep)+
  geom_hline(yintercept = 0.99, linetype = "dashed")+
  geom_point(aes(x = mean*100, y = pmax), shape = 21, 
             col = "white", fill = "black", stroke = 0.1)+
  theme_bw(base_size = 10)+
  xlab(label = expression(epsilon*" [%]"))+
  ylab(label = expression(pi))


ggsave(filename = "output/CDR3a_ClustIRR_volcano.pdf",
       plot = g_volcano,
       device = "pdf",
       width = 4,
       height = 4)

# CDR3 motif logos (only for TCR alpha chain)
ga <- get(load("output/ga.RData"))
meta <- ga$node_summary[duplicated(ga$node_summary$sample)==F, c("sample_id", "group")]

# simplify communities: if found in only one patient -> group under singletons
com <- ga$community_occupancy_matrix
i <- which(apply(com, MARGIN = 1, FUN = function(x) {
  sum(x)==1
}))
z <- com[i,]
z <- colSums(z)
nz <- com[-i,]
nz <- rbind(nz, z)
rm(z,i)

# Fig 1/2
i <- reshape2::melt(nz[unique(ep$community[which(ep$mean > 0 & ep$pmax>=pmax)]),]/colSums(nz))
i <- merge(x = i, y = meta, by.x = "Var2", by.y = "sample_id", all.x = TRUE)
j <- reshape2::melt(nz[unique(ep$community[which(ep$mean < 0 & ep$pmax>=pmax)]),]/colSums(nz))
j <- merge(x = j, y = meta, by.x = "Var2", by.y = "sample_id", all.x = TRUE)

g1 <- ggplot(data = i)+
  facet_wrap(~group, nrow = 1, scales = "free_x")+
  geom_point(aes(x = Var2, y = value*100, group = as.factor(Var1)),
             position = position_dodge(width = 0.5))+
  scale_y_log10("Freq [%]")+
  theme_bw(base_size = 10)
g2 <- ggplot(data = j)+
  facet_wrap(~group, nrow = 1, scales = "free_x")+
  geom_point(aes(x = Var2, y = value*100, group = as.factor(Var1)),
             position = position_dodge(width = 0.5))+
  scale_y_log10("Freq [%]")+
  theme_bw(base_size = 10)


com_i <- as.numeric(rownames(nz[unique(ep$community[which(ep$mean > 0 & ep$pmax>=pmax)]),]))
com_j <- as.numeric(rownames(nz[unique(ep$community[which(ep$mean < 0 & ep$pmax>=pmax)]),]))

w <- ga$node_summary[ga$node_summary$com %in% com_i,]
w$community <- as.character(w$community)
cs <- as.numeric(unique(w$community))
for(x in 1:length(cs)) {
  ov <- ClustIRR::get_community_feature_stats(node_summary = w, node_feature = "v", community_id = cs[x])
  oj <- ClustIRR::get_community_feature_stats(node_summary = w, node_feature = "j", community_id = cs[x])
  cat(paste0(gsub(pattern = "TCRB", replacement = '', x = ov$feature[1]), " ", 
             round(ov$feature_prop[1]*100, digits = 1), "%, ",
             gsub(pattern = "TCRB", replacement = '', x = oj$feature[1]), " ", 
             round(oj$feature_prop[1]*100, digits = 1), "%"), "\n")
}


w <- ga$node_summary[ga$node_summary$com %in% com_j,]
w$community <- as.character(w$community)
cs <- as.numeric(unique(w$community))
for(x in 1:length(cs)) {
  ov <- ClustIRR::get_community_feature_stats(node_summary = w, node_feature = "v", community_id = cs[x])
  oj <- ClustIRR::get_community_feature_stats(node_summary = w, node_feature = "j", community_id = cs[x])
  cat(paste0(gsub(pattern = "TCRB", replacement = '', x = ov$feature[1]), " ", 
             round(ov$feature_prop[1]*100, digits = 1), "%, ",
             gsub(pattern = "TCRB", replacement = '', x = oj$feature[1]), " ", 
             round(oj$feature_prop[1]*100, digits = 1), "%"), "\n")
}


g <- (g1/wrap_plots(lapply(X = com_i, FUN = get_cdr3_motifs, node_summary = ga$node_summary), ncol = 5)/
        g2/wrap_plots(lapply(X = com_j, FUN = get_cdr3_motifs, node_summary = ga$node_summary), ncol = 5))+
  patchwork::plot_layout(heights = c(1,4,1,0.4))
g

ggsave(filename = "output/Fig_CDR3a_ClustIRR.pdf",
       plot = g,
       device = "pdf",
       width = 7,
       height = 10)

