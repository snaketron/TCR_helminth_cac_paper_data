# Clonotype abundance diversity: Shannon's entropy + Normalized Rank Abundance Distributions

require(ggplot2)
require(reshape2)
require(patchwork)
require(MASS)
require(RADanalysis)
require(ggforce)

# get the data
mixcr_clones <- get(load("data/mixcr_clones.RData"))
mixcr_clones$sample_id <- paste0(mixcr_clones$experiment, "|", mixcr_clones$sample, "|", mixcr_clones$locus)

# get the V and J genes, deine clonotype: V+J+CDR3
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
meta <- mixcr_clones[,c("sample_id", "locus", "group")]
meta <- meta[duplicated(meta)==F, ]


#### Shannon entropy ####
hn <- function(x, d) {
  r <- d$readCount[d$sample_id == x]
  p <- r/sum(r)
  h <- -sum(p*log(p))
  return(h/log(length(p)))
  # return(h/length(p))
}

h <- function(x, d) {
  r <- d$readCount[d$sample_id == x]
  p <- r/sum(r)
  return(-sum(p*log(p)))
}
hs <- vapply(X = unique(mixcr_clones$sample_id), d = mixcr_clones, 
             FUN.VALUE = numeric(length = 1) , FUN = h)
hs <- data.frame(sample_id = names(hs), h = as.numeric(hs))
hs <- merge(x = hs, y = meta, by = "sample_id")
fig_h <- ggplot(data = hs)+
  geom_sina(aes(x = locus, y = h, col = group), size = 1)+
  theme_bw()+
  ylab(label = "Shannon's entropy (H)")+
  xlab(label = "Group")+
  theme(legend.position = "top")+
  scale_color_manual(name = "Group", 
                     values = c("CAC"="orange", "CAC+Hp"="steelblue"))

hs <- vapply(X = unique(mixcr_clones$sample_id), d = mixcr_clones, 
             FUN.VALUE = numeric(length = 1) , FUN = hn)
hs <- data.frame(sample_id = names(hs), h = as.numeric(hs))
hs <- merge(x = hs, y = meta, by = "sample_id")
fig_h <- ggplot(data = hs)+
  geom_sina(aes(x = locus, y = h, col = group), size = 1)+
  theme_bw()+
  ylab(label = "Shannon's entropy (H)")+
  xlab(label = "Group")+
  theme(legend.position = "top")+
  scale_color_manual(name = "Group", 
                     values = c("CAC"="orange", "CAC+Hp"="steelblue"))


wilcox.test(hs$h[hs$locus=="TRB" & hs$group=="CAC"], hs$h[hs$locus=="TRB" & hs$group!="CAC"])
wilcox.test(hs$h[hs$locus=="TRA" & hs$group=="CAC"], hs$h[hs$locus=="TRA" & hs$group!="CAC"])



#### NRADs ####

get_rads <- function(d) {
  d <- d[order(d$readCount, decreasing = T), ]
  o <- c()
  for(s in unique(d$sample_id)) {
    w <- d[d$sample_id == s, ]
    w$rank <- 1:nrow(w)
    o <- rbind(o, w)
  }
  return(o)
}

# get the RADs
mixcr_clones <- get_rads(d = mixcr_clones)
meta <- mixcr_clones[duplicated(mixcr_clones[,c("sample_id", "group", "locus")])==F,
                     c("sample_id", "group", "locus")]

y <- acast(mixcr_clones, formula = sample_id~rank, fill = 0, value.var = "readCount")

# that is the minimum rank
r <- apply(X = y, MARGIN = 1, FUN = function(x){sum(x!=0)})
min(r)


# Normalize the RADs -> NRADs
nrads <- RADnormalization_matrix(input = y, max_rank = 2875, average_over = 1000)
n <- melt(data = nrads$norm_matrix)
colnames(n) <- c("i", "rank", "nrad")
n$sample_id <- rownames(y)[n$i]
n <- merge(x = n, y = meta, by = "sample_id", all.x = TRUE)

# Plot
fig_nrads <- ggplot(data = n)+
  facet_wrap(facet = ~locus, scales = "free")+
  geom_line(aes(x=rank, y=nrad, col = group, group = sample_id))+
  scale_x_log10(name = "Rank")+
  scale_y_log10(name = "Normalized abundance")+
  theme_bw()+
  theme(legend.position = "top")+
  annotation_logticks(base = 10, sides = "lb")+
  scale_color_manual(name = "Group", 
                     values = c("CAC"="orange", "CAC+Hp"="steelblue"))


# Patch together
fig_diversity <- (fig_nrads|fig_h)+
  patchwork::plot_layout(widths = c(1,0.5))+
  patchwork::plot_annotation(tag_levels = "A")
fig_diversity

ggsave(filename = "outpu/diversity.pdf",
       plot = fig_diversity,
       device = "pdf",
       width = 7,
       height = 3)


