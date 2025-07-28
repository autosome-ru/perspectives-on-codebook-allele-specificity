library(ggseqlogo)
library(ggplot2)

args = commandArgs(TRUE)
filename = args[1]
out = args[2]

pfm <- read.table(filename, sep = "\t", stringsAsFactors = FALSE,
                   header = FALSE, skip = 3,row.names = 1)
colnames(pfm) <- c("A", "C", "G", "T")
pfm <- t(pfm)

colscheme <- make_col_scheme(chars = c("A", "C", "G", "T"), 
                             cols=c("#00CC00", "#0000CC", "#FFB302", "#CC0001"),
                             name="DNAalph")
plot <- ggplot() + 
  geom_logo(as.matrix(pfm), 
            namespace = colscheme$letter, 
            col_scheme = colscheme, 
            font="roboto_bold",
            method="prob") + 
  ylim(c(0, 1)) + 
  theme_logo() +
  theme(axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "null"),
        panel.spacing = unit(c(0, 0, 0, 0), "null"),
        panel.background = element_rect(fill = "transparent", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        legend.position = "none",
        axis.ticks.length = unit(0, "null"),
        axis.ticks.margin = unit(0, "null")) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,1)) 

png(paste0(out, ".pfm.png"), width=(1.5*ncol(pfm)/10), height=0.75, units='in', res=1000)
print(plot)
dev.off()

rev_pfm <- pfm[, ncol(pfm):1]
colnames(rev_pfm) <- 1:ncol(rev_pfm)
rev_pfm <- rev_pfm[c("T", "G", "C", "A"), ]
rownames(rev_pfm) <- c("A", "C", "G", "T")

plot <- ggplot() + 
  geom_logo(as.matrix(rev_pfm), 
            namespace = colscheme$letter, 
            col_scheme = colscheme, 
            font="roboto_bold",
            method="prob") + 
  ylim(c(0, 1)) + 
  theme_logo() +
  theme(axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "null"),
        panel.spacing = unit(c(0, 0, 0, 0), "null"),
        panel.background = element_rect(fill = "transparent", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        legend.position = "none",
        axis.ticks.length = unit(0, "null"),
        axis.ticks.margin = unit(0, "null")) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,1)) 

png(paste0(out, ".pfm.revcomp.png"), width=(1.5*ncol(pfm)/10), height=0.75, units='in', res=1000)
print(plot)
dev.off()