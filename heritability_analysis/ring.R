library(ggplot2)
library(broom)
library(tidyverse)
# Create test dat.alla
dat.all = data.frame(count=c(32.83, 8.90, 13.19, 45.07,
                             32.93, 8.88, 13.17, 45.02),
                 ring=c("A", "A","A","A", "B","B","B", "B"),
                 Category=c("Promoters","Exons", "Introns", "Intergenic regions",
                            "Promoters","Exons", "Introns", "Intergenic regions"))
dat.all$Category=factor(dat.all$Category, levels = c("Promoters","Exons", "Introns", "Intergenic regions"), ordered = T)

# compute fractions
#dat.all = dat.all[order(dat.all$count), ]
dat.all %<>% group_by(ring) %>% mutate(fraction = count / sum(count),
                                   ymax = cumsum(fraction),
                                   ymin = c(0,ymax[1:length(ymax)-1]))


# Add x limits
baseNum <- 4
#numCat <- length(unique(dat.all$ring))
dat.all$xmax <- as.numeric(dat.all$ring) + baseNum
dat.all$xmin = dat.all$xmax -1


# plot
p2.all = ggplot(dat.all, aes(fill=Category,
                     # alpha = ring,
                     ymax=ymax, 
                     ymin=ymin, 
                     xmax=xmax, 
                     xmin=xmin)) +
  geom_rect() +
  #geom_rect(colour="grey30") +
  coord_polar(theta="y") +
  xlim(c(0, 6)) +
  theme_bw() +
  theme(panel.grid=element_blank()) +
  theme(axis.text=element_blank()) +
  theme(axis.ticks=element_blank(),
        panel.border = element_blank()) +
  # labs(title="Customized ring plot") + 
  scale_fill_brewer(palette = "Set1") +
  scale_alpha_discrete(range = c(1,1))

p2.all

dat.hyper = data.frame(count=c(32.83, 8.90, 13.19, 45.07,
                               11.92, 27.62, 18.86, 41.60),
                       ring=c("A", "A","A","A", "B","B","B", "B"),
                       Category=c("Promoters","Exons", "Introns", "Intergenic regions",
                                  "Promoters","Exons", "Introns", "Intergenic regions"))
dat.hyper$Category=factor(dat.hyper$Category, levels = c("Promoters","Exons", "Introns", "Intergenic regions"), ordered = T)

# compute fractions
#dat.hyper = dat.hyper[order(dat.hyper$count), ]
dat.hyper %<>% group_by(ring) %>% mutate(fraction = count / sum(count),
                                         ymax = cumsum(fraction),
                                         ymin = c(0,ymax[1:length(ymax)-1]))


# Add x limits
baseNum <- 4
#numCat <- length(unique(dat.hyper$ring))
dat.hyper$xmax <- as.numeric(dat.hyper$ring) + baseNum
dat.hyper$xmin = dat.hyper$xmax -1


# plot
p2.hyper = ggplot(dat.hyper, aes(fill=Category,
                                 # alpha = ring,
                                 ymax=ymax, 
                                 ymin=ymin, 
                                 xmax=xmax, 
                                 xmin=xmin)) +
  geom_rect() +
  #geom_rect(colour="grey30") +
  coord_polar(theta="y") +
  xlim(c(0, 6)) +
  theme_bw() +
  theme(panel.grid=element_blank()) +
  theme(axis.text=element_blank()) +
  theme(axis.ticks=element_blank(),
        panel.border = element_blank()) +
  # labs(title="Customized ring plot") + 
  scale_fill_brewer(palette = "Set1") +
  scale_alpha_discrete(range = c(1,1))

p2.hyper

dat.hypo = data.frame(count=c(32.83, 8.90, 13.19, 45.07,
                              53.37, 2.39, 4.43, 39.81),
                      ring=c("A", "A","A","A", "B","B","B", "B"),
                      Category=c("Promoters","Exons", "Introns", "Intergenic regions",
                                 "Promoters","Exons", "Introns", "Intergenic regions"))
dat.hypo$Category=factor(dat.hypo$Category, levels = c("Promoters","Exons", "Introns", "Intergenic regions"), ordered = T)

# compute fractions
#dat.hypo = dat.hypo[order(dat.hypo$count), ]
dat.hypo %<>% group_by(ring) %>% mutate(fraction = count / sum(count),
                                        ymax = cumsum(fraction),
                                        ymin = c(0,ymax[1:length(ymax)-1]))


# Add x limits
baseNum <- 4
#numCat <- length(unique(dat.hypo$ring))
dat.hypo$xmax <- as.numeric(dat.hypo$ring) + baseNum
dat.hypo$xmin = dat.hypo$xmax -1


# plot
p2.hypo = ggplot(dat.hypo, aes(fill=Category,
                               # alpha = ring,
                               ymax=ymax, 
                               ymin=ymin, 
                               xmax=xmax, 
                               xmin=xmin)) +
  geom_rect() +
  #geom_rect(colour="grey30") +
  coord_polar(theta="y") +
  xlim(c(0, 6)) +
  theme_bw() +
  theme(panel.grid=element_blank()) +
  theme(axis.text=element_blank()) +
  theme(axis.ticks=element_blank(),
        panel.border = element_blank()) +
  # labs(title="Customized ring plot") + 
  scale_fill_brewer(palette = "Set1") +
  scale_alpha_discrete(range = c(1,1))

p2.hypo

dat.var = data.frame(count=c(32.83, 8.90, 13.19, 45.07,
                             9.03, 12.19, 24.56, 54.22),
                     ring=c("A", "A","A","A", "B","B","B", "B"),
                     Category=c("Promoters","Exons", "Introns", "Intergenic regions",
                                "Promoters","Exons", "Introns", "Intergenic regions"))
dat.var$Category=factor(dat.var$Category, levels = c("Promoters","Exons", "Introns", "Intergenic regions"), ordered = T)

# compute fractions
#dat.var = dat.var[order(dat.var$count), ]
dat.var %<>% group_by(ring) %>% mutate(fraction = count / sum(count),
                                       ymax = cumsum(fraction),
                                       ymin = c(0,ymax[1:length(ymax)-1]))


# Add x limits
baseNum <- 4
#numCat <- length(unique(dat.var$ring))
dat.var$xmax <- as.numeric(dat.var$ring) + baseNum
dat.var$xmin = dat.var$xmax -1


# plot
p2.var = ggplot(dat.var, aes(fill=Category,
                             # alpha = ring,
                             ymax=ymax, 
                             ymin=ymin, 
                             xmax=xmax, 
                             xmin=xmin)) +
  geom_rect() +
  #geom_rect(colour="grey30") +
  coord_polar(theta="y") +
  xlim(c(0, 6)) +
  theme_bw() +
  theme(panel.grid=element_blank()) +
  theme(axis.text=element_blank()) +
  theme(axis.ticks=element_blank(),
        panel.border = element_blank()) +
  # labs(title="Customized ring plot") + 
  scale_fill_brewer(palette = "Set1") +
  scale_alpha_discrete(range = c(1,1))

p2.var

# Promoters, exons, introns, intergenic regions: 
# Null distribution: 32.83, 8.90, 13.19, 45.07
# all: 32.93, 8.88, 13.17, 45.02 (p-value: 0.983, 0.9944, 0.9953, 0.992)
# hyper: 11.92, 27.62, 18.86, 41.60 (p-value: 1.192e-06 less, 6.329e-08 more, 0.1119, 0.4844)
# hypo: 53.37, 2.39, 4.43, 39.81 (p-value: 2.419e-05 more, 0.007331 less, 0.003167 less, 0.2883)
# variable: 9.03, 12.19, 24.56, 54.22 (p-value: 1.647e-08 less, 0.2715, 0.002224 more, 0.0668)

library(cowplot)
ring=plot_grid(p2.all, p2.hyper, p2.hypo, p2.var,
               labels = c("a", "b", "c", "d"),
               align = "vh")
ring
