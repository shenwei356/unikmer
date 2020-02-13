#!/usr/bin/env Rscript

library(scales)
library(ggplot2)
library(dplyr)
library(tidyr)
library(swr)
library(png)
library(cowplot)
library(ggthemes)

theme1 <- theme(
    text = element_text(size = 9),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 13),
    strip.background = element_rect(
        colour = "grey80",
        fill = "grey95",
        size = 0.5
    ),
    legend.position = "bottom",
    legend.spacing.x = unit(0.4, "cm"),
)

df <- read.csv("table.r.tsv", sep = "\t")
df$group <- factor(df$group, levels = unique(df$group), ordered = TRUE)
df$num <- factor(df$num, levels = unique(df$num), ordered = TRUE)

p <- ggplot(df, aes(x = k, y = value, color = group)) +
  # geom_line(aes(y=mean), color="#3e52a1") +
  geom_smooth(method = loess,
              se = FALSE,
              alpha = 0.2 ) +
  geom_point(aes(shape = group), size = 2) +
  scale_x_continuous(limits=c(10, 32), breaks = seq(1, 31, by = 4))+
  geom_vline(xintercept = seq(21,27,2), colour="grey80", linetype=2, size=0.2) + 
  scale_color_colorblind() +
  facet_grid(. ~ num, scales = "free_y") +
  xlab("K") +
  ylab("compression rate (%)") +
  ggtitle(NULL) +
  shenwei356.theme() +
  theme1

ggsave(
  p,
  file = "cr.jpg",
  width = 10,
  height = 4.5,
  dpi = 300
)

