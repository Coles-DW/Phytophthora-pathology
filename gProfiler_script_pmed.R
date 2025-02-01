##plot G:profiler result

##package
library(ggplot2)
library(tidyverse)
library(viridisLite)
library(viridis)
library(gridExtra)

#setwd

setwd("")

##input

gprofiler <- read.csv("Top 20 terms_122472_pmed.csv")

gprofiler <- gprofiler %>% 
  mutate_at(., vars(4:5), ~as.numeric(.)) %>% 
  mutate_at(., vars(source), ~as.factor(.)) %>% 
  arrange(intersection_size) %>%   
  arrange(source) %>%    
  #  mutate(term_name=factor(term_name, levels=term_name)) %>% 
  glimpse()

ggplot(gprofiler, aes(x = reorder(term_name, intersection_size), 
                      y = intersection_size,
                      color = adjusted_p_value)) +
  geom_segment(aes(xend=term_name, yend = 0), color = "darkgrey", size = 1.3) +
  geom_point(size = 4) +
  facet_grid(source~.,scales = "free",drop = T, space = "free") + #use this line when you have different types of GO terms (e.g. KEGG, GO:MF, GO:BP etc.)
  scale_color_viridis(begin = 0.95, end = 0.2)+
  coord_flip() +
  theme_bw() +
  theme(legend.position = "top", axis.text = element_text(size =13))

#Don't do the rest below
gprofiler.12 <- gprofiler %>% 
  filter(source =="12") %>% 
  glimpse()

g12 <- ggplot(gprofiler.12, aes(x = reorder(term_name, intersection_size), 
                                y = intersection_size,
                                color = adjusted_p_value)) +
  geom_segment(aes(xend=term_name, yend = 0), color = "darkgrey", size = 1.3) +
  geom_point(size = 4) +
  #  facet_grid(source~.,scales = "free",drop = T, space = "free") + #use this line when you have different types of GO terms (e.g. KEGG, GO:MF, GO:BP etc.)
  scale_color_viridis(begin = 0.95, end = 0.2)+
  ylim(c(0,60)) +
  coord_flip() +
  theme_bw() +
  theme(legend.position = "none", axis.text = element_text(size =13))

grid.arrange(all, g12, nrow = 2)

