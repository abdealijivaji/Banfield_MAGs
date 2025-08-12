rm(list = ls())

library(ggplot2)

setwd("Banfield_MAGs/")

cog_cat <- read.table("banfield_MAGs_curated.annotations.tsv", header = T, sep = "\t")

cog_sub <- cog_cat[,c(1,7)]
cog_sub <- cog_sub[!(cog_sub$COG_category == "-"), ]
#unique(cog_sub$COG_category)

cog_sub$MAG <- stringr::word(cog_sub$query, end = -2, sep = "_")

cog_sub <- data.frame(tidyr::separate_longer_position(cog_sub, COG_category, width = 1))

dict <- read.table("cog_dict.tsv", sep = "\t", header = T)


cog_sub$cog <- matchmaker::match_vec(cog_sub$COG_category, 
                      dictionary = dict)

  

bub <- ggplot(cog_sub, aes(MAG, cog)) + 
  #geom_count(col = "#F002F0") + 
  geom_count(aes(col = cog )) +
  labs(x= "MAG", y = "COG category") +
  theme_bw() +
  guides(col = "none") +
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1)) 
bub

ggsave("images/bubble_plot.svg", bub, device = "svg", width = 10, height = 6)

