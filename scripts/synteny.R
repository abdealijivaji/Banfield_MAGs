rm(list = ls())

library(gggenomes)
library(tidyverse)

seqs <- read_fai("data/synteny/cat_fasta_index")
seqs$bin_id <- c("GD2017_2_urea", "SSJD_D22_6_30_4", "SSJD_D22_6_30_4", "bjp_ig5181", 
                 "Faustovirus_ST1", "Pacmanvirus_A23", "DRTY6", "Pithovirus_LCDPAC01", 
                 "Pithovirus_LCDPAC01", "Pithovirus_LCDPAC01", "Pithovirus_LCDPAC01", 
                 "Pithovirus_LCDPAC01", "Pithovirus_LCDPAC01", "Pithovirus_LCDPAC01", 
                 "Pithovirus_LCDPAC01", "Pithovirus_LCDPAC01", "Pithovirus_LCDPAC01", 
                 "Pithovirus_LCDPAC01", "2019_SCN_bioreactor", "Hydrivirus", 
                 "SRVP_07252022_550855_L", "GD2017_2_strous", "Tupanvirus_soda_lake", 
                 "ALT_072018_38_16", "Organic_Lake_phycodnavirus_1", "AB_072018", 
                 "BML_coassembly", "Chrysochromulina_ericina_virus")

grp_list <- read.table("data/synteny/in_list", sep = "\t", col.names = c("group", "bin_id"))
seqs <- seqs %>% dplyr::left_join(grp_list, by = "bin_id")

#genes <- read_subfeats("data/synteny/cat_proteins.faa")
links <- read_paf("data/synteny/cat_minimap.paf", max_tags = 32)



p1 <- gggenomes(seqs = seqs, links = links) + 
  geom_seq() + geom_bin_label() +
  geom_link()
p1 %>%  flip_seqs(Tupanvirus_soda_lake)


# Just focusing on group 5 for now. Adding TIR annotation

tir_gp5_paf <- read_paf("data/synteny/grp_5_TIR_minimap.paf") %>%
  dplyr::filter(seq_id == seq_id2 & start < start2 & map_length > 99 & de < 0.1)
tir_gp5 <- dplyr::bind_rows(
  dplyr::select(tir_gp5_paf, seq_id=seq_id, start=start, end=end, de),
  dplyr::select(tir_gp5_paf, seq_id=seq_id2, start=start2, end=end2, de))

gp5 <- gggenomes(seqs = seqs[seqs$group == "grp_5",], links = links, feats = tir_gp5) + 
  geom_seq() + geom_bin_label() +
  geom_feat(colour = "darkred",  linewidth = 15, position = "identity") +
  geom_link()
gp5 %>%  flip_seqs(Tupanvirus_soda_lake)
