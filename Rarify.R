library(tidyverse)
library(vegan)

df <- read_csv('wf-metagenomics-counts-species.csv')%>%
  select(-MEY15_cntrl)


#df_prelim <- read_csv('wf-metagenomics-counts-species_prelim.csv')%>%
#  select(-MEY7_8_control)
#df2 <- merge(df, df_prelim, by = 'species', all = T)


df <- df%>%
  pivot_longer(
    cols = -c(species, total, superkingdom, kingdom, phylum, class, order,
              family, genus, tax),
    names_to = "SampleID",
    values_to = "value"
  )%>% # We will filter out unknown species
  filter(species != 'Unknown')

min_n_seqs <- df%>%
  group_by(SampleID)%>%
  summarize(n_seqs = sum(value))%>%
  summarize(min = min(n_seqs))%>%
  pull(min)

df2 <- df %>% 
  pivot_wider(names_from = "species", values_from = "value")

df2 <- df2%>%
  as.data.frame()

rownames(df2) <- df2$Name
df2 <- df2[,-1]
