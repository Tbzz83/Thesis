library (tidyverse)
library(vegan)

df <- read_csv('wf-metagenomics-counts-species.csv')%>%
  select(-c(MEY15_cntrl))%>%
  rename(name := !!names(.)[1])%>% #rename first col to "name"
  select(1:(which(names(.) == "total") - 1))%>%
  pivot_longer(!name, names_to = 'Group', values_to = 'value')%>%
  as_tibble()%>%
  group_by(Group)%>%
  mutate(total = sum(value))%>%
  #filter(total > 1800)
  group_by(name)%>%
  mutate(total = sum(value))%>% 
  filter(total !=0)%>%
  ungroup() %>%
  select(-total)%>%
  pivot_wider(Group, names_from = name)%>%
  as.data.frame()

rownames(df) <- df$Group
df <- df[, -1]
df <- as.matrix(df)

set.seed(19990430)
dist <- vegdist(df, method="bray")
nmds <- metaMDS(dist)

#cmdscale(df)

scores <- scores(nmds)%>%
  as_tibble(rownames = "SampleID")


scores <- scores%>%
  mutate(SampleID = ifelse(str_sub(SampleID, -2) == "_7", str_replace(SampleID, "_7$", "_6"), SampleID))%>%
  mutate(Timepoint = case_when(
    str_detect(SampleID, "_0") ~ "0",
    str_detect(SampleID, "_1") ~ "1",
    str_detect(SampleID, "_2") ~ "2",
    str_detect(SampleID, "_3") ~ "3",
    str_detect(SampleID, "_4") ~ "4",
    str_detect(SampleID, "_5") ~ "5",
    str_detect(SampleID, "_6") ~ "6",
    str_detect(SampleID, "_7") ~ "7",
    str_detect(SampleID, "_8") ~ "8",
    str_detect(SampleID, "_9") ~ "9"
  ), SampleID = str_replace(SampleID, "_[0-9]$", ""))
scores <- left_join(scores, metadata, by = 'SampleID')

# Adding age treatment column
columns_to_convert <- setdiff(names(scores), c("SampleID", "Treatment", "Sex", "Timepoint"))
scores <- scores %>%
  mutate_at(vars(columns_to_convert), as.numeric)%>%
  mutate(Age = ifelse(substr(SampleID, 3, 3) == "O", "Old", 
                      ifelse(substr(SampleID, 3, 3) == "Y", "young", NA)))

ggplot(data = scores, aes(x = NMDS1, y = NMDS2, color = Age ))+
  geom_point()

