library(tidyverse)
library(edgeR)

df <- read_csv('wf-metagenomics-counts-genus.csv')%>%
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
df <- df %>% filter(grepl("^MEO", Group))

df2 <- as_tibble(df$Group)
df2 <- df2 %>% rename(SampleID = value)
df2 <- df2 %>%
  mutate(SampleID = ifelse(str_sub(SampleID, -2) == "_7", str_replace(SampleID, "_7$", "_6"), SampleID))
df2 <- df2%>%
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
metadata <- read_csv("Metadata.csv")
df2 <-   left_join(df2, metadata, by = 'SampleID')
# Adding age treatment column
columns_to_convert <- setdiff(names(df2), c("SampleID", "Treatment", "Sex", "Timepoint"))
df2 <- df2 %>%
  mutate_at(vars(columns_to_convert), as.numeric)%>%
  mutate(Age = ifelse(substr(SampleID, 3, 3) == "O", "old", 
                      ifelse(substr(SampleID, 3, 3) == "Y", "young", NA)))

# df2 is our design dataframe, df is our count dataframe
Group <- factor(paste(df2$Treatment, df2$Timepoint, df2$Age, sep="."))
Age <- factor(df2$Age)
Treatment <- factor(df2$Treatment)
Timepoint <- factor(df2$Timepoint)
cbind(df2, Group=Group, Treatment=Treatment, Timepoint=Timepoint, Age=Age)

design <- model.matrix(~0+Group)
colnames(design) <- levels(Group)

# Need count matrix of species with Samples as column names, and Species as rownames
df3 <- read_csv('wf-metagenomics-counts-genus.csv')%>%
  select(-MEY15_cntrl)%>%
  select(-c(total, superkingdom, kingdom, phylum, class,
            order,family, tax))%>%
  rename_with(~str_replace(., "_7$", "_6"), ends_with("_7"))%>%
  as.data.frame()
row_names <- df3[, 1]  # Extract the values from the first column
df3 <- df3[, -1]  # Remove the first column from the data frame
rownames(df3) <- row_names
df3 <- select(df3, starts_with("MEO"))
df3 <- as.matrix(df3)

# Young only

d <- DGEList(counts = df3, group = Group)
keep <- filterByExpr(d, group=Group)
d <- d[keep, , keep.lib.sizes=FALSE]
d <- estimateDisp(d)
d <- calcNormFactors(d)


plotMDS(d, top = 1000, col = as.numeric(d$samples$group),
        main = "PCoA of genera between groups for old mice", pch=16)
legend("bottomright", col = as.numeric(d$samples$group), pch=c(15), legend=c("Exercise week 0","Exercise week 6",
                                                                             "Sedentary week 0", "Sedentary week 6"))
