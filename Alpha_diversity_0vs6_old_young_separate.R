library(tidyverse)
library(rstatix)
library(ggpubr)
library(ggbreak)


# Alpha diveristy calcualted based on species by epi2me
df <- read_csv('wf-metagenomics-diversity.csv')%>%
  select(-c(MEY15_cntrl))%>%
  select(-total)

metadata <- read_csv('Metadata.csv')

transposed_df <- as.data.frame(t(df))
colnames(transposed_df) <- transposed_df[1, ]
transposed_df <- transposed_df[-1, ]
transposed_df$SampleID <- rownames(transposed_df)
rownames(transposed_df) <- NULL

# Reorder the columns to move 'SampleID' to the left
transposed_df <- transposed_df[c("SampleID", setdiff(colnames(transposed_df), "SampleID"))]
df2 <- transposed_df

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

df2 <-   left_join(df2, metadata, by = 'SampleID')


# Adding age treatment column
columns_to_convert <- setdiff(names(df2), c("SampleID", "Treatment", "Sex", "Timepoint"))
df2 <- df2 %>%
  mutate_at(vars(columns_to_convert), as.numeric)%>%
  mutate(Age = ifelse(substr(SampleID, 3, 3) == "O", "Old", 
                      ifelse(substr(SampleID, 3, 3) == "Y", "young", NA)))
colnames(df2) <- gsub(" ", "_", colnames(df2))
colnames(df2) <- gsub("'", "", colnames(df2))

df2 <- df2%>%
  mutate(Simpsons_index_ranked_transformed = rank(df2$Simpsons_index))
# Grouping by Treatment and Timepoint 
df3 <- df2%>%group_by(Treatment, Timepoint, Age)%>%
  summarize(Total_counts = mean(Total_counts),
            Richness = mean(Richness),
            Shannon_Diversity_Index = mean(Shannon_diversity_index),
            Effective_N_Species = mean(Effective_number_of_species),
            Simpsons_Index = mean(Simpsons_index),
            Inverse_Simpson = mean(Inverse_Simpsons_index), 
            Pielou_evenness = mean(Pielous_evenness)
  )
# Creating new column TreatmentID which is a concatenation of Treatment and Timepoint using 'paste'
df3 <- df3%>%
  mutate(TreatmentID = paste(Treatment, Timepoint, sep = '_'))%>%
  mutate(TreatmentAge = paste(TreatmentID, Age, sep = '_'))

young.NC <- df2 %>%
  filter(Age == "young") %>%
  filter(Treatment == "NC")

young.NE <- df2 %>%
  filter(Age == "young") %>%
  filter(Treatment == "NE")

old.NC <- df2 %>%
  filter(Age == "Old") %>%
  filter(Treatment == "NC")

old.NE <- df2 %>%
  filter(Age == "Old") %>%
  filter(Treatment == "NE")
# ------------------------------

# Shannon diversity

diversity_metrics <- "Shannon_diversity_index"

# Young
# Shapiro-Wilk test for normality
shan_wilk <- young.NC%>%
  group_by(Timepoint, Treatment)%>%
  shapiro_test(Shannon_diversity_index)
shan_wilk
# For most groups data is not normally distributed

# QQ plot for residuals
ggqqplot(young.NC, "Shannon_diversity_index") + 
  facet_grid(Treatment ~ Timepoint, labeller = "label_both")


# 2 way repeat measure ANOVA
res.aov <- anova_test(
  data = young.NC, dv = diversity_metrics, wid = SampleID,
  within = Timepoint
)
get_anova_table(res.aov)
# No sig time effect

# Young.NE

# Shapiro-Wilk test for normality
shan_wilk <- young.NE%>%
  group_by(Timepoint, Treatment)%>%
  shapiro_test(Shannon_diversity_index)
shan_wilk
# For most groups data is not normally distributed

# QQ plot for residuals
ggqqplot(young.NE, "Shannon_diversity_index") + 
  facet_grid(Treatment ~ Timepoint, labeller = "label_both")


# 2 way repeat measure ANOVA
res.aov <- anova_test(
  data = young.NE, dv = diversity_metrics, wid = SampleID,
  within = Timepoint
)
get_anova_table(res.aov)
# No sig time effect

# Old.NC 

# Shapiro-Wilk test for normality
shan_wilk <- old.NC%>%
  group_by(Timepoint, Treatment)%>%
  shapiro_test(Shannon_diversity_index)
shan_wilk
# For most groups data is not normally distributed

# QQ plot for residuals
ggqqplot(old.NC, "Shannon_diversity_index") + 
  facet_grid(Treatment ~ Timepoint, labeller = "label_both")


# 2 way repeat measure ANOVA
res.aov <- anova_test(
  data = old.NC, dv = diversity_metrics, wid = SampleID,
  within = Timepoint
)
get_anova_table(res.aov)
# No sig time effect

# Old.NE

# Shapiro-Wilk test for normality
shan_wilk <- old.NE%>%
  group_by(Timepoint, Treatment)%>%
  shapiro_test(Shannon_diversity_index)
shan_wilk
# For most groups data is not normally distributed

# QQ plot for residuals
ggqqplot(old.NE, "Shannon_diversity_index") + 
  facet_grid(Treatment ~ Timepoint, labeller = "label_both")


# 2 way repeat measure ANOVA
res.aov <- anova_test(
  data = old.NE, dv = diversity_metrics, wid = SampleID,
  within = Timepoint
)
get_anova_table(res.aov)
# Significant timepoint effect

# Pairwise comparisons
pwc <- old.NE %>%
  pairwise_t_test(
    Shannon_diversity_index ~ Timepoint, paired = T, 
    p.adjust.method = "bonferroni"
  )
pwc

old.NE %>%
  group_by(Timepoint)%>%
  summarize(mean_shannon = mean(Shannon_diversity_index))

# -----------------------
# Shannon diversity

diversity_metrics <- "Shannon_diversity_index"

# Young
# Shapiro-Wilk test for normality
shan_wilk <- young.NC%>%
  group_by(Timepoint, Treatment)%>%
  shapiro_test(Shannon_diversity_index)
shan_wilk
# For most groups data is not normally distributed

# QQ plot for residuals
ggqqplot(young.NC, "Shannon_diversity_index") + 
  facet_grid(Treatment ~ Timepoint, labeller = "label_both")


# 2 way repeat measure ANOVA
res.aov <- anova_test(
  data = young.NC, dv = diversity_metrics, wid = SampleID,
  within = Timepoint
)
get_anova_table(res.aov)
# No sig time effect

# Young.NE

# Shapiro-Wilk test for normality
shan_wilk <- young.NE%>%
  group_by(Timepoint, Treatment)%>%
  shapiro_test(Shannon_diversity_index)
shan_wilk
# For most groups data is not normally distributed

# QQ plot for residuals
ggqqplot(young.NE, "Shannon_diversity_index") + 
  facet_grid(Treatment ~ Timepoint, labeller = "label_both")


# 2 way repeat measure ANOVA
res.aov <- anova_test(
  data = young.NE, dv = diversity_metrics, wid = SampleID,
  within = Timepoint
)
get_anova_table(res.aov)
# No sig time effect

# Old.NC 

# Shapiro-Wilk test for normality
shan_wilk <- old.NC%>%
  group_by(Timepoint, Treatment)%>%
  shapiro_test(Shannon_diversity_index)
shan_wilk
# For most groups data is not normally distributed

# QQ plot for residuals
ggqqplot(old.NC, diversity_metrics) + 
  facet_grid(Treatment ~ Timepoint, labeller = "label_both")


# 2 way repeat measure ANOVA
res.aov <- anova_test(
  data = old.NC, dv = diversity_metrics, wid = SampleID,
  within = Timepoint
)
get_anova_table(res.aov)
# No sig time effect

# Old.NE

# Shapiro-Wilk test for normality
shan_wilk <- old.NE%>%
  group_by(Timepoint, Treatment)%>%
  shapiro_test(Shannon_diversity_index)
shan_wilk
# For most groups data is not normally distributed

# QQ plot for residuals
ggqqplot(old.NE, "Shannon_diversity_index") + 
  facet_grid(Treatment ~ Timepoint, labeller = "label_both")


# 2 way repeat measure ANOVA
res.aov <- anova_test(
  data = old.NE, dv = diversity_metrics, wid = SampleID,
  within = Timepoint
)
get_anova_table(res.aov)


# -------------------

# Simpsons_rank transformed

diversity_metrics <- "Simpsons_index_ranked_transformed"

# Young
# Shapiro-Wilk test for normality
shan_wilk <- young.NC%>%
  group_by(Timepoint, Treatment)%>%
  shapiro_test(diversity_metrics)
shan_wilk
# For most groups data is not normally distributed

# QQ plot for residuals
ggqqplot(young.NC, diversity_metrics) + 
  facet_grid(Treatment ~ Timepoint, labeller = "label_both")


# 2 way repeat measure ANOVA
res.aov <- anova_test(
  data = young.NC, dv = diversity_metrics, wid = SampleID,
  within = Timepoint
)
get_anova_table(res.aov)
# No sig time effect

# Young.NE

# Shapiro-Wilk test for normality
shan_wilk <- young.NE%>%
  group_by(Timepoint, Treatment)%>%
  shapiro_test(diversity_metrics)
shan_wilk
# For most groups data is not normally distributed

# QQ plot for residuals
ggqqplot(young.NE, diversity_metrics) + 
  facet_grid(Treatment ~ Timepoint, labeller = "label_both")


# 2 way repeat measure ANOVA
res.aov <- anova_test(
  data = young.NE, dv = diversity_metrics, wid = SampleID,
  within = Timepoint
)
get_anova_table(res.aov)
# No sig time effect

# Old.NC 

# Shapiro-Wilk test for normality
shan_wilk <- old.NC%>%
  group_by(Timepoint, Treatment)%>%
  shapiro_test(diversity_metrics)
shan_wilk
# For most groups data is not normally distributed

# QQ plot for residuals
ggqqplot(old.NC, diversity_metrics) + 
  facet_grid(Treatment ~ Timepoint, labeller = "label_both")


# 2 way repeat measure ANOVA
res.aov <- anova_test(
  data = old.NC, dv = diversity_metrics, wid = SampleID,
  within = Timepoint
)
get_anova_table(res.aov)
# No sig time effect

# Old.NE

# Shapiro-Wilk test for normality
shan_wilk <- old.NE%>%
  group_by(Timepoint, Treatment)%>%
  shapiro_test(diversity_metrics)
shan_wilk
# For most groups data is not normally distributed

# QQ plot for residuals
ggqqplot(old.NE, diversity_metrics) + 
  facet_grid(Treatment ~ Timepoint, labeller = "label_both")


# 2 way repeat measure ANOVA
res.aov <- anova_test(
  data = old.NE, dv = diversity_metrics, wid = SampleID,
  within = Timepoint
)
get_anova_table(res.aov)


# -------------------------------
# Inverse Simpsons index
diversity_metrics <- "Inverse_Simpsons_index"

# Young
# Shapiro-Wilk test for normality
shan_wilk <- young.NC%>%
  group_by(Timepoint, Treatment)%>%
  shapiro_test(diversity_metrics)
shan_wilk
# For most groups data is not normally distributed

# QQ plot for residuals
ggqqplot(young.NC, diversity_metrics) + 
  facet_grid(Treatment ~ Timepoint, labeller = "label_both")


# 2 way repeat measure ANOVA
res.aov <- anova_test(
  data = young.NC, dv = diversity_metrics, wid = SampleID,
  within = Timepoint
)
get_anova_table(res.aov)
# No sig time effect

# Young.NE
young.NE <- young.NE %>%
  mutate(invsimp_ranked = rank(young.NE$Inverse_Simpsons_index))
diversity_metrics <- "invsimp_ranked"
# Shapiro-Wilk test for normality
shan_wilk <- young.NE%>%
  group_by(Timepoint, Treatment)%>%
  shapiro_test(diversity_metrics)
shan_wilk
# For most groups data is not normally distributed

# QQ plot for residuals
ggqqplot(young.NE, diversity_metrics) + 
  facet_grid(Treatment ~ Timepoint, labeller = "label_both")


# 2 way repeat measure ANOVA
res.aov <- anova_test(
  data = young.NE, dv = diversity_metrics, wid = SampleID,
  within = Timepoint
)
get_anova_table(res.aov)
# No sig time effect
diversity_metrics <- "Inverse_Simpsons_index"
# Old.NC 

# Shapiro-Wilk test for normality
shan_wilk <- old.NC%>%
  group_by(Timepoint, Treatment)%>%
  shapiro_test(diversity_metrics)
shan_wilk
# For most groups data is not normally distributed

# QQ plot for residuals
ggqqplot(old.NC, diversity_metrics) + 
  facet_grid(Treatment ~ Timepoint, labeller = "label_both")


# 2 way repeat measure ANOVA
res.aov <- anova_test(
  data = old.NC, dv = diversity_metrics, wid = SampleID,
  within = Timepoint
)
get_anova_table(res.aov)
# No sig time effect

# Old.NE

# Shapiro-Wilk test for normality
shan_wilk <- old.NE%>%
  group_by(Timepoint, Treatment)%>%
  shapiro_test(diversity_metrics)
shan_wilk
# For most groups data is not normally distributed

# QQ plot for residuals
ggqqplot(old.NE, diversity_metrics) + 
  facet_grid(Treatment ~ Timepoint, labeller = "label_both")


# 2 way repeat measure ANOVA
res.aov <- anova_test(
  data = old.NE, dv = diversity_metrics, wid = SampleID,
  within = Timepoint
)
get_anova_table(res.aov)

# -----------------
# Pielous Evenness

diversity_metrics <- "Pielous_evenness"

# Young
# Shapiro-Wilk test for normality
shan_wilk <- young.NC%>%
  group_by(Timepoint, Treatment)%>%
  shapiro_test(diversity_metrics)
shan_wilk
# For most groups data is not normally distributed

# QQ plot for residuals
ggqqplot(young.NC, diversity_metrics) + 
  facet_grid(Treatment ~ Timepoint, labeller = "label_both")


# 2 way repeat measure ANOVA
res.aov <- anova_test(
  data = young.NC, dv = diversity_metrics, wid = SampleID,
  within = Timepoint
)
get_anova_table(res.aov)
# No sig time effect

# Young.NE

# Shapiro-Wilk test for normality
shan_wilk <- young.NE%>%
  group_by(Timepoint, Treatment)%>%
  shapiro_test(diversity_metrics)
shan_wilk
# For most groups data is not normally distributed

# QQ plot for residuals
ggqqplot(young.NE, diversity_metrics) + 
  facet_grid(Treatment ~ Timepoint, labeller = "label_both")


# 2 way repeat measure ANOVA
res.aov <- anova_test(
  data = young.NE, dv = diversity_metrics, wid = SampleID,
  within = Timepoint
)
get_anova_table(res.aov)
# No sig time effect

# Old.NC 

# Shapiro-Wilk test for normality
shan_wilk <- old.NC%>%
  group_by(Timepoint, Treatment)%>%
  shapiro_test(diversity_metrics)
shan_wilk
# For most groups data is not normally distributed

# QQ plot for residuals
ggqqplot(old.NC, diversity_metrics) + 
  facet_grid(Treatment ~ Timepoint, labeller = "label_both")


# 2 way repeat measure ANOVA
res.aov <- anova_test(
  data = old.NC, dv = diversity_metrics, wid = SampleID,
  within = Timepoint
)
get_anova_table(res.aov)
# No sig time effect

# Old.NE

# Shapiro-Wilk test for normality
shan_wilk <- old.NE%>%
  group_by(Timepoint, Treatment)%>%
  shapiro_test(diversity_metrics)
shan_wilk
# For most groups data is not normally distributed

# QQ plot for residuals
ggqqplot(old.NE, diversity_metrics) + 
  facet_grid(Treatment ~ Timepoint, labeller = "label_both")


# 2 way repeat measure ANOVA
res.aov <- anova_test(
  data = old.NE, dv = diversity_metrics, wid = SampleID,
  within = Timepoint
)
get_anova_table(res.aov)
# Significant timepoint effect P = 0.024
old.NE %>%
  group_by(Timepoint)%>%
  summarize(mean(Pielous_evenness))

# ---------------------------------

# Richness

diversity_metrics <- "Richness"

# Young
# Shapiro-Wilk test for normality
shan_wilk <- young.NC%>%
  group_by(Timepoint, Treatment)%>%
  shapiro_test(diversity_metrics)
shan_wilk
# For most groups data is not normally distributed

# QQ plot for residuals
ggqqplot(young.NC, diversity_metrics) + 
  facet_grid(Treatment ~ Timepoint, labeller = "label_both")


# 2 way repeat measure ANOVA
res.aov <- anova_test(
  data = young.NC, dv = diversity_metrics, wid = SampleID,
  within = Timepoint
)
get_anova_table(res.aov)
# No sig time effect

# Young.NE

# Shapiro-Wilk test for normality
shan_wilk <- young.NE%>%
  group_by(Timepoint, Treatment)%>%
  shapiro_test(diversity_metrics)
shan_wilk
# For most groups data is not normally distributed

# QQ plot for residuals
ggqqplot(young.NE, diversity_metrics) + 
  facet_grid(Treatment ~ Timepoint, labeller = "label_both")


# 2 way repeat measure ANOVA
res.aov <- anova_test(
  data = young.NE, dv = diversity_metrics, wid = SampleID,
  within = Timepoint
)
get_anova_table(res.aov)
# No sig time effect

# Old.NC 

# Shapiro-Wilk test for normality
shan_wilk <- old.NC%>%
  group_by(Timepoint, Treatment)%>%
  shapiro_test(diversity_metrics)
shan_wilk
# For most groups data is not normally distributed

# QQ plot for residuals
ggqqplot(old.NC, diversity_metrics) + 
  facet_grid(Treatment ~ Timepoint, labeller = "label_both")


# 2 way repeat measure ANOVA
res.aov <- anova_test(
  data = old.NC, dv = diversity_metrics, wid = SampleID,
  within = Timepoint
)
get_anova_table(res.aov)
# No sig time effect

# Old.NE

# Shapiro-Wilk test for normality
shan_wilk <- old.NE%>%
  group_by(Timepoint, Treatment)%>%
  shapiro_test(diversity_metrics)
shan_wilk
# For most groups data is not normally distributed

# QQ plot for residuals
ggqqplot(old.NE, diversity_metrics) + 
  facet_grid(Treatment ~ Timepoint, labeller = "label_both")


# 2 way repeat measure ANOVA
res.aov <- anova_test(
  data = old.NE, dv = diversity_metrics, wid = SampleID,
  within = Timepoint
)
get_anova_table(res.aov)
# Not significant
old.NE %>%
  group_by(Timepoint)%>%
  summarize(mean(Pielous_evenness))

# ------------------------------------------
# Berger Parker Index

diversity_metrics <- "Berger_Parker_index"

# Young
# Shapiro-Wilk test for normality
shan_wilk <- young.NC%>%
  group_by(Timepoint, Treatment)%>%
  shapiro_test(diversity_metrics)
shan_wilk
# For most groups data is not normally distributed

# QQ plot for residuals
ggqqplot(young.NC, diversity_metrics) + 
  facet_grid(Treatment ~ Timepoint, labeller = "label_both")


# 2 way repeat measure ANOVA
res.aov <- anova_test(
  data = young.NC, dv = diversity_metrics, wid = SampleID,
  within = Timepoint
)
get_anova_table(res.aov)
# No sig time effect

# Young.NE

# Shapiro-Wilk test for normality
shan_wilk <- young.NE%>%
  group_by(Timepoint, Treatment)%>%
  shapiro_test(diversity_metrics)
shan_wilk
# For most groups data is not normally distributed

# QQ plot for residuals
ggqqplot(young.NE, diversity_metrics) + 
  facet_grid(Treatment ~ Timepoint, labeller = "label_both")


# 2 way repeat measure ANOVA
res.aov <- anova_test(
  data = young.NE, dv = diversity_metrics, wid = SampleID,
  within = Timepoint
)
get_anova_table(res.aov)
# No sig time effect

# Old.NC 

# Shapiro-Wilk test for normality
shan_wilk <- old.NC%>%
  group_by(Timepoint, Treatment)%>%
  shapiro_test(diversity_metrics)
shan_wilk
# For most groups data is not normally distributed

# QQ plot for residuals
ggqqplot(old.NC, diversity_metrics) + 
  facet_grid(Treatment ~ Timepoint, labeller = "label_both")


# 2 way repeat measure ANOVA
res.aov <- anova_test(
  data = old.NC, dv = diversity_metrics, wid = SampleID,
  within = Timepoint
)
get_anova_table(res.aov)
# No sig time effect

# Old.NE

# ---------------------------

# Fishers alpha
diversity_metrics <- "Fishers_alpha"

# Young
# Shapiro-Wilk test for normality
shan_wilk <- young.NC%>%
  group_by(Timepoint, Treatment)%>%
  shapiro_test(diversity_metrics)
shan_wilk
# For most groups data is not normally distributed

# QQ plot for residuals
ggqqplot(young.NC, diversity_metrics) + 
  facet_grid(Treatment ~ Timepoint, labeller = "label_both")


# 2 way repeat measure ANOVA
res.aov <- anova_test(
  data = young.NC, dv = diversity_metrics, wid = SampleID,
  within = Timepoint
)
get_anova_table(res.aov)
# No sig time effect

# Young.NE

# Shapiro-Wilk test for normality
shan_wilk <- young.NE%>%
  group_by(Timepoint, Treatment)%>%
  shapiro_test(diversity_metrics)
shan_wilk
# For most groups data is not normally distributed

# QQ plot for residuals
ggqqplot(young.NE, diversity_metrics) + 
  facet_grid(Treatment ~ Timepoint, labeller = "label_both")


# 2 way repeat measure ANOVA
res.aov <- anova_test(
  data = young.NE, dv = diversity_metrics, wid = SampleID,
  within = Timepoint
)
get_anova_table(res.aov)
# No sig time effect

# Old.NC 

# Shapiro-Wilk test for normality
shan_wilk <- old.NC%>%
  group_by(Timepoint, Treatment)%>%
  shapiro_test(diversity_metrics)
shan_wilk
# For most groups data is not normally distributed

# QQ plot for residuals
ggqqplot(old.NC, diversity_metrics) + 
  facet_grid(Treatment ~ Timepoint, labeller = "label_both")


# 2 way repeat measure ANOVA
res.aov <- anova_test(
  data = old.NC, dv = diversity_metrics, wid = SampleID,
  within = Timepoint
)
get_anova_table(res.aov)
# No sig time effect

# Old.NE

# Shapiro-Wilk test for normality
shan_wilk <- old.NE%>%
  group_by(Timepoint, Treatment)%>%
  shapiro_test(diversity_metrics)
shan_wilk
# For most groups data is not normally distributed

# QQ plot for residuals
ggqqplot(old.NE, diversity_metrics) + 
  facet_grid(Treatment ~ Timepoint, labeller = "label_both")


# 2 way repeat measure ANOVA
res.aov <- anova_test(
  data = old.NE, dv = diversity_metrics, wid = SampleID,
  within = Timepoint
)
get_anova_table(res.aov)
# Not significant
old.NE %>%
  group_by(Timepoint)%>%
  summarize(mean(Pielous_evenness))

# ---------------------------------------
# Effective number of species

diversity_metrics <- "Effective_number_of_species"

# Young
# Shapiro-Wilk test for normality
shan_wilk <- young.NC%>%
  group_by(Timepoint, Treatment)%>%
  shapiro_test(diversity_metrics)
shan_wilk
# For most groups data is not normally distributed

# QQ plot for residuals
ggqqplot(young.NC, diversity_metrics) + 
  facet_grid(Treatment ~ Timepoint, labeller = "label_both")


# 2 way repeat measure ANOVA
res.aov <- anova_test(
  data = young.NC, dv = diversity_metrics, wid = SampleID,
  within = Timepoint
)
get_anova_table(res.aov)
# No sig time effect

# Young.NE

# Shapiro-Wilk test for normality
shan_wilk <- young.NE%>%
  group_by(Timepoint, Treatment)%>%
  shapiro_test(diversity_metrics)
shan_wilk
# For most groups data is not normally distributed

# QQ plot for residuals
ggqqplot(young.NE, diversity_metrics) + 
  facet_grid(Treatment ~ Timepoint, labeller = "label_both")


# 2 way repeat measure ANOVA
res.aov <- anova_test(
  data = young.NE, dv = diversity_metrics, wid = SampleID,
  within = Timepoint
)
get_anova_table(res.aov)
# No sig time effect

# Old.NC 

# Shapiro-Wilk test for normality
shan_wilk <- old.NC%>%
  group_by(Timepoint, Treatment)%>%
  shapiro_test(diversity_metrics)
shan_wilk
# For most groups data is not normally distributed

# QQ plot for residuals
ggqqplot(old.NC, diversity_metrics) + 
  facet_grid(Treatment ~ Timepoint, labeller = "label_both")


# 2 way repeat measure ANOVA
res.aov <- anova_test(
  data = old.NC, dv = diversity_metrics, wid = SampleID,
  within = Timepoint
)
get_anova_table(res.aov)
# No sig time effect

# Old.NE

# Shapiro-Wilk test for normality
shan_wilk <- old.NE%>%
  group_by(Timepoint, Treatment)%>%
  shapiro_test(diversity_metrics)
shan_wilk
# For most groups data is not normally distributed

# QQ plot for residuals
ggqqplot(old.NE, diversity_metrics) + 
  facet_grid(Treatment ~ Timepoint, labeller = "label_both")


# 2 way repeat measure ANOVA
res.aov <- anova_test(
  data = old.NE, dv = diversity_metrics, wid = SampleID,
  within = Timepoint
)
get_anova_table(res.aov)
# Significant timepoint effect
old.NE %>%
  group_by(Timepoint)%>%
  summarize(mean(Effective_number_of_species))

# -----------------------------------------
# Total counts

diversity_metrics <- "Total_counts"

# Young
# Shapiro-Wilk test for normality
shan_wilk <- young.NC%>%
  group_by(Timepoint, Treatment)%>%
  shapiro_test(diversity_metrics)
shan_wilk
# For most groups data is not normally distributed

# QQ plot for residuals
ggqqplot(young.NC, diversity_metrics) + 
  facet_grid(Treatment ~ Timepoint, labeller = "label_both")


# 2 way repeat measure ANOVA
res.aov <- anova_test(
  data = young.NC, dv = diversity_metrics, wid = SampleID,
  within = Timepoint
)
get_anova_table(res.aov)
# No sig time effect

# Young.NE

# Shapiro-Wilk test for normality
shan_wilk <- young.NE%>%
  group_by(Timepoint, Treatment)%>%
  shapiro_test(diversity_metrics)
shan_wilk
# For most groups data is not normally distributed

# QQ plot for residuals
ggqqplot(young.NE, diversity_metrics) + 
  facet_grid(Treatment ~ Timepoint, labeller = "label_both")


# 2 way repeat measure ANOVA
res.aov <- anova_test(
  data = young.NE, dv = diversity_metrics, wid = SampleID,
  within = Timepoint
)
get_anova_table(res.aov)
# No sig time effect

# Old.NC 

# Shapiro-Wilk test for normality
shan_wilk <- old.NC%>%
  group_by(Timepoint, Treatment)%>%
  shapiro_test(diversity_metrics)
shan_wilk
# For most groups data is not normally distributed

# QQ plot for residuals
ggqqplot(old.NC, diversity_metrics) + 
  facet_grid(Treatment ~ Timepoint, labeller = "label_both")


# 2 way repeat measure ANOVA
res.aov <- anova_test(
  data = old.NC, dv = diversity_metrics, wid = SampleID,
  within = Timepoint
)
get_anova_table(res.aov)
# No sig time effect

# Old.NE

# Shapiro-Wilk test for normality
shan_wilk <- old.NE%>%
  group_by(Timepoint, Treatment)%>%
  shapiro_test(diversity_metrics)
shan_wilk
# For most groups data is not normally distributed

# QQ plot for residuals
ggqqplot(old.NE, diversity_metrics) + 
  facet_grid(Treatment ~ Timepoint, labeller = "label_both")


# 2 way repeat measure ANOVA
res.aov <- anova_test(
  data = old.NE, dv = diversity_metrics, wid = SampleID,
  within = Timepoint
)
get_anova_table(res.aov)
# Not significant
old.NE %>%
  group_by(Timepoint)%>%
  summarize(mean(Effective_number_of_species))


