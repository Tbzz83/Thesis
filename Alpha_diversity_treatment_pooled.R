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
                      ifelse(substr(SampleID, 3, 3) == "Y", "young", NA))) %>%
  filter(Timepoint == 0)
colnames(df2) <- gsub(" ", "_", colnames(df2))
colnames(df2) <- gsub("'", "", colnames(df2))


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



# ---------------------------------------------------------------------

# T-test of metrics young vs old week 0

# Check normality with shapiro
df2 %>%
  shapiro_test(Shannon_diversity_index)
# P = 0.0000254

df2 %>%
  shapiro_test(Berger_Parker_index)
# ^ sig

df2 %>%
  shapiro_test(Simpsons_index)
# ^ sig


df2 %>%
  shapiro_test(Inverse_Simpsons_index)
# ^ sig

df2 %>%
  shapiro_test(Fishers_alpha)

df2 %>%
  shapiro_test(Pielous_evenness)
# ^ sig

df2 %>%
  shapiro_test(Richness)
# ^ sig

df2 %>%
  shapiro_test(Effective_number_of_species)
# ^ Not sig

df2 %>%
  shapiro_test(Total_counts)
# ^ sig

# QQ plot of residuals
ggqqplot(df2, x ="Shannon_diversity_index")

ggqqplot(df2, x ="Simpsons_index")

ggqqplot(df2, x ="Inverse_Simpsons_index")

ggqqplot(df2, x ="Richness")

ggqqplot(df2, x ="Effective_number_of_species")

ggqqplot(df2, x ="Total_counts")

ggqqplot(df2, x ="Berger_Parker_index")

ggqqplot(df2, x ="Fishers_alpha")

# Several of these plots look slightly unstable so will use 
# non-parametric Wilcoxon test

shan <- wilcox.test(Shannon_diversity_index ~ Age, data = df2, 
                    exact = FALSE)
shan
# ^ not sig
simp <- wilcox.test(Simpsons_index ~ Age, data = df2, 
                    exact = FALSE)
simp
# ^ not sig

invsimp <- wilcox.test(Inverse_Simpsons_index ~ Age, data = df2, 
                       exact = FALSE)
invsimp
# ^ not sig

rich <- wilcox.test(Richness ~ Age, data = df2, 
                    exact = FALSE)
rich
# ^ not sig

pie <- wilcox.test(Pielous_evenness ~ Age, data = df2, 
                   exact = FALSE)
pie
# ^ not sig

f.alpha <- wilcox.test(Fishers_alpha ~ Age, data = df2, 
                       exact = FALSE)
f.alpha
# ^ not sig

b.park <- wilcox.test(Berger_Parker_index ~ Age, data = df2, 
                      exact = FALSE)
b.park
# ^ not sig

eff.sp <- wilcox.test(Effective_number_of_species ~ Age, data = df2, 
                      exact = FALSE)
eff.sp
# ^ not sig

counts <- wilcox.test(Total_counts ~ Age, data = df2, 
                      exact = FALSE)
counts
# ^ not sig

# --------------------------------------
# Boxplt
# Box plot of initial data
a <- ggboxplot(
  df2, x = 'Age', y = 'Shannon_diversity_index',
  color = 'Age', palette = 'jco',
  short.panel.labs = FALSE
)


b <- ggboxplot(
  df2, x = 'Age', y = 'Simpsons_index',
  color = 'Age', palette = 'jco',
  short.panel.labs = FALSE
)

c <- ggboxplot(
  df2, x = 'Age', y = 'Inverse_Simpsons_index',
  color = 'Age', palette = 'jco',
  short.panel.labs = FALSE
)

d <- ggboxplot(
  df2, x = 'Age', y = 'Richness',
  color = 'Age', palette = 'jco',
  short.panel.labs = FALSE
)

e <- ggboxplot(
  df2, x = 'Age', y = 'Pielous_evenness',
  color = 'Age', palette = 'jco',
  short.panel.labs = FALSE
)

f <- ggboxplot(
  df2, x = 'Age', y = 'Effective_number_of_species',
  color = 'Age', palette = 'jco',
  short.panel.labs = FALSE
)

g <- ggboxplot(
  df2, x = 'Age', y = 'Berger_Parker_index',
  color = 'Age', palette = 'jco',
  short.panel.labs = FALSE
)

h <- ggboxplot(
  df2, x = 'Age', y = 'Fishers_alpha',
  color = 'Age', palette = 'jco',
  short.panel.labs = FALSE
)

i <- ggboxplot(
  df2, x = 'Age', y = 'Total_counts',
  color = 'Age', palette = 'jco',
  short.panel.labs = FALSE
)
library(gridExtra)

grid.arrange(a,b,c,d,e,f,g,h, i)
