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

#prelim <- read_csv('prelim_diversity_to_join.csv')%>%
#  mutate_if(is.numeric, as.character)%>%
#  select(-Treatment)

#prelim <- left_join(prelim, metadata, by = 'SampleID')%>%
#  filter(Timepoint %in% c(0,7))%>%
#  mutate(Timepoint = ifelse(Timepoint == 7, 6, Timepoint))

###FIX####df3 <- bind_rows(df2, prelim)%>%
 ## filter(SampleID != 'MEY7control')%>%
 ## as_tibble()


# Adding age treatment column
columns_to_convert <- setdiff(names(df2), c("SampleID", "Treatment", "Sex", "Timepoint"))
df2 <- df2 %>%
  mutate_at(vars(columns_to_convert), as.numeric)%>%
  mutate(Age = ifelse(substr(SampleID, 3, 3) == "O", "Old", 
                      ifelse(substr(SampleID, 3, 3) == "Y", "young", NA)))

# Grouping by Treatment and Timepoint 
df3 <- df2%>%group_by(Treatment, Timepoint, Age)%>%
  summarize(Total_counts = mean(`Total counts`),
            Richness = mean(Richness),
            Shannon_Diversity_Index = mean(`Shannon diversity index`),
            Effective_N_Species = mean(`Effective number of species`),
            Simpsons_Index = mean(`Simpson's index`),
            Inverse_Simpson = mean(`Inverse Simpson's index`), 
            Pielou_evenness = mean(`Pielou's evenness`)
  )
# Creating new column TreatmentID which is a concatenation of Treatment and Timepoint using 'paste'
df3 <- df3%>%
  mutate(TreatmentID = paste(Treatment, Timepoint, sep = '_'))%>%
  mutate(TreatmentAge = paste(TreatmentID, Age, sep = '_'))

# Ordering df2 by Timepoint so can easily compare NC to NE
#df3$TreatmentID <- factor(df3$TreatmentID, levels = df3$TreatmentID[order(df3$Timepoint)])

# ----------------------------------------------------------------------------

# Shannon diversity

shan_summary <- df2%>%
  group_by(Timepoint, Treatment, Age)%>%
  get_summary_stats(`Shannon diversity index`, type = "mean_sd")
shan_summary

# Box plot of initial data
shan_bxp <- ggboxplot(
  df2, x = 'Treatment', y = 'Shannon diversity index',
  color = 'Timepoint', palette = 'jco',
  facet.by = 'Age', short.panel.labs = FALSE
)
shan_bxp

# Check assumptions
shan_outlier <- df2%>%
  group_by(Timepoint, Treatment, Age)%>%
  identify_outliers(`Shannon diversity index`)
shan_outlier
# No extreme outliers

# ------

# Have to convert column names to 'word_word' format for Rstatix to understand :(
colnames(df2) <- gsub(" ", "_", colnames(df2))
colnames(df2) <- gsub("'", "", colnames(df2))
 
# Shapiro-Wilk test for normality
shan_wilk <- df2%>%
  group_by(Timepoint, Treatment, Age)%>%
  shapiro_test(Shannon_diversity_index)
shan_wilk
# For most groups data is not normally distributed

# QQ plot for residuals
ggqqplot(df2, "Shannon_diversity_index") + 
  facet_grid(Treatment + Age ~ Timepoint, labeller = "label_both")

# Homogeneity of variance
shan_homog <- df2%>%
  group_by(Timepoint, Age)%>%
  levene_test(Shannon_diversity_index ~ Treatment)
shan_homog
# ^ All are homogenous



#2x2x2 Mixed Design ANOVA
res.aov <- anova_test(
  data = df2, dv = Shannon_diversity_index, wid = SampleID,
  between = c(Age, Treatment), within = Timepoint
)
get_anova_table(res.aov)
# * significance only in Treatment, no interaction effects

# Post-hoc simple main effect
treatment.effect.shan <- df2%>%
  group_by(Timepoint, Age) %>%
  anova_test(dv = Shannon_diversity_index, wid = SampleID, 
             between = Treatment)%>%
  get_anova_table()
treatment.effect.shan

# Pairwise comparisons
pwc.shan <- df2 %>%
  group_by(Timepoint, Age) %>%
  pairwise_t_test(
    Shannon_diversity_index ~ Treatment, paired = F, 
    p.adjust.method = "bonferroni"
  )
pwc.shan

# Bar plots of data
ShannonBar <- ggbarplot(df3, x = 'TreatmentAge', y = 'Shannon_Diversity_Index',
                        fill = 'Age', palette = 'jco')+scale_y_break(c(1, 3.8))
ShannonBar



# ----------------------------------------------------------------------------
# Simpson diversity
diversity_metrics <- "Simpsons_index"
## Looks like the data are non-normal with some outliers. Will try some 
## transformations 

df2$log_simpson <- log(df2$Simpsons_index)
df2$sqrt_simpson <- sqrt(df2$Simpsons_index)
diversity_metrics <- "Simpsons_index"

# Box plot of initial data
simp_bxp1 <- ggboxplot(
  df2, x = 'Treatment', y = 'Simpsons_index',
  color = 'Timepoint', palette = 'jco',
  facet.by = 'Age', short.panel.labs = FALSE
)
simp_bxp1


simp_summary <- df2%>%
  group_by(Timepoint, Treatment, Age)%>%
  get_summary_stats(diversity_metrics, type = "mean_sd")
simp_summary

# Check assumptions
simp_outlier <- df2%>%
  group_by(Timepoint, Treatment, Age)%>%
  identify_outliers(diversity_metrics)
simp_outlier
# ^ Are some extreme outliers

# Shapiro-Wilk test for normality
simp_wilk <- df2%>%
  group_by(Timepoint, Treatment, Age)%>%
  shapiro_test(diversity_metrics)
simp_wilk
# For most groups data is not normally distributed

# QQ plot for residuals
ggqqplot(df2, diversity_metrics) + 
  facet_grid(Treatment + Age ~ Timepoint, labeller = "label_both")

# Homogeneity of variance
simp_homog <- df2%>%
  group_by(Timepoint, Age)%>%
  levene_test(Simpsons_index ~ Treatment)
simp_homog
# ^ Are homogenous

## Looks like the data are non-normal with some outliers. Will try some 
## transformations 

#2x2x2 Mixed Design ANOVA
res.aov <- anova_test(
  data = df2, dv = diversity_metrics, wid = SampleID,
  between = c(Age, Treatment), within = Timepoint
)
get_anova_table(res.aov)
# * significance only in Treatment, no interaction effects

# Post-hoc simple main effect
treatment.effect.simp <- df2%>%
  group_by(Timepoint, Age) %>%
  anova_test(dv = diversity_metrics, wid = SampleID, 
             between = Treatment)%>%
  get_anova_table()
treatment.effect.simp

# Pairwise comparisons
pwc.simp <- df2 %>%
  group_by(Timepoint, Age) %>%
  pairwise_t_test(
    Simpsons_index ~ Treatment, paired = F, 
    p.adjust.method = "bonferroni"
  )
pwc.simp

# Doing an aligned ranks transformation ANOVA see https://rcompanion.org/handbook/F_16.html

df2 <- df2%>%
  mutate(Simpsons_index_ranked_transformed = rank(df2$Simpsons_index))
df2.young <- df2 %>%
  filter(Age == "young")
df2.old <- df2 %>%
  filter(Age == "Old")
# Simpsons anova for old
#2x2 Mixed Design ANOVA not looking at age as a factor now

# QQ plot for residuals
diversity_metrics = "Simpsons_index_ranked_transformed"

ggqqplot(df2.old, diversity_metrics) + 
  facet_grid(Treatment ~ Timepoint, labeller = "label_both")

res.aov <- anova_test(
  data = df2.old, dv = diversity_metrics, wid = SampleID,
  between = Treatment, within = Timepoint
)
get_anova_table(res.aov)
# No significance

# Box plot of initial data
simp_bxp <- ggboxplot(
  df2.old, x = 'Treatment', y = diversity_metrics,
  color = 'Timepoint', palette = 'jco', short.panel.labs = FALSE
) + 
  ggtitle("Old mice Simpsons index")
simp_bxp

# ---------
# Young mice rank simpsons

# QQ plot for residuals
diversity_metrics = "Simpsons_index_ranked_transformed"

ggqqplot(df2.young, diversity_metrics) + 
  facet_grid(Treatment ~ Timepoint, labeller = "label_both")

res.aov <- anova_test(
  data = df2.young, dv = diversity_metrics, wid = SampleID,
  between = Treatment, within = Timepoint
)
get_anova_table(res.aov)
# No significance

# Box plot of initial data
simp_bxp <- ggboxplot(
  df2.young, x = 'Treatment', y = diversity_metrics,
  color = 'Timepoint', palette = 'jco', short.panel.labs = FALSE
) + 
  ggtitle("Young mice Simpsons index")
simp_bxp





# Post-hoc simple main effect
treatment.effect.shan <- df2%>%
  group_by(Timepoint, Age) %>%
  anova_test(dv = Shannon_diversity_index, wid = SampleID, 
             between = Treatment)%>%
  get_anova_table()
treatment.effect.shan








SimpsonBar <- ggbarplot(df3, x = 'TreatmentAge', y = 'Simpsons_Index',
                        fill = 'Timepoint', palette = 'jco')+ scale_y_break(c(0.01, 0.915))
SimpsonBar
# ----------------------------------------------------------------------------
diversity_metrics <- "Inverse_Simpsons_index"

invsimp_summary <- df2%>%
  group_by(Timepoint, Treatment, Age)%>%
  get_summary_stats(diversity_metrics, type = "mean_sd")
invsimp_summary

# Box plot of initial data
invsimp_bxp <- ggboxplot(
  df2, x = 'Treatment', y = diversity_metrics,
  color = 'Timepoint', palette = 'jco',
  facet.by = 'Age', short.panel.labs = FALSE
)
invsimp_bxp

# Check assumptions
invsimp_outlier <- df2%>%
  group_by(Timepoint, Treatment, Age)%>%
  identify_outliers(diversity_metrics)
invsimp_outlier
# ^ Are some extreme outliers

# Shapiro-Wilk test for normality
invsimp_wilk <- df2%>%
  group_by(Timepoint, Treatment, Age)%>%
  shapiro_test(diversity_metrics)
invsimp_wilk
# For most groups data is not normally distributed

# QQ plot for residuals
ggqqplot(df2, diversity_metrics) + 
  facet_grid(Treatment + Age ~ Timepoint, labeller = "label_both")
## QQ plot better at higher sample size than shapiro wilk

# Homogeneity of variance
invsimp_homog <- df2%>%
  group_by(Timepoint, Age)%>%
  levene_test(Inverse_Simpsons_index ~ Treatment)
invsimp_homog
# ^ Are homogenous

## Looks like the data are normal
## transformations 

#2x2x2 Mixed Design ANOVA
invsimp.aov <- anova_test(
  data = df2, dv = diversity_metrics, wid = SampleID,
  between = c(Age, Treatment), within = Timepoint
)
get_anova_table(invsimp.aov)
# * significance only in Treatment, no interaction effects

# Post-hoc simple main effect
treatment.effect.invsimp <- df2%>%
  group_by(Timepoint, Age) %>%
  anova_test(dv = diversity_metrics, wid = SampleID, 
             between = Treatment)%>%
  get_anova_table()
treatment.effect.invsimp

# Pairwise comparisons
pwc.invsimp <- df2 %>%
  group_by(Timepoint, Age) %>%
  pairwise_t_test(
    Inverse_Simpsons_index ~ Treatment, paired = F, 
    p.adjust.method = "bonferroni"
  )
pwc.invsimp
# Adj p value not significant


InvSimpsonBar <- ggbarplot(df3, x = 'TreatmentID', y = 'Inverse_Simpson',
                           fill = 'Treatment', palette = 'jco')+scale_y_break(c(0.05, 1))
InvSimpsonBar
# ----------------------------------------------------------------------------
diversity_metrics <- "Pielous_evenness"


pie_summary <- df2%>%
  group_by(Timepoint, Treatment, Age)%>%
  get_summary_stats(diversity_metrics, type = "mean_sd")
pie_summary

# Box plot of initial data
pie_bxp <- ggboxplot(
  df2, x = 'Treatment', y = diversity_metrics,
  color = 'Timepoint', palette = 'jco',
  facet.by = 'Age', short.panel.labs = FALSE
)
pie_bxp

# Check assumptions
pie_outlier <- df2%>%
  group_by(Timepoint, Treatment, Age)%>%
  identify_outliers(diversity_metrics)
pie_outlier

# Shapiro-Wilk test for normality
pie_wilk <- df2%>%
  group_by(Timepoint, Treatment, Age)%>%
  shapiro_test(diversity_metrics)
pie_wilk
# For most groups data is not normally distributed

# QQ plot for residuals
ggqqplot(df2, diversity_metrics) + 
  facet_grid(Treatment + Age ~ Timepoint, labeller = "label_both")
## Residuals look good 

#2x2x2 Mixed Design ANOVA
pie.aov <- anova_test(
  data = df2, dv = diversity_metrics, wid = SampleID,
  between = c(Age, Treatment), within = Timepoint
)
get_anova_table(pie.aov)

# Post-hoc simple main effect
treatment.effect.pie <- df2%>%
  group_by(Timepoint, Age) %>%
  anova_test(dv = diversity_metrics, wid = SampleID, 
             between = Treatment)%>%
  get_anova_table()
treatment.effect.pie

# Pairwise comparisons
pwc.pie <- df2 %>%
  group_by(Timepoint, Age) %>%
  pairwise_t_test(
    Pielous_evenness ~ Treatment, paired = F, 
    p.adjust.method = "bonferroni"
  )
pwc.pie

EvennessBar <- ggbarplot(df3, x = 'TreatmentID', y = 'Pielou_evenness',
                         fill = 'Treatment', palette = 'jco')+scale_y_break(c(0.2, 0.65))
EvennessBar
# ----------------------------------------------------------------------------
diversity_metrics <- "Richness"


rich_summary <- df2%>%
  group_by(Timepoint, Treatment, Age)%>%
  get_summary_stats(diversity_metrics, type = "mean_sd")
rich_summary

# Box plot of initial data
rich_bxp <- ggboxplot(
  df2, x = 'Treatment', y = diversity_metrics,
  color = 'Timepoint', palette = 'jco',
  facet.by = 'Age', short.panel.labs = FALSE
)
rich_bxp

# Check assumptions
rich_outlier <- df2%>%
  group_by(Timepoint, Treatment, Age)%>%
  identify_outliers(diversity_metrics)
rich_outlier

# Shapiro-Wilk test for normality
rich_wilk <- df2%>%
  group_by(Timepoint, Treatment, Age)%>%
  shapiro_test(diversity_metrics)
rich_wilk
# For most groups data is not normally distributed

# QQ plot for residuals
ggqqplot(df2, diversity_metrics) + 
  facet_grid(Treatment + Age ~ Timepoint, labeller = "label_both")
# Residuals look fine
# Homogeneity of variance
#rich_homog <- df2%>%
#  group_by(Timepoint, Age)%>%
#  levene_test(diversity_metrics ~ Treatment)
#rich_homog



#2x2x2 Mixed Design ANOVA
rich.aov <- anova_test(
  data = df2, dv = diversity_metrics, wid = SampleID,
  between = c(Age, Treatment), within = Timepoint
)
get_anova_table(rich.aov)
# No significance


RichnessBar <- ggbarplot(df3, x = 'TreatmentID', y = 'Richness',
                         fill = 'Treatment', palette = 'jco')
RichnessBar
# ----------------------------------------------------------------------------
TotalCountsBar <- ggbarplot(df3, x = 'TreatmentID', y = 'Total_counts',
                            fill = 'Treatment', palette = 'jco')
TotalCountsBar
# ----------------------------------------------------------------------------

# ------------------------------------------
# Berger Parker index
diversity_metrics <- "Berger_Parker_index"

berg_summary <- df2%>%
  group_by(Timepoint, Treatment, Age)%>%
  get_summary_stats(diversity_metrics, type = "mean_sd")
berg_summary

# Box plot of initial data
berg_bxp <- ggboxplot(
  df2, x = 'Treatment', y = diversity_metrics,
  color = 'Timepoint', palette = 'jco',
  facet.by = 'Age', short.panel.labs = FALSE
)
berg_bxp

# Check assumptions
berg_outlier <- df2%>%
  group_by(Timepoint, Treatment, Age)%>%
  identify_outliers(diversity_metrics)
berg_outlier

# Shapiro-Wilk test for normality
berg_wilk <- df2%>%
  group_by(Timepoint, Treatment, Age)%>%
  shapiro_test(diversity_metrics)
berg_wilk
# For most groups data is not normally distributed

# QQ plot for residuals
ggqqplot(df2, diversity_metrics) + 
  facet_grid(Treatment + Age ~ Timepoint, labeller = "label_both")
## Residuals look good 

#2x2x2 Mixed Design ANOVA
berg.aov <- anova_test(
  data = df2, dv = diversity_metrics, wid = SampleID,
  between = c(Age, Treatment), within = Timepoint
)
get_anova_table(berg.aov)

# ---------------------------------------
# Effective species no.

diversity_metrics <- "Effective_number_of_species"
sp.no_summary <- df2%>%
  group_by(Timepoint, Treatment, Age)%>%
  get_summary_stats(diversity_metrics, type = "mean_sd")
sp.no_summary

# Box plot of initial data
sp.no_bxp <- ggboxplot(
  df2, x = 'Treatment', y = diversity_metrics,
  color = 'Timepoint', palette = 'jco',
  facet.by = 'Age', short.panel.labs = FALSE
)
sp.no_bxp

# Check assumptions
sp.no_outlier <- df2%>%
  group_by(Timepoint, Treatment, Age)%>%
  identify_outliers(diversity_metrics)
sp.no_outlier

# Shapiro-Wilk test for normality
sp.no_wilk <- df2%>%
  group_by(Timepoint, Treatment, Age)%>%
  shapiro_test(diversity_metrics)
sp.no_wilk
# For most groups data is not normally distributed

# QQ plot for residuals
ggqqplot(df2, diversity_metrics) + 
  facet_grid(Treatment + Age ~ Timepoint, labeller = "label_both")
## Residuals look good 

#2x2x2 Mixed Design ANOVA
sp.no.aov <- anova_test(
  data = df2, dv = diversity_metrics, wid = SampleID,
  between = c(Age, Treatment), within = Timepoint
)
get_anova_table(sp.no.aov)
# Significant main effect of Treatment AND interaction of age x treatment x timepoint
# -----
# Fisher's alpha

diversity_metrics <- "Fishers_alpha"
falpha_summary <- df2%>%
  group_by(Timepoint, Treatment, Age)%>%
  get_summary_stats(diversity_metrics, type = "mean_sd")
falpha_summary

# Box plot of initial data
falpha_bxp <- ggboxplot(
  df2, x = 'Treatment', y = diversity_metrics,
  color = 'Timepoint', palette = 'jco',
  facet.by = 'Age', short.panel.labs = FALSE
)
falpha_bxp

# Check assumptions
falpha_outlier <- df2%>%
  group_by(Timepoint, Treatment, Age)%>%
  identify_outliers(diversity_metrics)
falpha_outlier

# Shapiro-Wilk test for normality
falpha_wilk <- df2%>%
  group_by(Timepoint, Treatment, Age)%>%
  shapiro_test(diversity_metrics)
falpha_wilk
# For most groups data is not normally distributed

# QQ plot for residuals
ggqqplot(df2, diversity_metrics) + 
  facet_grid(Treatment + Age ~ Timepoint, labeller = "label_both")
## Residuals look good 

#2x2x2 Mixed Design ANOVA
falpha.aov <- anova_test(
  data = df2, dv = diversity_metrics, wid = SampleID,
  between = c(Age, Treatment), within = Timepoint
)
get_anova_table(falpha.aov)
# Significant effect in Age x Timepoint

# -----
# total counts 

diversity_metrics <- "Total_counts"
counts_summary <- df2%>%
  group_by(Timepoint, Treatment, Age)%>%
  get_summary_stats(diversity_metrics, type = "mean_sd")
counts_summary

# Box plot of initial data
counts_bxp <- ggboxplot(
  df2, x = 'Treatment', y = diversity_metrics,
  color = 'Timepoint', palette = 'jco',
  facet.by = 'Age', short.panel.labs = FALSE
)
counts_bxp

# Check assumptions
counts_outlier <- df2%>%
  group_by(Timepoint, Treatment, Age)%>%
  identify_outliers(diversity_metrics)
counts_outlier

# Shapiro-Wilk test for normality
counts_wilk <- df2%>%
  group_by(Timepoint, Treatment, Age)%>%
  shapiro_test(diversity_metrics)
counts_wilk
# For most groups data is not normally distributed

# QQ plot for residuals
ggqqplot(df2, diversity_metrics) + 
  facet_grid(Treatment + Age ~ Timepoint, labeller = "label_both")
## Residuals look good 

#2x2x2 Mixed Design ANOVA
counts.aov <- anova_test(
  data = df2, dv = diversity_metrics, wid = SampleID,
  between = c(Age, Treatment), within = Timepoint
)
get_anova_table(counts.aov)
# No significance





library(gridExtra)

grid.arrange(shan_bxp, simp_bxp1, invsimp_bxp, rich_bxp, pie_bxp, sp.no_bxp,
             berg_bxp, falpha_bxp, counts_bxp)

