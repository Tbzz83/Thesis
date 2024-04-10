library(tidyverse)
library(edgeR)

df <- read_csv('wf-metagenomics-counts-phylum.csv')%>%
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
df3 <- read_csv('wf-metagenomics-counts-phylum.csv')%>%
  select(-MEY15_cntrl)%>%
  select(-c(total, superkingdom, kingdom,tax))%>%
  rename_with(~str_replace(., "_7$", "_6"), ends_with("_7"))%>%
  as.data.frame()

row_names <- df3[, 1]  # Extract the values from the first column
df3 <- df3[, -1]  # Remove the first column from the data frame
rownames(df3) <- row_names
df3 <- as.matrix(df3)

d <- DGEList(counts = df3, group = Group)
keep <- filterByExpr(d, group=Group)
d <- d[keep, , keep.lib.sizes=FALSE]
d <- estimateDisp(d)
d <- calcNormFactors(d)

# MDS plot
cpm <- cpm(d)%>%
  as.data.frame()
logcpm <- cpm(d, log=TRUE)



plotMDS(d, top = 1000, col = as.numeric(d$samples$group), pch = 16,
        main = "PCoA of phyla between groups")
fit <- glmQLFit(d, design)
tr <- glmTreat(fit, coef = 2, lfc = 1)
topTags(tr)

my.contrasts <- makeContrasts(
  # Within young comparisons
  NE.0.youngvsNC.0.young = NC.0.young-NE.0.young,
  NE.6.youngvsNC.6.young = NC.6.young-NE.6.young,
  NE.0.youngvsNE.6.young = NE.6.young-NE.0.young,
  NC.0.youngvsNC.6.young = NC.6.young-NC.0.young,
  
  # Within old comparisons
  NE.0.oldvsNC.0.old = NC.0.old-NE.0.old,
  NE.6.oldvsNC.6.old = NC.6.old-NE.6.old,
  NE.0.oldvsNE.6.old = NE.6.old-NE.0.old,
  NC.0.oldvsNC.6.old = NC.6.old-NC.0.old,
  
  # Between age comparisons
  NE.0.youngvsNE.0.old = NE.0.young-NE.0.old,
  NE.6.youngvsNE.6.old = NE.6.young-NE.6.old,
  
  NC.0.youngvsNC.0.old = NC.0.young-NC.0.old,
  NC.6.youngvsNC.6.old = NC.6.young-NC.6.old,
  
  # Between age overall
  
  
  
  levels = design
  
)

# Can look at total counts per million by group
group_cpm <- cpmByGroup(d)


## All plots are with regard to second group in title. So if downregulate, it
## means it was expressed less in the second group

qlf <- glmQLFTest(fit, contrast=my.contrasts[,"NE.0.youngvsNC.0.young"])
test_NE.0.youngvsNC.0.young <- qlf$table%>%
  mutate_all(as.numeric)%>%
  filter(PValue < 0.05)
plotMD(qlf)
# Filter for features with significant adj p-value (FDR < 0.05)
#res <- as.data.frame(topTags(qlf, n = 3, sort.by = "p.value"))%>%
#  filter(FDR < 0.05)
#text(res$logCPM, res$logFC, labels = rownames(res), col = "black", cex = 0.5,
#     pos = 3)
is.de <- decideTestsDGE(qlf)
summary(is.de)

# Plot species by comparison
#res2 = res%>%
#  rownames_to_column(var = "tax")
y = as.data.frame(group_cpm)%>%
  select(c("NE.0.young","NC.0.young"))%>%
  rownames_to_column(var="tax")
y = inner_join(y,res2,by='tax')
y <- pivot_longer(y, cols = 2:3,
                  names_to = "group",
                  values_to = "value")
y$comparison <- 'NE.0.youngvsNC.0.young'
# all.comp will be master dataframe for bar charts
all.comp <- y
#ggplot(all.comp, aes(x = group, y = value, fill = tax)) +
  # Implement a grouped bar chart
#  geom_bar(position = "dodge", stat = "identity")

qlf <- glmQLFTest(fit, contrast=my.contrasts[,"NE.6.oldvsNC.6.old"])
test_NE.6.oldvsNC.6.old <- qlf$table%>%
  mutate_all(as.numeric)%>%
  filter(PValue < 0.05)
is.de <- decideTestsDGE(qlf)
topTags(qlf)
plotMD(qlf)
# Filter for features with significant adj p-value (FDR < 0.05)
res <- as.data.frame(topTags(qlf, sort.by = "p.value", n=4))%>%
  filter(FDR < 0.05)
#text(res$logCPM, res$logFC, labels = rownames(res), col = "black", cex = 0.5,
#     pos = 3)
is.de <- decideTestsDGE(qlf)
summary(is.de)
# Plot species by comparison
res2 = res%>%
  rownames_to_column(var = "tax")
y = as.data.frame(group_cpm)%>%
  select(c("NE.6.old","NC.6.old"))%>%
  rownames_to_column(var="tax")
y = inner_join(y,res2,by='tax')
y <- pivot_longer(y, cols = 2:3,
                  names_to = "group",
                  values_to = "value")
y$comparison <- 'NE.6.oldvsNC.6.old'
all.comp <- rbind(all.comp, y)

qlf <- glmQLFTest(fit, contrast=my.contrasts[,"NE.6.youngvsNC.6.young"])
test_NE.6.youngvsNC.6.young <- qlf$table%>%
  mutate_all(as.numeric)%>%
  filter(PValue < 0.05)
is.de <- decideTestsDGE(qlf)
summary(is.de)
topTags(qlf)
plotMD(qlf)
# Filter for features with significant adj p-value (FDR < 0.05)
res <- as.data.frame(topTags(qlf))%>%
 filter(FDR < 0.05)
#text(res$logCPM, res$logFC, labels = rownames(res), col = "black", cex = 0.5,
#       pos = 3)
is.de <- decideTestsDGE(qlf)
summary(is.de)
# ^ No DE species
# Plot species by comparison
res2 = res%>%
  rownames_to_column(var = "tax")
y = as.data.frame(group_cpm)%>%
  select(c("NE.6.young","NC.6.young"))%>%
  rownames_to_column(var="tax")
y = inner_join(y,res2,by='tax')
y <- pivot_longer(y, cols = 2:3,
                  names_to = "group",
                  values_to = "value")
y$comparison <- 'NE.6.youngvsNC.6.young'
all.comp <- rbind(all.comp, y)

qlf <- glmQLFTest(fit, contrast=my.contrasts[,"NE.0.youngvsNE.6.young"])
test_NE.0.youngvsNE.6.young <- qlf$table%>%
  mutate_all(as.numeric)%>%
  filter(PValue < 0.05)
is.de <- decideTestsDGE(qlf)
topTags(qlf)
plotMD(qlf)
# Filter for features with significant adj p-value (FDR < 0.05)
res <- as.data.frame(topTags(qlf, n=2, sort.by='p.value'))%>%
  filter(FDR < 0.05)
text(res$logCPM, res$logFC, labels = rownames(res), col = "black", cex = 0.5,
     pos = 3)
is.de <- decideTestsDGE(qlf)
summary(is.de)

# Plot species by comparison
res2 = res%>%
  rownames_to_column(var = "tax")
y = as.data.frame(group_cpm)%>%
  select(c("NE.0.young","NE.6.young"))%>%
  rownames_to_column(var="tax")
y = inner_join(y,res2,by='tax')
y <- pivot_longer(y, cols = 2:3,
                  names_to = "group",
                  values_to = "value")
y$comparison <- 'NE.0.youngvsNE.6.young'
all.comp <- rbind(all.comp, y)

qlf <- glmQLFTest(fit, contrast=my.contrasts[,"NC.0.youngvsNC.6.young"])
test_NC.0.youngvsNC.6.young <- qlf$table%>%
  mutate_all(as.numeric)%>%
  filter(PValue < 0.05)
is.de <- decideTestsDGE(qlf)
topTags(qlf)
plotMD(qlf)
# Filter for features with significant adj p-value (FDR < 0.05)
res <- as.data.frame(topTags(qlf, n = 2))%>%
  filter(FDR < 0.05)
#text(res$logCPM, res$logFC, labels = rownames(res), col = "black", cex = 0.5,
#     pos = 3)
is.de <- decideTestsDGE(qlf)
summary(is.de)

# Plot species by comparison
res2 = res%>%
  rownames_to_column(var = "tax")
y = as.data.frame(group_cpm)%>%
  select(c("NC.0.young","NC.6.young"))%>%
  rownames_to_column(var="tax")
y = inner_join(y,res2,by='tax')
y <- pivot_longer(y, cols = 2:3,
                  names_to = "group",
                  values_to = "value")
y$comparison <- 'NC.0.youngvsNC.6.young'
all.comp <- rbind(all.comp, y)

qlf <- glmQLFTest(fit, contrast=my.contrasts[,"NE.0.oldvsNC.0.old"])
test_NE.0.oldvsNC.0.old <- qlf$table%>%
  mutate_all(as.numeric)%>%
  filter(PValue < 0.05)
is.de <- decideTestsDGE(qlf)
summary(is.de)
topTags(qlf)
plotMD(qlf)
# Filter for features with significant adj p-value (FDR < 0.05)
res <- as.data.frame(topTags(qlf, n=1))%>%
  filter(FDR < 0.05)
text(res$logCPM, res$logFC, labels = rownames(res), col = "black", cex = 0.5,
     pos = 3)
is.de <- decideTestsDGE(qlf)
summary(is.de)

# Plot species by comparison
res2 = res%>%
  rownames_to_column(var = "tax")
y = as.data.frame(group_cpm)%>%
  select(c("NE.0.old","NC.0.old"))%>%
  rownames_to_column(var="tax")
y = inner_join(y,res2,by='tax')
y <- pivot_longer(y, cols = 2:3,
                  names_to = "group",
                  values_to = "value")
y$comparison <- 'NE.0.oldvsNC.0.old'
all.comp <- rbind(all.comp, y)

qlf <- glmQLFTest(fit, contrast=my.contrasts[,"NE.0.oldvsNE.6.old"])
test_NE.0.oldvsNE.6.old <- qlf$table%>%
  mutate_all(as.numeric)%>%
  filter(PValue < 0.05)
is.de <- decideTestsDGE(qlf)
topTags(qlf)
plotMD(qlf)
# Filter for features with significant adj p-value (FDR < 0.05)
res <- as.data.frame(topTags(qlf, n=12))%>%
  filter(FDR < 0.05)
#text(res$logCPM, res$logFC, labels = rownames(res), col = "black", cex = 0.5,
#     pos = 3)
is.de <- decideTestsDGE(qlf)
summary(is.de)

# Plot species by comparison
res2 = res%>%
  rownames_to_column(var = "tax")
y = as.data.frame(group_cpm)%>%
  select(c("NE.0.old","NE.6.old"))%>%
  rownames_to_column(var="tax")
y = inner_join(y,res2,by='tax')
y <- pivot_longer(y, cols = 2:3,
                  names_to = "group",
                  values_to = "value")
y$comparison <- 'NE.0.oldvsNE.6.old'
all.comp <- rbind(all.comp, y)


qlf <- glmQLFTest(fit, contrast=my.contrasts[,"NC.0.oldvsNC.6.old"])
test_NC.0.oldvsNC.6.old <- qlf$table%>%
  mutate_all(as.numeric)%>%
  filter(PValue < 0.05)
is.de <- decideTestsDGE(qlf)
topTags(qlf)
plotMD(qlf)
# Filter for features with significant adj p-value (FDR < 0.05)
res <- as.data.frame(topTags(qlf, n=1))%>%
  filter(FDR < 0.05)
#  text(res$logCPM, res$logFC, labels = rownames(res), col = "black", cex = 0.5,
#       pos = 3)
is.de <- decideTestsDGE(qlf)
summary(is.de)

# Plot species by comparison
res2 = res%>%
  rownames_to_column(var = "tax")
y = as.data.frame(group_cpm)%>%
  select(c("NC.0.old","NC.6.old"))%>%
  rownames_to_column(var="tax")
y = inner_join(y,res2,by='tax')
y <- pivot_longer(y, cols = 2:3,
                  names_to = "group",
                  values_to = "value")
y$comparison <- 'NC.0.oldvsNC.6.old'
all.comp <- rbind(all.comp, y)

qlf <- glmQLFTest(fit, contrast=my.contrasts[,"NE.0.youngvsNE.0.old"])
test_NE.0.youngvsNE.0.old <- qlf$table%>%
  mutate_all(as.numeric)%>%
  filter(PValue < 0.05)
is.de <- decideTestsDGE(qlf)
topTags(qlf, n=100)
plotMD(qlf)
# Filter for features with significant adj p-value (FDR < 0.05)
res <- as.data.frame(topTags(qlf, n=84))%>%
  filter(FDR < 0.05)
#text(res$logCPM, res$logFC, labels = rownames(res), col = "black", cex = 0.5,
#     pos = 3)
is.de <- decideTestsDGE(qlf)
summary(is.de)

# Plot species by comparison
res2 = res%>%
  rownames_to_column(var = "tax")
y = as.data.frame(group_cpm)%>%
  select(c("NE.0.young","NE.0.old"))%>%
  rownames_to_column(var="tax")
y = inner_join(y,res2,by='tax')
y <- pivot_longer(y, cols = 2:3,
                  names_to = "group",
                  values_to = "value")
y$comparison <- 'NE.0.youngvsNE.0.old'
all.comp <- rbind(all.comp, y)

qlf <- glmQLFTest(fit, contrast=my.contrasts[,"NE.6.youngvsNE.6.old"])
test_NE.6.youngvsNE.6.old <- qlf$table%>%
  mutate_all(as.numeric)%>%
  filter(PValue < 0.05)
is.de <- decideTestsDGE(qlf)
topTags(qlf)
plotMD(qlf)
# Filter for features with significant adj p-value (FDR < 0.05)
res <- as.data.frame(topTags(qlf, n=1))%>%
  filter(FDR < 0.05)
text(res$logCPM, res$logFC, labels = rownames(res), col = "black", cex = 0.5,
     pos = 3)
is.de <- decideTestsDGE(qlf)
summary(is.de)

# Plot species by comparison
res2 = res%>%
  rownames_to_column(var = "tax")
y = as.data.frame(group_cpm)%>%
  select(c("NE.6.young","NE.6.old"))%>%
  rownames_to_column(var="tax")
y = inner_join(y,res2,by='tax')
y <- pivot_longer(y, cols = 2:3,
                  names_to = "group",
                  values_to = "value")
y$comparison <- 'NE.6.youngvsNE.6.old'
all.comp <- rbind(all.comp, y)

qlf <- glmQLFTest(fit, contrast=my.contrasts[,"NC.0.youngvsNC.0.old"])
test_NC.0.youngvsNC.0.old <- qlf$table%>%
  mutate_all(as.numeric)%>%
  filter(PValue < 0.05)
is.de <- decideTestsDGE(qlf)
topTags(qlf)
plotMD(qlf)
# Filter for features with significant adj p-value (FDR < 0.05)
res <- as.data.frame(topTags(qlf, n = 28))%>%
  filter(FDR < 0.05)
text(res$logCPM, res$logFC, labels = rownames(res), col = "black", cex = 0.5,
     pos = 3)
is.de <- decideTestsDGE(qlf)
summary(is.de)

# Plot species by comparison
res2 = res%>%
  rownames_to_column(var = "tax")
y = as.data.frame(group_cpm)%>%
  select(c("NC.0.young","NC.0.old"))%>%
  rownames_to_column(var="tax")
y = inner_join(y,res2,by='tax')
y <- pivot_longer(y, cols = 2:3,
                  names_to = "group",
                  values_to = "value")
y$comparison <- 'NC.0.youngvsNC.0.old'
all.comp <- rbind(all.comp, y)

qlf <- glmQLFTest(fit, contrast=my.contrasts[,"NC.6.youngvsNC.6.old"])
test_NC.6.youngvsNC.6.old <- qlf$table%>%
  mutate_all(as.numeric)%>%
  filter(PValue < 0.05)
is.de <- decideTestsDGE(qlf)
topTags(qlf)
plotMD(qlf)
# Filter for features with significant adj p-value (FDR < 0.05)
res <- as.data.frame(topTags(qlf, n=6))%>%
  filter(FDR < 0.05)
#text(res$logCPM, res$logFC, labels = rownames(res), col = "black", cex = 0.5,
#     pos = 3)
is.de <- decideTestsDGE(qlf)
summary(is.de)

# Plot species by comparison
res2 = res%>%
  rownames_to_column(var = "tax")
y = as.data.frame(group_cpm)%>%
  select(c("NC.6.young","NC.6.old"))%>%
  rownames_to_column(var="tax")
y = inner_join(y,res2,by='tax')
y <- pivot_longer(y, cols = 2:3,
                  names_to = "group",
                  values_to = "value")
y$comparison <- 'NC.6.youngvsNC.6.old'
all.comp <- rbind(all.comp, y)

# -----------------------------------------------

library(gridExtra)



bar.plots <- all.comp%>%
  filter(comparison %in% c("NC.0.youngvsNC.6.young",
                           "NC.0.oldvsNC.0.old",
                           "NE.0.youngvsNE.6.young",
                           "NE.0.oldvsNE.6.old"))
# Plotting all comparisons 
ggplot(bar.plots, aes(x = group, y = value, fill = tax)) +
  # Implement a grouped bar chart
  geom_bar(position = "dodge", stat = "identity") + 
  facet_wrap(~comparison, scales = "free") + 
  theme_bw() + 
  ylab("counts per million (cpm)") + 
  guides(fill = guide_legend(ncol = 1))
# ----------------------------------------------
# Comparing the change in ex old to ex young and contrl old to control young







# ----------------------------------------------
# Below pools all the treatment groups together mainly just to look 
# At differences between ages initially at the 0 timepoint



# df2 is our design dataframe, df is our count dataframe
df4 <- df2
Group2 <- factor(paste(df4$Age,df4$Timepoint, sep="."))
Age2 <- factor(df2$Age)
Timepoint2 <- factor(df2$Timepoint)
cbind(df4, Group=Group, Treatment=Treatment, Timepoint=Timepoint, Age=Age)

design2 <- model.matrix(~0+Group2)
colnames(design2) <- levels(Group2)

d2 <- DGEList(counts = df3, group = Group2)
keep <- filterByExpr(d2, group=Group2)
d2 <- d2[keep, , keep.lib.sizes=FALSE]
d2 <- estimateDisp(d2)
d2 <- calcNormFactors(d2)

# MDS plot
cpm <- cpm(d2)%>%
  as.data.frame()
logcpm <- cpm(d2, log=TRUE)



plotMDS(d2, top = 1000, col = as.numeric(d2$samples$group))
fit2 <- glmQLFit(d2, design2)
tr2 <- glmTreat(fit2, coef = 2, lfc = 1)
topTags(tr2)

my.contrasts2 <- makeContrasts(
  
  old.0vsyoung.0 = old.0 - young.0,
  old.6vsyoung.6 = old.6 - young.6,
  
  
  
  
  levels = design2
  
)

qlf <- glmQLFTest(fit2, contrast=my.contrasts2[,"old.0vsyoung.0"])
is.de <- decideTestsDGE(qlf)
topTags(qlf)
plotMD(qlf)
# Filter for features with significant adj p-value (FDR < 0.05)
res <- as.data.frame(topTags(qlf, n=2))%>%
  filter(FDR < 0.05)
text(res$logCPM, res$logFC, labels = rownames(res), col = "black", cex = 0.5,
     pos = 3)
is.de <- decideTestsDGE(qlf)
summary(is.de)


qlf <- glmQLFTest(fit2, contrast=my.contrasts2[,"old.6vsyoung.6"])
is.de <- decideTestsDGE(qlf)
topTags(qlf)
plotMD(qlf)
summary(is.de)
# Filter for features with significant adj p-value (FDR < 0.05)
res <- as.data.frame(topTags(qlf, n=4))%>%
  filter(FDR < 0.05)
text(res$logCPM, res$logFC, labels = rownames(res), col = "black", cex = 0.5,
     pos = 3)
is.de <- decideTestsDGE(qlf)
summary(is.de)

# ---------------------------------
# Relative abundance plot



rel.abund.table <- cpm
rel.abund.table$phylum <- rownames(rel.abund.table)
rel.abund.table <- rel.abund.table[, c(ncol(rel.abund.table), 1:(ncol(rel.abund.table)-1))]
# Restore default row names
row.names(rel.abund.table) <- NULL
rel.abund.table <- rel.abund.table%>%
  pivot_longer(cols = -phylum, names_to = 'SampleID', values_to = 'cpm')

rel.abund.table <- rel.abund.table %>%
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
  ), SampleID = str_replace(SampleID, "_[0-9]$", ""))%>%
  left_join(metadata, by = "SampleID") %>%
  mutate(Age = ifelse(substr(SampleID, 3, 3) == "O", "old", 
                      ifelse(substr(SampleID, 3, 3) == "Y", "young", NA)))%>%
  mutate(Treatment_time = paste(Age, Treatment, Timepoint, sep = "_"))

x <- rel.abund.table %>%
  group_by(Treatment_time, phylum)%>%
  summarize(total_cpm = mean(cpm))%>%
  ungroup()%>%
  group_by(Treatment_time) %>%
  mutate(abundance = total_cpm/sum(total_cpm)*100)%>%
  ungroup()

ggplot(data = x,aes(x = Treatment_time, y = abundance, fill = phylum)) + 
  geom_bar(stat = "identity")+ theme_minimal()  +
  theme(axis.text.x = element_text(size = 12))


# -------------------------------------------------
# Trimming all.comp to use for excel heatmap

heat.df <- all.comp %>%
  select(c('tax', 'logFC', 'comparison'))%>%
  filter(comparison %in% c('NE.0.youngvsNE.6.young',
                           'NC.0.youngvsNC.6.young',
                           'NE.0.oldvsNE.6.old',
                           'NC.0.oldvsNC.6.old'))%>%
  distinct(logFC, .keep_all=T)%>%
  pivot_wider(names_from = comparison, values_from = logFC)%>%
  as.data.frame()

heat.df[is.na(heat.df)] <- 0
ylab <- heat.df[,1]
rownames(heat.df) <- heat.df[,1]
heat.df <- heat.df[,-1]%>%
  as.matrix()

heatmap(heat.df, labRow = ylab,
        Colv = NA,
        Rowv = NA,
        las = 2)
legend("topleft", legend = c("Low", "Medium", "High"), fill = heat.colors(4))




