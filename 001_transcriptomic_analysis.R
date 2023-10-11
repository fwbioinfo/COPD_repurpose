
# 0. preperation ----
library(data.table)
library(tidyr) 
library(magrittr) 
library(GEOquery)
library(dplyr)

library(limma)
library(hgug4112a.db)
library(annotate)

library(ggpubr)
library(ggplot2)
library(ggvenn)

library(conflicted)
conflicts_prefer(dplyr::select); conflicts_prefer(dplyr::filter)
conflicts_prefer(base::intersect); conflicts_prefer(base::setdiff)

extractVar = function(keyword = '', df){
  ### used for columns not aligned properly, as in some GEO metadata
  ### search the whole data table and extract keywords per row
  ### output is a vector in original row order.
  
  ## test
  #keyword = 'gold'
  #df <- gse1.raw %>% select(`characteristics_ch1.2`:characteristics_ch1.11)
  
  result <- apply(df, 1, function(x) x[grep(keyword, x)] %>% unname())
  ## make result into vector
  result[sapply(result, function(x) identical(character(0), x))] <- NA
  result <- unlist(result)
}


# save.image(file = 'tmp.RData')
# load(file = 'tmp.RData')
#load(file = 'tmp-Aug-07-2023.RData')

###### 1. load GEO data by getGEO ##############
#Sys.setenv(VROOM_CONNECTION_SIZE = 5000072)
gse <- getGEO("GSE47460", GSEMatrix = TRUE)

### expression data ########
dim(gse[[1]]); dim(gse[[2]])
expr1 <- exprs(gse[[1]])
expr2 <- exprs(gse[[2]])
sum(row.names(expr1) != row.names(expr1))
expr.raw <- cbind(expr1, expr2)
expr.raw[1:3,1:3]
rm(expr1); rm(expr2)

### meta data ###########
show(gse)
gse1.raw <- pData(phenoData(gse[[1]]))

gse1 <- gse1.raw %>%
  mutate(sampleid = title,
         geo_accession = geo_accession,
         batch = 'GPL14550',
         disease = `disease state:ch1`,
         sex = characteristics_ch1.1,
         age = characteristics_ch1.2,
         gold = '',
         smoke = '') %>%
  dplyr::select(sampleid, geo_accession, batch, disease, sex, age, gold, smoke)

gse1$gold <- extractVar(keyword = 'gold', df = gse1.raw %>% select(`characteristics_ch1.2`:characteristics_ch1.11))
gse1$smoke <- extractVar(keyword = 'smoker', df = gse1.raw %>% select(`characteristics_ch1.2`:characteristics_ch1.11))


gse2.raw <- pData(phenoData(gse[[2]]))

gse2 <- gse2.raw %>%
  mutate(sampleid = title,
         geo_accession = geo_accession,
         batch = 'GPL6480',
         disease = `disease state:ch1`,
         sex = characteristics_ch1.1,
         age = characteristics_ch1.2,
         gold = '',
         smoke = '') %>%
  select(sampleid, geo_accession, batch, disease, sex, age, gold, smoke)

gse2$gold <- extractVar(keyword = 'gold', df = gse2.raw %>% select(`characteristics_ch1.2`:characteristics_ch1.11))
gse2$smoke <- extractVar(keyword = 'smoker', df = gse2.raw %>% select(`characteristics_ch1.2`:characteristics_ch1.11))



meta.raw <- rbind(gse1, gse2)
meta <- meta.raw %>% 
  mutate(Sample_ID = sapply(strsplit(sampleid, "_"), function(x) x[1]) ) %>% 
  ## remove ILD
  filter(disease != "Interstitial lung disease") %>%
  filter(!grepl( 'ILD', sampleid)) %>%
  ## removing no smoking status
  filter(!is.na(smoke)) %>%
  ## re-label control
  mutate(gold = ifelse(disease == 'Control', 'Control', gold)) %>%
  ## remove gold0 in COPD
  filter(! (disease == "Chronic Obstructive Lung Disease" & gold == 'gold stage: 0-At Risk' )) %>%
  ## remove never smoker in control
  filter(! (gold == 'Control' & smoke == 'smoker?: 3-Never')) %>%
  #### value format
  mutate( age = gsub('age: ','', age, fixed = T),
          age = as.numeric(age)) %>%
  mutate(sex = gsub("Sex: ", "", sex, fixed = T),
         sex = factor(sex),
         smoke = gsub("smoker?: ", '', smoke, fixed = T),
         smoke = factor(smoke),
         gold = gsub("gold stage: ", '', gold, fixed = T),
         gold = factor(gold)) %>%
  ### rename
  dplyr::rename(goldgroup = gold,
         Gender = sex,
         Age = age)
dim(meta)

### demographic info #######
meta %>% count(disease, goldgroup)
## Sex
meta %>% count(goldgroup, Gender) %>% tidyr::pivot_wider(., names_from = 'Gender', values_from = 'n')
fisher.test(table(meta[,c('Gender','goldgroup')]))
reporttools::pairwise.fisher.test(meta$Gender, meta$goldgroup, p.adjust.method = 'BH')


meta %>% count(goldgroup, smoke)  %>% tidyr::pivot_wider(., names_from = 'smoke', values_from = 'n', values_fill = 0)
test <- meta %>% count(goldgroup, smoke)  %>% tidyr::pivot_wider(., names_from = 'smoke', values_from = 'n', values_fill = 0)
fisher.test(table(meta[,c('smoke','goldgroup')]))
#fmsb::pairwise.fisher.test(test$`1-Current`,test$`2-Ever (>100)`, test$`3-Never`, p.adjust.method = 'BH')
reporttools::pairwise.fisher.test(meta$smoke, meta$goldgroup, p.adjust.method = 'BH')
#fisher.test(table(meta[,c('smoke','goldgroup')])) %>% fisher.test
#table(meta[,c('smoke','goldgroup')])  %>% RVAideMemoire::fisher.multcomp(., p.method = 'BH')

c(
fisher.test(table(meta[meta$goldgroup %in% c('2-Moderate COPD','4-Very Severe COPD'),c('smoke','goldgroup')]))$p.value,
fisher.test(table(meta[meta$goldgroup %in% c('1-Mild COPD','4-Very Severe COPD'),c('smoke','goldgroup')]))$p.value,
fisher.test(table(meta[meta$goldgroup %in% c('3-Severe COPD','4-Very Severe COPD'),c('smoke','goldgroup')]))$p.value,
fisher.test(table(meta[meta$goldgroup %in% c('Control','4-Very Severe COPD'),c('smoke','goldgroup')]))$p.value,

fisher.test(table(meta[meta$goldgroup %in% c('1-Mild COPD','3-Severe COPD'),c('smoke','goldgroup')]))$p.value,
fisher.test(table(meta[meta$goldgroup %in% c('2-Moderate COPD','3-Severe COPD'),c('smoke','goldgroup')]))$p.value,
fisher.test(table(meta[meta$goldgroup %in% c('Control','3-Severe COPD'),c('smoke','goldgroup')]))$p.value,

fisher.test(table(meta[meta$goldgroup %in% c('1-Mild COPD','2-Moderate COPD'),c('smoke','goldgroup')]))$p.value,
fisher.test(table(meta[meta$goldgroup %in% c('2-Moderate COPD','Control'),c('smoke','goldgroup')]))$p.value,

fisher.test(table(meta[meta$goldgroup %in% c('1-Mild COPD','Control'),c('smoke','goldgroup')]))$p.value

) %>% p.adjust(., method = 'BH')

meta %>% group_by(goldgroup) %>%
  summarise(n = n(),
            mean = mean(Age, na.rm = T),
            sd = sd(Age, na.rm = T),
            se = sd(Age, na.rm = T)/sqrt(n))
meta %>% compare_means(Age~goldgroup, data = ., method = 't.test', ref.group = 'Control')
aov(Age~goldgroup, data = meta) %>% summary()
aov(Age~goldgroup, data = meta) %>% TukeyHSD()
pairwise.t.test(meta$Age, meta$goldgroup, p.adjust.method = 'BH')

# 2. Differentially expressed genes -----
### limma on ex smoker ########
set.seed(seed = 2023)

df_limma <- meta %>% 
  filter(smoke == '2-Ever (>100)')
dim(meta); dim(df_limma)

expr <- expr.raw[,df_limma$geo_accession]

design <- model.matrix(~ 0 + goldgroup + Age + Gender + batch, data=df_limma)


head(design)
colnames(design) <- c("GOLD1","GOLD2",'GOLD3','GOLD4','Healthy_Control',
                      'Age','gender', 'batch')
head(design);dim(design)

# fit the linear model to the filtered expression set
fit <- lmFit(expr, design)

# GOLD 4
contrast.matrix <- makeContrasts(GOLD4_Control = GOLD4 - Healthy_Control,levels=design)
disease_fits <- contrasts.fit(fit, contrast.matrix)
disease_ebFit <- eBayes(disease_fits)

nrow(topTable(disease_ebFit,  number=Inf, p.value = 0.05,lfc=log2(2)))

sigG4g <- topTable(disease_ebFit,  number=Inf)
sigG4g$probe <- row.names(sigG4g)


# GOLD 3
contrast.matrix <- makeContrasts(GOLD3_Control = GOLD3 - Healthy_Control,levels=design)
disease_fits <- contrasts.fit(fit, contrast.matrix)
disease_ebFit <- eBayes(disease_fits)
sigG3g <- topTable(disease_ebFit,  number=Inf)
sigG3g$probe <- row.names(sigG3g)


# GOLD 2
contrast.matrix <- makeContrasts(GOLD2_Control = GOLD2 - Healthy_Control,levels=design)
disease_fits <- contrasts.fit(fit, contrast.matrix)
disease_ebFit <- eBayes(disease_fits)
sigG2g <- topTable(disease_ebFit,  number=Inf)
sigG2g$probe <- row.names(sigG2g)


# GOLD 1
contrast.matrix <- makeContrasts(GOLD1_Control = GOLD1 - Healthy_Control,levels=design)
disease_fits <- contrasts.fit(fit, contrast.matrix)
disease_ebFit <- eBayes(disease_fits)
sigG1g <- topTable(disease_ebFit,  number=Inf) # 16
sigG1g$probe <- row.names(sigG1g)

## sig genes
sigG4 <- sigG4g %>% filter(adj.P.Val < 0.05 & abs(logFC) > log2(2))  #168
sigG3 <- sigG3g %>% filter(adj.P.Val < 0.05 & abs(logFC) > log2(2)) #83
sigG2 <- sigG2g %>% filter(adj.P.Val < 0.05 & abs(logFC) > log2(2)) # 9
sigG1 <- sigG1g %>% filter(adj.P.Val < 0.05 & abs(logFC) > log2(2)) # 5

sigG1 %>% mutate(updown=sign(logFC)) %>% count(updown)
sigG2 %>% mutate(updown=sign(logFC)) %>% count(updown)
sigG3 %>% mutate(updown=sign(logFC)) %>% count(updown)
sigG4 %>% mutate(updown=sign(logFC)) %>% count(updown)

## explore intersection
ggvenn(list(g4=row.names(sigG4),g3=row.names(sigG3),g2=row.names(sigG2),g1=row.names(sigG1)),
               show_percentage = F)
ggvenn(list(g4=row.names(sigG4),
            g3=row.names(sigG3),
            #g1=row.names(sigG1) ,
            g2=row.names(sigG2)
            ),
       show_percentage = F)
ggvenn(list(g4=row.names(sigG4),g3=row.names(sigG3),g2=row.names(sigG2)), show_percentage = F)
ggvenn(list(g4=row.names(sigG4),g3=row.names(sigG3),g1=row.names(sigG1)), show_percentage = F)
ggvenn(list(g2=row.names(sigG2),g3=row.names(sigG3),g1=row.names(sigG1)), show_percentage = F)
ggvenn(list(g2=row.names(sigG2),g3=row.names(sigG3),g4=row.names(sigG4)), show_percentage = F)
intersect(row.names(sigG1), row.names(sigG3))



#### Alluvial plot: genes across stages ##########
#install.packages("ggalluvial")
library(ggplot2)
library(ggalluvial)

sigG4g$probe <- row.names(sigG4g)
sigG3g$probe <- row.names(sigG3g)
sigG2g$probe <- row.names(sigG2g)
sigG1g$probe <- row.names(sigG1g)

df_alv <- rbind(sigG4g,sigG3g,sigG2g,sigG1g)
dim(df_alv)
df_alv$gene_count <- 1
df_alv$GOLD_stage <- c(rep('GOLD IV',15261),rep('GOLD III',15261),rep('GOLD II',15261),rep('GOLD I',15261))
head(df_alv)
df_alv <- df_alv %>%
  mutate(dechange = case_when(adj.P.Val<0.05 & df_alv$logFC > log2(2) ~ 'Up',
                              adj.P.Val<0.05 & df_alv$logFC < -log2(2) ~ 'Down',
                              T ~ 'Not DEG'))

table(df_alv[1:15261,]$dechange)
str(df_alv)
df_alv$dechange <- as.factor(df_alv$dechange)
head(df_alv)

### excluding genes not changing (prep 1st table changed genes)
nc4 <- subset(df_alv,GOLD_stage =='GOLD IV' & dechange == 'Not DEG')$probe
nc3 <- subset(df_alv,GOLD_stage =='GOLD III' & dechange == 'Not DEG')$probe
nc2 <- subset(df_alv,GOLD_stage =='GOLD II' & dechange == 'Not DEG')$probe
nc1 <- subset(df_alv,GOLD_stage =='GOLD I' & dechange == 'Not DEG')$probe

nc1234 <- intersect(nc1,intersect(nc2,intersect(nc3,nc4)))
length(nc1234) 
head(nc1234)

df_alv_changed <- subset(df_alv, ! probe  %in% nc1234)
dim(df_alv_changed)


table(df_alv_changed[,c("GOLD_stage",'dechange')])


# order: up, NDE, down
str(df_alv)
df_alv_changed$dechange <- factor(df_alv_changed$dechange,levels=c('Not DEG','Up','Down'))
levels(df_alv_changed$dechange)
table(df_alv_changed[,c("GOLD_stage",'dechange')])

ggplot(df_alv_changed,
       aes(x = GOLD_stage, stratum = dechange, alluvium = probe,
           y = gene_count,
           fill = dechange, label = dechange)) +
  scale_x_discrete(expand = c(.1, .1)) +
  geom_flow(stat = "alluvium", lode.guidance = "frontback") +
  #scale_fill_brewer(palette = "RdYlBu") +
  geom_stratum(alpha = .5) +
  theme_classic() +
  geom_text(stat = "stratum", size = 10) +
  theme(legend.position = "none", 
        axis.text.x = element_text(size=20), axis.text.y=element_text(size=20),
        axis.title.x=element_text(size=20), axis.title.y=element_text(size=20)) +
  labs(x='')
ggsave('output/acrossStage_FC2.png', height = 5, width = 7)

ggplot(df_alv_changed,
       aes(x = GOLD_stage, stratum = dechange, alluvium = probe,
           y = gene_count,
           fill = dechange, label = dechange)) +
  scale_x_discrete(expand = c(.1, .1)) +
  geom_flow(stat = "alluvium", lode.guidance = "frontback") +
  geom_stratum(alpha = .5) +
  theme_classic() +
  #geom_text(stat = "stratum", size = 10) +
  scale_fill_manual( values = c('grey','Salmon','deepskyblue')) +
  theme(legend.position = "none", 
        axis.text.x = element_text(size=16), axis.text.y=element_text(size=16),
        axis.title.x=element_text(size=16), axis.title.y=element_text(size=16)) +
  labs(x='')
ggsave('output/acrossStage_FC2_nolabel.png', height = 5.5, width = 7)


### intersect gene names
#gene.symbols <- getSYMBOL(row.names(sigG1g),"hgug4112a.db")
ref.genes <- pData(featureData(gse[[1]]))

intersect(row.names(sigG1), row.names(sigG3))
mygenes <- intersect(row.names(sigG2), intersect(row.names(sigG4), row.names(sigG3)))

sigG1234 <- bind_rows(sigG1 %>% mutate(group = 'g1'),
                      sigG2 %>% mutate(group = 'g2')) %>%
  bind_rows(., sigG3 %>% mutate(group = 'g3')) %>%
  bind_rows(., sigG4 %>% mutate(group = 'g4'))
row.names(sigG1234) <- NULL
sigG1234 %>% filter(probe %in% mygenes) %>% arrange(probe) %>% 
  mutate(updown = sign(logFC)) %>% count(probe, updown)



sigG1234g <- bind_rows(sigG1g %>% mutate(group = 'g1'),
                      sigG2g %>% mutate(group = 'g2')) %>%
  bind_rows(., sigG3g %>% mutate(group = 'g3')) %>%
  bind_rows(., sigG4g %>% mutate(group = 'g4')) %>%
  left_join(., ref.genes %>% select(probe = ID, GENE_SYMBOL))

row.names(sigG1234g) <- NULL

sigG1234g %>% filter(probe %in% mygenes) %>% arrange(probe) %>%
  mutate(FC = ifelse(adj.P.Val>0.05, 'adj.p>0.05', 2^logFC )) %>% View()

### saving sig genes file
output <- sigG1234 %>%
  left_join(., ref.genes %>% select(probe = ID, GENE_SYMBOL, EntrezID = GENE))
write.csv(output, './COPD_transcriptomics_sig_genes.csv', row.names = F)

########## pathway  #######

# psig1.2G4 <- row.names(subset(sigG4g,adj.P.Val <0.05 & abs(logFC) > log2(1.5))) #737
# psig1.2G3 <- row.names(subset(sigG3g,adj.P.Val <0.05 & abs(logFC) > log2(1.5))) #307
# psig1.2G2 <- row.names(subset(sigG2g,adj.P.Val <0.05 & abs(logFC) > log2(1.5))) #60
# psig1.2G1 <- row.names(subset(sigG1g,adj.P.Val <0.05 & abs(logFC) > log2(1.5))) #5


g4 <- gene.symbols[row.names(subset(sigG4g,adj.P.Val <0.05 & abs(logFC) > log2(2)))]
g3 <- gene.symbols[row.names(subset(sigG3g,adj.P.Val <0.05 & abs(logFC) > log2(2)))]
g2 <- gene.symbols[row.names(subset(sigG2g,adj.P.Val <0.05 & abs(logFC) > log2(2)))]
g1 <- gene.symbols[row.names(subset(sigG1g,adj.P.Val <0.05 & abs(logFC) > log2(2)))]


library(clusterProfiler)
library(enrichplot)

g4_entrezid <- AnnotationDbi::select(org.Hs.eg.db, keys=g4, columns=c("SYMBOL", "ENTREZID"),keytype='SYMBOL')
gog4 <- enrichGO(g4_entrezid$ENTREZID, OrgDb = "org.Hs.eg.db", ont="BP", readable=TRUE)

#write.csv(g4_entrezid,'output/g4_723entrezid.csv',row.names=F)

# g3_entrezid <- AnnotationDbi::select(org.Hs.eg.db, keys=g3, columns=c("SYMBOL", "ENTREZID"),keytype='SYMBOL')
# gog3 <- enrichGO(g3_entrezid$ENTREZID, OrgDb = "org.Hs.eg.db", ont="BP", readable=TRUE)
# 
# g2_entrezid <- AnnotationDbi::select(org.Hs.eg.db, keys=g2, columns=c("SYMBOL", "ENTREZID"),keytype='SYMBOL')
# gog2 <- enrichGO(g2_entrezid$ENTREZID, OrgDb = "org.Hs.eg.db", ont="BP", readable=TRUE)
# 
# g1_entrezid <- AnnotationDbi::select(org.Hs.eg.db, keys=g1, columns=c("SYMBOL", "ENTREZID"),keytype='SYMBOL')
# gog1 <- enrichGO(g1_entrezid$ENTREZID, OrgDb = "org.Hs.eg.db", ont="BP", readable=TRUE)
# 
# 
# ggvenn(list(gog4$ID,gog3$ID,gog2$ID,gog1$ID) %>% setNames(., c(4:1)))

dim(gog4)
gog4
summary(gog4) %>% filter(p.adjust<0.05) %>% nrow()

barplot(gog4, showCategory=20, xlab='gene count') +
 scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 90)) +
  #theme(legend.position = 'top') +
  NULL
ggsave('output/g4.gopathway.png', width = 6, height = 4)

#write.table(gog4,'157genes_clusterProfilerEnrichment.txt',row.names=F,quote=F,sep='\t')


### volcano plot #################
# Make a basic volcano plot
results <- sigG4g %>%
  left_join(., ref.genes %>% select(probe = ID, gene = GENE_SYMBOL)) %>%
  mutate(DEG = case_when(logFC>=1& adj.P.Val<0.05 ~ 'Up',
                         logFC<=-1& adj.P.Val<0.05 ~ 'Down',
                         T ~ 'Not Sig')) %>%
  mutate(p.rank = dense_rank(adj.P.Val),
         plot.label = ifelse(p.rank < 30 & DEG != 'Not Sig', gene, NA))

# with(results, plot(logFC, -log10(adj.P.Val),
#                    xlab="log2 Fold Change", ylab="-log10(adjusted p-value)",main="GOLD 4 vs Control",
#                    pch=20,xlim=c(-2.5,2.5), cex.lab=1.5,cex.axis=1.5,cex=1.5))
# abline(h=-log10(0.05),v=c(-1,1),lty=2)
# text(-2,-log10(0.05) + 0.3 ,labels = 'adjusted p = 0.05')
# 
# # Add colored points: red if padj<0.05, orange of log2FC>1, green if both)
# with(subset(results, logFC>=1& adj.P.Val<0.05), points(logFC, -log10(adj.P.Val), pch=20, col="red"))
# with(subset(results, logFC<=-1& adj.P.Val<0.05), points(logFC, -log10(adj.P.Val), pch=20, col="blue"))
# 
# 
# # Label points with the textxy function from the calibrate plot
# library(calibrate)
# dev.off()
# 
# with(subset(results, -log10(adj.P.Val)>8.5 & abs(logFC)>=1.38), textxy(logFC, -log10(adj.P.Val), labs=gene, cex=1,offset=0.5))
# 
# subset(results, -log10(adj.P.Val)>8.5 & abs(logFC)>=1.38)$gene

ggplot(data=results, aes(x=logFC, y=-log10(adj.P.Val), col=DEG, label=plot.label)) +
  geom_point() + 
  theme_classic() +
  geom_text_repel(show.legend = FALSE) +
  scale_color_manual(values=c("blue", "black", "red")) +
  geom_vline(xintercept=c(-1, 1), col="grey", linetype = 'dashed') +
  geom_hline(yintercept=-log10(0.05), col="grey", linetype = 'dashed') +
  theme(legend.position = 'top')
ggsave('output/volcano.g4.png', width = 5, height = 5)
