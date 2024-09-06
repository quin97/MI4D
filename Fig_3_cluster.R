library(tidyverse)
library(magrittr)
library(phyloseq)
library(ghibli)
library(dendextend)
library(ComplexHeatmap)
library(caret)
library(microViz)
library(vegan)
library(wesanderson)

theme_set(
  theme_classic(base_size = 20)
)
sample_imputed<-read.csv("samplemat_pbmc_scaled_imputed.csv")%>%
  filter(!ID%in%c("MI4D022", "MI4D028", "MI4D044"))%>%
  mutate(ID2 = if_else(grepl("PRE",ID),paste0("MI4D0",as.numeric(str_extract(ID,"\\d{3}?"))+52),ID),
         female = ifelse(Sex=="F",1,0),
         biomass= ifelse(counts >=median(counts),"High","Low"))%>%
  column_to_rownames("ID2")

ps.mi4d.asv<-readRDS("ps.mi4d2.rds")
ps.mi4d.asv.otu<-ps.mi4d.asv%>%otu_table()%>%data.frame()%>%as.matrix()
dimnames(ps.mi4d.asv.otu)[[2]]<-rownames(sample_imputed)
ps.mi4d.asv<-phyloseq(sample_data(sample_imputed),otu_table(ps.mi4d.asv.otu,taxa_are_rows=T),tax_table(ps.mi4d.asv))

mi4d.comp.core <- filter_taxa(ps.mi4d.asv, function(x) sum(x > 0) > nsamples(ps.mi4d.asv)*.1 & mean(x/sample_sums(ps.mi4d.asv),na.rm = T) > .0005, prune = T)


# hierarchical clustering
mi4d.comp.core.tax<-data.frame(tax_table(mi4d.comp.core))%>%mutate(Species = case_when(is.na(Genus) & is.na(Species) ~ paste(Family, "Genus"),
                                                                                       !(is.na(Genus)) & is.na(Species) ~ paste(Genus, "Sp."),
                                                                                       TRUE ~ paste(Genus, Species)))
species_name<-mi4d.comp.core.tax$Species

mi4d.comp.core.rel<-mi4d.comp.core%>%tax_transform("clr")%>%otu_table()%>%t()%>%data.frame()
colnames(mi4d.comp.core.rel)<-species_name
mi4d.sample.dist <- vegdist(mi4d.comp.core.rel,method ="euclidean")
mi4d.sample.hclust<-hclust(mi4d.sample.dist, method = "ward.D2")
mi4d.sample.dend<-as.dendrogram(mi4d.sample.hclust)

# clustering of ASVs
mi4d.asv.dist <-vegdist(t(mi4d.comp.core.rel),method ="euclidean")
mi4d.asv.hclust<-hclust(mi4d.asv.dist, method = "ward.D2")
mi4d.asv.dend<-as.dendrogram(mi4d.asv.hclust)

meta <- data.frame(phyloseq::sample_data(mi4d.comp.core)) %>%rownames_to_column("ID2")

mi4d.sample.hclust3<-cutree(mi4d.sample.hclust,h=80)
mi4d.sample.hclust3.2<-as.numeric(gsub(0,1,gsub(1, 2,gsub(2,0,mi4d.sample.hclust3)))) # switch 1 and 2 label for color-coding
table(mi4d.sample.hclust3)

meta$cluster<-paste0("cluster",mi4d.sample.hclust3)
sample_data(mi4d.comp.core)<-meta%>%column_to_rownames("ID2")


col_ha = columnAnnotation(Biomass = meta$biomass,Cluster = meta$cluster,
                          col = list(Cluster = c("cluster1"= "#F0D77B", "cluster2" = "#B4DAE5","cluster3"="#AE93BE"),
                                     Biomass = c("High" = "#2E9FDF", "Low" = "#E7B800")))
# fig dim: 20 15
mi4d.comp.core.rel%>%t()%>%# apply(.,2,function(x) x = x/sum(x))%>%
  Heatmap(name = "ASV abundance\n (clr-transformed)",
          cluster_columns=mi4d.sample.dend,cluster_rows = mi4d.asv.dend, 
          column_title = "sample ID",
          top_annotation = col_ha,
          row_title = "ASVs", column_title_side = "bottom",
          row_names_gp = gpar(fontsize = 6),
          col = rev(hcl.colors(10, "RdYlBu")))

# clustering at Genus level
# fig dim: 12 12
mi4d.comp.core%>% tax_fix(unknowns = c("bacterium", "caccae", "finegoldii", "hominis", "intestinalis", "massiliensis", "sanguinis", "stercoris"))%>%
  tax_agg(.,rank = "Genus")%>%
  tax_transform("clr")%>%otu_table()%>%data.frame()%>%
  Heatmap(name = "Genus abundance\n (clr-transformed)",
          cluster_columns=mi4d.sample.dend,
          column_title = "sample ID",
          top_annotation = col_ha,
          row_title = "Genera", column_title_side = "bottom",
          row_names_gp = gpar(fontsize = 10),
          col = rev(hcl.colors(10, "RdYlBu")))





# plot on dbRDA 
source("Fig1_2_Mb_formal.R")
rda.axes$cluster<-meta$cluster
rda.axes%>%ggplot(aes(x= CAP1, y = CAP2))+geom_point(aes(color = cluster,shape = biomass),size = 3.5)+
  stat_ellipse(aes(color = cluster),linetype = 2,type = "t", level = .95)+
  scale_color_manual(values = ghibli_palette("LaputaMedium")[7:5])+
  labs(x = "dbRDA1(5.7%)",y = "dbRDA2(4.0%)")+scale_x_continuous(limits = c(-.4,.4))


# plot composition in clusters
mi4d.comp.core.phy<-mi4d.comp.core%>%tax_fix()%>%tax_agg(.,rank = "Phylum")
core.phy.abundance<-data.frame(otu_table(mi4d.comp.core.phy))
  
core.phy.abundance_long<-core.phy.abundance%>%rownames_to_column("Phylum")%>%
  pivot_longer(.,!Phylum,names_to = "ID2",values_to = "read")%>%
  left_join(meta)

core.phy.abundance_long%>%
  group_by(cluster,Phylum)%>%summarise(n = sum(read))%>%
  group_by(cluster)%>%mutate(freq = n/sum(n))%>%
  ggplot(aes(x = cluster, y = freq, fill= Phylum))+
  geom_bar(stat = "identity")+
  scale_fill_ghibli_d(labels = c("Actinomycetota","Bacteroidota","Desulfobacterota","Bacillota","Pseudomonadota","Verrucomicrobiota"),
                      "SpiritedMedium",direction = -1)+
  labs(x = "cluster", y = "Relative abundance")

proteo_abund<-core.phy.abundance_long%>%filter(Phylum == "Proteobacteria")
proteo_abund%>%ggplot(.,aes(x =cluster, y =read, color = cluster))+geom_boxplot(outlier.color = NA)+geom_point(aes(shape = Sex),size = 3,position = position_jitterdodge(.5))+
  labs(x = "Hierarchical Cluster", y ="Pseudomonadota(Proteobacteria)\n relative abundance")+scale_color_manual(values = ghibli_palette("LaputaMedium")[7:5], name = NULL)

anova_test<-aov(read ~ cluster,data = proteo_abund)
summary(anova_test)
TukeyHSD(anova_test)


proteo<-mi4d.comp.core%>%subset_taxa( Phylum == "Proteobacteria")
proteo_long<-proteo%>%tax_fix()%>%tax_agg(.,rank = "Genus")%>%
  otu_table()%>%data.frame()%>%rownames_to_column("ASV")%>%
  pivot_longer(.,!ASV,names_to = "ID2",values_to = "read")%>%
  left_join(meta)
proteo_long%>%group_by(ASV)%>%summarise(n = sum(read),p = sum(read>0))
summary(proteo_lm<-glm(read~Calprotectin_fecal*Sex+BMI+Age+Calprotectin_plasma, family = Gamma(link = "log"),
                       data = proteo_long%>%
                         filter(ASV =="Escherichia-Shigella"& read!=0)))

proteo_long%>%filter(ASV =="Escherichia-Shigella"& read!=0)%>%
  ggplot(.,aes(x = Calprotectin_fecal, y = read,color= Sex, fill = Sex))+
  geom_smooth(method = "glm", method.args = list(family = "Gamma"))+geom_point(size = 3.5)+
  labs(x = "Fecal calprotectin (normalized)", y ="Escherichia reads (log transformed)")+
  scale_fill_manual(values = wes_palette("GrandBudapest2")[3:4])+
  scale_color_manual(values = wes_palette("GrandBudapest2")[3:4])+
  scale_y_continuous(trans = "log",labels = function(x) format(round(x,2)))

proteo_mod2<-glm(read~Calprotectin_fecal*Sex, family = Gamma(link = "log"),
                 data = proteo_long%>%
                   filter(ASV =="Escherichia-Shigella"& read!=0))
proteo_long%>%filter(ASV =="Escherichia-Shigella"& read!=0)%>%
  mutate(mod2.residuals=proteo_mod2$residuals)%>%
  ggplot(.,aes(x = Calprotectin_plasma, y = mod2.residuals))+
  geom_smooth(method = "lm")+geom_point(size = 3.5)+scale_y_continuous(limits = c(-2.5,6))+
  labs(x = "Plasma calprotectin (normalized)", y ="Escherichia-Shigella Genus read \n(Residuals from \nFecal Calprotectin & Sex)")


summary(lm(Calprotectin_fecal~Calprotectin_plasma+Age+Sex+BMI, data = meta))
meta%>%
  ggplot(aes(x = Calprotectin_fecal, y=Calprotectin_plasma))+
  geom_point(size = 3)+
  labs(x = "Fecal calprotectin (normalized)", y = "Plasma calprotectin (normalized)")



# cluster association with metadata ---------------------------------------

summary(mod1<-lm(HOMA.IR ~ BMI+Sex*cluster+Age,data = meta)) 
summary(mod2<-lm(Triglycerides ~ BMI+Sex*cluster+Age,data = meta)) 

meta%>%
  ggplot(.,aes(x =cluster, y = HOMA.IR, color = cluster,shape = Sex))+geom_boxplot(outlier.color = NA)+geom_point(size = 3,position = position_jitterdodge(.5))+
  # scale_x_discrete(breaks = c("cluster1","cluster3","cluster2"))+
  labs(x = "hierarchical Cluster", y ="HOMA-IR (normalized)")+scale_color_manual(values = ghibli_palette("LaputaMedium")[7:5], name = NULL)


meta%>%
  ggplot(.,aes(x =cluster, y = Triglycerides, color = cluster,shape = Sex))+geom_boxplot(outlier.color = NA)+geom_point(size = 3,position = position_jitterdodge(.5))+
  labs(x = "hierarchical Cluster", y ="Triglycerides (normalized)")+scale_color_manual(values = ghibli_palette("LaputaMedium")[7:5], name = NULL)



library(nnet)
meta$cluster2<-relevel(as.factor(meta$cluster), ref = 3)
multi_mo <- multinom(cluster2 ~ counts, data = meta,model=TRUE)
summary(multi_mo)
exp(coef(multi_mo))
(z <- summary(multi_mo)$coefficients/summary(multi_mo)$standard.errors)
(p <- (1 - pnorm(abs(z), 0, 1)) * 2)
# new_data<-data.frame(counts_raw = seq(min(meta$counts_raw),max(meta$counts_raw),length.out = 100))
new_data<-data.frame(counts = seq(min(meta$counts),max(meta$count),length.out = 100))
pred<-cbind(counts = new_data,predict(multi_mo,new_data,type="probs"))%>%data.frame()
pred%>%pivot_longer(cols = !counts, names_to = "cluster",values_to = "prob")%>%
  ggplot(aes(x = counts, y = prob,fill = cluster))+geom_area()+# coord_trans(x = "log10")+scale_x_continuous(breaks = 10^seq(10,12,0.5),labels = function(x) format(x, scientific = TRUE))+
  scale_fill_manual(values = ghibli_palette("LaputaMedium")[7:5],guide="none")+
  labs(x="Bacterial counts in dry stool (normalized)", y = "Probability of cluster")



# confirmatory factor analysis  --------------------------------------------

require(lavaan)

meta2<-meta%>%mutate(cluster1 = ifelse(cluster == "cluster1", 1,0),
                     cluster2 = ifelse(cluster == "cluster2", 1,0),
                     female = female+1)%>%
  select(-cluster)

require(Hmisc)
analytes<-meta2%>%dplyr::select(ELA2, Calprotectin_plasma, MPO, LBP, CD14)
res <- rcorr(as.matrix(analytes))
res2<-getCov(res)

model1 <- '
# efa block 
efa("efa1")*meta1+
efa("efa1")*meta2 =~ Triglycerides+LDL+HOMA.IR+CRPhs+Age+BMI
# efa block 2
efa("efa2")*inflammation =~ ELA2+CD14+LBP

# efa block 3
efa("efa3")*gut =~ cluster1+cluster2+Calprotectin_fecal+water_content


# regressions
meta1+meta2=~ gut+inflammation+counts*female

# residual correlations
Calprotectin_plasma ~~ ELA2+MPO
LBP ~~ CRPhs
# ALT ~~ AST
Cholesterol ~~ Triglycerides+LDL
HOMA.IR ~~ WBISI+IGI
counts ~~ water_content+Calprotectin_fecal
gut ~~ inflammation+counts*female
'

set.seed(42)
fit1 <- sem(model = model1, data = meta2, rotation = "geomin")
summary(fit1,fit.measures = TRUE)
# saveRDS(fit1,"20240901_SEM_fit1.rds")

