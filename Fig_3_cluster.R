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
counts_raw<-data.frame(counts_raw = sample_data(ps.mi4d.asv)$counts)
counts_raw$ID2<-rownames(sample_imputed)
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
# Supp 3A dim: 20 15
mi4d.comp.core.rel%>%t()%>%
  Heatmap(name = "ASV abundance\n (clr-transformed)",
          cluster_columns=mi4d.sample.dend,cluster_rows = mi4d.asv.dend, 
          column_title = "sample ID",
          top_annotation = col_ha,
          row_title = "ASVs", column_title_side = "bottom",
          row_names_gp = gpar(fontsize = 6),
          col = rev(hcl.colors(10, "RdYlBu")))

# clustering at Genus level
# Fig 3A dim: 12 12
mi4d.comp.core.genus<-mi4d.comp.core%>% tax_fix(unknowns = c("bacterium", "caccae", "finegoldii", "hominis", "intestinalis", "massiliensis", "sanguinis", "stercoris"))%>%
  tax_agg(.,rank = "Genus")

mi4d.comp.core.genus%>%
  tax_transform("clr")%>%otu_table()%>%data.frame()%>%
  Heatmap(name = "Genus abundance\n (clr-transformed)",
          cluster_columns=mi4d.sample.dend,
          column_title = "sample ID",
          top_annotation = col_ha,
          row_title = "Genera", column_title_side = "bottom",
          row_names_gp = gpar(fontsize = 10),
          col = rev(hcl.colors(10, "RdYlBu")))




# Supp 3B dbRDA 
source("Fig1_2_Mb_formal.R") # disregard warnings
rda.axes$cluster<-meta$cluster
table_sex_cluster<-with(meta,table(Sex,cluster))
rda.axes%>%ggplot(aes(x= CAP1, y = CAP2))+geom_point(aes(color = cluster,shape = Sex),size = 3.5)+
  stat_ellipse(aes(color = cluster),linetype = 2,type = "t", level = .95)+
  scale_color_manual(values = ghibli_palette("LaputaMedium")[7:5])+
  labs(x = "dbRDA1(5.7%)",y = "dbRDA2(4.0%)")+scale_x_continuous(limits = c(-.4,.4)) +
  annotation_custom(gridExtra::tableGrob(table_sex_cluster), xmin=.05, xmax=.4, ymin=.38, ymax=.45)


# Supp 3C absolute abundance of Pseudomonadota
mi4d.comp.core.phy<-mi4d.comp.core%>%tax_fix()%>%tax_agg(.,rank = "Phylum")
core.phy.abundance<-data.frame(otu_table(mi4d.comp.core.phy))
 
core.phy.abundance_freq<-apply(core.phy.abundance,2,function(x) x/sum(x))
proteo_abund_abs<-core.phy.abundance_freq%>%data.frame()%>%
  rownames_to_column("Phylum")%>%
  pivot_longer(.,!Phylum,names_to = "ID2",values_to = "read")%>%
  filter(Phylum == "Proteobacteria")%>%merge(counts_raw)%>%left_join(meta)%>%
  mutate(absabun = read*counts_raw)
summary(aov(absabun~cluster,data = proteo_abund_abs))
proteo_abund_abs%>%
  ggplot(.,aes(x =cluster, y =absabun, color = cluster))+geom_boxplot(outlier.color = NA)+geom_point(aes(shape = Sex),size = 3,position = position_jitterdodge(.5))+
  labs(x = "Hierarchical Cluster", y ="Pseudomonodota(Proteobacteria)\n absolute abundance")+scale_color_manual(values = ghibli_palette("LaputaMedium")[7:5], name = NULL)



# Fig 3C
core.phy.abundance_long<-core.phy.abundance%>%rownames_to_column("Phylum")%>%
  pivot_longer(.,!Phylum,names_to = "ID2",values_to = "read")%>%
  left_join(meta)
proteo_abund<-core.phy.abundance_long%>%filter(Phylum == "Proteobacteria")
proteo_abund%>%ggplot(.,aes(x =cluster, y =read, color = cluster))+geom_boxplot(outlier.color = NA)+geom_point(aes(shape = Sex),size = 3,position = position_jitterdodge(.5))+
  labs(x = "Hierarchical Cluster", y ="Pseudomonadota(Proteobacteria)\n relative abundance")+scale_color_manual(values = ghibli_palette("LaputaMedium")[7:5], name = NULL)

anova_test<-aov(read ~ cluster,data = proteo_abund)
summary(anova_test)
TukeyHSD(anova_test)

# for reviewer
mi4d.comp.core.genus.abs<-
  mi4d.comp.core.genus%>%transform_sample_counts(., function(x) x /sum(x))%>%
  otu_table()%>%data.frame()%>%
  rownames_to_column("Genus")%>%
  pivot_longer(.,!Genus,names_to = "ID2",values_to = "read")%>%merge(counts_raw)%>%left_join(meta)%>%
  mutate(absabun = read*counts_raw)
   
mi4d.comp.core.genus.abs%>%filter(Genus == "Escherichia-Shigella")%>%
  ggplot(.,aes(x =cluster, y =log1p(absabun), color = cluster))+geom_boxplot(outlier.color = NA)+
  geom_point(aes(shape = Sex),size = 3,position = position_jitterdodge(.5))+
  labs(x = "Hierarchical Cluster", y ="Escherichia-Shigella\n absolute abundance (log-transformed)")+scale_color_manual(values = ghibli_palette("LaputaMedium")[7:5], name = NULL)
summary(aov(absabun~cluster,data = mi4d.comp.core.genus.abs%>%filter(Genus == "Escherichia-Shigella")))

genera<-mi4d.comp.core%>%tax_fix()%>%tax_agg(.,rank = "Genus")%>%
  otu_table()
genera_long<-genera%>%data.frame()%>%rownames_to_column("ASV")%>%
  pivot_longer(.,!ASV,names_to = "ID2",values_to = "read")%>%
  left_join(meta)

genera_long_sum<-genera_long%>%group_by(ASV,cluster)%>%
  summarise(n = sum(read),p = sum(read>0))%>%
  mutate(prevalence = case_when(cluster == "cluster1"~ p/13,
                                cluster == "cluster2"~ p/15,
                                cluster == "cluster3"~ p/25))
# test for prevalence ass'n with cluster
prevalent_ASV<-genera_long%>%group_by(ASV)%>%
  summarise(n = sum(read),p = sum(read>0))%>%filter(p>=40)

rare_genera<-genera_long%>%filter(read==0)%>%group_by(ASV,cluster)%>%tally()%>%mutate(prevalence_threshold = case_when(cluster == "cluster1"~ floor(13*.1),
                                                                                                                       cluster == "cluster2"~ floor(15*.1),
                                                                                                                       cluster == "cluster3"~ floor(25*.1)))%>%
  filter(n>prevalence_threshold)

rare_genera_full<-rare_genera%>%select(-prevalence_threshold)%>%
  full_join(crossing(ASV =genera_long_sum$ASV, cluster = genera_long_sum$cluster))
# chi-square test for whether missing genera is uniformly distributed across cluster
rare_genera_full[is.na(rare_genera_full$n),"n"]<-0
p_vals<-rare_genera_full%>%filter(!(ASV%in%prevalent_ASV$ASV))%>%
  split(.$ASV)%>%map_dbl(~chisq.test(x=.$n, p = c(13,15,25), rescale.p = TRUE)$p.value)%>%data.frame(p=.)
chisq_result<-p_vals%>%mutate(p.adj = p.adjust(p,n=32))%>%filter(p.adj<.05)%>%
  rownames_to_column("ASV")%>%
  mutate(cluster = "cluster2",p.adj2 = format(signif(p.adj, 3),scientific = TRUE))

# Supp Fig 3E 15*15?
genera_long_sum%>%filter(!(ASV%in%prevalent_ASV$ASV))%>%
  ggplot(aes(x = cluster, y =ASV))+geom_tile(aes(fill = prevalence,alpha = prevalence<.5))+
  scale_fill_viridis_c()+scale_alpha_manual(values =c(.5,1),guide = "none")+
  geom_text(data = chisq_result, size = 3,
            color = "black",aes(label = paste("p.adj=",p.adj2)))+labs(y = "Genus", x = "")

# test for abundance ass'n with cluster
filtered_data<-genera_long%>%filter(read>0 & ASV%in%prevalent_ASV$ASV)%>%
  nest_by(ASV)
custom_control <- glm.control(maxit = 200, epsilon = 1e-8, trace = TRUE)
glm_test<-filtered_data%>%
  mutate(glm_test = list(
    tryCatch(
      glm(read ~ Age + cluster + Sex+BMI, family = Gamma(link = "log"), data = data, control = custom_control),
      error = function(e) NA
    )
  ))
glm_test_result<-glm_test%>%summarise(broom::tidy(glm_test))%>%
  filter(grepl("cluster",term))%>%mutate(p.adj = p.adjust(p.value,n=38))%>%filter(p.adj<.05)%>%
  mutate(cluster = gsub("clustercluster","cluster",term),p.adj2 = format(signif(p.adj, 3),scientific = TRUE),
         symbol = case_when(p.adj < .001 ~ "***", p.adj < .01 ~ "**", p.adj < .05 ~ "*"))

# Supp Fig 3D dim 12 * 20
genera_long%>%filter(ASV%in%prevalent_ASV$ASV)%>%select(ASV,ID2,read, cluster)%>%
  # left_join(glm_test_result)%>%
  ggplot(aes(y = ASV,x = log1p(read), color = cluster,group = interaction(ASV,cluster)))+
  geom_boxplot()+geom_point(position = position_jitterdodge())+scale_color_manual(values = ghibli_palette("LaputaMedium")[7:5], name = NULL)+
  geom_text(data = glm_test_result%>%filter(cluster=="cluster2"), size = 7,fontface = "bold",
            x = 6,color = ghibli_palette("LaputaDark")[6],aes(label = symbol))+
  geom_text(data = glm_test_result%>%filter(cluster=="cluster3"), size = 7,fontface = "bold",
            x = 6,vjust = -.5,color = ghibli_palette("LaputaDark")[5],aes(label = symbol))+
  labs(x = "Relative abundance \n(log1p-transformed)", y = "")


summary(proteo_lm<-glm(read~Calprotectin_fecal*Sex+BMI+Age+Calprotectin_plasma, family = Gamma(link = "log"),
                       data = genera_long%>%
                         filter(ASV =="Escherichia-Shigella"& read!=0)))

# Fig 3D
genera_long%>%filter(ASV =="Escherichia-Shigella"& read!=0)%>%
  ggplot(.,aes(x = Calprotectin_fecal, y = read,color= Sex, fill = Sex))+
  geom_smooth(method = "glm", method.args = list(family = "Gamma"))+geom_point(size = 3.5)+
  labs(x = "Fecal calprotectin (normalized)", y ="Escherichia reads (log transformed)")+
  scale_fill_manual(values = wes_palette("GrandBudapest2")[3:4])+
  scale_color_manual(values = wes_palette("GrandBudapest2")[3:4])+
  scale_y_continuous(trans = "log",labels = function(x) format(round(x,2)))

# Supp 3G
proteo_mod2<-glm(read~Calprotectin_fecal*Sex, family = Gamma(link = "log"),
                 data = genera_long%>%
                   filter(ASV =="Escherichia-Shigella"& read!=0))
genera_long%>%filter(ASV =="Escherichia-Shigella"& read!=0)%>%
  mutate(mod2.residuals=proteo_mod2$residuals)%>%
  ggplot(.,aes(x = Calprotectin_plasma, y = mod2.residuals))+
  geom_smooth(method = "lm")+geom_point(size = 3.5)+scale_y_continuous(limits = c(-2.5,6))+
  labs(x = "Plasma calprotectin (normalized)", y ="Escherichia-Shigella Genus read \n(Residuals from \nFecal Calprotectin & Sex)")

# Supp 3F
summary(lm(Calprotectin_fecal~Calprotectin_plasma+Age+Sex+BMI, data = meta))
meta%>%
  ggplot(aes(x = Calprotectin_fecal, y=Calprotectin_plasma))+
  geom_point(size = 3)+
  labs(x = "Fecal calprotectin (normalized)", y = "Plasma calprotectin (normalized)")



# cluster association with metadata ---------------------------------------

summary(mod1<-lm(HOMA.IR ~ BMI+Sex*cluster+Age,data = meta)) 
summary(mod2<-lm(Triglycerides ~ BMI+Sex*cluster+Age,data = meta)) 

# Fig 3E
meta%>%
  ggplot(.,aes(x =cluster, y = HOMA.IR, color = cluster,shape = Sex))+geom_boxplot(outlier.color = NA)+geom_point(size = 3,position = position_jitterdodge(.5))+
  # scale_x_discrete(breaks = c("cluster1","cluster3","cluster2"))+
  labs(x = "hierarchical Cluster", y ="HOMA-IR (normalized)")+scale_color_manual(values = ghibli_palette("LaputaMedium")[7:5], name = NULL)

# Fig 3F
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
new_data<-data.frame(counts = seq(min(meta$counts),max(meta$count),length.out = 100))
pred<-cbind(counts = new_data,predict(multi_mo,new_data,type="probs"))%>%data.frame()

# Fig 3B
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

# set.seed(42)
# fit1 <- sem(model = model1, data = meta2, rotation = "geomin")
fit1 <-readRDS("20240901_SEM_fit1.rds")
summary(fit1,fit.measures = TRUE)
# saveRDS(fit1,"20240901_SEM_fit1.rds")

