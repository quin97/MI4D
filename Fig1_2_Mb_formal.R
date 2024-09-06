suppressPackageStartupMessages(library(microViz))
library(tidyverse)
library(magrittr)
library(phyloseq)
library(ComplexHeatmap)
library(caret)
library(ghibli)
library(wesanderson)
library(vegan)
theme_set(
  theme_classic(base_size = 25)
)

tax.levels<-c("Kingdom", "Phylum", "Class", "Order","Family", "Genus","Species")

ps.mi4d.asv<-readRDS("ps.mi4d2.rds")

ps.mi4d.asv%>%tax_fix()%>%
  phyloseq_validate()

ps.mi4d.asv%<>%
  ps_mutate(
    biomass = if_else(counts >= 5.9e10, true = "High", false = "Low"),
    female = if_else(Sex == "F", true = 1, false = -1),
    ID2 = if_else(grepl("PRE",ID),paste0("MI4D0",as.numeric(str_extract(ID,"\\d{3}?"))+52),ID)
  )

# histogram of mean relative abundance
mean.abund<-ps.mi4d.asv%>%transform_sample_counts(., function(x) x / sum(x)*100)%>% # normalize total read of each individual to 1
  otu_table()%>%data.frame()
mean.ra<-mean.abund%>%reframe(ASV = rownames(mean.abund),mean.ra =rowMeans(.))%>%mutate(mean.ra.log10 = log10(mean.ra))
summary(mean.ra)
ggplot(mean.ra,aes(x = mean.ra.log10))+geom_histogram(fill = "slateblue",bins = 20)+
  geom_vline(xintercept = quantile(mean.ra$mean.ra.log10,c(0.25,.5,.75)),linetype = "dashed")+
  labs(x = "Mean relative abundance in samples\n(log10-transformed)", y = "Number of ASVs")

# rare fraction curve 15*6.18
source("Fig_1_rarefaction_fun.r")
rareplot <- ggrare(ps.mi4d.asv, step = 100, label = "Sample", color = "ID2", 
                   plot = FALSE, title = "MI4D Rarefaction Curve", parallel = TRUE, se = FALSE) 
(rareplot = rareplot + labs(y = "Number of unique ASVs")+ scale_color_viridis_d()+
  guides(color = guide_legend(override.aes = list(size=4))))

adiv<-estimate_richness(ps.mi4d.asv, measures = c("Shannon","Chao1"))
adiv$ID<-sample_data(ps.mi4d.asv)$ID
adiv2<-sample_data(ps.mi4d.asv)%>%data.frame()%>%merge(adiv)%>%column_to_rownames("ID2")
test<-ps.mi4d.asv%>%otu_table()%>%data.frame()%>%as.matrix()%>%t()
ps.mi4d.asv.tax<-tax_table(ps.mi4d.asv)
ps.mi4d.asv.tax2<-ps.mi4d.asv.tax%>%data.frame()%>%mutate(Phylum = 
                                             case_when(Phylum == "Firmicutes" ~ "Bacillota",
                                              Phylum == "Actinobacteriota" ~ "Actinomycetota",
                                              Phylum == "Proteobacteria" ~ "Pseudomonadota",
                                              Phylum == "Desulfobacterota" ~ "Pseudomonadota",
                                              Phylum == "Cyanobacteria" ~ "Cyanobacteriota",
                                              Phylum == "Campylobacterota" ~ "Pseudomonadota",
                                              Phylum == "Euryarchaeota" ~ "Methanobacteriota",
                                              is.na(Phylum)  ~ "Bacteria Kingdom",
                                              TRUE ~ as.character(Phylum)))
ps.mi4d.asv.tax2<-tax_table(ps.mi4d.asv.tax2)
colnames(ps.mi4d.asv.tax2)<-tax.levels
taxa_names(ps.mi4d.asv.tax2)<-taxa_names(ps.mi4d.asv.tax)
dimnames(test)[[1]]<-rownames(adiv2)
ps.mi4d.asv2<-phyloseq(sample_data(adiv2),otu_table(test,taxa_are_rows = F),ps.mi4d.asv.tax2)

# bar plot, absolute abundance
# for relative abundance, comment out 'tax_transform_for_plot = "identity"'
adiv2%<>%dplyr::arrange(counts)
sample_order<-rownames(adiv2)
comp.otu<-ps.mi4d.asv2%>%transform_sample_counts(function(x){x/sum(x)})%>%otu_table()%>%data.frame()%>%as.matrix()
abs.otu<-sample_data(ps.mi4d.asv2)[rownames(comp.otu),]$counts*comp.otu
rownames(abs.otu)<-rownames(comp.otu)
ps.mi4d.asv.abs<-ps.mi4d.asv2
otu_table(ps.mi4d.asv.abs)<-otu_table(abs.otu,taxa_are_rows = F)
# Fig dim inches: 15, 7
ps.mi4d.asv.abs%>%
  tax_fix()%>%
  comp_barplot(
    sample_order = sample_order,
    tax_level = "Phylum",n_taxa = 15, 
    palette = colorRampPalette(c(ghibli_palette("LaputaMedium")[7:3],ghibli_palette("KikiMedium")[4:7],wes_palette("Darjeeling2")))(16),
    other_name = "Other",# tax_transform_for_plot = "identity",# for absolute abundance
    merge_other = F, bar_outline_colour = "black" 
  ) +
   # scale_y_continuous(breaks = c(3, 50, 100, 200, 400)*10**9)+
  labs(x = NULL, y = "Relative abundance") + # "Bacterial counts/g dry stool"
  theme(# axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


adiv_ppc<-preProcess(adiv2%>%select(-Sex,-se.chao1),method = c("center", "scale", "YeoJohnson"))
adiv_ppc$method
adiv_scaled<-predict(adiv_ppc,adiv2)


# histograms for bacteria counts and water content
# Fig dim inches: 10, 6.18
p<-adiv2%>%mutate(log_biomass = log10(counts))%>%
  ggplot(aes(x = log_biomass))+geom_density()
pg <- ggplot_build(p)
pg_data<-data.frame(pg$data[[1]],stringsAsFactors = F)
ggplot(data=pg_data,aes(x,y))+geom_line()+geom_area(aes(fill = x<log10(5.92e+10)), alpha = 0.5)+
  scale_fill_manual(values = c("#2E9FDF","#E7B800"),name = "Biomass",labels = c("High","Low"))+geom_vline(xintercept = log10(5.92e+10),linetype = "dashed",size = 1)+
  labs(x="Bacterial counts in dry stool \n(log10-transformed)", y ="Proportion of participants")

fit_shannon<-summary(lm(Shannon~counts+Sex+Age+BMI,data = adiv_scaled))
fit_chao1<-summary(lm(Chao1~counts+Sex+Age+BMI,data = adiv_scaled))
fit_cal_fecal<-summary(lm(water_content~Calprotectin_fecal+counts+Sex+Age+BMI,data = adiv_scaled))

# Fig dim inches: 9.2, 5.8
adiv_scaled%>%ggplot(.,aes(x=counts,y=Shannon))+
  geom_point(aes(color=Sex),size =3)+
  geom_smooth(method = "lm")+
  scale_color_manual(values = c("lightcoral","lightslateblue"))+
  labs(x = "Bacterial counts in dry stool (normalized)", y = "Shannon index (normalized)") 
# Fig dim inches: 8.5, 5.1
ggplot(adiv_scaled%>%filter(!is.na(Sex)),aes(x=counts,y=Chao1))+
  geom_point(size =3,aes(color=Sex))+
  geom_smooth(method = "lm")+
  scale_color_manual(values = c("lightcoral","lightslateblue"))+
  labs(x = "Bacterial counts in dry stool (normalized)", y = "Chao1 index (normalized)") 


adiv_scaled2<-apply(otu_table(ps.mi4d.asv2),1,function(x) sum(x>0))%>%data.frame(asvcount = .)%>%merge(adiv_scaled, by = "row.names")
fit_asvcount<-summary(lm(asvcount~counts+Sex+Age+BMI,data = adiv_scaled2)) # p: counts = 8.13e-05, SexM = 0.059,Age = 0.048 
ggplot(adiv_scaled2,aes(x=counts,y=asvcount))+geom_point(aes(color = Sex),size = 3)+
  geom_smooth(method = "lm")+
  scale_color_manual(values = c("lightcoral","lightslateblue"))+
  labs(x = "Bacterial counts in dry stool (normalized)", y = "Number of unique ASVs")

 
# ORDINATION Plot
## only choose taxa >  0.1 prevalance, i.e. at least 5 samples and mean relative abundance > 0.0005 (first quartile of relative abundance)
ps.mi4d2<- filter_taxa(ps.mi4d.asv2, function(x) sum(x > 0) > nsamples(ps.mi4d.asv2)*.1 & mean(x/sample_sums(ps.mi4d.asv2),na.rm = T) > .0005, prune = T)
length(get_taxa_unique(ps.mi4d2, "Genus"))

sample_imputed<-read.csv("samplemat_pbmc_scaled_imputed.csv")%>%
  filter(!ID%in%c("MI4D022", "MI4D028", "MI4D044"))%>%
  mutate(ID2 = if_else(grepl("PRE",ID),paste0("MI4D0",as.numeric(str_extract(ID,"\\d{3}?"))+52),ID),
         female = ifelse(Sex=="F",1,-1),
         biomass= ifelse(counts >=median(counts),"High","Low"))%>%
  column_to_rownames("ID2")

sample_data(ps.mi4d2)<-sample_data(sample_imputed)


bray_dists <- ps.mi4d2 %>%
  tax_fix()%>%
  dist_calc("bray") 
pcoa<-bray_dists%>% 
  ord_calc("PCoA")

bray_perm <- bray_dists %>%
  dist_permanova(
    seed = 42, # set.seed to ensure reproducibility of random process
    n_processes = 1, n_perms = 6006, 
    variables = "counts+Sex+Age+BMI+Calprotectin_fecal" 
  )

perm_get(bray_perm) %>% as.data.frame() # record R2 and P

perm2 <- bray_dists %>%
  dist_permanova(variables = c("counts","female","Age","BMI","Calprotectin_fecal"), seed = 42) 

bray_dists_ord<-perm2 %>% 
  ord_calc(method = "CAP",constraints = c("counts","female","Age","BMI","Calprotectin_fecal"))

# add species score to ordination object
vegan::sppscores(bray_dists_ord@ord)<-ps.mi4d2%>%otu_table()
mi4d.mds.score<-data.frame(bray_dists_ord@ord$CCA$biplot)

rda.axes<-data.frame(bray_dists_ord@ord$CCA$u)%>%
  merge(adiv_scaled,by="row.names")
rda.axes%>%ggplot(aes(x= CAP1, y = CAP2))+geom_point(aes(color = biomass),size = 3.5)+
  stat_ellipse(aes(color = biomass),linetype = 2,type = "t", level = .95)+
  scale_color_manual(values = c("#2E9FDF","#E7B800"),name = "Biomass")+
  labs(x = "dbRDA1(5.7%)",y = "dbRDA2(4.0%)")+scale_x_continuous(limits = c(-.4,.4))
# Fig  6.18*6.18
rda.axes%>%ggplot(aes(x= biomass, y = CAP1,color = biomass))+
  geom_boxplot(outlier.colour = NA)+
  geom_point(position = position_jitter(),size = 3.5, alpha = .75)+
  scale_color_manual(values = c("#2E9FDF","#E7B800"),name = "Biomass")+
  labs(y = "dbRDA1(5.7%)",x = "")

rda.axes%>%ggplot(aes(x= CAP1, y = CAP2))+geom_point(aes(color = Sex),size = 3.5)+
  stat_ellipse(aes(color = Sex),linetype = 2,type = "t")+
  scale_color_manual(values = c("lightcoral","lightslateblue"))+
  labs(x = "dbRDA1(5.7%)",y = "dbRDA2(4.0%)")+scale_x_continuous(limits = c(-.4,.4))
# Fig  6.18*6.18
rda.axes%>%ggplot(aes(x= Sex, y = CAP2,color = Sex))+
  geom_boxplot(outlier.colour = NA)+
  geom_point(position = position_jitter(),size = 3.5, alpha = .75)+
  scale_color_manual(values = c("lightcoral","lightslateblue"),name = "Sex")+
  labs(y = "dbRDA2(4.0%)",x = "")


# top taxa contributing to axis1
mi4d.dbrda.asv<-data.frame(bray_dists_ord@ord$CCA$v)%>%arrange(desc(abs(CAP1)))
mi4d.dbrda.asv2<-data.frame(bray_dists_ord@ord$CCA$v)%>%arrange(desc(abs(CAP2)))

top20.dbrda1<-ps.mi4d2%>%
  prune_taxa(rownames(mi4d.dbrda.asv)[1:20],.)

top20.dbrda1.tax<-data.frame(tax_table(top20.dbrda1))%>%
  mutate(Species = case_when(is.na(Genus) & is.na(Species) ~ paste(Family, "Genus"),
                             !(is.na(Genus)) & is.na(Species) ~ paste(Genus, "Sp."),
                             TRUE ~ paste(Genus, Species)))
species_name<-top20.dbrda1.tax$Species
top20.dbrda1.tax<-top20.dbrda1.tax[match(rownames(mi4d.dbrda.asv)[1:20],rownames(top20.dbrda1.tax)),]
top20.dbrda1.asv<-top20.dbrda1%>%tax_transform("clr")%>%otu_table()
# write.table(data.frame(top20.dbrda1.asv),"20230109_Mb_top20.dbrda1.txt",sep="\t",quote = F)
# write.table(top20.dbrda1.tax,"20230109_Mb_top20.dbrda1_tax.txt",sep="\t",quote = F)
# colnames(top20.dbrda1.asv)<-species_name

top20.dbrda1.tax["ASV_0063",]$Species<-"Eubacterium\n coprostanoligenes Genus"
adiv_scaled2%<>%column_to_rownames("Row.names")%>%rownames(top20.dbrda1.asv)

top20.dbrda1.test<-cbind(top20.dbrda1.asv,adiv_scaled2[,c("biomass","Sex")])
top20.dbrda1.scaled.test_result<- top20.dbrda1.test%>%
  gather(key = ASV, value = abund, -biomass,-Sex) %>%
  split(.$ASV)%>%
  map(~compare.2.vectors(.[.$biomass=="Low",]$abund,.[.$biomass=="High",]$abund,tests = "nonparametric"))%>%
  walk(print)

top20.dbrda1.scaled.test_result_nonpara<-top20.dbrda1.scaled.test_result%>%map_df(~.$nonparametric$p)
top20.dbrda1.scaled.test_result_nonpara<-top20.dbrda1.scaled.test_result_nonpara%>%t()%>%data.frame()%>%
  rename_at(vars(colnames(.)),~c("stats_Wilcoxon","permutation","coin_Wilcoxon","median"))
top20.dbrda1.scaled.test_result_nonpara$ASV <-rownames(top20.dbrda1.scaled.test_result_nonpara)
top20.dbrda1.scaled.test_result_nonpara$p.adj<-p.adjust(top20.dbrda1.scaled.test_result_nonpara$permutation,method = "BH")
top20.dbrda1.scaled.test_result_nonpara%<>%mutate(symbol = case_when(p.adj < .05 ~ "p < 0.05",p.adj < .1 ~ "p < 0.1"))

top20.dbrda1.test.long<-top20.dbrda1.test%>%pivot_longer(cols = !c(biomass,Sex), names_to = "ASV", values_to = "ASV_read")%>%
  left_join(top20.dbrda1.scaled.test_result_nonpara)

# 12 * 12
top20.dbrda1.test.long%>%ggplot(aes(x= ASV_read, y = ASV,color =biomass))+
  geom_boxplot(outlier.shape = NA)+geom_point(position=position_jitterdodge())+
  scale_y_discrete(limits = rownames(top20.dbrda1.tax),labels =  top20.dbrda1.tax$Species)+
  scale_color_manual(name = "Biomass",values = c("#2E9FDF","#E7B800"))+
  geom_text(data = top20.dbrda1.test.long%>%group_by(ASV)%>%slice_max(order_by = ASV_read), size = 5,
            x = 7,color = "black",aes(label = symbol))+labs(y = "",x="ASV reletive abundance\n(center-log transformed)")

top20.dbrda1.test.long.sign<-top20.dbrda1.test.long%>%filter(!is.na(symbol))
top20.dbrda1.tax.sign<-top20.dbrda1.tax[unique(top20.dbrda1.test.long.sign$ASV),]
top20.dbrda1.test.long.sign%>%ggplot(aes(y= ASV_read, x = ASV,color =biomass))+
  geom_boxplot(outlier.shape = NA)+geom_point(position=position_jitterdodge())+
  scale_x_discrete(limits = rownames(top20.dbrda1.tax.sign),labels =  top20.dbrda1.tax.sign$Species)+
  scale_color_manual(name = "Biomass",values = c("#2E9FDF","#E7B800"))+
  labs(x = "",y="ASV reletive abundance\n(center-log transformed)")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# dbrda2 ------------------------------------------------------------------

top20.dbrda2<-ps.mi4d2%>%
  prune_taxa(rownames(mi4d.dbrda.asv2)[1:20],.)

top20.dbrda2.tax<-data.frame(tax_table(top20.dbrda2))%>%
  mutate(Species = case_when(is.na(Genus) & is.na(Species) ~ paste(Family, "Genus"),
                             !(is.na(Genus)) & is.na(Species) ~ paste(Genus, "Sp."),
                             TRUE ~ paste(Genus, Species))
           )
species_name2<-top20.dbrda2.tax$Species
top20.dbrda2.tax<-top20.dbrda2.tax[match(rownames(mi4d.dbrda.asv2)[1:20],rownames(top20.dbrda2.tax)),]
top20.dbrda2.asv<-top20.dbrda2%>%tax_transform("clr")%>%otu_table()
# write.table(data.frame(top20.dbrda2.asv),"20230109_Mb_top20.dbrda2.txt",sep="\t",quote = F)
# write.table(top20.dbrda2.tax,"20230109_Mb_top20.dbrda2_tax.txt",sep="\t",quote = F)
# colnames(top20.dbrda2.asv)<-species_name2

top20.dbrda2.test<-cbind(top20.dbrda2.asv,adiv_scaled2[,c("biomass","Sex")])


top20.dbrda2.scaled.test_result<- top20.dbrda2.test%>%
  gather(key = ASV, value = abund, -biomass,-Sex) %>%
  split(.$ASV)%>%
  map(~compare.2.vectors(.[.$Sex == "F",]$abund,.[.$Sex == "M",]$abund))%>%
  walk(print)

top20.dbrda2.scaled.test_result_nonpara<-top20.dbrda2.scaled.test_result%>%map_df(~.$nonparametric$p)
top20.dbrda2.scaled.test_result_nonpara<-top20.dbrda2.scaled.test_result_nonpara%>%t()%>%data.frame()%>%
  rename_at(vars(colnames(.)),~c("stats_Wilcoxon","permutation","coin_Wilcoxon","median"))
top20.dbrda2.scaled.test_result_nonpara$ASV <-rownames(top20.dbrda2.scaled.test_result_nonpara)
top20.dbrda2.scaled.test_result_nonpara$p.adj<-p.adjust(top20.dbrda2.scaled.test_result_nonpara$permutation,method = "fdr")
top20.dbrda2.scaled.test_result_nonpara%<>%mutate(symbol = case_when(p.adj < .05 ~ "p < 0.05",p.adj < .1 ~ "p < 0.1"))

top20.dbrda2.scaled.test_result_nonpara<-top20.dbrda2.scaled.test_result_nonpara[rownames(top20.dbrda2.tax),]



top20.dbrda2.test.long<-top20.dbrda2.test%>%pivot_longer(cols = !c(biomass,Sex), names_to = "ASV", values_to = "ASV_read")%>%
  left_join(top20.dbrda2.scaled.test_result_nonpara)
# 12*12
top20.dbrda2.test.long%>%ggplot(aes(x= ASV_read, y = ASV,color =Sex))+
  geom_boxplot(outlier.shape = NA)+geom_point(position=position_jitterdodge())+
  scale_y_discrete(limits = rownames(top20.dbrda2.tax),labels =  top20.dbrda2.tax$Species)+
  scale_color_manual(values = c("lightcoral","lightslateblue"))+
  geom_text(data = top20.dbrda2.test.long%>%group_by(ASV)%>%slice_max(order_by = ASV_read), size = 5,
            x = 6,color = "black",aes(label = symbol))+labs(y = "",x="ASV reletive abundance\n(center-log transformed)")




