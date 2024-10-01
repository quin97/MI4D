# Required Packages -------------------------------------------------------
packages <- c("BiocManager", 
              "plyr", 
              "dplyr",
              "data.table",
              "stringr", 
              "reshape2", 
              "ggplot2", 
               "phyloseq",
              "microViz",
              "tidyr",
              "tibble",
              "RColorBrewer",
              "readxl",
              "ghibli",
              "wesanderson")

# Install regular packages if needed
#install.packages(setdiff(packages, rownames(installed.packages())))  

# Load packages -----------------------------------------------------------
lapply(packages, require, character.only = TRUE)
Sys.setenv(PATH = paste(Sys.getenv("PATH"), "./16S_ribosomal_RNA/", sep = .Platform$path.sep))

# Data Analysis -----------------------------------------------------------
# Setup -------------------------------------------------------------------
# Declare file prefixes to attach

errPoolName = "SEQ291_298_2023_11_02"


# Read in tables ----------------------------------------------------------

asv = read.table("MI4D Chemostat pseudo_pooled ASV_Table run291_298_SILVA_taxonomy_2023_11_02.txt", header = TRUE, row.names = NULL, sep = "")

# Renames ASV_Num to ASV_Number (consistency)
asv <- asv %>%
  dplyr::rename(OTU_Number = OTU_Num)

# Creates extra ASV_Number_Original column for ease of tracking
asv <- asv %>%
  dplyr::mutate(OTU_Number_Original = OTU_Number, .after = "OTU_Number") 

# Select desired columns, if needed
asv = asv[,!grepl("(?<=[HF|LF])+(D2.|D3.|D00)|PREMI4D|H2O|ZYMO|Control|Media", names(asv), perl = TRUE)]

# Zero padding sample names
colnames(asv)[grep("MI4D", colnames(asv))] = gsub("MI4D", "MI4D0", colnames(asv)[grep("MI4D", colnames(asv))])

# This just reorganizes s.t. inocula are at the beginning of the table
asv = asv %>% 
  relocate(grep("(MI4D002LFStool)", names(asv), value = T, perl = T), .before = (contains("MI4D002LFD01"))) %>%
  relocate(grep("(MI4D007HFStool)", names(asv), value = T, perl = T), .before = (contains("MI4D007HFD01"))) %>%
  relocate(grep("(MI4D016HFStool)", names(asv), value = T, perl = T), .before = (contains("MI4D016HFD01"))) %>%
  relocate(grep("(MI4D046LFStool)", names(asv), value = T, perl = T), .before = (contains("MI4D046LFD01"))) %>%
  relocate(grep("(MI4D008LFStool)", names(asv), value = T, perl = T), .before = (contains("MI4D008LFD01"))) %>%
  relocate(grep("(MI4D011HFStool)", names(asv), value = T, perl = T), .before = (contains("MI4D011HFD01"))) %>%
  relocate(grep("(MI4D021HFStool)", names(asv), value = T, perl = T), .before = (contains("MI4D021HFD01"))) %>%
  relocate(grep("(MI4D036LFStool)", names(asv), value = T, perl = T), .before = (contains("MI4D036LFD01"))) 
  
# Rename ASV sample stool to D00
names(asv)
names(asv) = gsub("Stool", "D00", names(asv))

# Sort columns by number of reads
asv = asv [order(rowSums(asv[15:ncol(asv)]),decreasing = TRUE),] 
asv = asv [order(asv$dummy_BLASTphylum),] 
asv$OTU_Number = 1:nrow(asv)
asv$OTU_Number = sprintf("%04d", asv$OTU_Number)

# Phyloseq ----------------------------------------------------------------

# Plate 1 and 2 merged meta data
meta.asv <- as.data.table(grep("OTU|BLAST|sequences|per_id|bit|evalue", names(asv), value = T, invert = "TRUE"))
meta.asv = meta.asv %>% 
  dplyr::rename(Sample = V1) 
meta.asv = tidyr::extract (meta.asv, Sample, into = "Day", regex = "(?<=[HL]F|PREMI4D)+(.*)", remove = F)
meta.asv$Day = gsub ("D", "", meta.asv$Day)
meta.asv$Day = gsub ("Stool", "00", meta.asv$Day)
meta.asv = separate (meta.asv, Sample, into = "Patient", sep = 8, remove = F)
meta.asv = extract(meta.asv, Patient, into = "Biomass", regex = "(?<=MI4D...)+(.)", remove = F)
meta.asv$Patient = gsub("H", "", meta.asv$Patient)
meta.asv$Patient = gsub("L", "", meta.asv$Patient)
meta.asv$Biomass = gsub("H", "High", meta.asv$Biomass)
meta.asv$Biomass = gsub("L", "Low", meta.asv$Biomass)
meta.asv$Sample2 = meta.asv$Sample
meta.asv$Sample2 = gsub("MI4D0", "MI4D", meta.asv$Sample2)
meta.asv$Patient = gsub ("MI4D016", "MI4D013", meta.asv$Patient)
meta.asv$Day = gsub ("^0", "", meta.asv$Day)
meta.asv$Day = factor(meta.asv$Day, levels = c("0", "1", "2", "4", "6", "8", "10", "12", "14", "16", "18"))
meta.asv$Patient = factor(meta.asv$Patient, levels = c("MI4D002", "MI4D008", "MI4D036", "MI4D046", "MI4D007", "MI4D011", "MI4D013", "MI4D021"))
meta.asv$Biomass = factor(meta.asv$Biomass, levels = c("High", "Low"))
meta.asv$Biomass2 = meta.asv$Biomass


# Setup OTU tables 
otu = asv[,!grepl("OTU|BLAST|sequences|per_id|bit|evalue", names(asv))]

# Setup taxa ID table
tax = asv[,grepl("BLAST", names(asv))]
names(tax) = gsub("dummy_BLAST", "", names(tax))
names(tax) = str_to_title(names(tax))

# Set up metadata tables 
meta.asv <- meta.asv %>%
  tibble::column_to_rownames("Sample")

# Set up all OTU and taxonomic tables as matrices, metadata table can be left as is
otu.mat <- as.matrix(otu)
tax.mat <- as.matrix(tax)

# Import OTU, taxonomic matrices and metadata table into Phyloseq object file
OTU = otu_table(otu.mat, taxa_are_rows = TRUE)
TAX = tax_table(tax.mat)
tax2<-tax%>%mutate(Phylum = case_when(Phylum == "Firmicutes" ~ "Bacillota",
                                      Phylum == "Proteobacteria" ~ "Pseudomonadota",
                                      Phylum == "Actinobacteriota" ~ "Actinomycetota",
                                      Phylum == "Desulfobacterota" ~ "Pseudomonadota",
                                      Phylum == "Euryarchaeota" ~ "Methanobacteriota",
                                      Phylum == "Cyanobacteria" ~ "Cyanobacteriota",
                                      Phylum == "Campylobacterota" ~ "Pseudomonadota",
                                      Phylum == "Chloroflexi" ~ "Chloroflexota",
                                      Phylum == "Dependentiae" ~ "Babelota",
                                      is.na(Phylum)  ~ "Bacteria Kingdom",
                                      TRUE ~ as.character(Phylum)))
tax2<-tax_table(tax2)
colnames(tax2)<-rank_names(TAX)
taxa_names(tax2)<-taxa_names(TAX)

meta = sample_data(meta.asv)
ps = phyloseq(OTU, tax2, meta)

# Check Phyloseq Objects
ps

# Theme Settings ----------------------------------------------------------
quin.theme = theme_classic(base_size = 25)
quin.theme.text =  theme(strip.background = element_rect(colour="white", fill="white"),
                         axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5))

# Set manual colour scale
palette = colorRampPalette(c(ghibli_palette("LaputaMedium")[7:1],ghibli_palette("KikiMedium")[4:7],wes_palette("Darjeeling2")))(18)

# Reorder phyla to be nicer colours
phyla.list = ps %>% 
  tax_fix() %>% 
  tax_sort(by = sum) 
phyla.list =  unique(tax_table(phyla.list)[,grep("Phylum", colnames(tax_table(phyla.list)))])
names(palette) = phyla.list


# Alpha diversity --------------------------------------------
ps.fix = subset_samples(ps, Day != "0")
ps.fix = ps.fix %>% tax_fix()

# Fig 1F 12.5*6.18
alpha_chao <- plot_richness(ps.fix, x = "Day", color = "Biomass", measures = c("Chao1")) + 
  labs(x = "Timepoint (day of bioreactor culture)", y = "Chao1 index") +  stat_summary(fun = mean, geom = "line", aes(group = Biomass, color = Biomass, linetype = Biomass), linewidth = 1) +
  stat_summary(fun = mean, geom = "point", aes(group = Biomass, color = Biomass), size = 3) + 
 stat_summary(fun.data = mean_se, geom = "errorbar", aes(group = Biomass, linetype = Biomass), alpha = 0.7, position = position_dodge(width = .2)) + 
  scale_color_manual(values =  c("#2E9FDF", "#E7B800"))

alpha_chao$layers <- alpha_chao$layers[c(-1, -2)]

alpha_chao + quin.theme + quin.theme.text

# Shannon
alpha_shannon <- plot_richness(ps.fix, x = "Day", color = "Biomass", measures = c("shannon")) + 
  xlab("Day") +  stat_summary(fun = mean, geom = "line", aes(group = Biomass, color = Biomass, linetype = Biomass), linewidth  = 1) +
  stat_summary(fun = mean, geom = "point", aes(group = Biomass, color = Biomass), size = 3) + 
    stat_summary(fun.data = mean_se, geom = "errorbar", aes(group = Biomass, color = Biomass), alpha = 0.7, position = position_dodge(width = .2)) +  
  scale_color_manual(values =  c("#2E9FDF", "#E7B800"))

alpha_shannon$layers <- alpha_shannon$layers[c(-1)]

alpha_shannon + quin.theme + quin.theme.text

# Export
tiff(eval(paste0("./",errPoolName,"_chao.tiff")), width = 1920, height = 1080)
alpha_chao + quin.theme + quin.theme.text
dev.off()




# NMLE testing on alpha diversity
require(nlme)

df = estimate_richness(ps.fix, measures = "Chao1")
df2 = sample_data(ps.fix)
df3 = merge(df, df2, by = 'row.names')

fit.chao = df3 %>% 
  mutate(Day = factor(Day, levels = rev(unique(Day)))) %>%
    lme(Chao1 ~ Day + Biomass, random = ~1 | Row.names, data = .) 
round(summary(fit.chao)$tTable, 2)


ps%>%
  tax_fix()%>%
  comp_barplot(
    x = "Day",
    tax_level = "Phylum",n_taxa = 14, 
    sample_order = "-1",
    palette = colorRampPalette(c(ghibli_palette("LaputaMedium")[7:3],ghibli_palette("KikiMedium")[4:7],wes_palette("Darjeeling2")))(16),
    other_name = "Other",
    facet_by = "Patient",
    nrow = 2,
    merge_other = F, bar_outline_colour = "black" 
  ) +
  labs(x = NULL, y = "Relative abundance") + 
  theme(# axis.text.y = element_blank(), axis.ticks.y = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


