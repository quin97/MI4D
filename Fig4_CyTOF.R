library(readxl)
library(tidyverse)
library(afex)

theme_set(
  theme_classic(base_size = 25)
)
mi4d_met<-read.table("samplemat_pbmc.txt",header = T,sep = "\t",as.is = T)%>%select(-c(calories_reported,starts_with("OGTT"),notes))
mi4d_sex<-read_xlsx("MI4D_PREMI4D_sex.xlsx",col_types = c("text","logical","logical","date","numeric","text","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","text"))

#diffcyt masks dplyr
mi4d_met_full<-mi4d_met%>%left_join(mi4d_sex)%>%select(-`Date of Study Visit`,-contains("YN"))%>%filter(!ID%in%c("PREMI4D002","PREMI4D012","PREMI4D013","PREMI4D015","PREMI4D016"))




mi4d.md<-mi4d_met_full%>%dplyr::filter(ID%in%c("MI4D042", "MI4D002", "MI4D051", "MI4D025", "MI4D008", "MI4D024", "MI4D021", "MI4D011", "MI4D007", "MI4D019"))%>%
  transmute(ID = ID,
            Age = Age,
            Sex = Sex,
            BMI = BMI,
            Biomass = factor(ifelse(counts > 5.9E10,"High","Low"),levels = c("High","Low")))


mi4d.md$Sex <- factor(mi4d.md$Sex, levels = c("F", "M"))
mi4d.md$ID <- factor(mi4d.md$ID, levels = c("MI4D042", "MI4D002", "MI4D051", "MI4D025", "MI4D008", "MI4D024", "MI4D021", "MI4D011", "MI4D007", "MI4D019"))
mi4d.md$cluster<-factor(paste0("cluster",c(1,3,3,2,2,2,3,3,1,3)))

mi4d.cells<-mi4d.md%>%# filter(!ID =="MI4D042")%>%
  mutate(neutrophil = c(5.52, 3.1, 27.79, 4.32, 5.74, 2.2, 28.77, 49.5, 62.69, 29.33))

compare.2.vectors(mi4d.cells[mi4d.cells$Biomass=="Low",]$neutrophil,mi4d.cells[mi4d.cells$Biomass=="High",]$neutrophil,tests = "nonparametric")
ggplot(mi4d.cells,aes(x = Biomass,y =neutrophil,color = Biomass))+
  stat_summary(fun.data = mean_se, geom = "pointrange", fun.args = list(mult = 2),position = position_dodge(.7),size = .05)+
  stat_summary(fun.y = mean, geom = "crossbar",position = position_dodge(.7),linewidth = 1) + 
  geom_point(size = 3,position = position_jitterdodge(.5))+
  scale_color_manual(name = "Biomass\n group",values = c("#2E9FDF","#E7B800"))+labs(x = "",y = expression(atop(paste("% CD45"^"lo", "CD66b"^"+", "in"),"Live single cells")))
# "% CD45lo CD66b+ in \nLive single cells"



# ######## permutation test for var
# group_var <- mi4d.cells %>% 
#   group_by(Biomass) %>% 
#   summarise(vars = var(neutrophil))
# 
# diff_var <- group_var %>%
#   summarise(test_stat = diff(vars))
# 
# # Simulation st uo\p
# set.seed(42)
# repetitions <- 1000
# simulated_values <- rep(NA, repetitions)
# 
# # Run our simulations
# for(i in 1:repetitions){
#   simdata <-  mi4d.cells %>% 
#     mutate(Biomass = sample(Biomass))
# 
# 
#   sim_value <- simdata %>% 
#     group_by(Biomass) %>%
#     summarise(vars= var(neutrophil)) %>%
#     summarise(value = diff(vars))
#   
#   # Store your simulated statistics
#   simulated_values[i] <- as.numeric(sim_value)
# }
# 
# # Make into a tibble so it is easy to work with
# sim <- tibble(var_diff = simulated_values)
# 
# # Plot our estimated sampling distribution
# sim %>% ggplot(aes(x=var_diff)) + 
#   geom_histogram(binwidth=1, color="black", fill="gray") +
#   geom_vline(xintercept = abs(diff_var$test_stat), colour = "red") +
#   geom_vline(xintercept = -abs(diff_var$test_stat), colour = "red")
# 
# 
# # Calculate p-value
# num_more_extreme <- sim %>% 
#   filter(abs(var_diff) >= abs(diff_var$test_stat)) %>% 
#   summarise(n())
# 
# p_value <- as.numeric(num_more_extreme / repetitions)



# 12 * 6.18
cells2<-read_csv("./20221122_CyTOF_Experiement_prep/% major Lineages in CD45hi CD66b- Live Singlets % of CD45hi CD66b-s.csv")
ggplot(cells2,aes(x = Population,y =Frequency,color = Biomass))+
  stat_summary(fun.data = mean_se, geom = "pointrange", fun.args = list(mult = 2),position = position_dodge(.7),size = .05)+
  stat_summary(fun.y = mean, geom = "crossbar",position = position_dodge(.7),linewidth = 1) + 
  geom_point(size = 3,position = position_jitterdodge(.2))+
  scale_color_manual(values = c("#2E9FDF","#E7B800"))+
  scale_x_discrete(limits=c("TCRab cells","TCRgd cells","B cells","NK/ILC", "myeloid cells"), labels = c("TCRαβ cells","TCRγδ cells","B cells","NK+ILC", "Myeloid cells"))+
  labs(x = "",y = expression(paste("% CD45"^"hi", "CD66b"^"-")))+
  theme(axis.text.x = element_text(angle = -30, vjust = .1, hjust=0.3))

cells2stats<-cells2%>%split(.$Population)%>%map(~compare.2.vectors(.[.$Biomass=="Low",]$Frequency,.[.$Biomass=="High",]$Frequency,tests = "nonparametric"))%>%
  map_df(~.$nonparametric$p)


cells3<-read_csv("./20221122_CyTOF_Experiement_prep/% T cell subsets in Total T cells % of T cells.csv")

ggplot(cells3,aes(x = Population,y =Frequency,color = Biomass))+
  stat_summary(fun.data = mean_se, geom = "pointrange", position = position_dodge(1),fun.args = list(mult = 2),size = .05)+
  stat_summary(fun.y = mean, geom = "crossbar",position = position_dodge(1),linewidth = 1) + 
  geom_point(size = 3,position = position_jitterdodge(.2))+
  scale_color_manual(values = c("#2E9FDF","#E7B800"))+
  scale_x_discrete(limits=c("Naive CD4","Memory CD4","Naive CD8","Memory CD8", "DN T cells","CD8a- MAIT","CD8a+ MAIT"),labels = c("Naive CD4","Memory CD4","Naive CD8","Memory CD8", "DN T cells",expression(paste("CD8a"^"-"," MAIT")),expression(paste("CD8a"^"+"," MAIT"))))+ 
  labs(x = "",y = "% of T cells")+
  theme( axis.text.x = element_text(angle = -30, vjust = .1, hjust=0.3))
cells3stats<-cells3%>%split(.$Population)%>%map(~compare.2.vectors(.[.$Biomass=="Low",]$Frequency,.[.$Biomass=="High",]$Frequency,tests = "nonparametric"))%>%
  map_df(~.$nonparametric$p)

cells4<-read_csv("./20221122_CyTOF_Experiement_prep/% IgDIgM subsets in B cells % of B cells.csv")

ggplot(cells4,aes(x = Population,y =Frequency,color = Biomass))+
  stat_summary(fun.data = mean_se, geom = "pointrange", position = position_dodge(1),fun.args = list(mult = 2),size = .05)+
  stat_summary(fun.y = mean, geom = "crossbar",position = position_dodge(1),linewidth = 1) + 
  geom_point(size = 3,position = position_jitterdodge(.2))+
  scale_color_manual(values = c("#2E9FDF","#E7B800"),guide = "none")+
  labs(x = "",y = "% of B cells")+
  scale_x_discrete(limits=c("IgDhi B cells","IgMhi B cellls","IgD- IgM- B cells"),labels=c("IgDhi B cells","IgMhi B cellls",expression(paste("IgD"^"-", "IgM"^"-", "B cells"))))+
  theme( axis.text.x = element_text(angle = -30, vjust = .1, hjust=0.3),plot.margin=unit(c(1,5,1,1), 'cm'))
cells4stats<-cells4%>%split(.$Population)%>%map(~compare.2.vectors(.[.$Biomass=="Low",]$Frequency,.[.$Biomass=="High",]$Frequency,tests = "nonparametric"))%>%
  map_df(~.$nonparametric$p)

cells5<-read_csv("./20221122_CyTOF_Experiement_prep/% INKILC Mono DC in non-TB cells % of non-TB cells.csv")

ggplot(cells5,aes(x = Population,y =Frequency,color = Biomass))+
  stat_summary(fun.data = mean_se, geom = "pointrange", position = position_dodge(1),fun.args = list(mult = 2),size = .05)+
  stat_summary(fun.y = mean, geom = "crossbar",position = position_dodge(1),linewidth = 1) + 
  geom_point(size = 3,position = position_jitterdodge(.2))+
  scale_color_manual(values = c("#2E9FDF","#E7B800"))+
  labs(x = "",y = "% of non-T/B cells")+
  scale_x_discrete(limits=c("CD16+ NK","CD16lo CD56hi NK","CD14+ monocytes","CD16+ monocytes","cDC1","cDC2","pDC"),
                   labels = c(expression(paste("CD16"^"+", "NK")),expression(paste("CD16"^"lo", "CD56"^"hi", "NK")),
                              expression(paste("CD14"^"+", "monocytes")),expression(paste("CD16"^"+", "monocytes")),"cDC1","cDC2","pDC"))+
  theme( axis.text.x = element_text(angle = -30, vjust = .1, hjust=0.3))
cells5stats<-cells5%>%split(.$Population)%>%map(~compare.2.vectors(.[.$Biomass=="Low",]$Frequency,.[.$Biomass=="High",]$Frequency,tests = "nonparametric"))%>%
  map_df(~.$nonparametric$p)

cells6<-read_csv("./20221122_CyTOF_Experiement_prep/Paper Figure Markers on CD66b+ CD45lo  cells vs CD14+ monocytes Medians.csv")

cells6$Marker<-unlist(lapply(str_split(cells6$Marker,"-"),"[",1))

cells6%>%
  ggplot(.,aes(x = Population,y =Median,color = Biomass))+
  stat_summary(fun.data = mean_se, geom = "pointrange", fun.args = list(mult = 2),position = position_dodge(.7),size = .05)+
  stat_summary(fun.y = mean, geom = "crossbar",position = position_dodge(.7),linewidth = 1) + 
  geom_point(size = 3,position = position_jitterdodge(.2))+
  scale_color_manual(values = c("#2E9FDF","#E7B800"))+
  scale_x_discrete(limits= c("CD14+ monocytes","CD16+ monocytes","CD45lo CD66b+"),
                   labels = c(expression(paste("CD14"^"+","monocytes")),expression(paste("CD16"^"+","monocytes")),expression(paste("CD45"^"lo", "CD66b"^"+")))
                   # labels=c(expression(atop("CD14"^"+", "monocytes")), expression(atop("CD16"^"+", "monocytes")),
                           #  expression(paste("CD45"^"lo", "CD66b"^"+")))
                   )+
  labs(x = "",y = paste("MMI"))+facet_wrap(~Marker, ncol = 3,dir = "v",scales = "free_y")+
  theme(axis.text.x = element_text(angle = -30,hjust = -.1, vjust = 0.5),legend.justification = c(-0.1, -1))
  



