library(tidyverse)
library(magrittr)
library(readxl)
library(ghibli)

theme_set(
  theme_classic(base_size = 25)
)

chemo.counts<-read_xlsx("MI4D_chemostat_bacterial_counts_AG7_241_asof20220411.xlsx")
chemo.counts%<>%rename(ID =`Participant ID`, Time = `Timepoint (day of chemostat culture)`,Counts = `[total cells]/mL (media subtracted)`)%>%
  mutate(Biomass = ifelse(Group == "low","Low","High"))
chemo.counts$ID<-gsub("PREMI4D016","MI4D013",chemo.counts$ID)
days<-unique(chemo.counts$Time)
chemo.counts$`Date of chemostat run`

# chemo.counts%>%filter(Time > 0 &Time <=18 & `Date of chemostat run`!="03/2021")%>%
#   ggplot(.,aes(x = Time, y = Counts,color = ID,shape = Biomass, group =ID))+
#   geom_point(size =3)+
#   scale_color_manual(values =c(ghibli_palette("LaputaMedium")[7:4],ghibli_palette("KikiMedium")[4:7]))+
#   geom_line(size = 1,aes(linetype = Biomass))+labs(x = "Timepoint (day of chemostat culture)",y = "Bacteria counts per mL")+
#   scale_x_continuous(breaks = days[days<=14])


counts.toplot<-chemo.counts%>%filter(Time > 0 &Time <=18 & `Date of chemostat run`!="03/2021")%>%
  group_by(Biomass,Time)%>%
  summarise(counts.mean=mean(Counts,na.rm = T),
            counts.sd=sd(Counts,na.rm = T))%>%
  mutate(upper = counts.mean + counts.sd,
         lower = counts.mean - counts.sd)
# Fig 1G 12*6.18
ggplot(counts.toplot, aes(x =as.factor(Time),y =counts.mean,linetype = Biomass,color = Biomass,group = Biomass))+
  geom_point(size = 3)+
  geom_errorbar(aes(ymin = lower , ymax = upper),alpha = .7,position = position_dodge(width = .2))+
  geom_line(size =1)+
  labs(y="Bacterial counts/mL",x="Timepoint (day of bioreactor culture)",title = "")+
  scale_color_manual(values =c("#2E9FDF","#E7B800"))


require(nlme)
fit<-chemo.counts%>%
  filter(Time > 0 &Time <=18)%>%
  mutate(Time = factor(Time, levels = rev(unique(Time))))%>%
  lme(Counts~ Time+Biomass,random = ~ 1|ID,data=.)
round(summary(fit)$tTable,2)
write.csv(round(summary(fit)$tTable,2), "Fig1G_stats.csv",quote = F)




# compare fecal count to https://doi.org/10.1038/s41467-021-27098-7
qmp<-read_xlsx("manuscript2_revision/41467_2021_27098_MOESM3_ESM.xlsx",sheet = "S1-3",
               col_types = c("numeric","numeric","numeric","text","text",rep("numeric",31)))
unique(qmp$Day_Number)
unique(qmp$ID_Number)
fit2<-qmp%>%filter(!is.na(Cell_count_per_gram)&!is.na(Day_Number))%>%
  mutate
  # mutate(Day_Number = factor(Day_Number))%>%
  lme(Cell_count_per_gram~ Day_Number+Enterotype_nr,random = ~ Day_Number|ID_Number,data=.)
round(summary(fit2)$tTable, 3)

ggplot(qmp%>%filter(!is.na(Cell_count_per_gram)&!is.na(Day_Number)), 
       aes(x = Day_Number,y =Cell_count_per_gram,color = ID_Number,group = ID_Number))+
  geom_point(size = 3)+
  # geom_errorbar(aes(ymin = lower , ymax = upper),alpha = .7,position = position_dodge(width = .2))+
  geom_line(size =1)+
  labs(y="Fecal bacterial biomass",x="Day",title = "")

  
qmp%>%filter(!is.na(Cell_count_per_gram)&!is.na(Day_Number)&!is.na(Enterotype_nr))%>%
  group_by(Enterotype_nr,Day_Number)%>%
  summarise(counts.mean=mean(Cell_count_per_gram,na.rm = T),
            counts.sd=sd(Cell_count_per_gram,na.rm = T))%>%
  mutate(upper = counts.mean + counts.sd,
         lower = counts.mean - counts.sd)%>%
  ggplot(aes(x = Day_Number,y =counts.mean,
             color = Enterotype_nr,linetype = Enterotype_nr,
           group = Enterotype_nr))+
  geom_point(size = 3)+
  geom_errorbar(aes(ymin = lower , ymax = upper),alpha = .7,position = position_dodge(width = .2))+
  geom_line(size =1)+
  labs(y="Fecal bacterial biomass",x="Day",title = "")


