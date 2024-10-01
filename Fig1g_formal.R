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

