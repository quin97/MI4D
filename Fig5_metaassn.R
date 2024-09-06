library(tidyverse)
library(readxl)
library(caret)
library(wesanderson)
library(magrittr)

theme_set(
  theme_classic(base_size = 25)
)

mi4d_met<-read.table("samplemat_pbmc.txt",header = T,sep = "\t",as.is = T)%>%select(-c(calories_reported,starts_with("OGTT"),notes))
mi4d_sex<-read_xlsx("MI4D_PREMI4D_sex.xlsx",col_types = c("text","logical","logical","date","numeric","text","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","text"))



mi4d_met_full<-mi4d_met%>%left_join(mi4d_sex)%>%select(-`Date of Study Visit`,-contains("YN"))%>%filter(!ID%in%c("PREMI4D002","PREMI4D012","PREMI4D013","PREMI4D015","PREMI4D016"))
mi4d_met<-mi4d_met_full%>%filter(!is.na(Sex))%>%filter(!grepl("PRE",ID))


mi4d_ppc<-preProcess(mi4d_met,
                     method = c("center", "scale", "YeoJohnson"))
mi4d_ppc$yj

mi4d_scaled<-predict(mi4d_ppc,mi4d_met)


# correlation plot --------------------------------------------------------------

analytes<-mi4d_scaled%>%dplyr::select(Elastase=ELA2, Calprotectin_plasma, Myeloperoxidase=MPO, LPS_binding_Protein = LBP, sCD14=CD14)
res<-cor(analytes,use = "na.or.complete")# "pairwise.complete.obs"
require(Hmisc)
res2 <- rcorr(as.matrix(analytes))
require(corrplot)
corrplot(res, method="color", col= colorRampPalette(hcl.colors(100,"RdYlBu"))(10),  
         type="upper", # order="hclust", 
         addCoef.col = NULL, # Add coefficient of correlation
         tl.col="black", tl.srt=45, #Text label color and rotation
         # Combine with significance
         p.mat = res2$P, sig.level = 0.05, insig = "label_sig",  pch = "*",pch.col = "black",pch.cex = 5,
         # hide correlation coefficient on the principal diagonal
         diag=T
)

lipids<-mi4d_scaled%>%dplyr::select(Cholesterol, HDL, LDL, Triglycerides)
res<-cor(lipids,use = "na.or.complete")# "pairwise.complete.obs"
require(Hmisc)
res2 <- rcorr(as.matrix(lipids))
require(corrplot)
corrplot(res, method="color", col= colorRampPalette(hcl.colors(100,"RdYlBu"))(10),  
         type="upper", # order="hclust", 
         addCoef.col = NULL, # Add coefficient of correlation
         tl.col="black", tl.srt=45, #Text label color and rotation
         # Combine with significance
         p.mat = res2$P, sig.level = 0.05, insig = "label_sig",  pch = "*",pch.col = "black",pch.cex = 5,
         # hide correlation coefficient on the principal diagonal
         diag=T
)

glycemic<-mi4d_scaled%>%dplyr::select(Fasting.glucose,Fasting.insulin,HOMA.IR,IGI,WBISI)
res<-cor(glycemic,use = "na.or.complete")# "pairwise.complete.obs"
require(Hmisc)
res2 <- rcorr(as.matrix(glycemic))
require(corrplot)
corrplot(res, method="color", col= colorRampPalette(hcl.colors(100,"RdYlBu"))(10),  
         type="upper", # order="hclust", 
         addCoef.col = NULL, # Add coefficient of correlation
         tl.col="black", tl.srt=45, #Text label color and rotation
         # Combine with significance
         p.mat = res2$P, sig.level = 0.05, insig = "label_sig",  pch = "*",pch.col = "black",pch.cex = 5,
         # hide correlation coefficient on the principal diagonal
         diag=T
)

# Elastase ----------------------------------------------------------------
# baseline is F for Sex
summary(lm(Triglycerides~Age+BMI+Sex+ELA2,mi4d_scaled))
summary(count.ela2<-lm(counts~Age+BMI+Sex*ELA2,mi4d_scaled))

mi4d_scaled%>%
  ggplot(aes(x = ELA2, y=counts,color = Sex,fill = Sex))+# geom_histogram(aes(x = ELA2),bins = 10)
  geom_point(size = 3)+geom_smooth(method = "lm")+
  scale_color_manual(values = wes_palette("GrandBudapest2")[3:4])+
  scale_fill_manual(values = wes_palette("GrandBudapest2")[3:4])+
  labs(x = "Plasma neutrophil elastase (normalized)", y = "Bacterial counts in dry stool \n(normalized)")

mi4d_scaled%>%
  ggplot(aes(x = ELA2, y=Triglycerides))+
  geom_point(size = 3,aes(color=Sex))+geom_smooth(method="lm",se = T)+
  scale_color_manual(values = c("lightcoral","lightslateblue"))+
  labs(x = "Plasma neutrophil elastase (normalized)", y = "Triglycerides (normalized)")


# Calprotectin (Plasma) ------------------------------------------------------------

summary(lm(Triglycerides~Age+BMI+Sex+Calprotectin_plasma,mi4d_scaled))
summary(lm(counts~Age+BMI+Sex*Calprotectin_plasma,mi4d_scaled))


mi4d_scaled%>%
  ggplot(aes(x = Calprotectin_plasma, y=counts,color = Sex,fill = Sex))+
  geom_point(size = 3)+geom_smooth(method = "lm")+
  scale_color_manual(values = wes_palette("GrandBudapest2")[3:4])+
  scale_fill_manual(values = wes_palette("GrandBudapest2")[3:4])+
  labs(x = "Plasma calprotectin (normalized)", y = "Bacterial counts in dry stool \n(normalized)")

mi4d_scaled%>%
  ggplot(aes(x = Calprotectin_plasma, y=Triglycerides))+
  geom_point(size = 3,aes(color = Sex))+geom_smooth(method="lm",se = T)+
  scale_color_manual(values = c("lightcoral","lightslateblue"))+
  labs(x = "Plasma calprotectin (normalized)", y = "Triglycerides (normalized)")


# MPO ------------------------------------------------------

summary(lm(Triglycerides~Age+BMI+Sex+MPO,mi4d_scaled))
summary(lm(counts~Age+BMI+Sex*MPO,mi4d_scaled))

summary(mpo.ela<-lm(MPO~ELA2+Age+Sex+BMI,mi4d_scaled))
summary(lm(MPO~Calprotectin_plasma+Age+Sex+BMI,mi4d_scaled))
summary(lm(ELA2~Calprotectin_plasma+Age+Sex+BMI,mi4d_scaled))

mi4d_scaled%>%
  ggplot(aes(x = MPO, y=counts,color = Sex,fill = Sex))+
  geom_point(size = 3)+geom_smooth(method = "lm")+
  scale_color_manual(values = wes_palette("GrandBudapest2")[3:4])+
  scale_fill_manual(values = wes_palette("GrandBudapest2")[3:4])+
  labs(x = "Plasma myeloperoxidase (normalized)", y = "Bacterial counts in dry stool\n (normalized)")

mi4d_scaled%>%
  ggplot(aes(x = MPO, y=Triglycerides))+
  geom_point(size = 3,aes(color = Sex))+geom_smooth(method="lm",se = T)+
  scale_color_manual(values = c("lightcoral","lightslateblue"))+
  labs(x = "Plasma myeloperoxidase (normalized)", y = "Triglycerides (normalized)")


# Calprotectin fecal ------------------------------------------------------
summary(cal<-lm(Calprotectin_plasma~Calprotectin_fecal+Sex+Age+BMI,mi4d_scaled))
toprint<-round(summary(cal)$coefficients,3)%>%matrix(ncol = 4)

mi4d_scaled%>%
  ggplot(aes(x = Calprotectin_fecal, y=Calprotectin_plasma,color = Sex,fill = Sex))+
  geom_point(size = 3)+
  scale_color_manual(values = c("lightcoral","lightslateblue"))+
  scale_fill_manual(values = c("lightcoral","lightslateblue"))+
  labs(y = "Plasma calprotectin (normalized)", x = "Fecal calprotectin (normalized)")

# CD14 ------------------------------------------------------

summary(cd14.lbp<-lm(CD14~LBP+Age+Sex+BMI,mi4d_scaled)) 

summary(lm(counts~Age+BMI+Sex*CD14,mi4d_scaled))


mi4d_scaled%>%
  ggplot(aes(x = CD14, y=LBP))+
  geom_point(size = 3,aes(color = Sex))+geom_smooth(method="lm",se = T)+
  scale_color_manual(values = c("lightcoral","lightslateblue"))+
  labs(x = "Plasma soluble CD14 (normalized)", y = "Plasma LPS-binding protein \n(normalized)")

# LBP ------------------------------------------------------
summary(lbp.crp<-lm(CRPhs~Age+BMI+Sex+LBP,mi4d_scaled))
toprint<-round(summary(lbp.crp)$coefficients,3)%>%matrix(ncol = 4)


mi4d_scaled%>%
  ggplot(aes(x = LBP, y=CRPhs,color = Sex,fill = Sex))+
  geom_point(size = 3)+# geom_smooth(method = "lm")+
  geom_abline(size = 1.5,slope = coef(lbp.crp)[5], intercept =coef(lbp.crp)[1]+coef(lbp.crp)[4],color = "lightslateblue")+
  geom_abline(size = 1.5,slope = coef(lbp.crp)[5], intercept =coef(lbp.crp)[1],color = "lightcoral")+
  scale_color_manual(values =c("lightcoral","lightslateblue"))+
  labs(x = "Plasma LPS-binding protein (normalized)", y = "Plasma C-reactive protein\n (normalized)")


