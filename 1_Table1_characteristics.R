library(tidyverse)
library(readxl)
library(caret)
library(ghibli)
library(magrittr)

theme_set(
  theme_classic(base_size = 25)
)

mi4d_met<-read.table("samplemat_pbmc.txt",header = T,sep = "\t",as.is = T)%>%select(-c(calories_reported,starts_with("OGTT"),notes))%>%filter(!ID%in%c("PREMI4D002","PREMI4D012","PREMI4D013","PREMI4D015","PREMI4D016"))
mi4d_sex<-read_xlsx("MI4D_PREMI4D_sex.xlsx",col_types = c("text","logical","logical","date","numeric","text","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","text"))
mi4d_met%<>%left_join(mi4d_sex)%>%select(-`Date of Study Visit`,-contains("YN"))%>%mutate(biomass = if_else(counts > 5.9e10, true = "High", false = "Low"))

meta.selected<-mi4d_met%>%
  select(ID, Sex,Age,Tanner=Self_Assessed_Tanner_Stage,Weight_kg, BMI, `Waist Circumference_cm`,Fasting.glucose,Fasting.insulin, IGI,HOMA.IR,WBISI,HbA1c, 
         CRPhs,Cholesterol, Triglycerides,LDL,HDL)%>%
  group_by(Sex)%>%
  summarise(n.sex = sum(!is.na(Sex)),
            n.age = sum(!is.na(Age)),
            med.age = median(Age,na.rm = T),
            min.age = min(Age),
            max.age = max(Age),
            n.tanner = sum(!is.na(Tanner)),
            p.tanner3 = sum(Tanner <=3,na.rm = T),
            p.tanner5 = sum(Tanner >3,na.rm = T),
            p.tanner.na = sum(is.na(Tanner))/n(),
            n.weight = sum(!is.na(Weight_kg)),
            med.weight = median(Weight_kg,na.rm = T),
            min.weight = min(Weight_kg),
            max.weight = max(Weight_kg),
            n.bmi = sum(!is.na(BMI)),
            med.bmi = median(BMI,na.rm = T),
            min.bmi = min(BMI),
            max.bmi = max(BMI),
            n.waist= sum(!is.na(`Waist Circumference_cm`)),
            med.waist = median(`Waist Circumference_cm`,na.rm = T),
            min.waist = min(`Waist Circumference_cm`,na.rm = T),
            max.waist = max(`Waist Circumference_cm`,na.rm = T),
            n.gluc = sum(!is.na(Fasting.glucose)),
            med.gluc = median(Fasting.glucose,na.rm = T),
            min.gluc = min(Fasting.glucose,na.rm = T),
            max.gluc = max(Fasting.glucose,na.rm = T),
            n.ins = sum(!is.na(Fasting.insulin)),
            med.ins= median(Fasting.insulin,na.rm = T),
            min.ins = min(Fasting.insulin,na.rm = T),
            max.ins = max(Fasting.insulin,na.rm = T),
            n.igi = sum(!is.na(IGI)),
            med.igi= median(IGI,na.rm = T),
            min.igi = min(IGI,na.rm = T),
            max.igi = max(IGI,na.rm = T),
            n.homa = sum(!is.na(HOMA.IR)),
            med.homa = median(HOMA.IR,na.rm = T),
            min.homa = min(HOMA.IR,na.rm = T),
            max.homa = max(HOMA.IR,na.rm = T),
            n.wb = sum(!is.na(WBISI)),
            med.wb = median(WBISI,na.rm = T),
            min.wb = min(WBISI,na.rm = T),
            max.wb = max(WBISI,na.rm = T),
            n.hb= sum(!is.na(HbA1c)),
            med.hb = median(HbA1c,na.rm = T),
            min.hb = min(HbA1c,na.rm = T),
            max.hb= max(HbA1c,na.rm = T),
            n.crp = sum(!is.na(CRPhs)),
            med.crp = median(CRPhs,na.rm = T),
            min.crp = min(CRPhs,na.rm = T),
            max.crp= max(CRPhs,na.rm = T),
            n.chol = sum(!is.na(Cholesterol)),
            med.chol = median(Cholesterol,na.rm = T),
            min.chol = min(Cholesterol,na.rm = T),
            max.chol= max(Cholesterol,na.rm = T),
            n.trig = sum(!is.na(Triglycerides)),
            med.trig = median(Triglycerides,na.rm = T),
            min.trig = min(Triglycerides,na.rm = T),
            max.trig= max(Triglycerides,na.rm = T),
            n.ldl = sum(!is.na(LDL)),
            med.ldl = median(LDL,na.rm = T),
            min.ldl = min(LDL,na.rm = T),
            max.ldl= max(LDL,na.rm = T),
            n.hdl = sum(!is.na(HDL)),
            med.hdl = median(HDL,na.rm = T),
            min.hdl = min(HDL,na.rm = T),
            max.hdl= max(HDL,na.rm = T))



wc.test<-mi4d_met%>%
  select(ID, Sex,Age,Weight_kg, BMI, `Waist Circumference_cm`,Fasting.glucose, Fasting.insulin,IGI,HOMA.IR,WBISI,HbA1c,
         CRPhs,Cholesterol, Triglycerides,LDL,HDL)%>%
  pivot_longer(cols = !c(ID,Sex))%>%
  split(.$name)%>%
  map(~wilcox.test(value ~ Sex, data = .))%>%
  map_dbl(~.$p.value)

mi4d_met%>%select(ID, Sex,Tanner=Self_Assessed_Tanner_Stage)%>%column_to_rownames("ID")%>%transmute(Tanner.cat = ifelse(Tanner > 3, "Mature","Immature"),Sex = Sex)%>%table()%>%fisher.test()

