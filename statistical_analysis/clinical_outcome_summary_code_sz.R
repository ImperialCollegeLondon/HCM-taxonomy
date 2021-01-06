#script to run survival analysis - by Sean L Zheng

#load in required packages
library(data.table)
library(tidyverse)
library("survival")
library("survminer")
library("lessR")
library(RColorBrewer)
library(gtsummary)
library(cmprsk)

####Load in survival dataset####
sa <- fread('/Path_to_survival_dataset')

#filter out those that have withdrawn consent
sa <- sa %>% filter(hasWithdrawn==FALSE)

#if running incident cases then include line below to remove individuals who have had an event preceding recruitment
#if not then bypass
sa <- sa %>% filter(is.na(pre_death_mace))

#####fit survival model####
fit <- survfit(Surv(time_death_mace, event_death_mace) ~ geno, data=sa)

#produce CH survival curve
ggsurv<- ggsurvplot(fit, 
                    ylim=c(0,1), 
                    xlab = "Years",
                    ylab="Probability freedom from events", 
                    legend='bottom', 
                    legend.title = 'Genotype', font.legend=c(14),
                    legend.labs = c('SARC-NEG', 'SARC-IND', 'SARC-P/LP'),
                    palette=c("#BDBDBD", "#FFA500",'#A52A2A'),
                    risk.table = TRUE, risk.table.title = NULL, risk.table.height=0.20,
                    linetype='strata',
                    fontsize=4, font.x=14, font.y=14, font.tickslab=12,
                    break.time.by = 10,
                    break.y.by=0.1,
                    fun='cumhaz',
                    title=NULL)
ggsurv$plot <- ggsurv$plot + theme(legend.text=element_text(size=14))
ggsurv

#####Cox PH - compare between pairs of arms ####
res.cox <- coxph(Surv(sa$time_death_mace, event_death_mace) ~ 
                   gene_positive  + gene_vus + sexgenetic, 
                 data=sa)
print(res.cox)
coxph(Surv(sa$time_death_mace, event_death_mace) ~ 
        gene_positive + gene_vus + sexgenetic, 
      data=sa) %>% gtsummary::tbl_regression(exp=TRUE)

#if plotting - then separate into 1 vs 1 comparison and then combine later
sa1 <- sa %>% filter(geno!=3)
geno_df <- with(sa1, data.frame(geno=c(1,2), sexgenetic = c(1,1), agerecruit=rep(mean(sa1$agerecruit, na.rm=TRUE), 2)))
res.cox <- coxph(Surv(sa1$time_death_mace, sa1$event_death_mace) ~ geno + agerecruit + sexgenetic, data=sa1)
fit <- survfit(res.cox, data=sa1, newdata=geno_df)
plots <- list()
plots[[1]] <- ggsurvplot(fit, conf.int = TRUE, 
                         legend.labs=c("Sarc-ve", "Sarc Other"), 
                         ylim=c(0,1), 
                         xlim=c(0,80),
                         xlab = "Years", ylab="Freedom from outcome probability", 
                         legend='bottom', 
                         legend.title = 'Genotype',
                         palette=c("#BDBDBD", "#FFA500"),
                         risk.table = FALSE, tables.height=0.2,
                         linetype = 'strata',
                         ggtheme=theme_bw(),
                         fontsize=3,
                         censor=FALSE,
                         break.time.by = 20,
                         break.y.by=0.2,
                         title=NULL,
                         size=1)

sa2 <- sa %>% filter(geno!=2)
geno_df <- with(sa2, data.frame(geno=c(1,3), sexgenetic = c(1,1), agerecruit=rep(mean(sa2$agerecruit, na.rm=TRUE), 2)))
res.cox <- coxph(Surv(sa2$time_death_mace, sa2$event_death_mace) ~ geno + agerecruit + sexgenetic, data=sa2)
fit1 <- survfit(res.cox, data=sa2, newdata=geno_df)
plots[[2]] <- ggsurvplot(fit1, conf.int = TRUE, 
                         legend.labs=c("Sarc-ve", "Sarc+ve"), 
                         ylim=c(0,1), 
                         xlim=c(0,80),
                         xlab = "Years", ylab="Freedom from outcome probability", 
                         legend='bottom', 
                         legend.title = 'Genotype',
                         palette=c("#BDBDBD", "#A52A2A"),
                         risk.table = FALSE, tables.height=0.2,
                         linetype = 'strata',
                         ggtheme=theme_bw(),
                         fontsize=3,
                         censor=FALSE,
                         break.time.by = 20,
                         break.y.by=0.2,
                         title=NULL,
                         size=1)

plots1 <- arrange_ggsurvplots(plots, print=TRUE, nrow=1, ncol=2)

####Competing risk analysis####
#labels
sa$cr_death_mace <- factor(sa$cr_death_mace, levels = c(0,1,2), labels = c('Censored', 'MACE', 'Death'))

sa1 <- sa %>% filter(sexgenetic==0)

#Event free survival - regression
res <- coxph(formula = Surv(sa1$time_death_mace, cr_death_mace =='MACE') ~ gene_positive + gene_vus, data=sa1, ties = c("efron","breslow","exact")[1])
print(res)
coxph(Surv(sa1$time_death_mace, cr_death_mace =='MACE') ~ gene_positive + gene_vus, data=sa1, ties    = c("efron","breslow","exact")[1]) %>% gtsummary::tbl_regression(exp=TRUE)

#####Test for assumption of proportional hazards#####
#Using the Schoenfeld residuals
fit <- coxph(Surv(sa$time_death_mace, sa$event_death_mace) ~ geno + sexgenetic, data=sa)
test <- cox.zph(fit)
print(test)
plot(test)

#####Sex stratified analysis####
sa$sexgenetic1 <- ifelse(sa$sexgenetic==0, 2, 1)
res.cox <- coxph(Surv(sa$time_death_mace, event_death_mace) ~ 
                   (gene_positive * sexgenetic1) + (gene_vus * sexgenetic1), 
                 data=sa)
print(res.cox)
coxph(Surv(sa$time_death_mace, event_death_mace) ~ 
        (gene_positive * sexgenetic1) + (gene_vus * sexgenetic1), 
      data=sa) %>% gtsummary::tbl_regression(exp=TRUE)

fit2<- survfit(Surv(sa$time_death_mace, event_death_mace) ~ geno + sexgenetic, data=sa)
ggsurv <- ggsurvplot(fit2, 
                     ylim=c(0,0.4), 
                     ggtheme=theme_bw(),
                     xlab = "Years from recruitment",
                     ylab="Probability freedom from events", 
                     fun='cumhaz',
                     legend='bottom', 
                     legend.title = 'Genotype', font.legend=c(14),
                     legend.labs = c('SARC-NEG (F)', 'SARC-NEG (M)','SARC-IND (F)', 'SARC-IND (M)','SARC-P/LP (F)','SARC-P/LP (M)'),
                     palette=c("#BDBDBD", "#BDBDBE","#FFA500","#FFA501",'#A52A2A','#A52A2B'),
                     #palette=c("#BDBDBD", "#FFA500",'#A52A2A'),
                     risk.table = TRUE, risk.table.title = NULL, 
                     risk.table.height = 0.3,
                     linetype=c('solid', 'dotted', 'solid', 'dotted','solid', 'dotted'),
                     censor=FALSE,
                     fontsize=4, font.x=14, font.y=14, font.tickslab=14)
ggsurv$plot <- ggsurv$plot + theme(legend.text=element_text(size=14))
ggsurv

####Wall thickness####
#note if running FD then use same dataset but change covariate to fd

wt <- fread('/Path_to_survival_dataset_in_participants_with_cmr')
wt <- wt %>% filter(hasWithdrawn==FALSE)
fit <- survfit(Surv(time_death_mace, event_death_mace) ~ wt, data=wt)

ggsurv<- ggsurvplot(fit, 
                    ylim=c(0,1), 
                    xlim=c(0,80),
                    xlab = "Years",
                    ylab="Probability freedom from events", 
                    legend='bottom', 
                    legend.title = 'LV max wall thickness (mm)', font.legend=c(14),
                    legend.labs = c('<13', '13 to 15', '>15'),
                    palette='Blues',
                    risk.table = TRUE, risk.table.title = NULL, risk.table.height=0.20,
                    linetype='strata',
                    fontsize=4, font.x=14, font.y=14, font.tickslab=12,
                    break.time.by = 10,
                    break.y.by=0.1,
                    fun='cumhaz',
                    title=NULL)
ggsurv$plot <- ggsurv$plot + theme(legend.text=element_text(size=14))
ggsurv

#Cox PH - change between using categorical and continuous variable for LV wall thickness
res.cox <- coxph(Surv(time_death_mace, event_death_mace) ~ 
                   medium + thick
                 #  WT_Max..mm.
                 + autoSysBP + sexgenetic, data=wt)
print(res.cox)
coxph(Surv(time_death_mace, event_death_mace) ~ 
        medium + thick
      #    WT_Max..mm.
      + autoSysBP +  sexgenetic, data=wt) %>% gtsummary::tbl_regression(exp=TRUE)

