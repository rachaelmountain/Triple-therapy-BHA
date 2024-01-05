# load required packages

#library(tidyverse)
library(dplyr)
library(ggplot2)
library(remotes)
install_github("rachaelmountain/epicR", ref="triple-therapy-BHA")
library(epicR)
library(tidyverse)
library(grid)
library(magrittr)
library(scales)




# common input values
time_horizon <- 20 # change for sensitivity analysis
discount_rate <- 0.015 # change for sensitivity analysis
closed_cohort <- 1
med_adherence <- 0.7 # change for sensitivity analysis
smoking_adherence <- 0.7

GOLD_min <- 1 # change for subgroup analysis
dyspnea <- 0 # change for subgroup analysis (1=symptomatic)



# settings
record_mode <- c(record_mode_none=0,
                 record_mode_agent=1,
                 record_mode_event=2,
                 record_mode_some_event=3)

settings <- get_default_settings()
settings$record_mode <- record_mode["record_mode_none"]
settings$random_number_agent_refill <- 1
settings$runif_buffer_size <- time_horizon*70 
settings$rexp_buffer_size <- time_horizon*140 
settings$rnorm_buffer_size <- time_horizon*8 + 10 
settings$rgamma_buffer_size <- time_horizon + 2 

settings$n_base_agents <- 10000000





# function
BHA_simul <- function(early_triple, GOLD_min=1, dyspnea=0){
  
  early_triple=1; GOLD_min=1; dyspnea=0;
  
  init_session(settings=settings)
  input <- get_input()
  
  input$values$global_parameters$time_horizon <- time_horizon
  input$values$global_parameters$discount_qaly <- discount_rate
  input$values$global_parameters$closed_cohort <- closed_cohort
  
  input$values$medication$medication_adherence <- med_adherence
  input$values$smoking$smoking_cessation_adherence <- smoking_adherence
  
  input$values$medication$triple_therapy_early_initiation <- early_triple
  input$values$medication$triple_therapy_early_initiation_criteria[1:2] <- c(GOLD_min, dyspnea)
  
  #input$values$utility$pneumonia_dutil <- -0.118 # change for subgroup analysis
  
  #set.seed(445)
  run(input=input$values)
  
  output <- Cget_output()
  output_ex <- Cget_output_ex()
  
  terminate_session()
  
  #output$n_agents/output$n_cohort # about 0.4% fit analysis criteria
  
  res <- data.frame(agents=output$n_agents,
                    cohort=output$n_cohort,
                    exacs=sum(output_ex$n_exac_by_ctime_severity),
                    exacs_sev=sum(output_ex$n_exac_by_ctime_severity[,c(3,4)]),
                    death_exacs=sum(output_ex$n_exac_death_by_ctime_severity),
                    death_other=output$n_deaths-sum(output_ex$n_exac_death_by_ctime_severity),
                    pneumonia=sum(output_ex$n_pneumonia_by_ctime),
                    treat_yrs_dual=output_ex$medication_time_by_class[2]-output_ex$medication_time_by_class[4],
                    treat_yrs_triple=output_ex$medication_time_by_class[4],
                    life_yrs=output$cumul_time,
                    qalys=output$total_qaly
  )
  
  exacs <- as.data.frame(output_ex$n_exac_by_gold_severity) %>% 
    rename(GOLD1=V1,GOLD2=V2,GOLD3=V3,GOLD4=V4) %>% 
    mutate(severity=c("mild","mod","sev","vsev"), .before=GOLD1)

  
  return(list(res=res,exacs=exacs))
  
}

base <- BHA_simul(early_triple=0, GOLD_min=2, dyspnea=0)
early <- BHA_simul(early_triple=1, GOLD_min=2, dyspnea=0)



# save results - change file name accordingly
# saveRDS(list(base=results_base,early=results_early), "BHA_results.Rda")






# Results -----------------------------------------------------------------


# results function

results <- function(dat, analysis){
  
  results <- (rbind(dat$base$res/dat$base$res$cohort*100,
                   dat$early$res/dat$early$res$cohort*100,
                   dat$early$res/dat$early$res$cohort*100-dat$base$res/dat$base$res$cohort*100))
  
  results <- results %>% 
    mutate(analysis=paste(analysis), .before=agents) %>% 
    mutate(group=c("base","early","diff"), .before=agents) %>% 
    pivot_longer(cols=agents:qalys, names_to="outcome")
  
  return(results)
  
}



exac_qaly <- matrix(c(0.0225, 0.0155, 0.0488, 0.0488,
                      0.0225, 0.0155, 0.0488, 0.0488,
                      0.0728,	0.0683,	0.0655,	0.0655,
                      0.0728,	0.0683,	0.0655,	0.0655),
                    nrow=4,ncol=4,byrow=T)



# primary analysis results

dat <- readRDS("BHA_results.Rda") 

res_main <- results(dat=dat, analysis="Primary")





# Demographics ------------------------------------------------------------

dem <- readRDS("BHA_cohort_dem.Rda")

prop.table(dem$sex)*100 # 54.7% female

age_tab <- cbind(40:100,dem$age[which(dem$age>0)])
sum(age_tab[,1]*age_tab[,2])/sum(age_tab[,2]) # 74.0 mean age

round(prop.table(dem$gold[-1])*100,1)




# Subgroups plot ----------------------------------------------------------

dat <- readRDS("BHA_results.Rda")
datSG1 <- readRDS("BHA_results_SG_GOLD.Rda")
datSG2 <- readRDS("BHA_results_SG_dysp.Rda")


results <- rbind(dat$early$res/dat$early$res$cohort*100-dat$base$res/dat$base$res$cohort*100,
                 datSG1$early$res/datSG1$early$res$cohort*100-datSG1$base$res/datSG1$base$res$cohort*100,
                 datSG2$early$res/datSG2$early$res$cohort*100-datSG2$base$res/datSG2$base$res$cohort*100)


results <- results %>% 
  mutate(scenario=c("GOLD C/D","GOLD 2","GOLD D"),.before=agents) %>% 
  pivot_longer(cols=agents:qalys,names_to = "outcome", values_to = "value")

plotdat <- results %>% 
  filter(outcome!="agents" & outcome!="cohort" & outcome!="death_other") %>% 
  mutate(scenario=fct_relevel(scenario,c("GOLD C/D","GOLD 2","GOLD D"))) %>% 
  mutate(outcome=fct_recode(outcome,`Total exacerbations`="exacs",
                            `Severe/very severe exacerbations`="exacs_sev",
                            `Pneumonia events`="pneumonia",
                            `Mortality due to exacerbations`="death_exacs",
                            `Dual therapy treatment years`="treat_yrs_dual",
                            `Triple therapy treatment years`="treat_yrs_triple",
                            `Life years`="life_yrs",
                            `Total QALYs`="qalys"
                            )) %>% 
  mutate(outcome=fct_relevel(outcome,c("Total exacerbations","Severe/very severe exacerbations",
                                       "Mortality due to exacerbations","Pneumonia events",                   
                                       "Dual therapy treatment years","Triple therapy treatment years",
                                       "Life years","Total QALYs")))

ggplot(plotdat, aes(y=value, fill=scenario, x=scenario)) + 
  geom_bar(position="dodge", stat="identity", width=0.75) +
  scale_fill_manual("Subgroup", values=c("skyblue","thistle","lightgreen"),
                    labels=c("Primary analysis - GOLD C/D",
                             "Subgroup analysis - GOLD \u2265 2","Subgroup analysis - GOLD D")) +
  xlab("") + ylab("Difference in outcome measures for early initiation relative to standard care") +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  facet_wrap(~outcome, scales="free") +
  theme(legend.position = c(0.83,0.18), legend.title=element_blank(),
        legend.text = element_text(size=11),
        legend.spacing.y = unit(0.5, 'cm')) +
  guides(fill = guide_legend(byrow = TRUE))

ggsave("subgroups.tiff",width=5000,height=3750,units="px",dpi=600)



ggsave("subgroups.pdf",width=5000,height=3750,units="px",dpi=600)






# Sensitivity analysis results --------------------------------------------

SATElow <- readRDS("BHA_results_TElow.Rda")
SATEhigh <- readRDS("BHA_results_TEhigh.Rda")
SARlow <- readRDS("BHA_results_Rlow.Rda")
SARhigh <- readRDS("BHA_results_Rhigh.Rda")
SAUlow <- readRDS("BHA_results_Ulow.Rda")
SAUhigh <- readRDS("BHA_results_Uhigh.Rda")
SAMA30 <- readRDS("BHA_results_MA30.Rda")
SAMA50 <- readRDS("BHA_results_MA50.Rda")
SADR0 <- readRDS("BHA_results_D0.Rda")
SADR3 <- readRDS("BHA_results_D3.Rda")
SATH5 <- readRDS("BHA_results_TH5.Rda")
SATH35 <- readRDS("BHA_results_TH35.Rda")

SA_results <- function(dat){
  
  results <- rbind(dat$base$res$qalys/dat$base$res$cohort*100,
                     dat$early$res$qalys/dat$early$res$cohort*100,
                     dat$early$res$qalys/dat$early$res$cohort*100-dat$base$res$qalys/dat$base$res$cohort*100)
  
  results <- as.data.frame(t(results))
  names(results) <- c("base","early","diff")
  
  return(results)
  
}

SA_restab <- rbind(
  SA_results(dat=SATElow),
  SA_results(dat=SATEhigh),
  SA_results(dat=SARlow),
  SA_results(dat=SARhigh),
  SA_results(dat=SAUlow),
  SA_results(dat=SAUhigh),
  SA_results(dat=SAMA30),
  SA_results(dat=SAMA50),
  SA_results(dat=SADR0),
  SA_results(dat=SADR3),
  SA_results(dat=SATH5),
  SA_results(dat=SATH35)
)


# original value of output
base.value <- 4.751033

# width of columns in plot (value between 0 and 1)
width <- 0.75

order.parameters <- SA_restab %>% 
  mutate(parameter=rep(c("Treatment effectiveness",
                         "Pneumonia risk","Pneumonia utility",
                         "Medication adherence","Discount rate","Time horizon"),each=2)) %>% 
  mutate(type=rep(c("Low","High"),6)) %>% 
  mutate(type=factor(type,levels=c("Low","High"))) %>% 
  filter(parameter!="MA" & parameter!="DR") %>% 
  select(-base) %>% select(-early) %>% 
  pivot_wider(names_from=type,values_from=diff) %>% 
  rowwise() %>% 
  mutate(UL_diff=max(abs(High-Low),abs(Low-base.value),abs(High-base.value))) %>% 
  arrange(UL_diff) %>%
  mutate(parameter=factor(x=parameter, levels=parameter)) %>%
  select(parameter) %>% unlist() %>% levels()

df <- SA_restab %>% 
  mutate(parameter=rep(c("Treatment effectiveness",
                         "Pneumonia risk","Pneumonia utility",
                         "Medication adherence","Discount rate","Time horizon"),each=2)) %>% 
  mutate(type=rep(c("Low","High"),6)) %>% 
  mutate(type=factor(type,levels=c("Low","High"))) %>% 
  filter(parameter!="MA" & parameter!="DR") %>% 
  select(-base) %>% select(-early) %>% 
  mutate(parameter=factor(parameter, levels=order.parameters),
         ymin=pmin(diff, base.value),
         ymax=pmax(diff, base.value),
         xmin=as.numeric(parameter)-width/2,
         xmax=as.numeric(parameter)+width/2)


# create plot
# (use scale_x_continuous to change labels in y axis to name of parameters)
p <- ggplot() + 
  geom_rect(data = df, 
            aes(ymax=ymax, ymin=ymin, xmax=xmax, xmin=xmin, fill=type)) +
  scale_fill_manual(values=c("skyblue","lightgreen")) +
  theme_bw() + 
  ylab("Net QALY") +
  theme(axis.title.y=element_blank(),
        axis.text=element_text(size=11),
        legend.position = 'bottom',
        legend.title = element_blank()) + 
  geom_hline(yintercept = base.value) +
  scale_x_continuous(breaks = c(1:length(order.parameters)),
                     labels = c(label_wrap(20)("Pneumonia risk (1.25 / 1.87)"),
                                label_wrap(25)("Pneumonia disutility\n(0.0081 / 0.12)"),
                                label_wrap(15)("Discount rate \n(0% / 3%)"),
                                label_wrap(25)("Medication adherence\n(30% / 50%)"),
                                label_wrap(20)("Time horizon\n(5-year / 35-year)"),
                                label_wrap(30)("Triple therapy effectiveness\n(0.35 / 0.48)"))) +
  scale_y_continuous(breaks=c(2,3,4,5,6,7)) +
  coord_flip()


p





# define a list (don't use `c(...)` here) of desired y-axis labels, starting with the
# bottom-most label in your plot & work up from there
desired.labels <- list(label_wrap(20)("Pneumonia risk (1.25 / 1.87)"),
                       label_wrap(25)("Pneumonia disutility\n(0.0081 / 0.12)"),
                       label_wrap(15)("Discount rate \n(0% / 3%)"),
                       label_wrap(25)("Medication adherence\n(30% / 50%)"),
                       label_wrap(20)("Time horizon\n(5-year / 35-year)"),
                       label_wrap(30)("Triple therapy effectiveness\n(0.35 / 0.48)"))



# convert to grob object
gp <- ggplotGrob(p)

# locate label grob in the left side y-axis
old.label <- gp$grobs[[grep("axis-l", gp$layout$name)]]$children[["axis"]]$grobs[[1]]$children[[1]]

# define each label as its own text grob, replacing the values with those from
# our list of desired y-axis labels
new.label <- lapply(seq_along(old.label$label),
                    function(i) textGrob(label = desired.labels[[i]],
                                         x = old.label$x[i], y = old.label$y[i],
                                         just = old.label$just, hjust = old.label$hjust,
                                         vjust = old.label$vjust, rot = old.label$rot,
                                         check.overlap = old.label$check.overlap,
                                         gp = old.label$gp))

# remove the old label
gp$grobs[[grep("axis-l", gp$layout$name)]]$children[["axis"]]$grobs[[1]] %<>%
  removeGrob(.$children[[1]]$name)

# add new labels
for(i in seq_along(new.label)) {
  gp$grobs[[grep("axis-l", gp$layout$name)]]$children[["axis"]]$grobs[[1]] %<>%
    addGrob(new.label[[i]])
}

# check result
grid.draw(gp)


ggsave("sensitivity.tiff",width=5000,height=3000,units="px",dpi=600)

ggsave("sensitivity.pdf",width=5000,height=3000,units="px",dpi=600)

