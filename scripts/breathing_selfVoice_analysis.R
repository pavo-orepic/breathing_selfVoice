rm(list=ls(all=TRUE)); graphics.off() # clear out workspace and close opened figures

# Load libraries
library(ggplot2); 
library(afex); afex_options(es_aov='pes');afex_options(correction_aov='GG') # anova
library(readr)
library(lme4)
library(lmerTest)
library(gridExtra)
library(dplyr)
library(nlme)
library(magrittr)
library(modelfree)
library(MASS)      
library(sjPlot)
library(simr)
library(Rmisc)

theme_set(theme_bw()) # set ggplot theme to b&w
rootdir='D:\\switchdrive\\Projects\\PhD\\respiration\\github\\pavo2020'
setwd(rootdir)


computeCohenD <- function(n, mean1, sd1, mean2, sd2) {
  denom = (((n-1)*sd1^2 + (n-1)*sd2^2) /(n+n-2))^0.5
  d = (mean1 - mean2)/denom 
  
  return(d)
}




#######################################
#######################################
##                                   ##
##                                   ##
##           RESPIRATION             ##
##                                   ##
##                                   ##
#######################################
#######################################

# ---------
# load data
# ---------

a<-read.csv(file='breathing_selfVoice_data.csv',header=FALSE,sep=',') 
names(a)<-c('subject', 'task', 'condition', 'correct', 'phase', 'reactionTime', 'respiration', 'taskValue', 'response')

a$subject = as.factor(a$subject-2)
a$task = ifelse(a$task==1,'loudness','voice')
a$condition = ifelse(a$condition==1,'ASYNC','SYNC')

a$correctFactor = ifelse(a$correct==1,'yes','no')
a$respiration = ifelse(a$respiration==1,'inspiration','expiration')

a$sine = sin(a$phase)

a$condition = as.factor(a$condition)
a$task = as.factor(a$task)
a$correctFactor = as.factor(a$correctFactor)
a$respiration = as.factor(a$respiration)
a$comtask[a$task=="loudness"] = transform(a[a$task=="loudness",],id=as.numeric(factor(taskValue)))$id
a$comtask[a$task=="voice"] = transform(a[a$task=="voice",],id=as.numeric(factor(taskValue)))$id  

a$responseTime = a$reactionTime
a$response = ifelse(a$response==1,0,1) # in the MATLAB script 1 was left button (reference), 3 the right one (answer)

males = c('1', '2', '12', '17', '18', '23', '24', '27', '28')
a$gender = ifelse(is.element(a$subject,males), 'male', 'female')
a$gender = as.factor(a$gender)


a = as.data.frame(a)

a_voice = a[a$task=='voice',]
a_loudness = a[a$task=='loudness',]


# ----------------
# outliers removal
# ----------------
# Trials with reaction times greater or smaller than two interquartile ranges from the median for each subject
# were considered as outliers and excluded. 

# loudness dataset
trial_l = 1
while (trial_l != 0) {
  a_loudness <- dplyr::left_join(a_loudness, a_loudness %>% dplyr::group_by(subject) %>%
                                   dplyr::summarise(limsup = median(responseTime)+2*IQR(responseTime),liminf = median(responseTime)-2*IQR(responseTime)))
  a_loudness$rm = ifelse(a_loudness$responseTime<a_loudness$liminf | a_loudness$responseTime>a_loudness$limsup,1,0)
  trial_l = mean(a_loudness$rm) # trial exclusion
  
  a_loudness = a_loudness %>% filter(rm==0) %>% dplyr::select(-rm,-limsup,-liminf)
}

# voice dataset
trial_v = 1
while (trial_v != 0) {
  a_voice <- dplyr::left_join(a_voice, a_voice %>% dplyr::group_by(subject) %>%
                                dplyr::summarise(limsup = median(responseTime)+2*IQR(responseTime),liminf = median(responseTime)-2*IQR(responseTime)))
  a_voice$rm = ifelse(a_voice$responseTime<a_voice$liminf | a_voice$responseTime>a_voice$limsup,1,0)
  trial_v = mean(a_voice$rm) # trial exclusion
  
  a_voice = a_voice %>% filter(rm==0) %>% dplyr::select(-rm,-limsup,-liminf)
}




# --------------
# GLMER response
# --------------

# By-subject random slopes for the main effects were added following model selection based on maximum likelihood. 
m_voice = glmer(response ~  condition * respiration * comtask  + (comtask|subject), data=a_voice,family = binomial)
m_loudness = glmer(response ~  condition * respiration * comtask  + (1|subject), data=a_loudness,family = binomial)

summary(m_voice)
summary(m_loudness)


m_voice_async = glmer(response ~  respiration * comtask  + (comtask|subject), data=a_voice[a_voice$condition=='ASYNC',],family = binomial)
m_voice_sync = glmer(response ~  respiration * comtask  + (comtask|subject), data=a_voice[a_voice$condition=='SYNC',],family = binomial)
summary(m_voice_async)
summary(m_voice_sync)



# --------------
# gender effects
# --------------

m_voice_gen = glmer(response ~  condition * respiration * comtask + gender * respiration   + (1|subject), data=a_voice,family = binomial)
summary(m_voice_gen)
anova(m_voice_gen, m_voice)



# ------------------
# LMER response time
# ------------------

mrt_voice = lmer(responseTime ~ respiration  * condition * poly(comtask,2) + (condition|subject), data=a_voice, REML = FALSE)
mrt_loudness = lmer(responseTime ~ respiration  * condition * poly(comtask,2) + (condition+respiration|subject), data=a_loudness, REML = FALSE)

summary(mrt_voice)
summary(mrt_loudness)



# --------
# quickpsy
# --------

library(quickpsy) 


fit <- quickpsy(a_voice[a_voice$condition=='ASYNC',],comtask,response,within=.(respiration),random = .(subject),fun=logistic_fun, 
                prob = .5,guess=T,lapses=T,parini = list(c(1,6),c(0.05,6),c(0,.3),c(0, .3)),bootstrap='none') # parini: thresh,slope,guess,lapse

var = as.data.frame(fit$par)
aov_car(par ~  parn * respiration  + Error(subject|parn*respiration), data=var)

n = 28

# PSE
ttest_pse = t.test(var[var$parn=='p1' & var$respiration=='inspiration',]$par, var[var$parn=='p1' & var$respiration=='expiration',]$par, paired=T)
pse_p = ttest_pse["p.value"][[1]]
d_pse = computeCohenD(n, 
                      summary(var[var$parn=='p1' & var$respiration=='inspiration',]$par)[[4]],
                      sd(var[var$parn=='p1' & var$respiration=='inspiration',]$par),
                      summary(var[var$parn=='p1' & var$respiration=='expiration',]$par)[[4]],
                      sd(var[var$parn=='p1' & var$respiration=='expiration',]$par))
PSE = pwr::pwr.t.test(n, d_pse)


# slope
ttest_slope = t.test(var[var$parn=='p2' & var$respiration=='inspiration',]$par, var[var$parn=='p2' & var$respiration=='expiration',]$par, paired=T)
slope_p = ttest_slope["p.value"][[1]]
d_slope = computeCohenD(n, 
                        summary(var[var$parn=='p2' & var$respiration=='inspiration',]$par)[[4]],
                        sd(var[var$parn=='p2' & var$respiration=='inspiration',]$par),
                        summary(var[var$parn=='p2' & var$respiration=='expiration',]$par)[[4]],
                        sd(var[var$parn=='p2' & var$respiration=='expiration',]$par))
slope = pwr::pwr.t.test(n, d_slope)


# left asymptote
ttest_la = t.test(var[var$parn=='p3' & var$respiration=='inspiration',]$par, var[var$parn=='p3' & var$respiration=='expiration',]$par, paired=T)
la_p = ttest_la["p.value"][[1]]
d_la = computeCohenD(n, 
                     summary(var[var$parn=='p3' & var$respiration=='inspiration',]$par)[[4]],
                     sd(var[var$parn=='p3' & var$respiration=='inspiration',]$par),
                     summary(var[var$parn=='p3' & var$respiration=='expiration',]$par)[[4]],
                     sd(var[var$parn=='p3' & var$respiration=='expiration',]$par))
la = pwr::pwr.t.test(n, d_la)


# right asymptote
ttest_ra = t.test(var[var$parn=='p4' & var$respiration=='inspiration',]$par, var[var$parn=='p4' & var$respiration=='expiration',]$par, paired=T)
ra_p = ttest_ra["p.value"][[1]]
d_ra = computeCohenD(n, 
                     summary(var[var$parn=='p4' & var$respiration=='inspiration',]$par)[[4]],
                     sd(var[var$parn=='p4' & var$respiration=='inspiration',]$par),
                     summary(var[var$parn=='p4' & var$respiration=='expiration',]$par)[[4]],
                     sd(var[var$parn=='p4' & var$respiration=='expiration',]$par))
ra = pwr::pwr.t.test(n, d_ra) 



# -----------------------------------
# BRMS bayes factor for null effects
# -----------------------------------

library(brms)
chains=4
iter=10000
warmup = 2000


# prior: slope from self-other task
priors <- c(set_prior("normal(-1.02, 1.56)", class = "b", coef= "respirationinspiration"),
            set_prior("normal(0.52, 0.54)", class = "b", coef= "comtask"),
            set_prior("normal(0.24, 0.4)", class = "b", coef= "respirationinspiration:comtask"))

# loudness task
m_l_brm = brm(response ~ condition * respiration * comtask + (1 | subject), data = a_loudness,
         prior=priors,family = brms::bernoulli(), sample_prior=T,
         chains=chains,iter=iter,warmup = warmup, cores=chains,save_all_pars = TRUE,
         file='loudness_brms')
h <- c("respiration*comtask" = "respirationinspiration:comtask  = 0")
(hyp_loudness=hypothesis(m_l_brm, h))


# prior with mean zero
priors_zero <- c(set_prior("normal(-1.02, 1.56)", class = "b", coef= "respirationinspiration"),
            set_prior("normal(0.52, 0.54)", class = "b", coef= "comtask"),
            set_prior("normal(0, 0.4)", class = "b", coef= "respirationinspiration:comtask"))

m_l_zero = brm(response ~ condition * respiration * comtask + (1 | subject), data = a_loudness,
              prior=priors_zero,family = brms::bernoulli(), sample_prior=T,
              chains=chains,iter=iter,warmup = warmup, cores=chains,save_all_pars = TRUE,
              file='loudness_brms_zero')
(hyp_loudness=hypothesis(m_l_zero, h))

# synchronous condition
m2_sync = brm(response ~  respiration * comtask  + (comtask|subject), data=a_voice[a_voice$condition=='SYNC',]
              ,prior=priors,family = brms::bernoulli(), sample_prior=T,
              chains=chains,iter=iter,warmup = warmup, cores=chains,save_all_pars = TRUE,
              file='sync_with_priors')
(hyp_sync=hypothesis(m2_sync, h))


# ------
# sjplot
# ------

avgdata_voice = aggregate(response ~ subject + comtask + respiration,a_voice,mean)
pt_voice=plot_model(m_voice,type = "pred",terms=c('comtask','respiration'),title='', alpha=.1)
voice_plot = pt_voice + 
  geom_hline(yintercept=.5, linetype='dashed',alpha=.3) +
  xlab('\nStimulus') +  ylab(expression(paste("Ratio of ", italic("self"), " responses"))) +
  ylim(0, 1) +
  ggtitle('Self-other') +
  scale_x_discrete(name ="Self-voice (%)", limits=c("15","30","45", "55","70","85")) +
  scale_color_manual(name="Respiration", values=c('#F8766D', '#00BFC4'),
                     labels=c("Expiration", "Inspiration")) +
  scale_fill_manual(name="Respiration", values=c('#F8766D', '#00BFC4'),
                    labels=c("Expiration", "Inspiration")) +
  theme(axis.text.x = element_text(colour="grey20",size=13,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="grey20",size=13,angle=0,hjust=1,vjust=0,face="plain"),
        axis.title.x = element_text(colour="grey20",size=14,angle=0,face="plain", margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=14,angle=90,face="plain", margin = margin(t = 0, r = 10, b = 0, l = 0)),
        strip.text = element_text(size=14),
        plot.title = element_text(colour="grey20",size=16,angle=0,hjust=.5,vjust=.5,face="plain"),
        #panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position="right",
        legend.text = element_text(size=11),
        legend.title=element_text(size=12),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = 'black'))    + 
  geom_line(size=1.1, linetype="solid") +
  stat_summary(data=avgdata_voice ,aes(x = comtask, y = response, color = respiration),fun.data=mean_cl_normal,geom='point',inherit.aes = F)

voice_plot

avgdata_voice_async = aggregate(response ~ subject + comtask + respiration,a_voice[a_voice$condition=='ASYNC',],mean)
pt_async=plot_model(m_voice_async,type = "pred",terms=c('comtask','respiration'),title='', alpha=.1)
async_plot = pt_async + geom_hline(yintercept=.5, linetype='dashed',alpha=.3) +
  xlab('\nStimulus') +  ylab(expression(paste("Ratio of ", italic("self"), " responses"))) +
  ylim(0, 1) +
  ggtitle('Asynchronous') +
  scale_x_discrete(name ="Self-voice (%)", limits=c("15","30","45", "55","70","85")) +
  scale_color_manual(name="Respiration", values=c('#F8766D', '#00BFC4'),
                     labels=c("expiration", "inspiration")) +
  theme(axis.text.x = element_text(colour="grey20",size=13,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="grey20",size=13,angle=0,hjust=1,vjust=0,face="plain"),
        axis.title.x = element_text(colour="grey20",size=14,angle=0,face="plain", margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=14,angle=90,face="plain", margin = margin(t = 0, r = 10, b = 0, l = 0)),
        strip.text = element_text(size=14),
        plot.title = element_text(colour="grey20",size=16,angle=0,hjust=.5,vjust=.5,face="plain"),
        #panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position="none",
        legend.text = element_text(size=11),
        legend.title=element_text(size=12),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = 'black'))    + 
  geom_line(size=1.1, linetype="solid") +
  stat_summary(data=avgdata_voice_async ,aes(x = comtask, y = response, color = respiration),fun.data=mean_cl_normal,geom='point',inherit.aes = F)


avgdata_voice_sync = aggregate(response ~ subject + comtask + respiration,a_voice[a_voice$condition=='SYNC',],mean)
pt_sync=plot_model(m_voice_sync,type = "pred",terms=c('comtask','respiration'),title='', alpha=.1)
sync_plot = pt_sync + geom_hline(yintercept=.5, linetype='dashed',alpha=.3) +
  #geom_hline(yintercept=.5, linetype='dashed',alpha=.3) +
  xlab('\nStimulus') +  ylab(expression(paste("Ratio of ", italic("self"), " responses"))) +
  #ylab('') +
  ylim(0, 1) +
  ggtitle('Synchronous') +
  scale_x_discrete(name ="Self-voice (%)", limits=c("15","30","45", "55","70","85")) +
  scale_color_manual(name="Respiration", values=c('#F8766D', '#00BFC4'),
                     labels=c("expiration", "inspiration")) +
  theme(axis.text.x = element_text(colour="grey20",size=13,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="grey20",size=13,angle=0,hjust=1,vjust=0,face="plain"),
        axis.title.x = element_text(colour="grey20",size=14,angle=0,face="plain", margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=14,angle=90,face="plain", margin = margin(t = 0, r = 10, b = 0, l = 0)),
        strip.text = element_text(size=14),
        plot.title = element_text(colour="grey20",size=16,angle=0,hjust=.5,vjust=.5,face="plain"),
        #panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position="none",
        legend.text = element_text(size=11),
        legend.title=element_text(size=12),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = 'black'))    + 
  geom_line(size=1.1, linetype="solid") +
  stat_summary(data=avgdata_voice_sync ,aes(x = comtask, y = response, color = respiration),fun.data=mean_cl_normal,geom='point',inherit.aes = F)


grid.arrange(async_plot, sync_plot, nrow=1)

avgdata_loudness = aggregate(response ~ subject + comtask + respiration,a_loudness,mean)
pt_loud=plot_model(m_loudness,type = "pred",terms=c('comtask','respiration'),title='', alpha=.1)
loud_plot = pt_loud + geom_hline(yintercept=.5, linetype='dashed',alpha=.3) +
  #geom_hline(yintercept=.5, linetype='dashed',alpha=.3) +
  xlab('\nStimulus') +  ylab(expression(paste("Ratio of ", italic("loud"), " responses"))) +
  ylim(0, 1) +
  ggtitle('Loudness') +
  scale_x_discrete(name ="Sound intensity (dBFS)", limits=c("-14","-13","-12.5", "-11.5","-11","-10")) +
  scale_color_manual(name="Respiration", values=c('#F8766D', '#00BFC4'),
                     labels=c("Expiration", "Inspiration")) +
  theme(axis.text.x = element_text(colour="grey20",size=13,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="grey20",size=13,angle=0,hjust=1,vjust=0,face="plain"),
        axis.title.x = element_text(colour="grey20",size=14,angle=0,face="plain", margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=14,angle=90,face="plain", margin = margin(t = 0, r = 10, b = 0, l = 0)),
        strip.text = element_text(size=14),
        plot.title = element_text(colour="grey20",size=16,angle=0,hjust=.5,vjust=.5,face="plain"),
        #panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position="right",
        legend.text = element_text(size=11),
        legend.title=element_text(size=12),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = 'black'))    + 
  geom_line(size=1.1, linetype="solid") +
  stat_summary(data=avgdata_loudness ,aes(x = comtask, y = response, color = respiration),fun.data=mean_cl_normal,geom='point',inherit.aes = F)

loud_plot




#####################################
#####################################
##                                 ##
##                                 ##
##         QUESTIONNAIRE           ##
##                                 ##
##                                 ##
#####################################
#####################################

a_q = c()
a_q<-read.csv(file='q_long.csv',header=FALSE,sep=',')
names(a_q)<-c('Subjects', 'Condition', 'ConditionOrder', 'Questions', 'Ratings')

a_q$Subjects = as.factor(a_q$Subjects-2) #because they were from 3 to 32
a_q$Condition = as.factor(a_q$Condition)
a_q$Questions = as.factor(a_q$Questions)


#keep only questionnaire questions
questionsOut = c('q4_loc', 'hand','pq_1', 'pq_2')


a_q = a_q[!is.element(a_q$Questions,questionsOut),]

a_q_all = a_q
# subj 9 and 4 are missing physiology data (device error)
a_q = a_q[a_q$Subjects!='4' & a_q$Subjects!='9',]


# ------------
# mixed models
# ------------
# They reduce to t-tests. This is just a quick check 

# significant questions
summary(lmer(Ratings~Condition+(1|Subjects),data=a_q[a_q$Questions=='q1',], REML = FALSE)) # self-touch
summary(lmer(Ratings~Condition+(1|Subjects),data=a_q[a_q$Questions=='q2',], REML = FALSE)) # other-touch

# other questions
summary(lmer(Ratings~Condition+(1|Subjects),data=a_q[a_q$Questions=='q3',], REML = FALSE)) # control
summary(lmer(Ratings~Condition+(1|Subjects),data=a_q[a_q$Questions=='q4',], REML = FALSE)) # presence hallucination
summary(lmer(Ratings~Condition+(1|Subjects),data=a_q[a_q$Questions=='q5',], REML = FALSE)) # friend-voice
summary(lmer(Ratings~Condition+(1|Subjects),data=a_q[a_q$Questions=='q6',], REML = FALSE)) # friend-presence
summary(lmer(Ratings~Condition+(1|Subjects),data=a_q[a_q$Questions=='q7',], REML = FALSE)) # other-presence
summary(lmer(Ratings~Condition+(1|Subjects),data=a_q[a_q$Questions=='q8',], REML = FALSE)) # other-voice
summary(lmer(Ratings~Condition+(1|Subjects),data=a_q[a_q$Questions=='q9',], REML = FALSE)) # bodily sensations


# -----------------------------
# t-test: effect size and power
# -----------------------------
# They are one-tailed, because the effects are known from before (Blanke et al., 2014; Salomon et al., 2020)

n = 28

# self touch
avg_st = aggregate(Ratings ~ Subjects + Condition ,a_q[a_q$Questions=='q1',],mean)
ttest_st = t.test(avg_st[avg_st$Condition=='S',]$Ratings, avg_st[avg_st$Condition=='A',]$Ratings, paired=T, alternative = "greater")
st_p = ttest_st["p.value"][[1]]
d_st = computeCohenD(n, 
                      summary(a_q[a_q$Questions=='q1' & a_q$Condition=='S',]$Ratings)[[4]],
                      sd(a_q[a_q$Questions=='q1' & a_q$Condition=='S',]$Ratings),
                      summary(a_q[a_q$Questions=='q1' & a_q$Condition=='A',]$Ratings)[[4]],
                      sd(a_q[a_q$Questions=='q1' & a_q$Condition=='A',]$Ratings))
self_touch = pwr::pwr.t.test(n, d_st, alternative="greater") 


# other touch (somatic passivity)
avg_ot = aggregate(Ratings ~ Subjects + Condition ,a_q[a_q$Questions=='q2',],mean)
ttest_ot = t.test(avg_ot[avg_ot$Condition=='A',]$Ratings, avg_ot[avg_ot$Condition=='S',]$Ratings, paired=T, alternative = "greater")
ot_p = ttest_ot["p.value"][[1]]/2
d_ot = computeCohenD(n, 
                     summary(a_q[a_q$Questions=='q2' & a_q$Condition=='S',]$Ratings)[[4]],
                     sd(a_q[a_q$Questions=='q2' & a_q$Condition=='S',]$Ratings),
                     summary(a_q[a_q$Questions=='q2' & a_q$Condition=='A',]$Ratings)[[4]],
                     sd(a_q[a_q$Questions=='q2' & a_q$Condition=='A',]$Ratings))
other_touch = pwr::pwr.t.test(n, d_ot, alternative="greater") 



# -----------------------------------
# summaries for significant questions
# -----------------------------------

CI(a_q$Ratings[a_q$Questions=='q1' & a_q$Condition=='A'])
CI(a_q$Ratings[a_q$Questions=='q1' & a_q$Condition=='S'])

CI(a_q$Ratings[a_q$Questions=='q2' & a_q$Condition=='A'])
CI(a_q$Ratings[a_q$Questions=='q2' & a_q$Condition=='S'])



# ---------------------------------------------------
# identify the subjects who experienced the illusions
# ---------------------------------------------------

tbo_1 = c()           # touch by other (somatic passivity)
st_1 = c()            # self-touch

tbo_pos = c()         # positive touch by other
st_pos = c()          # positive self touch

sub_que_1 = c()       # to check if it is sequentially taken

i = 1
for (s in unique(a_q$Subjects)) {
    aqs = a_q[a_q$Subjects==s,]
    sub_que_1 = c(sub_que_1, s)
    
    tbo_1[i] = aqs$Ratings[aqs$Questions=='q2' & aqs$Condition=='A'] - aqs$Ratings[aqs$Questions=='q2' & aqs$Condition=='S']
    st_1[i] = aqs$Ratings[aqs$Questions=='q1' & aqs$Condition=='A'] - aqs$Ratings[aqs$Questions=='q1' & aqs$Condition=='S']

    if (tbo_1[i] > 0) { tbo_pos = c(tbo_pos, s) }
    if (st_1[i] < 0) { st_pos = c(st_pos, s) }      # self-touch is actually always higher in sync
    
    i = i+ 1
}



# ---------------------------
# binary covariates to glmers
# ---------------------------

# it is better without condition, because 'condition' is incorporated in passivity already (as Async - Sync)

a_voice$tbo = as.factor(ifelse(is.element(a_voice$subject,tbo_pos),'Passivity+', 'Passivity-'))
a_voice$st = as.factor(ifelse(is.element(a_voice$subject,st_pos),'ST+', 'ST-'))
a_voice = a_voice[a_voice$subject!='4',] # he/she had physiology recordings in only 1 condition - not balances


# Passivity

m_tbo = glmer(response ~   tbo  * comtask * respiration + (1|subject), data=a_voice,family = binomial)
summary(m_tbo)

m_tbo_plus = glmer(response ~ comtask * respiration + (comtask|subject), data=a_voice[a_voice$tbo=='Passivity+',],family = binomial)
summary(m_tbo_plus)

m_tbo_minus = glmer(response ~ comtask * respiration + (1|subject), data=a_voice[a_voice$tbo=='Passivity-',],family = binomial)
summary(m_tbo_minus)

avgdata_voice_tbo_plus = aggregate(response ~ subject + comtask + respiration ,a_voice[a_voice$tbo=='Passivity+',],mean)
pt_voice_tbo_plus=plot_model(m_tbo_plus,type = "pred",terms=c('comtask','respiration'),title='', alpha=.1)
plot_voice_tbo_plus = pt_voice_tbo_plus + 
  geom_hline(yintercept=.5, linetype='dashed',alpha=.3) +
  xlab('\nStimulus') +  ylab(expression(paste("Ratio of ", italic("self"), " responses"))) +
  ylim(0, 1) +
  ggtitle('Passivity+') +
  scale_x_discrete(name ="Self-voice (%)", limits=c("15","30","45", "55","70","85")) +
  scale_color_manual(name="Respiration", values=c('#F8766D', '#00BFC4'),
                     labels=c("Expiration", "Inspiration")) +
  theme(axis.text.x = element_text(colour="grey20",size=13,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="grey20",size=13,angle=0,hjust=1,vjust=0,face="plain"),
        axis.title.x = element_text(colour="grey20",size=14,angle=0,face="plain", margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=14,angle=90,face="plain", margin = margin(t = 0, r = 10, b = 0, l = 0)),
        strip.text = element_text(size=14),
        plot.title = element_text(colour="grey20",size=16,angle=0,hjust=.5,vjust=.5,face="plain"),
        #panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position="none",
        legend.text = element_text(size=11),
        legend.title=element_text(size=12),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = 'black'))    + 
  geom_line(size=1.1, linetype='solid') +
  stat_summary(data=avgdata_voice_tbo_plus ,aes(x = comtask, y = response, color = respiration),fun.data=mean_cl_normal,geom='point',inherit.aes = F)


avgdata_voice_tbo_minus = aggregate(response ~ subject + comtask + respiration ,a_voice[a_voice$tbo=='Passivity-',],mean)
pt_voice_tbo_minus=plot_model(m_tbo_minus,type = "pred",terms=c('comtask','respiration'),title='', alpha=.1)
plot_voice_tbo_minus = pt_voice_tbo_minus + 
  geom_hline(yintercept=.5, linetype='dashed',alpha=.3) +
  xlab('\nStimulus') +  ylab(expression(paste("Ratio of ", italic("self"), " responses"))) +
  ylim(0, 1) +
  ggtitle('Passivity-') +
  scale_x_discrete(name ="Self-voice (%)", limits=c("15","30","45", "55","70","85")) +
  scale_color_manual(name="Respiration", values=c('#F8766D', '#00BFC4'),
                     labels=c("Expiration", "Inspiration")) +
  theme(axis.text.x = element_text(colour="grey20",size=13,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="grey20",size=13,angle=0,hjust=1,vjust=0,face="plain"),
        axis.title.x = element_text(colour="grey20",size=14,angle=0,face="plain", margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=14,angle=90,face="plain", margin = margin(t = 0, r = 10, b = 0, l = 0)),
        strip.text = element_text(size=14),
        plot.title = element_text(colour="grey20",size=16,angle=0,hjust=.5,vjust=.5,face="plain"),
        #panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position="none",
        legend.text = element_text(size=11),
        legend.title=element_text(size=12),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = 'black'))    + 
  geom_line(size=1.1, linetype='solid') +
  stat_summary(data=avgdata_voice_tbo_minus ,aes(x = comtask, y = response, color = respiration),fun.data=mean_cl_normal,geom='point',inherit.aes = F)

grid.arrange(plot_voice_tbo_plus, plot_voice_tbo_minus, nrow=1)


# Self-touch

m_st = glmer(response ~   st  * comtask * respiration  + (respiration|subject), data=a_voice,family = binomial)
summary(m_st)




#################################
#################################
##                             ##
##                             ##
##         HEARTBEAT           ##
##                             ##
##                             ##
#################################
#################################

a<-read.csv(file='heartbeat_selfVoice_data.csv',header=FALSE,sep=',') 
names(a)<-c('subject', 'task', 'condition', 'taskValue', 'response', 'phase', 'responseTime', 'correct', 'binPhase')

a$subject = as.factor(a$subject)
a$task = as.factor(a$task)
a$condition = as.factor(a$condition)
a$binPhase = as.factor(a$binPhase)

a$comtask[a$task=="loudness"] = transform(a[a$task=="loudness",],id=as.numeric(factor(taskValue)))$id
a$comtask[a$task=="voice"] = transform(a[a$task=="voice",],id=as.numeric(factor(taskValue)))$id 

a = a[a$binPhase!='other',] # keep only systole & diastole
a_voice = a[a$task=='voice',]
a_loudness = a[a$task=='loudness',]


# --------------
# outlier trials
# --------------

# loudness dataset
trial_l = 1
while (trial_l != 0) {
  a_loudness <- dplyr::left_join(a_loudness, a_loudness %>% dplyr::group_by(subject) %>%
                                   dplyr::summarise(limsup = median(responseTime)+2*IQR(responseTime),liminf = median(responseTime)-2*IQR(responseTime)))
  a_loudness$rm = ifelse(a_loudness$responseTime<a_loudness$liminf | a_loudness$responseTime>a_loudness$limsup,1,0)
  trial_l = mean(a_loudness$rm) # trial exclusion
  
  a_loudness = a_loudness %>% filter(rm==0) %>% dplyr::select(-rm,-limsup,-liminf)
}

# voice dataset
trial_v = 1
while (trial_v != 0) {
  a_voice <- dplyr::left_join(a_voice, a_voice %>% dplyr::group_by(subject) %>%
                                dplyr::summarise(limsup = median(responseTime)+2*IQR(responseTime),liminf = median(responseTime)-2*IQR(responseTime)))
  a_voice$rm = ifelse(a_voice$responseTime<a_voice$liminf | a_voice$responseTime>a_voice$limsup,1,0)
  trial_v = mean(a_voice$rm) # trial exclusion
  
  a_voice = a_voice %>% filter(rm==0) %>% dplyr::select(-rm,-limsup,-liminf)
}


# ------------
# mixed models
# ------------

m_voice = glmer(response ~  condition * binPhase * comtask  + (condition|subject), data=a_voice,family = binomial)
m_loudness = glmer(response ~  condition * binPhase * comtask  + (binPhase|subject), data=a_loudness,family = binomial) 
summary(m_voice)
summary(m_loudness)

m_rt_v = lmer(responseTime ~ binPhase  * condition * poly(comtask,2) + (condition|subject), data=a_voice, REML = FALSE)
m_rt_l = lmer(responseTime ~ binPhase  * condition * poly(comtask,2) + (condition|subject), data=a_loudness, REML = FALSE)
summary(m_rt_v)
summary(m_rt_l)


# ----------------------------------------
# bayes factor for null effect of loudness
# ----------------------------------------
library(brms)
chains=4
iter=10000
warmup = 2000

priors_heart <- c(set_prior("normal(-1.02, 1.56)", class = "b", coef= "binPhasesystole"),
            set_prior("normal(0.52, 0.54)", class = "b", coef= "comtask"),
            set_prior("normal(0.24, 0.4)", class = "b", coef= "binPhasesystole:comtask"))

# voice

m_heart_voice = brm(response ~  condition * binPhase * comtask  + (condition|subject), data=a_voice,
         prior=priors_heart,family = brms::bernoulli(), sample_prior=T,
         chains=chains,iter=iter,warmup = warmup, cores=chains,save_all_pars = TRUE,
         file='heartbeat_voice')
h <- c("binPhase*comtask" = "binPhasesystole:comtask  = 0")
(hyp_heart_voice=hypothesis(m_heart_voice, h))


# loudness

m_heart_loudness = brm(response ~  condition * binPhase * comtask  + (binPhase|subject), data=a_loudness,
                    prior=priors_heart,family = brms::bernoulli(), sample_prior=T,
                    chains=chains,iter=iter,warmup = warmup, cores=chains,save_all_pars = TRUE,
                    file='heartbeat_loudness')
h <- c("binPhase*comtask" = "binPhasesystole:comtask  = 0")
(hyp_heart_loudness=hypothesis(m_heart_loudness, h))




##################################
##################################
##                              ##
##                              ##
##         PHYSIOLOGY           ##
##                              ##
##                              ##
##################################
##################################


breathingRate<-read.csv(file='breathingRate.csv',sep=',') 
breathingRate$subject = as.factor(breathingRate$subject)
breathingRate$condition = as.factor(breathingRate$condition)
breathingRate$task = as.factor(breathingRate$task)

breathingRateVar<-read.csv(file='breathingRateVar.csv',sep=',') 
breathingRateVar$subject = as.factor(breathingRateVar$subject)
breathingRateVar$condition = as.factor(breathingRateVar$condition)
breathingRateVar$task = as.factor(breathingRateVar$task)

heartrate<-read.csv(file='heartrate.csv',sep=',') 
heartrate$subject = as.factor(heartrate$subject)
heartrate$condition = as.factor(heartrate$condition)
heartrate$task = as.factor(heartrate$task)

hrv<-read.csv(file='hrv.csv',sep=',') 
hrv$subject = as.factor(hrv$subject)
hrv$condition = as.factor(hrv$condition)
hrv$task = as.factor(hrv$task)


summary(aov_car(rate ~   condition*task   + Error(subject| condition*task ), data=breathingRate))
summary(aov_car(rateVar ~   condition*task   + Error(subject| condition*task ), data=breathingRateVar))
summary(aov_car(heartrate ~   condition*task   + Error(subject| condition*task ), data=heartrate))
summary(aov_car(hrv ~   condition*task   + Error(subject| condition*task ), data=hrv))


# -------------
# CI: condition
# -------------

# mean(breathingRate[breathingRate$condition=='Asynchronous',]$rate)
# sd(breathingRate[breathingRate$condition=='Asynchronous',]$rate)
# mean(breathingRate[breathingRate$condition=='Synchronous',]$rate)
# sd(breathingRate[breathingRate$condition=='Synchronous',]$rate)
CI(breathingRate[breathingRate$condition=='Asynchronous',]$rate)
CI(breathingRate[breathingRate$condition=='Synchronous',]$rate)

# mean(breathingRateVar[breathingRateVar$condition=='Asynchronous',]$rateVar)
# sd(breathingRateVar[breathingRateVar$condition=='Asynchronous',]$rateVar)
# mean(breathingRateVar[breathingRateVar$condition=='Synchronous',]$rateVar)
# sd(breathingRateVar[breathingRateVar$condition=='Synchronous',]$rateVar)
CI(breathingRateVar[breathingRateVar$condition=='Asynchronous',]$rateVar)
CI(breathingRateVar[breathingRateVar$condition=='Synchronous',]$rateVar)


# mean(heartrate[heartrate$condition=='Asynchronous',]$heartrate)
# sd(heartrate[heartrate$condition=='Asynchronous',]$heartrate)
# mean(heartrate[heartrate$condition=='Synchronous',]$heartrate)
# sd(heartrate[heartrate$condition=='Synchronous',]$heartrate)
CI(heartrate[heartrate$condition=='Asynchronous',]$heartrate)
CI(heartrate[heartrate$condition=='Synchronous',]$heartrate)


# mean(hrv[hrv$condition=='Asynchronous',]$hrv)
# sd(hrv[hrv$condition=='Asynchronous',]$hrv)
# mean(hrv[hrv$condition=='Synchronous',]$hrv)
# sd(hrv[hrv$condition=='Synchronous',]$hrv)
CI(hrv[hrv$condition=='Asynchronous',]$hrv)
CI(hrv[hrv$condition=='Synchronous',]$hrv)


# --------
# CI: task
# --------


CI(breathingRate[breathingRate$task=='loudness',]$rate)
CI(breathingRate[breathingRate$task=='voice',]$rate)


CI(breathingRateVar[breathingRateVar$task=='loudness',]$rateVar)
CI(breathingRateVar[breathingRateVar$task=='voice',]$rateVar)


CI(heartrate[heartrate$task=='loudness',]$heartrate)
CI(heartrate[heartrate$task=='voice',]$heartrate)


CI(hrv[hrv$task=='loudness',]$hrv)
CI(hrv[hrv$task=='voice',]$hrv)

