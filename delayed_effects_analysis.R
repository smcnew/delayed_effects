# Analysis for McNew et al. 2019(20?) Oecologia


#Packages
library(lmerTest)
library(lme4)
library(MASS)
library(plotrix)
library(smatr)
library(reshape2)
library(ggplot2)
library(dplyr)
library(sjPlot)
library(car)
library(glmmTMB)


# Data --------------------------------------------------------------------

# Data from nests in 2015, with just repeat nests from round 2
nests2015M <- read.csv("./CSVs/2015study2.csv")
nests2015M$round <- as.factor(nests2015M$round)
nests2015M$dead <- nests2015M$nestlings-nests2015M$fledglings
nests2015M$hatchdate <- as.Date(nests2015M$hatchdate, format= "%m/%d/%y")
nests2015M$ndvi.date <- as.Date(nests2015M$ndvi.date, format= "%m/%d/%y")
nests2015M$failure.date <- as.Date(nests2015M$failure.date, format= "%d-%b-%y")
nests2015M$renest.date <- as.Date(nests2015M$renest.date, format= "%m/%d/%y")
nests2015M$deadeggs <- nests2015M$eggs- nests2015M$nestlings
head(nests2015M)

#Parent conditions data from both years Each parent has a row. No duplicated parents.
parents <- read.csv("./CSVs/parents1516.csv")
parents$hatch.date <- as.Date(parents$hatch.date, format= "%m/%d/%y")
parents$Date.second.nest <- as.Date(parents$Date.second.nest, format= "%m/%d/%y")
parents$year <- as.factor(parents$year)

#Physiology for just 2015, including just 2nd round nests that were repeats
#Replicate nest, 1 = 1 nest studied, 2 = 2 nests studied.
nestlings2015 <- read.csv("./CSVs/nestlings2015M.csv")
nestlings2015$nest.round <- paste(nestlings2015$nest,nestlings2015$round,sep=".")
nestlings2015$round <- as.factor(nestlings2015$round)

time2015 <- read.csv("./CSVs/behavior-2015.csv")
time2015[is.na(time2015)] <- 0
time2015$standing <- time2015$standing.at.rim + time2015$standing.in.nest
time2015$cleaning <- time2015$sanitation + time2015$allopreening
time2015$present <- 100-time2015$out.of.sight
time2015 <- time2015[time2015$first.rnd==1,]
time2015$round <- as.factor(time2015$round)


#
# Reproductive success  ---------------------------------------------------

glmmTMB(intensity ~ treatment + round + (1|id),
        ziformula = ~treatment, data=nests2015M, family="nbinom1") %>% summary()

glmer.nb(intensity ~ treatment + round + (1|id),
         data=nests2015M) %>% summary()

aggregate(intensity ~ round, data=nests2015M[nests2015M$treatment=="W",], std.error)

#Eggs
cmp.eggs <- glmmTMB(eggs ~ round * treatment + (1|id), data=nests2015M, family='compois')
tab_model(cmp.eggs, transform="exp")
aggregate(nestlings ~ round + treatment, data=nests2015M, mean)


#Brood size
glm.cmp(nestlings ~ treatment * round, data=nests2015M)
lmer(nestlings ~ treatment*round+(1|id), data=nests2015M) %>%
  tab_model() #Nestlings: P increased but W did not in round 2

glmer(cbind(nestlings, deadeggs) ~ treatment*round + (1|id),
      family="binomial", data=nests2015M) %>% tab_model(transform=NULL)

#Fledging success
glmer.nb(fledglings ~ treatment*round+(1|id), data=nests2015M) %>%
  tab_model(transform="exp") #Fledgies: 2 round both treatments increase, interaction ns

glmer(cbind(fledglings,dead) ~ treatment*round+(1|id), data=nests2015M,
      family="binomial") %>% tab_model(transform=NULL) #% fledging success improves both treat

aggregate(fledglings ~ treatment + round, nests2015M, mean)

#Nestling mass
#dplyr::select(nests2015M, nest.id, nestlings, eggs) %>% rename(nest.round = nest.id) %>%
#  merge(., nestlings2015, all.y=T) -> nestlings2015 #adding clutch size though I don't think this is necessary

lmer(mass ~ treatment * round +  (1|nest),
     data=nestlings2015[nestlings2015$visit==1,])  %>% tab_model()

#
# Parent condition and survival -----------------------------------------------------------------

head(nests2015M)
aggregate(hatchdate ~ treatment, nests2015M[nests2015M$round==1,], median)
aggregate(hatchdate ~ treatment, nests2015M[nests2015M$round==1,], min)
aggregate(hatchdate ~ treatment, nests2015M[nests2015M$round==1,], max)

#probability of renesting

summary(lm(time.to.renest~treatment, data=nests2015M)) ##NS
glm(as.factor(repeat.) ~ treatment,
    data = nests2015M[nests2015M$round ==1,], family="binomial") %>% tab_model() #NS

#Scaled mass of parents
lmer(scaledmass ~ nestage + sex + treatment + nest1egg +(1|year) , data=parents) %>%
  tab_model() #Full model

lmer(scaledmass ~ nestage +(1|year) , data=parents) %>%
  tab_model() #Best model


#immunology from both years combined
lmer(phil.215cal ~  sex  + nest1egg + (1|year), data=parents) %>% tab_model() #bestish model
lmer(phil.215cal ~  sex  + nest1egg + treatment  + nestage +
       (1|year), data=parents) %>% tab_model()


##Survival
#Full model
glm(present.2016 ~  sex + treatment + nest1nestling + hatch.date,
    data=parents[parents$year==2015,], family="binomial") %>% tab_model(transform=NULL)
#Minimal model
glm(present.2016 ~  sex,
    data=parents[parents$year==2015,], family="binomial") %>% tab_model(transform=NULL)




#
# Behavior ----------------------------------------------------------------

mean(nests2015M$nestlings)


lmer(feeding.nestlings ~ treatment * round + (1|nest), time2015) %>% tab_model()
lmer(feeding.nestlings ~ treatment * round + nestlings +  (1|nest), time2015) %>% tab_model()
lmer(sitting ~ treatment * round + (1|nest), time2015) %>% tab_model() #nestlings NS
lmer(cleaning ~ treatment * round   + (1|nest), time2015) %>% tab_model()#nestlings NS
lmer(standing ~ treatment * round + (1|nest), time2015) %>% tab_model() #nestlings NS
lmer(present ~ treatment * round + (1|nest), time2015) %>% tab_model() #nestlings NS


# Plots -------------------------------------------------------------------
par(mfrow=c(1,1), mar=c(5,5,4,3))
xlab1 <- "Nesting Attempt"
s.cex=2
yupper=.1
ylower=.3
lwd=2
l.cex=2
a.cex=2
cols<-c("gray21", "gray82")
error.bar <- function(x, mean, err, lwd) {
  arrows(x0=x,y0=mean+err, y1=mean-err,
         length=0.10, angle=90, code=3, lwd=lwd)
}


dotplot2<- function(data, mean1, colorindex) {
  meanwrap <- paste(substitute(mean1))
  mean <- with(data, get(meanwrap))
  errorwrap<-paste(substitute(mean1), "err", sep=".")
  error <-with(data, get(errorwrap))
  plot(NA,NA, xlim=c(.5,length(mean)/2+1), ylim=c((min(mean)-ylower*min(mean)),
                                                  (max(mean)+yupper*max(mean))), xaxt="n", xlab= xlab1,
       cex.axis=a.cex, cex.lab=l.cex, ylab=ylab1, main=meanwrap, bty="l")
  axis(1, at = c(1:(length(mean)/2)), labels = c("1st", "2nd"), cex.axis=l.cex)
  arrows(x0 = c(sort(rep(1:(length(mean)/2),2)))[1:2],  y0 = mean[1:2],
         x1 = c(sort(rep(1:(length(mean)/2),2)))[3:4],  y1 = mean[3:4], length=0,
         lwd=lwd)
  error.bar(c(sort(rep(1:(length(mean)/2),2))), mean, error, lwd)
  points(c(sort(rep(1:(length(mean)/2),2))), mean, bg = cols[colorindex],
         pch=21, cex=s.cex)


  # legend("topright", legend = c("Permethrin", "Water"), pch=21, cex=.8,
  #pt.bg = cols[colorindex])
}

#Aggregate data for plots
massagg<-aggregate(mass~treatment+round, data=nestlings2015[nestlings2015$visit==1,], mean)
massagg$mass.err <- aggregate(mass~treatment*round, data=nestlings2015[nestlings2015$visit==1,], std.error)[,3]

nests2015agg<- aggregate(cbind(eggs, nestlings, fledglings)~treatment+round, data = nests2015M, mean)
nests2015err<- aggregate(cbind(eggs, nestlings, fledglings)~treatment+round, data = nests2015M, std.error)
colnames(nests2015err) <- paste(colnames(nests2015agg),"err", sep=".")
nests2015agg <- cbind(nests2015agg, nests2015err[,-c(1:2)])
nests2015agg <- nests2015agg[,c(2,1,3:8)]

#Figure 1
pdf("oecologia_fig1.pdf", useDingbats = F)
par(mfrow=c(2,2), mar=c(5,5,4,4))
ylab1 = "Clutch size"
dotplot2(nests2015agg, eggs, nests2015agg$treatment)
ylab1 = "Brood size"
dotplot2(nests2015agg, nestlings, nests2015agg$treatment)
ylab1 = "Mass (g)"
dotplot2(massagg, mass, massagg$treatment)
ylab1 = "Fledglings"
dotplot2(nests2015agg, fledglings, nests2015agg$treatment)
dev.off()


#AGGREGATE BEHAVIOR DATA
time2015agg <- aggregate(cbind(feeding.nestlings, sitting, cleaning, standing) ~
                           treatment + round, data=time2015, mean)
time2015agg.err <- aggregate(cbind(feeding.nestlings, sitting, cleaning, standing) ~
                               treatment + round, data=time2015, std.error)
colnames(time2015agg.err) <- paste(colnames(time2015agg), "err", sep=".")
time2015agg <- cbind(time2015agg, time2015agg.err[,-c(1,2)])
time2015agg <- time2015agg[,c(2,1,3:10)]

#Figure 2
yupper=.2
ylower=.5
l.cex = 1.8
a.cex = 1.8
pdf("oecologia_fig2.pdf", useDingbats = F)
par(mfrow=c(2,2), mar=c(5,5,4,4))
ylab1 = "Percent time"
dotplot2(time2015agg, feeding.nestlings, nests2015agg$treatment)
dotplot2(time2015agg, sitting, nests2015agg$treatment)
dotplot2(time2015agg, cleaning, nests2015agg$treatment)
dotplot2(time2015agg, standing, nests2015agg$treatment)
dev.off()


### figure 3 parent mass
pdf("Oecologia_fig3.pdf", useDingbats = F)
par(mfrow=c(1,1), mar=c(5,5,4,4))
plot(scaledmass~nestage, data=parents, bg=cols[treatment],  cex=1.5, ylab="Scaled mass of parents (g)", xlab="Age of clutch (days)",
     cex.lab=2, cex.axis=1.5, bty="l", pch = 21)
legend("topright", legend=c("Fumigated", "Sham-fumigated"), pch=21, pt.bg=cols, cex=1.5)
abline(lm(scaledmass~nestage, data=parents), lwd=2)
dev.off()


