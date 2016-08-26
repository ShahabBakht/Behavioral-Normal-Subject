library(lme4)  # load library
library(arm)  # convenience functions for regression in R
library(R.matlab)
library(ggplot2)
library(ggthemes)

numPoints <- 5
numSubjects <- 9
numSamples <- numPoints * numSubjects


tble <- readMat("D:\\Analysis\\Behavioral-Normal-Subject\\For R\\DATA_2.mat")
DATA.2 <- tble$DATA.2
X.2 <- DATA.2[,2,]
Y.2 <- DATA.2[,1,]
X.2.c <- matrix(t(X.2),nrow = numSamples,ncol = 1)
Y.2.c <- matrix(t(Y.2),nrow = numSamples,ncol = 1)
C.2 <- vector()
for (i in 1:numSubjects){C.2 <- c(C.2,array(i,c(numPoints,1)))}
# C.2 <- c(matrix(1,5,1),matrix(2,5,1),matrix(3,5,1),matrix(4,5,1),+
#             matrix(5,5,1),matrix(6,5,1),matrix(7,5,1),matrix(8,5,1),matrix(9,5,1),matrix(10,5,1))
C.2.c <- matrix(C.2,nrow = numSamples,ncol = 1)
  
tble <- readMat("D:\\Analysis\\Behavioral-Normal-Subject\\For R\\DATA_6.mat")
DATA.6 <- tble$DATA.6
X.6 <- DATA.6[,2,]
Y.6 <- DATA.6[,1,]
X.6.c <- matrix(t(X.6),nrow = numSamples,ncol = 1)
Y.6.c <- matrix(t(Y.6),nrow = numSamples,ncol = 1)
C.6 <- vector()
for (i in 1:numSubjects){C.6 <- c(C.6,array(i,c(numPoints,1)))}
C.6.c <- matrix(C.6,nrow = numSamples,ncol = 1)

tble <- readMat("D:\\Analysis\\Behavioral-Normal-Subject\\For R\\DATA_20.mat")
DATA.20 <- tble$DATA.20
X.20 <- DATA.20[,2,]
Y.20 <- DATA.20[,1,]
X.20.c <- matrix(t(X.20),nrow = numSamples,ncol = 1)
Y.20.c <- matrix(t(Y.20),nrow = numSamples,ncol = 1)
C.20 <- vector()
for (i in 1:numSubjects){C.20 <- c(C.20,array(i,c(numPoints,1)))}
C.20.c <- matrix(C.20,nrow = numSamples,ncol = 1)


####################################################################################################################
## Linear modeling

GoFolr <- matrix(0,3,1)
GoFvsi <- matrix(0,3,1)
GoFvi <- matrix(0,3,1)
Coeff <- matrix(0,3,1)

### linear fit for 2 degrees
df1 <- data.frame(X.2.c,C.2.c,Y.2.c)

# ordinary linear regression
MLexamp.2 <- glm(Y.2.c ~ X.2.c, data = df1)
# GoFolr[1] <- AIC(MLexamp1.1)
GoFolr[1] <- BIC(MLexamp.2)

# varying slope and intercept linear regression with lmer library
REexamp1.2 <- lmer(Y.2.c ~ X.2.c + (1+X.2.c|C.2.c), data = df1) 
# GoFvsi[1] <- AIC(REexamp1.1)
GoFvsi[1] <- BIC(REexamp1.2)
betas <-  fixef(REexamp1.2)
Coeff[1]<-betas[2]

# varying intercept linear regression with lmer library
REexamp2.2 <- lmer(Y.2.c ~ X.2.c + (1|C.2.c), data = df1)
# GoFvi[1] <- AIC(REexamp1.2)
GoFvi[1] <- BIC(REexamp2.2)
CI.2 <- confint(REexamp2.2,"X.2.c",level = 0.8)

####################################################################################################################
### linear fit for 6 degrees
df1 <- data.frame(X.6.c,C.6.c,Y.6.c)

# ordinary linear regression
MLexamp.6 <- glm(Y.6.c ~ X.6.c, data = df1)
# GoFolr[2] <- AIC(MLexamp1.1)
GoFolr[2] <- BIC(MLexamp.6)

# varying slope and intercept linear regression with lmer library
REexamp1.6 <- lmer(Y.6.c ~ X.6.c + (1+X.6.c|C.6.c), data = df1) 
# GoFvsi[2] <- AIC(REexamp1.1)
GoFvsi[2] <- BIC(REexamp1.6)
betas <-  fixef(REexamp1.6)
Coeff[2]<-betas[2]

# varying intercept linear regression with lmer library
REexamp2.6 <- lmer(Y.6.c ~ X.6.c + (1|C.6.c), data = df1)
# GoFvi[2] <- AIC(REexamp1.2)
GoFvi[2] <- BIC(REexamp2.6)
CI.6 <- confint(REexamp2.6,"X.6.c",level = 0.8)

####################################################################################################################
### linear fit for 20 degrees
df1 <- data.frame(X.20.c,C.20.c,Y.20.c)

# ordinary linear regression
MLexamp.20 <- glm(Y.20.c ~ X.20.c, data = df1)
# GoFolr[3] <- AIC(MLexamp1.1)
GoFolr[3] <- BIC(MLexamp.20)

# varying slope and intercept linear regression with lmer library
REexamp1.20 <- lmer(Y.20.c ~ X.20.c + (1+X.20.c|C.20.c), data = df1) 
GoFvsi[3] <- AIC(REexamp1.20)
GoFvsi[3] <- BIC(REexamp1.20)
betas <-  fixef(REexamp1.20)
Coeff[3]<-betas[2]


# varying intercept linear regression with lmer library
REexamp2.20 <- lmer(Y.20.c ~ X.20.c + (1|C.20.c), data = df1)
# GoFvi[3] <- AIC(REexamp1.2)
GoFvi[3] <- BIC(REexamp2.20)
CI.20 <- confint(REexamp2.20,"X.20.c",level = 0.8)



####################################################################################################################
### Compare models 
# anova(REexamp1.2)
# anova(REexamp1.6)
# anova(REexamp1.20)
### plot the results

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

x2 <- X.2.c
y2 <- Y.2.c
d2 <- data.frame(x2, y2)
x6 <- X.6.c
y6 <- Y.6.c
d6 <- data.frame(x6, y6)
x20 <- X.20.c
y20 <- Y.20.c
d20 <- data.frame(x20, y20)


p <- ggplot(d2, aes(x2, y2)) + geom_point(size=2.5,colour = cbbPalette[1]) + geom_rangeframe() + theme_tufte() +
  xlab("Target Velocity") + ylab("Eye Velocity") + 
  theme(axis.title.x = element_text(vjust=-0.5), axis.title.y = element_text(vjust=1.5)) 
p <- p + geom_point(data=d6,aes(x6, y6) ,size=2.5,colour = cbbPalette[4])  
p <- p + geom_point(data=d20,aes(x20, y20) ,size=2.5,colour = cbbPalette[7]) 
print(p)

p <- ggplot(d2, aes(x2, y2)) + geom_point(size=2,colour = cbbPalette[1]) + geom_rug(colour = cbbPalette[1]) + theme_tufte(ticks=F) + 
  xlab("Target Velocity") + ylab("Eye Velocity") + 
  theme(axis.title.x = element_text(vjust=-0.5), axis.title.y = element_text(vjust=1))
p <- p + geom_point(data=d20,aes(x6, y6),size=2,colour = cbbPalette[4]) + geom_rug(data=d20,aes(x6, y6),colour = cbbPalette[4])
p <- p + geom_point(data=d20,aes(x20, y20),size=2,colour = cbbPalette[7]) + geom_rug(data=d20,aes(x20, y20),colour = cbbPalette[7])
# to add another layer for the lines estimated from the fixed effects

betas <- fixef(REexamp1.2)
predicty.2 <- betas[1]+betas[2]*(5:20)
xtest <- 5:20
testdata <- data.frame(xtest,predicty.2)
p <- p + geom_line(aes(x=xtest,y = predicty.2),data=testdata,colour = cbbPalette[1])


betas <- fixef(REexamp1.6)
predicty.6 <- betas[1]+betas[2]*(5:20)
xtest <- 5:20
testdata <- data.frame(xtest,predicty.6)
p <- p + geom_line(aes(x=xtest,y = predicty.6),data=testdata,colour = cbbPalette[4])

betas <- fixef(REexamp1.20)
predicty.20 <- betas[1]+betas[2]*(5:20)
xtest <- 5:20
testdata <- data.frame(xtest,predicty.20)
p <- p + geom_line(aes(x=xtest,y = predicty.20),data=testdata,colour = cbbPalette[7])
p <- p + annotate("text", x = 30, y = 19, adj=1,  family="serif",parse = TRUE,size = 8,
                  label = c("italic(y[i] == alpha[ij]%.%x[i] + beta[j] +e[i])"))


print(p)
ggsave(file="D:\\Analysis\\Behavioral-Normal-Subject\\From R\\all_subjects_data.eps", fonts=c("serif", "Palatino"))



# x <- c(x2,x6,x20)
# y <- c(y2,y6,y20)
# group <- c(matrix("2",1,45),matrix("6",1,45),matrix("20",1,45))
# dataall <- data.frame(x, y, group)
# bp <- ggplot(data=dataall, aes(x=x, y=y, fill=group)) + geom_rug(aes(colour=group)) + geom_point(aes(colour=group),size=2) + xlab("Target Velocity") + ylab("Eye Velocity") + theme_tufte() 
# bp <- bp + geom_line(aes(x = x, y = c(predict(REexamp2.2),predict(REexamp2.6),predict(REexamp2.20))),size=1)
# print(bp)

### plot coefficients


df <- data.frame(
  diameter = (c(2, 6, 20)),
  group = factor(c(1, 1, 1)),
  slope = Coeff
)
limits <- aes(ymax = c(CI.2[2],CI.6[2],CI.20[2]), ymin=c(CI.2[1],CI.6[1],CI.20[1]))
dodge <- position_dodge(width=0.9)
p <- ggplot(df, aes(y=slope, x=diameter))
p <- p + geom_point(size=4) 
p <- p + geom_line(aes(group=group))
p <- p + geom_errorbar(limits, position=dodge, width=0.1)
p <- p + theme_tufte()
print(p)
ggsave(file="D:\\Analysis\\Behavioral-Normal-Subject\\From R\\all_slopes.eps", fonts=c("serif", "Palatino"))

####################################################################################################################
### Fit model to each subject
all.coeff <- matrix(0,numSubjects,3)
conditionList <- c(2,6,20)
CI <- array(0,c(2,numSubjects,3))
for (condcount in 1:3){
  thisCond <- conditionList[condcount]
  tble <- readMat(paste("D:\\Analysis\\Behavioral-Normal-Subject\\For R\\DATA_",thisCond,".mat",sep = ""))
  DATA <- tble$DATA
  for (subcount in 1:numSubjects){
    
    X <- DATA[subcount,2,]
    Y <- DATA[subcount,1,]
    X <-matrix(X,numPoints,1)
    Y <-matrix(Y,numPoints,1)
    df <- data.frame(X,Y)
    MLexamp <- glm(Y ~ X, data = df)
    CI[,subcount,condcount] <- confint(MLexamp,"X")
    all.coeff[subcount,condcount] <- MLexamp$coefficients[2]
    
  }
  
}


### plot one selected subject
# plot regression line



whichSubject <- 4
df <- data.frame(
  diameter = (c(2, 6, 20)),
  group = factor(c(1, 1, 1)),
  slope = all.coeff[whichSubject,]
)
limits <- aes(ymax = c(CI[2,whichSubject,1],CI[2,whichSubject,2],CI[2,whichSubject,3]), ymin=c(CI[1,whichSubject,1],CI[1,whichSubject,2],CI[1,whichSubject,3]))
dodge <- position_dodge(width=0.9)
p <- ggplot(df, aes(y=slope, x=diameter,colour = c("red","green","black")))
p <- p + geom_point(aes(y=slope, x=diameter),data = df, size=5,colour=c(cbbPalette[1],cbbPalette[1],cbbPalette[1])) 
p <- p + geom_line(aes(y=slope, x=diameter,group=group),data = df,colour="black",size=1)
p <- p + geom_errorbar(limits, position=dodge, width=0.1,colour=c(cbbPalette[1],cbbPalette[1],cbbPalette[1])) + scale_x_continuous(breaks = c(2,6,20))
p <- p + theme_tufte()


print(p)


### plot all coefficients
for (subcount in 1:numSubjects){
  df <- data.frame(
    diameter = (c(2, 6, 20)),
    group = factor(c(1, 1, 1)),
    slope = all.coeff[subcount,]
  )
    limits <- aes(ymax = c(CI[2,subcount,1],CI[2,subcount,2],CI[2,subcount,3]), ymin=c(CI[1,subcount,1],CI[1,subcount,2],CI[1,subcount,3]))
    dodge <- position_dodge(width=0.9)
    if (subcount < 2){p <- ggplot(df, aes(y=slope, x=diameter))+ scale_x_continuous(breaks = c(2,6,20))}
    p <- p + geom_point(aes(y=slope, x=diameter),data = df, size=4,colour=cbPalette[1]) 
    p <- p + geom_line(aes(y=slope, x=diameter,group=group),data = df,colour=cbPalette[1])
    #p <- p + geom_errorbar(limits, position=dodge, width=0.1)
    p <- p + theme_tufte()
  
}
# add the mixed-model result on top
df <- data.frame(
  diameter = (c(2, 6, 20)),
  group = factor(c(1, 1, 1)),
  slope = Coeff
)
limits <- aes(ymax = c(CI.2[2],CI.6[2],CI.20[2]), ymin=c(CI.2[1],CI.6[1],CI.20[1]))
dodge <- position_dodge(width=0.9)
p <- p + geom_point(aes(y=slope, x=diameter),data = df,size=6,colour=cbbPalette[1]) 
p <- p + geom_line(aes(y=slope, x=diameter,group=group),size = 1,data = df,colour=cbbPalette[1])
p <- p + geom_errorbar(limits, position=dodge, width=0.1,colour=cbbPalette[1])
p <- p + theme_tufte()
print(p)
ggsave(file="D:\\Analysis\\Behavioral-Normal-Subject\\From R\\all_subjects_slopes.eps", fonts=c("serif", "Palatino"))

#####
# plot quantiles and regression for one subject

alldata_sb <- readMat("D:\\Analysis\\Behavioral-Normal-Subject\\Final Clean Data\\SVC_ls.mat")
SVC <- alldata_sb$SVC
SVCoI <- c(SVC[3,],SVC[4,])
TV <- alldata_sb$TV
TVoI <- c(TV[3,],TV[4,])
dfsvc <- data.frame(SVCoI,TVoI)
p <- ggplot(dfsvc, aes(TVoI, SVCoI)) + geom_point(size=1.5,colour = cbPalette[1],na.rm=TRUE) +  
  #theme_tufte() +
  xlab("Target Velocity") + ylab("Eye Velocity")



whichSubject <-4
whichCond <- 2 # 6 degrees
thisCond <- conditionList[whichCond]
tble <- readMat(paste("D:\\Analysis\\Behavioral-Normal-Subject\\For R\\DATA_",thisCond,".mat",sep = ""))
DATA <- tble$DATA
X <- DATA[whichSubject,2,]
Y <- DATA[whichSubject,1,]
X <-matrix(X,numPoints,1)
Y <-matrix(Y,numPoints,1)
df <- data.frame(X,Y)
MLexamp <- glm(Y ~ X, data = df)
xtest <- 5:20
ytest <- MLexamp$coefficients[1] + MLexamp$coefficients[2] * xtest
dftest <- data.frame(xtest,ytest)
#p <- ggplot(df, aes(X, Y))
p <- p + geom_point(aes(x=X, y=Y),data=df,size = 6,colour = cbbPalette[4]) +
  xlab("Target Velocity") + ylab("Eye Velocity") + 
  theme(axis.title.x = element_text(vjust=-0.5), 
        axis.title.y = element_text(vjust=1.5),
        panel.background = element_rect(fill = "white")
        #axis.line.y = element_line(colour = "black", size = .5),
        #axis.line.x = element_line(colour = "black", size = .5)
  )
p <- p + geom_line(aes(x=xtest,y = ytest),data=dftest,colour = cbbPalette[4],linetype = 2,size = 1) 
print(p)
ggsave(file="D:\\Analysis\\Behavioral-Normal-Subject\\From R\\linearfit_ls_6degrees.eps", fonts=c("serif", "Palatino"))

#####
# statistics of the best fit (here mixed effect with only varying intercept: REexami.j)

# to get a combination of fixed and random effects:
COEFs.mixed.2 <- coef(REexamp2.2)
COEFs.mixed.6 <- coef(REexamp2.6)
COEFs.mixed.20 <- coef(REexamp2.20)

# to get only the fixed effects:
COEFs.fixed.2 <- coef(summary(REexamp2.2))
COEFs.fixed.6 <- coef(summary(REexamp2.6))
COEFs.fixed.20 <- coef(summary(REexamp2.20))

# to get only the random effects:
COEFs.random.2 <- ranef(REexamp2.2)
COEFs.random.6 <- ranef(REexamp2.6)
COEFs.random.20 <- ranef(REexamp2.20)

#####
# anova for mixed models fit with and with condition effect

X.c <- c(X.2.c,X.6.c,X.20.c)
Y.c <- c(Y.2.c,Y.6.c,Y.20.c)
C.c <- c(C.2.c,C.6.c,C.20.c)
COND.c <- c(matrix(2,1,45),matrix(6,1,45),matrix(20,1,45))

df <- data.frame(X.c,C.c,Y.c,COND.c)
REexamp.noCONDeffect <- lmer(Y.c ~ X.c + (1+X.c|C.c), data = df1) 
REexamp.withCONDeffect <- lmer(Y.c ~ X.c + (1+X.c|C.c) + (1+X.c|COND.c), data = df1) 
REexamp.withCONDeffect2 <- lmer(Y.c ~ X.c + (1+X.c|C.c) + (1+X.c|COND.c) + (1+X.c|COND.c:C.c) , data = df1)
REexamp.withCONDeffect2 <- lmer(Y.c ~ X.c + (1+X.c|C.c) + (1+X.c|COND.c) + (1+X.c|COND.c:C.c) , data = df1)

REexamp.noCONDeffect <- lmer(Y.c ~ X.c + (1|C.c), data = df1) 
REexamp.withCONDeffect <- lmer(Y.c ~ X.c + (1|C.c) + (1|COND.c), data = df1) 
REexamp.withCONDeffect2 <- lmer(Y.c ~ X.c + (1|C.c) + (1|COND.c) + (1|COND.c:C.c) , data = df1)
REexamp.withCONDeffect2 <- lmer(Y.c ~ X.c + (1|C.c) + (1|COND.c) + (1|COND.c:C.c) , data = df1)
