#FINAL PROJECT MATH E-23C
#JEAN-SEBASTIEN PAUL
#EXPLORING COVID-19 DEATHRATES


#NOTES:

#Correlation heatmap code learnt from:
#http://www.sthda.com/english/wiki/ggplot2-quick-correlation-matrix-heatmap-r-software-and-data-visualization


#PACKAGES:

#install.packages("plyr")
library(plyr)
#install.packages("fitdistrplus")
library(fitdistrplus)
#install.packages("utils")
library(utils)
#install.packages("httr")
library(httr)
#install.packages("TSstudio")
library(TSstudio)
#install.packages("tidyverse")
library(tidyverse)
# install.packages("chron")
library(chron)
#install.packages("expss")
library(expss)
#install("ggplot2")
library(ggplot2)
#install.packages("reshape2")
library(reshape2)
#install.packages("stats4")
library(stats4)
#install.packages("RColorBrewer")
library(RColorBrewer)
#install.packages("latticeExtra")
library(latticeExtra)
#install.packages("Hmisc")     
library(Hmisc)
#install.packages("sem")
library(sem)
#install.packages("rgl")
library(rgl)
#install.packages("multcomp")
library(multcomp)
#install.packages("leaps")
library(leaps)
#install.packages("aplpack")
library(aplpack)
#install.packages("Rcmdr")
library(Rcmdr)
#install.packages("MASS")
library(MASS)
#install.packages("car")
library(car)
#install.packages("quantmod")
library(quantmod)
#install.packages("nnet")
library(nnet)
#install.packages("neuralnet")
library(neuralnet)
#install.packages("glmnet")
library(glmnet)
#install.packages("miscTools")
library(miscTools)


#SETTING SEED FOR CONSISTENT RESULTS
set.seed(3.141592)



#FUNCTIONS
#(REQ: Professional Looking Software Engineering - functions)

#formatJHU:
#Formats JHU COVID 19 data so that it is formatted like
#Europa Open Data
formatJHU=function(JHUdata,title){
  JHUdata = arrange(JHUdata, Country.Region)
  
  track=1
  Temp1=subset(JHUdata,,-c(Province.State,Lat,Long))
  Temp2=subset(Temp1,,-c(Country.Region)) #using for col summing
  Temp1=Temp1[0,] #going to fill this in
  k=1
  attach(JHUdata)
  for(i in 2:nrow(JHUdata)){
    if(Country.Region[i]!=Country.Region[i-1]){
      Temp1[k,1]=Country.Region[i-1]
      Temp1[k,2:ncol(Temp1)]=colSums(Temp2[track:(i-1),])
      k=k+1
      track=i
    }
  }
  JHUdata=Temp1
  colnames(JHUdata)[1]="countriesAndTerritories"
  
  detach(JHUdata)
  
  start=as.Date(colnames(JHUdata[2]),format="X%m.%d.%y")
  latest=as.Date(colnames(JHUdata[ncol(JHUdata)]),format="X%m.%d.%y")
  seq(start:latest)
  Date=rep(seq(start,latest,by="days"),nrow(JHUdata))
  
  countries=c()
  for(i in 1:nrow(JHUdata)){
    countries=c(countries, rep(toString(JHUdata$countriesAndTerritories[i]),(ncol(JHUdata)-1)))
  }
  results=c()
  for(i in 1:nrow(JHUdata)){
    results=c(results, as.integer(JHUdata[i,2:ncol(JHUdata)]))
  }
  
  data <- data.frame(countries, results, Date)
  colnames(data)[2]=title
  colnames(data)[1]="countriesAndTerritories"
  return(data)
}

#dccum:
#determines cumulative deaths and cases for data 
#that follows a structure as in 'geographicdata'
#and data2
dccum = function(datain){
  data=datain
  N=nrow(data)
  track=data$countriesAndTerritories[N]
  dsum=0
  csum=0
  #exploiting how data is arranged by time
  for(k in 1:N){
    if(track==data$countriesAndTerritories[k]){
      dsum=dsum+geographicdata$deaths[k]
      csum=csum+geographicdata$cases[k]
    }
    else{
      track=data$countriesAndTerritories[k]
      dsum=data$deaths[k]
      csum=data$cases[k]
    }
    data$deaths2[k]=dsum
    data$cases2[k]=csum
  }
  return(data)
}

#data2Setup:
#Sets up data2 once predictors, testing and geographicdata
#combined. Also for use when reading data in as a .csv
data2Setup = function(data){
  data <- subset(data, select = -c(Tests, Test.Pop,Quarantine,Schools,Restrictions,Total.Recovered,Total.Deaths,Total.Infected,popData2018))
  
  #(lose rows with NA Values in data)
  data2=drop_na(data)
  data2$days=data2$Date-data2$Date[which.min(data2$Date)]
  
  #Categorical Variables:
  #Rich beign 1 defined as being above mean GDP/capita. 
  #Crime being 1 defined as being above mean crime
  #index (indicates more crime)
  #(REQ: At least 2 categorical or logical columns)
  
  data2$Rich=as.numeric(data2$GDP.2018/data2$Population.2020>(mean(data2$GDP.2018/data2$Population.2020))*1+0)
  data2$Crime=as.numeric(data2$Crime.Index>(mean(data2$Crime.Index))*1+0)
  return(data2)
}

#dateCol:
#Takes a dataframe of the same format as Geographic Data and data2
#(sharing time columns) and generates a column of dates for these
dateCol = function(data){
  data$Date= as.Date(paste(data$year,data$month,data$day,sep="-"),"%Y-%m-%d")
  return(data)
}

#ChiSq:
#Returns chi squared value
ChiSq <-function(Obs,Exp){
  sum((Obs-Exp)^2/Exp)
}

#get_upper_tri:
#Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

#%.%:
#Dot product function
"%.%" <- function(x,y) sum(x*y)

#linearm:
#Performs linear regression
linearm = function(predictors,dependent){
  v1=rep(1,nrow(predictors))
  A <- cbind(v1,predictors)
  B <- t(A)%*%A; B
  P <- A%*%solve(B)%*%t(A)
  return(solve(B)%*%t(A)%*%dependent)
}

#ggbar:
#Returns a ggplot which would plot a side by side barplot
#of the 2 inputted vectors, deaths and model. X axis noted 
#as days.
ggbar = function(a,b){
  df=as.data.frame(cbind(a,b))
  df$days=seq(1:nrow(df))
  colnames(df)[1]="Real Total Deaths"
  colnames(df)[2]="Model"
  df = melt(df, id.vars=c("days"))
  
  return(ggplot(df, aes(days, value, fill=variable)) 
         + geom_bar(stat='Identity',position=position_dodge())
         +ggtitle("Modelling the time series data for total COVID-19 deaths"))
}

#dropvif1:
#For use in model building section. Takes linear model
#Finds largest vif. If greater than 10, it drops the relevant
#variable. Returns data with dropped variable
dropvif1=function(fit,train.data){
  maxvif=which.max(vif(fit))
  if(vif(fit)[[maxvif[[1]] ]]>10){
    print(paste("Dropped",colnames(train.data)[maxvif[[1]]]))
    train.data <- subset(train.data, select=-c(maxvif[[1]]))
  }
  return(train.data)
}

#normalize:
#Normalizes data for neural net models
normalize=function(x){
  return((x-min(x))/(max(x) -min(x) ))
}


#DATA SCRAPING
#NOTE: A CSV file has also been provided. This is here if you
#are at all interested at how I collected and formatted that data

#We will have 2 main dataframes:
#geographicdata-in depth data in recoveries, cases and deaths
#data2-data with predictors and testing data, but less in depth
#in recoveries, cases and deaths, due to lack of data in some 
#of those predictors and rows with NA values being dropped


# #Get data from internet
# JHURecovs=read.csv("https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_recovered_global.csv")
# predictors=read.csv("https://raw.githubusercontent.com/jspaul2003/nCoV2019/master/more_data/covid19_by_country.csv")
# testing=read.csv("https://raw.githubusercontent.com/jspaul2003/nCoV2019/master/testing.csv")
# GET("https://opendata.ecdc.europa.eu/covid19/casedistribution/csv", authenticate(":", ":", type="ntlm"), write_disk(tf <- tempfile(fileext = ".csv")))
# geographicdata=read.csv(tf)
# 
# #Formatting this data
# JHURecovs1=formatJHU(JHURecovs,"recoveries")
# predictors$countriesAndTerritories=predictors$Country
# testing$countriesAndTerritories=testing$Entity
# testing$Date=as.Date(testing$Date, format= "%b %d, %y")
# geographicdata=dateCol(geographicdata)
# geographicdata=geographicdata%>%right_join(JHURecovs1, by=c("countriesAndTerritories","Date"))
# geographicdata=drop_na(geographicdata)
# 
# #Making Our Big Datasets
# #(REQ: A dataframe, At least 2 numeric columns, A data set with lots of columns,
# #allowing comparison of many different variables.)
# geographicdata$days=geographicdata$Date-geographicdata$Date[which.min(geographicdata$Date)]
# 
# #determine cumulative deaths and cases for geographic data-> use this
# #to find active cases, deathrate
# geographicdata=dccum(geographicdata)
# geographicdata$active=geographicdata$cases2-geographicdata$deaths2-geographicdata$recoveries
# geographicdata$deathrate=geographicdata$deaths2/(geographicdata$cases2+1*(geographicdata$cases2==0))
# 
# 
# data2=geographicdata %>% inner_join(testing, by=c("countriesAndTerritories","Date"))
# data2=predictors %>% right_join(data2, by=c("countriesAndTerritories"))
# data2=data2Setup(data2)

#READ FROM CSV (ALTERNATIVE TO DATASCRAPING)
#(REQ: A dataframe, At least 2 numeric columns, A data set with lots of columns, 
#allowing comparison of many different variables.)

data=read.csv("data.csv")
data=dateCol(data)
geographicdata=data[,which(colnames(data)=="countriesAndTerritories"):ncol(data)]
geographicdata$days=geographicdata$Date-geographicdata$Date[which.min(geographicdata$Date)]
geographicdata <- subset(geographicdata, select = -c(Entity, Code, Total.tests))
data2=data2Setup(data)
data2 <- subset(data2, select = -c(X))

#(REQ: At least 20 rows, A data set so large it can be used as a population from 
#which samples are taken)

nrow(geographicdata)
nrow(data2)


#ANALYSIS

#PART 1:
#Exploration regarding Distribution of Variables related to deathrate from COVID-19
#using probability distributions.

#I)
#Can we model the time series data of new deaths world wide as binomial distribution?
#We will look at the interval between days 1 to 106. We will use the more in depth 
#in deathrate, geographicdata dataframe. It will contain over 50 countries.
#(REQ: A barplot,Nicely labeled graphics using ggplot, with good use of color, line 
#styles...)

#max days
max(geographicdata$days)

world=geographicdata[1:max(geographicdata$days),]
world$deathrate=1

a=subset(geographicdata,geographicdata$days==1)
logic=geographicdata$days==max(geographicdata$days)&geographicdata$countriesAndTerritories%in%a$countriesAndTerritories
a=subset(geographicdata,logic)

#number of countries at days 1 and last day
nrow(a)

for(i in 1:(max(geographicdata$days))){
  index=which(geographicdata$days==i&(geographicdata$countriesAndTerritories%in%a$countriesAndTerritories))
  index2=which(geographicdata$days==i)
  
  world$Date[i]=i+min(geographicdata$Date)
  world$days[i]=i
  
  world$deathrate[i]=sum(geographicdata$deaths2[index2])/(sum(geographicdata$cases2[index2])+1*(sum(geographicdata$cases2[index2])==0))
  world$deaths2[i]=sum(geographicdata$deaths2[index])
  world$cases2[i]=sum(geographicdata$cases2[index])
  world$deaths[i]=sum(geographicdata$deaths[index])
  world$cases[i]=sum(geographicdata$cases[index])
  world$active[i]=sum(geographicdata$active[index])
}
world$days=as.numeric(world$days)

#Could we plot a binomial distribution for deaths
model=world$active*0.17*pbinom(1:nrow(world),(round(world$active/1000)),mean(world$deathrate))
ggbar(world$deaths2,model)+xlab("Days since 2020-01-22")+ylab("Total COVID-19 Deaths")

#Seems like a very weak model.

#How about a poisson model?
fitdistr(world$deaths2, "poisson")
model=world$active*0.17*ppois(1:nrow(world), 17.83948)
ggbar(world$deaths2,model)+xlab("Days since 2020-01-22")+ylab("Total COVID-19 Deaths")

#Seems just as weak as the binomial model.

#II)
#How are death rates distributed: do they converge at a certain value or do they vary wildly?
#(REQ: A histogram, A probability density function overlaid on a histogram, A p-value or other
#statistic based on a distribution function, Appropriate use of R functions for a probability 
#distribution other than binomial, normal, or chi-square, Nicely labeled graphics using ggplot, 
#with good use of color, line styles...)

toplot=subset(geographicdata, deaths2!=0)
p=ggplot(toplot, aes(x=deathrate))  +
  geom_histogram(aes(y=..density..), binwidth=0.01, colour="steelblue", fill="white",bins="fd") +
  ylab("Density")+xlab("Death Rate")+ggtitle("Density Histogram of death rate")
p

#Lets try modelling with a gamma function
fitdistr(geographicdata$deathrate[which(geographicdata$deaths2!=0)], "gamma")
p+stat_function(fun=dgamma, args=list(shape=1.06011822, rate=21.99330112),col="red")+xlim(0,1)
#Seems good

#Lets test if we can indeed model this with a gamma function


#create bins by breaking data into deciles
bins=qgamma(0.1*(0:10),1.06011822 ,21.99330112)

binstuff=cut(geographicdata$deathrate[which(geographicdata$deaths2!=0)], breaks=bins,labels=F); binstuff

obs=as.vector(table(binstuff))

exp=rep(sum(obs)/10,10)
chisq=ChiSq(obs,exp); chisq

#P-value
#10 Categories, imposed 2 parameters, set total of actual and expected equal. Therefore we
#lose 3 dfs to get 7dfs

pval=pchisq(chisq,df=7,lower.tail = F); pval

#We got a p-value of 1.651946e-37. It is thus extremely improbably that we observe a test 
#statistic this extreme from the relevant chi square distribution, so under the current 
#evidence we strongly reject the null hypothesis at the 0.05 level of significance that the 
#data follows a gamma distribution.


#PART 2:
#Exploring Correlations and Independences within our dataframes to Deathrate
#These could be use later on for modelling.


#I)
#Is there a significant difference in death rates between rich countries 
#and poor countries?
#We will check via a permutation test
#(REQ: A permutation test, Comparison of analysis by classical methods (chi-square, 
#CLT) and simulation methods, An example where permutation tests or other computational 
#techniques clearly work better than classical methods Nicely labeled graphics using 
#ggplot, with good use of color, line styles...)

#we will look only at the latest data from the latest date
data2$Date[which.max(data2$Date)]

rich=which(data2$Rich==1&data2$Date==data2$Date[which.max(data2$Date)])
poor=which(data2$Rich==0&data2$Date==data2$Date[which.max(data2$Date)])
deathRateDif=mean(data2$deathrate[rich])-mean(data2$deathrate[poor]); deathRateDif
#Unsuprisingly Higher, but not by much. Likely to be insignificant.

temp=data2[which(data2$Date==data2$Date[which.max(data2$Date)]),]

N=100000; diff=numeric(N)
for(i in 1:N){
  samp=sample(nrow(temp),sum(temp$Rich==0)) 
  drsamp=mean(temp$deathrate[samp])
  drother=mean(temp$deathrate[-samp])
  diff[i]=drsamp-drother
}

p=ggplot() + 
  aes(diff)+ 
  geom_histogram(binwidth=0.001,bins=100, colour=rgb(0,0.8,0.107,1), fill=rgb(0,0.62,0.107,1))+
  ylab("Count")+
  xlab("Differences in death rates between rich and poor countries")+
  ggtitle("Histogram of differences in death rates between rich and poor countries")

p+geom_vline(xintercept =deathRateDif,col="blue")

#Does not look significant, reaffirming what we thought before

pv.1t=(sum(diff>=deathRateDif)+1)/(N+1)
pv.2t=2*pv.1t;pv.2t

#Not significant at 0.05 level of significance somewhat suprisingly-
#There is a 58.88741% chance of discrepency by chance. 
#This is probably because poorer countries do not have the resources
#to test everyone and thus do not know exactly how some might have died.
#Moreover, it is possible that these countries are not tracking all deaths
#in them.

#How about using a t test?
#We will use a two-sample t-test to check whether there is evidence against the 
#null hypothesis that two population means are equal
t.test(data2$deathrate[rich],data2$deathrate[poor])

#same conclusion at the 0.05 level of significance, albeit weaker with
#p-value 0.04360523, with us rejecting the null hypothesis of no mean 
#difference

#However the weakness could be explained by whether the deathrate follows 
#a normal distribution. (Which it does not as shown in I)). We will plot
#again to confirm, and run a shapiro wilk test

p=ggplot(data2, aes(x=deathrate))  +
  geom_histogram(aes(y=..density..), binwidth=0.01, colour="steelblue", fill="white",bins="fd")
p+stat_function(fun=dnorm, args=list(mean=mean(data2$deathrate), sd=sd(data2$deathrate)),col="red")+xlim(0,1)
shapiro.test(data2$deathrate)

#we receive p-value 2.2e-16 in the shapiro wilk test for normality
#meaning that we reject the null hypothesis of no significant difference
#with the normal distribution, and conclude that under the current 
#evidence, death rate does not follow a normal distribution

#As such the t test is fundamentally flawed. And the permutation test's
#result of a strong rejection is given much more weight. This method is 
#more appropriate than the classical method.

#II)
#Are crime and death rates are independent?
#hdr is a logical variable which is 1 for countries in the upper half
#of death rates, and 0 for countries in the lower half.
#We will check via a chi-squared test on hdr and crime.
#(REQ: A contingency table, Analysis of a contingency table)

data2$hdr=0+(data2$deathrate>median(data2$deathrate))
lookingat=data2[which(data2$deaths2>0),]

#cro from expss generates a contingency table
cro(data2$hdr,data2$Crime)  

#Run Chi-square test
Observed=table(data2$hdr,data2$Crime)
Expected <- outer(rowSums(Observed), colSums(Observed))/sum(Observed)
chisq <- ChiSq(Observed,Expected); chisq
pval= 1 - pchisq(chisq,1); pval
#We strongly reject null hypothosis of independence at the 0.05 level of 
#significance, with a P-value so small that the computer rounded it to 0, 
#(much less than 0.05). This looks promising for modelling later on.


#III)
#Exploring Correlation with all our other data variables
#We will carry this out using a correlation heatmap
#(REQ: A graphical display different from those in the textbook 
#or in the class scripts, Appropriate use of correlation, Nicely labeled 
#graphics using ggplot, with good use of color, line styles...)

#extract numeric columns, ensure everything is numeric
temp <- subset(data2, select = -c(Country,countriesAndTerritories, geoId, countryterritoryCode, continentExp, Date, day, month, year, hdr, Entity, Code, dateRep))
sapply(temp, is.numeric)
temp$days=as.numeric(temp$days)

cormat <- round(cor(temp),2)


# Melt the correlation matrix
upper_tri <- get_upper_tri(cormat)
melted_cormat <- melt(upper_tri, na.rm = TRUE)
# Heatmap
ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 1, 
                                   size = 10, hjust = 1))+
  xlab("")+ylab("")+ggtitle("Correlation heatmap for COVID-19 data")
coord_fixed()

#Death rate actually isnt strongly correlated with much at all, 
#with the other variables. Interestingly, we dont see high degrees of 
#similarity. Perhaps this indicates that deathrate is relatively constant.


#PART 3:
#Modelling deathrate and deathtoll.

#I)
#Modelling deathrate via logistic regression.
#lets train some models on a random sample of 30% of the data
#and then test on the rest of the sample.
#then we can compare fitted R^2 to determine best model
#(REQ: Logistic Regression)

temp <- subset(data2, select = -c(Country,countriesAndTerritories, geoId, countryterritoryCode, continentExp, Date, day, month, year, hdr, Entity, Code, dateRep))

#would be unfair to feed death data here
temp2=subset(temp,select=-c(deaths,deaths2,active))

boxplot(temp2)
#there seems to be major spread of outliers and general spread in
#GDP.2018. We will log(GDP.2018) to improve this.
temp2$GDP.2018=I(log(temp2$GDP))
boxplot(temp2)

#there seems to be major spread of outliers and general spread in
#Population.2020, Total.tests. We will log these to improve this.
temp2$Population.2020=I(log(temp2$Population.2020))
temp2$Total.tests=I(log(temp2$Total.tests))
boxplot(temp2)

#there seems to be major spread of outliers and general spread in
#Cases2, recoveries. We cannot log these as they have 0 values.
#lets try squarerooting these

temp2$recoveries=I(sqrt(temp2$recoveries))
temp2$cases2=I(sqrt(temp2$cases2))
boxplot(temp2)

#there seems to be major spread of outliers and general spread in
#cases. We cannot log this as it has 0 values. lets try 
#squarerooting these
temp2$cases=I(sqrt(temp2$cases))
boxplot(temp2)


#80/20 split for training and testing data
sample <-sample.int(n =nrow(temp2),size =floor(.80*nrow(temp2)), replace = F)
train.data <- temp2[sample, ]
test.data <-temp2[-sample,]
train.data=drop_na(train.data)

fit=glm(deathrate~.,data=train.data,family="binomial")

#Check for multicollinearity
vif(fit)

#Deal with multicollinearity
length=1
while(length!=ncol(train.data)){
  length=ncol(train.data)
  train.data=dropvif1(fit,train.data)
  fit=glm(deathrate~.,data=train.data,family="binomial")
}

vif(fit)

#Male.lung seems quite high; although not 10 its very close to it. I will drop
#this as such
train.data=subset(train.data,select = -c(Male.Lung))
fit=glm(deathrate~.,data=train.data,family="binomial")

vif(fit)
#Seems good-we now have our first model
fit; summary(fit)
#R2

#model 2
fit4=glm(deathrate~.^2,data=train.data,family="binomial"); summary(fit4)


#model 5, LASSO Model
X = model.matrix(deathrate~.,train.data )[,-1]
Y=train.data$deathrate
cv=cv.glmnet(X,Y,alpha =1)
model=glmnet(X,Y,alpha =1, lambda=cv$lambda.min)

#model 6, Ridge Model
cv=cv.glmnet(X,Y,alpha =0)
model2=glmnet(X,Y,alpha =0, lambda=cv$lambda.min)

#model 7, neural net
#we need to normalize the data
train.data2=subset(train.data,select=-c(days)) #otherwise get errors when normalizing
ntrain.data=as.data.frame(lapply(train.data2,normalize))
ntrain.data=drop_na(ntrain.data)
nn=neuralnet(deathrate ~ ., data=ntrain.data, hidden=c(16,25), linear.output =T, threshold=0.1)

#model 8, the constant
fit8=mean(train.data$deathrate)

#Now Lets test accuracy!
#We will look at SSE
#we also need to remove the columns in testing data (for testing lasso and ridge)
#that we removed in training
test.data=test.data[which(colnames(test.data)%in%colnames(train.data))]

sse1=sum(test.data$deathrate-predict(fit,new=test.data))^2
sse2=sum(test.data$deathrate-predict(fit2,new=test.data))^2


X.test = model.matrix(deathrate~.,test.data )
fits = X.test%*%coef(model)
sse5=sum(test.data$deathrate-fits)^2

fits2 =X.test%*%coef(model2)
sse6=sum(test.data$deathrate-fits2)^2

nn.results <- compute(nn, test.data)
sse7=sum(test.data$deathrate-nn.results$net.result)^2

sse8=sum(test.data$deathrat-fit8)^2/nrow(test.data)


c(mean(sse1),mean(sse2),mean(sse5),mean(sse6),mean(sse7),mean(sse8))/(nrow(test.data))

#Suprsingly the constant model was the most accurate with sse 1.336694e-06! Deathrate 
#is best approximated by a constant


#II)
#Can we do a better job modelling the total deaths?
#(REQ: Use of linear regression)

temp2=subset(temp,select=-c(deaths2,deathrate,active))

#we can run same code to deal with spread of outliers
temp2$GDP.2018=I(log(temp2$GDP))
temp2$Population.2020=I(log(temp2$Population.2020))
temp2$Total.tests=I(log(temp2$Total.tests))
temp2$recoveries=I(sqrt(temp2$recoveries))
temp2$cases2=I(sqrt(temp2$cases2))
temp2$cases=I(sqrt(temp2$cases))

boxplot(temp2)

#80/20 split for training and testing data
sample <-sample.int(n =nrow(temp2),size =floor(.80*nrow(temp2)), replace = F)
train.data <- temp2[sample, ]
test.data <-temp2[-sample,]
train.data=drop_na(train.data)

fit=lm(deaths~.,data=train.data)
vif(fit)

#deal with multicollinearity
length=1
while(length!=ncol(train.data)){
  length=ncol(train.data)
  train.data=dropvif1(fit,train.data)
  fit=lm(deaths~.,data=train.data)
}

vif(fit)

#looking good, we have our first model
summary(fit)
plot(fit)
#diagnostic plots dont look promising, residuals vs fitted not totally
#random, standardized residuals dont seem to follow normal distribution well

#model 2
fit2=stepwise(fit,direction='forward/backward',criterion='BIC',trace='false'); summary(fit2)
plot(fit2)
#diagnostic plots dont look promising, residuals vs fitted not totally
#random, standardized residuals dont seem to follow normal distribution well

#model 3
fit3=stepwise(fit,direction='forward/backward',criterion='AIC',trace='false'); summary(fit3)
plot(fit3)
#diagnostic plots dont look promising, residuals vs fitted not totally
#random, standardized residuals dont seem to follow normal distribution well

#model 4
fit4=lm(deaths~.^2,data=train.data)
plot(fit4)
#diagnostic plots seem better, residuals vs fitted not totally
#random at lower fitted values, but well spaced as this gets higher,
#standardized residuals dont seem to follow normal distribution well

#model 5
fit5=stepwise(fit4,direction='forward/backward',criterion='AIC',trace='false'); summary(fit3)
plot(fit5)
#diagnostic plots dont look promising, residuals vs fitted not totally
#random, standardized residuals dont seem to follow normal distribution well

#model 6
fit6=stepwise(fit4,direction='forward/backward',criterion='AIC',trace='false'); summary(fit3)
plot(fit6)
#diagnostic plots seem better, residuals vs fitted not totally
#random at lower fitted values, but well spaced as this gets higher,
#standardized residuals dont seem to follow normal distribution well

#model 7, LASSO Model
X = model.matrix(deaths~.,train.data )[,-1]
Y=train.data$deaths
cv=cv.glmnet(X,Y,alpha =1)
model=glmnet(X,Y,alpha =1, lambda=cv$lambda.min)

#model 8, Ridge Model
cv=cv.glmnet(X,Y,alpha =0)
model2=glmnet(X,Y,alpha =0, lambda=cv$lambda.min)

#model 9, neural net
#we need to normalize the data
ntrain.data=as.data.frame(lapply(train.data,normalize))
nn=neuralnet(deaths ~ .-cases , data=ntrain.data, hidden=c(16,25), linear.output =T, threshold=0.1)
#dropping cases or else get:
#Error: the error derivative contains a NA; varify that the derivative function 
#does not divide by 0 (e.g. cross entropy)

#which shouldnt be possible as probability never 0 implying this column led
#R to rounding some really small number

#Checking our models
#we also need to remove the columns in testing data (for testing lasso and ridge)
#that we removed in training
test.data=test.data[which(colnames(test.data)%in%colnames(train.data))]



sse1=sum(test.data$deaths-predict(fit,new=test.data))^2
sse2=sum(test.data$deaths-predict(fit2,new=test.data))^2
sse3=sum(test.data$deaths-predict(fit3,new=test.data))^2
sse4=sum(test.data$deaths-predict(fit4,new=test.data))^2
sse5=sum(test.data$deaths-predict(fit5,new=test.data))^2
sse6=sum(test.data$deaths-predict(fit6,new=test.data))^2

X.test = model.matrix(deaths~.,test.data )
fits = X.test%*%coef(model)
sse7=sum(test.data$deaths-fits)^2

fits2 =X.test%*%coef(model2)
sse8=sum(test.data$deaths-fits2)^2

nn.results <- compute(nn, test.data)
sse9=sum(test.data$deaths-nn.results$net.result)^2

c(mean(sse1),mean(sse2),mean(sse3),mean(sse4),mean(sse5),mean(sse6),mean(sse7),mean(sse8),mean(sse9))/(nrow(test.data))


#Our fifth model seems best



