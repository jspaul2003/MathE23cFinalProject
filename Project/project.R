#FINAL PROJECT MATH E-23C
#JEAN-SEBASTIEN PAUL
#EXPLORING COVID-19 DEATHRATES



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

#data2Setup
#Sets up data2 once predictors, testing and geographicdata
#combined. Also for use when reading data in as a .csv
data2Setup = function(data){
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

#ChiSq
#Returns chi squared value
ChiSq <-function(Obs,Exp){
  sum((Obs-Exp)^2/Exp)
}



#DATA SCRAPING
#NOTE: A CSV file has also been provided. This is here if you
#are at all interested at how I collected and formatted that data

#We will have 2 main dataframes:
#geographicdata-in depth data in recoveries, cases and deaths
#data2-data with predictors and testing data, but less in depth
#in recoveries, cases and deaths, due to lack of data in some 
#of those predictors and rows with NA values being dropped


#Get data from internet
JHURecovs=read.csv("https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_recovered_global.csv")
predictors=read.csv("https://raw.githubusercontent.com/jspaul2003/nCoV2019/master/more_data/covid19_by_country.csv")
testing=read.csv("https://raw.githubusercontent.com/jspaul2003/nCoV2019/master/testing.csv")
GET("https://opendata.ecdc.europa.eu/covid19/casedistribution/csv", authenticate(":", ":", type="ntlm"), write_disk(tf <- tempfile(fileext = ".csv")))
geographicdata=read.csv(tf)

#Formatting this data
JHURecovs1=formatJHU(JHURecovs,"recoveries")
predictors$countriesAndTerritories=predictors$Country
testing$countriesAndTerritories=testing$Entity
testing$Date=as.Date(testing$Date, format= "%b %d, %y")
geographicdata=dateCol(geographicdata)
geographicdata=geographicdata%>%right_join(JHURecovs1, by=c("countriesAndTerritories","Date"))
geographicdata=drop_na(geographicdata)

#Making Our Big Datasets
#(REQ: A dataframe, At least 2 numeric columns, A data set with lots of columns, 
#allowing comparison of many different variables.)
geographicdata$days=geographicdata$Date-geographicdata$Date[which.min(geographicdata$Date)]

#determine cumulative deaths and cases for geographic data-> use this
#to find active cases, deathrate
geographicdata=dccum(geographicdata)
geographicdata$active=geographicdata$cases2-geographicdata$deaths2-geographicdata$recoveries
geographicdata$deathrate=geographicdata$deaths2/(geographicdata$cases2+1*(geographicdata$cases2==0))


data2=geographicdata %>% inner_join(testing, by=c("countriesAndTerritories","Date"))
data2=predictors %>% right_join(data2, by=c("countriesAndTerritories"))
data2 <- subset(data2, select = -c(Tests, Test.Pop,Quarantine,Schools,Restrictions,Total.Recovered,Total.Deaths,Total.Infected,popData2018))
data2=data2Setup(data2)

#READ FROM CSV ALTERNATIVE (TO DATASCRAPING)
#Note that REQs of dataframes are all specified above in the datascraping 
#code comments

#data=read.csv("data.csv")
#data=dateCol(data)
#geographicdata=data[,which(colnames(data)=="countriesAndTerritories"):ncol(data)]
#data2=data2Setup(data)
#(REQ: A dataframe, At least 2 numeric columns, A data set with lots of columns, 
#allowing comparison of many different variables.)

#(REQ: At least 20 rows)
nrow(geographicdata)
nrow(data2)


#ANALYSIS

#PART 1:
#Exploration regarding Distribution of Variables related to deathrate from COVID-19
#using probability distributions.

#I)
#Can we model the time series data of new deaths world wide as binomial distribution?
#We will look at the interval between days 1 to 122. We will use the more in depth 
#in deathratel, geographicdata dataframe. It will contain 67 countries always having data 
#between 1 to 122.
#(REQ: A barplot)

world=geographicdata[c(seq(1:max(geographicdata$days)),max(geographicdata$days)),]
world$deathrate=0
a=geographicdata[which(geographicdata$days==1),]
a=geographicdata[which(geographicdata$days==max(geographicdata$days)&geographicdata$countriesAndTerritories%in%a$countriesAndTerritories),]
for(i in 1:(max(data2$days))){
  world$Date[i]=i+min(geographicdata2$Date)
  world$days[i]=i
  
  world$deathrate[i]=sum(geographicdata$deaths2[which(geographicdata$days==i)])/(sum(geographicdata$cases2[which(geographicdata$days==i)])+1*(sum(geographicdata$cases2[which(geographicdata$days==i)])==0))
  world$deaths2[i]=sum(geographicdata$deaths2[which(geographicdata$days==i&(geographicdata$countriesAndTerritories%in%a$countriesAndTerritories))])
  world$cases2[i]=sum(geographicdata$cases2[which(geographicdata$days==i&(geographicdata$countriesAndTerritories%in%a$countriesAndTerritories))])
  world$deaths[i]=sum(geographicdata$deaths[which(geographicdata$days==i&(geographicdata$countriesAndTerritories%in%a$countriesAndTerritories))])
  world$cases[i]=sum(geographicdata$cases[which(geographicdata$days==i&(geographicdata$countriesAndTerritories%in%a$countriesAndTerritories))])
  world$active[i]=sum(geographicdata$active[which(geographicdata$days==i&(geographicdata$countriesAndTerritories%in%a$countriesAndTerritories))])
}

#Could we plot a binomial distribution for deaths
world$dr2=world$deaths/world$active
barplot(rbind((world$deaths2),world$active/6*pbinom(1:nrow(world),(round(world$active/1000)),mean(world$deathrate))), beside = TRUE, col = c("red", "blue")) #no

#Seems like a very weak model.


#II)
#How are death rates distributed: do they converge at a certain value or do they vary wildly?
#(REQ: A histogram, A probability density function overlaid on a histogram, A p-value or other
#statistic based on a distribution function, Appropriate use of R functions for a probability 
#distribution other than binomial, normal, or chi-square.)

hist(geographicdata$deathrate[which(geographicdata$deaths2!=0)],prob=T,breaks="fd")

#Lets try modelling with a gamma function
curve( dgamma(x,0.45,13)    ,add=T,col="red")
#Seems relatively well approximated

#Lets test if we can indeed model with a gamma function

#create bins by breaking data into deciles
bins=qgamma(0.1*(1:10),0.45,13)

binstuff=cut(geographicdata$deathrate[which(geographicdata$deaths2!=0)], breaks=bins,labels=F); binstuff
obs=as.vector(table(binstuff))

exp=rep(sum(obs)/10,10)
chisq=sum((obs-exp)^2/exp); chisq

#p-value
#10 Categories, imposed 2 parameters, set total of actual and expected equal. Therefore we
#lose 3 dfs to get 7dfs

pval=pchisq(chisq,df=7,lower.tail = F); pval





#PART 2:
#Exploring Correlations and Independences within our dataframes to Deathrate
#These could be use later on for modelling.


#I)
#Is there a significant difference in death rates between rich countries 
#and poor countries?
#We will check via a permutation test
#(REQ: A permutation test)

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
hist(diff,breaks=100, col=rgb(0,0.62,0.107,1))
abline(v=deathRateDif,col="blue") 

#Does not look significant, reaffirming what we thought before

pv.1t=(sum(diff>=deathRateDif)+1)/(N+1)
pv.2t=2*pv.1t;pv.2t

#Not significant at 0.05 level of significance somewhat suprisingly-
#There is a 58.88741% chance of discrepency by chance. 
#This is probably because poorer countries do not have the resources
#to test everyone and thus do not know exactly how some might have died.
#Moreover, it is possible that these countries are not tracking all deaths
#in them.


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
chisq <- sum((Observed - Expected)^2/Expected); chisq
pval= 1 - pchisq(chisq,1); pval
#We strongly reject null hypothosis of independence at the 0.05 level of significance, with a 
#P-value of 3.116538e-09, much less than 0.05.
#This looks promising for modelling later on.

























