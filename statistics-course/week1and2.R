### This is the code for Week 1: Getting started

#### Video
dat <- read.csv("femaleMiceWeights.csv")
str(dat)
dat[12,2]
weight <- dat$Bodyweight
weight[11]
length(dat)
dat

#### Exercizes
mean(dat$Bodyweight[13:24])
?sample
sample(13:24, 1)
set.seed(1)
sample(13:24, 1)
dat[21,2]

#### Video
install.packages("dplyr")

dat <- read.csv("femaleMiceWeights.csv")
library(dplyr)

# filter out the mouses with chow diet, find their bodyweight an make a numerical vector
controls <- filter(dat, Diet=="chow")
controls <- select(controls, Bodyweight)
unlist(controls)

# same but faster and more clear
controls <- filter(dat, Diet=="chow") %>% select(Bodyweight) %>% unlist

### Exercises
install.packages("downloader")
library(downloader)
url="https://raw.githubusercontent.com/genomicsclass/dagdata/master/inst/extdata/msleep_ggplot2.csv"
filename <- basename(url)
download(url,filename)
dat <- read.csv("msleep_ggplot2.csv")
class(dat)
primates <- filter(dat, order == "Primates")
nrow(primates)
class(primates)
select(primates, "sleep_total") %>% class 
select(primates, "sleep_total") %>% unlist %>% mean
?summarize
filter(dat, order == "Primates") %>% summarize(mean(sleep_total))

### Excersize
load("skew.RData")
dim(dat)
par(mfrow = c(3,3))
for (i in 1:9){
  qqnorm(dat[,i])
}

### Excersize
head(InsectSprays)
boxplot(InsectSprays)
str(InsectSprays)
InsectSprays

?split
data <- split(InsectSprays, InsectSprays$spray)
select(data, "count")

select(data$F, "count") %>% unlist %>% mean

plot(nym.2002)
data(nym.2002, package="UsingR")
par(mfrow = c(1,4))
males <- filter(nym.2002, gender=="Male") %>% select(time) %>% unlist
females <- filter(nym.2002, gender=="Female") %>% select(time) %>% unlist
plot(females_time, xlab = "#", ylab = "time")
points(females_time, col = "green")
points(males_time, col = "red")
boxplot(females_time, males_time)
hist(females,xlim=c(range( nym.2002$time)))
hist(males,xlim=c(range( nym.2002$time)))
?hist





### This is the code for Week 2: Random variables and probability distribution

### Video
library(dplyr)  # required for %>%

dat <- read.csv("femaleMiceWeights.csv")  #load the data

control <- filter(dat, Diet=="chow") %>% select(Bodyweight) %>% unlist  # separte measurements in a list
treatment <- filter(dat, Diet=="hf") %>% select(Bodyweight) %>% unlist

mean( treatment)
mean( control)
obs <- mean( treatment) - mean( control)  # random variables --> is it chance that I obtained this value? 

population <-read.csv("femaleControlsPopulation.csv")
population <- unlist (population)

mean( sample( population, 12))

control <- sample( population, 12)
treatment <- sample( population, 12)
mean( treatment) - mean( control)

n <- 10000  # record all the differences that we see under this null hypothesis
nulls <- vector("numeric", n)
for(i in 1:n){
  control <- sample( population, 12)
  treatment <- sample( population, 12)
  nulls[i] <- mean( treatment) - mean( control)
}

max(nulls)
hist(nulls)  # picture of the null distribution under the null hypothesis

sum( nulls > obs)/n
mean( nulls > obs)  # this is the same as the line above
mean( abs( nulls) > obs)  # this is the p-value --> difference in means under null hypothesis

library(rafalib)
mypar()

qqnorm(nulls)
qqline(nulls)


# p-value: what is the probability that an outcome from the null distribution is bigger 
# than what we observed when the null hypothesis is true


### Exercizes
mean_pop <- mean(population)
set.seed(5)
mean_sample <- mean(sample(population, 5))
mean_dif <- abs(mean_pop-mean_sample)

# make averages5
set.seed(1)
n <- 1000
averages5 <- vector("numeric", n)
for(i in 1:n){
  sample5 <- sample( population, 5)
  averages5[i] <- mean( sample5)
}
hist( averages5)
mean_population <- mean( population)
(sum( averages5 < mean_population - 1) + sum( averages5 > mean_population + 1)) / 1000
mean( abs( averages5 - mean(population)) > 1)

# make averages50
set.seed(1)
n <- 1000
averages50 <- vector("numeric", n)
for(i in 1:n){
  sample50 <- sample( population, 50)
  averages50[i] <- mean( sample50)
}
hist( averages50)
mean( abs( averages50 - mean(population)) > 1)

install.packages("gapminder")
library(gapminder)
library(dplyr)
data(gapminder)
head(gapminder)
x <- filter(gapminder, year == 1952) %>% select(lifeExp) %>% unlist  # life expectancies of each country for 1952

# cumulative distribution fucntion: the function F(a) for any a, 
# which tells you the proportion of the values which are less than or equal to a

# mean(x <= a): calculates the number of values in x which are less than or equal to a,
# divided by the total number of values in x, or: the proportion of values <= a

mean(x <= 40)  # proportion of countries in 1952 with life expectancy <= 40
mean(x <= 60) - mean(x <= 40)

prop = function(q){  # where x is life expectancy of countries in 1952
  mean(x <= q)       # where q is life expectancy
}
prop(40)

qs = seq(from=min(x), to=max(x), length=20)
props = sapply(qs, prop)

props = sapply(qs, function(q) mean(x <= q))  # same as above, but with anonymous function

plot(props, qs) # self-made function for cummulative distrubution function
plot(ecdf(x))   # build-in function for cummulative distribution function

### Video
# Normal distribution --> summary is mean & standard deviation --> proportion of any interval
# mu = average/mean
# sigma = standard deviation
# sigma^2 == variance == std^2  --> distance between each individual and the mean
# quantile-quantile plot --> check the distrubution of the data 
# plot percintiles for data v.s. normal distrubution
# standard units (Z) --> mean 0 and std 1, differenc from the mean, no more units

### Exercize
par(mfrow = c(1,2))
hist(averages5)
hist(averages50)
mean(averages50 <= 25) - mean(averages50 <= 23)

#to find the proportion of observations below a cutoff x: pnorm(x, mu, sigma)
pnorm(25, 23.9, 0.43) - pnorm(23, 23.9, 0.43)

### Video
# random variable --> if we take another sample, these values will change
# How much do they change? 
# How far away are they from the population averages?

### Exercize
dat <- read.csv("mice_pheno.csv")
dat <- na.omit( dat)  # remove lines that contain missing values
library(rafalib)

# difference in bodyweight depending on diet for males
x <- filter(dat, Diet == "chow") %>% filter(Sex == "M") %>% select(Bodyweight) %>% unlist
mean(x)
popsd(x)  # determine the standard deviation of the population

set.seed(1)
X <- sample(x, 25)
mean(X)

y <- filter(dat, Diet == "hf") %>% filter(Sex == "M") %>% select(Bodyweight) %>% unlist
mean(y)
popsd(y)

set.seed(1)
Y <- sample(y, 25)
mean(Y)

abs_dif_male <- abs((mean(y) - mean(x)) - (mean(Y) - mean(X)))


# difference in bodyweight depending on diet for females
x <- filter(dat, Diet == "chow") %>% filter(Sex == "F") %>% select(Bodyweight) %>% unlist
mean(x)
popsd(x)  # determine the standard deviation of the population

set.seed(2)
X <- sample(x, 25)
mean(X)

y <- filter(dat, Diet == "hf") %>% filter(Sex == "F") %>% select(Bodyweight) %>% unlist
mean(y)
popsd(y)

set.seed(2)
Y <- sample(y, 25)
mean(Y)

abs_dif_female <- abs((mean(y) - mean(x)) - (mean(Y) - mean(X)))
abs_dif_female


### Video
# Centrolal Limit Theorem: sample average follows a normal distribution
# sigmaX: population standard deviation --> average distance to population mean
# std (of )sample average) = std(population)/sample sizer --> larger sample size = smaller spread

### Exercizes
pnorm(1)-pnorm(-1)  # proportion of numbers with 1 std away from the list's average
pnorm(2)-pnorm(-2)  # proportion of numbers with 2 std away from the list's average
pnorm(3)-pnorm(-3)  # proportion of numbers with 3 std away from the list's average

y <- filter(dat, Diet == "chow" & Sex == "M") %>% select(Bodyweight) %>% unlist
mean <- mean(y)
sd = popsd(y)
mean = mean(y)
hist(y)

prop = function(q){  # where x is life expectancy of countries in 1952
  mean(y <= q)       # where q is life expectancy
}

# proportion of the mice within 1 std away from the average weight
proportion <- 1 - prop(mean-sd) - (1 - prop(mean+sd))

z <- (y - mean(y))/popsd(y)  # translate to standard units (Z)
mean( abs(z) <= 1)           # same answer!
mean( abs(z) <= 2)
mean( abs(z) <= 3)

mypar(2,2)
y <- filter(dat, Sex=="M" & Diet=="chow") %>% select(Bodyweight) %>% unlist
z <- ( y - mean(y) ) / popsd(y)
qqnorm(z);abline(0,1)
y <- filter(dat, Sex=="F" & Diet=="chow") %>% select(Bodyweight) %>% unlist
z <- ( y - mean(y) ) / popsd(y)
qqnorm(z);abline(0,1)
y <- filter(dat, Sex=="M" & Diet=="hf") %>% select(Bodyweight) %>% unlist
z <- ( y - mean(y) ) / popsd(y)
qqnorm(z);abline(0,1)
y <- filter(dat, Sex=="F" & Diet=="hf") %>% select(Bodyweight) %>% unlist
z <- ( y - mean(y) ) / popsd(y)
qqnorm(z);abline(0,1)

y <- filter(dat, Sex=="M" & Diet=="chow") %>% select(Bodyweight) %>% unlist
set.seed(1)
avgs <- replicate(10000, mean( sample(y, 25)))
mypar(1,2)
hist(avgs)
qqnorm(avgs)
qqline(avgs)
mean(avgs)
popsd(avgs)


### Video
# How do we use the Central Limit Theorem to obtain p-values to obtain confidence intervals
# and to perform statistical interference in general? 

library(dplyr)
# This is the sample data 
dat <- read.csv("femaleMiceWeights.csv")
control <- filter(dat, Diet == "chow") %>% select(Bodyweight) %>% unlist
treatment <- filter(dat, Diet == "hf") %>% select(Bodyweight) %>% unlist

N <- length(treatment)
obs <- mean(treatment) - mean(control)
se <- sqrt(var(treatment))/N + var(control)/N
tstat <- obs / se

# pnorm: what proportuion of normally distrubuted data with means 0 and std 1 are lover than 
# whatever value you put here --> only used sample data
2*(1 - pnomr(tstat)) # p-value --> two tails


# This is the data of the entire population
population <-read.csv("femaleControlsPopulation.csv")
population <- unlist (population)

n <- 10000  # record all the differences that we see under this null hypothesis
nulls <- vector("numeric", n)
for(i in 1:n){
  control <- sample( population, N)
  treatment <- sample( population, N)
  se <- sqrt(var(treatment))/N + var(control)/N  # standard error estimates
  nulls[i] <- (mean( treatment) - mean( control))/se
}

library(rafalib)
mypar()
qqnorm(nulls)
abline(0,1)


# Without population data

library(dplyr)
# This is the sample data 
dat <- read.csv("femaleMiceWeights.csv")
control <- filter(dat, Diet == "chow") %>% select(Bodyweight) %>% unlist
treatment <- filter(dat, Diet == "hf") %>% select(Bodyweight) %>% unlist

t.test(treatment, control)  # this will compute mean and std, assuming t-approximation

# check if it is somewhat similar to normal distribution
qqnorm(control)
qqline(control)

qqnorm(treatment)
qqline(treatment)


### Exercize
n <- 100
x <- sample(1:6, n, replace = TRUE)
mean(x == 6)

mypar(2,2)
set.seed(1)
n <- 100
sides <- 6
p <- 1/sides

zs <- replicate(10000, {
  x <- sample(1:sides, n, replace = TRUE)
  (mean(x == 6) - p) / sqrt(p*(1-p)/n)
})

qqnorm(zs)
abline(0,1)
mean(abs(zs) > 2)


ps <- c(0.5,0.5,0.01,0.01)
ns <- c(5,30,30,100)
library(rafalib)
mypar(4,2)
for(i in 1:4){
  p <- ps[i]
  sides <- 1/p
  n <- ns[i]
  zs <- replicate(10000,{
    x <- sample(1:sides,n,replace=TRUE)
    (mean(x==1) - p) / sqrt(p*(1-p)/n)
  }) 
  hist(zs,nclass=7)
  qqnorm(zs)
  abline(0,1)
}

# A major difference with binary data, for which we know the variance is  p(1−p) , 
# is that with quantitative data we need to estimate the population standard deviation

dat <- read.csv("femaleMiceWeights.csv")
X <- filter(dat, Diet == "chow") %>% select(Bodyweight) %>% unlist
Y <- filter(dat, Diet == "hf")   %>% select(Bodyweight) %>% unlist
mean(X)
sd(X)
# How well does mean(X) (average of the sample) approximates mu_X (average of the population)?
# mean(X) follows a normal distribution with mean mu_X and std std_X/sqrt(sample size)

# The probability that our estimate mean(X) is off by more than 2 grams from mu_X ?!
2 * ( 1-pnorm( 2/sd(X) * sqrt(12) ) )

SE <- sqrt(var(Y)/12 + var(X)/12) # standard error of (mean(X) - mean(Y))
t <- ( mean(Y) - mean(X) )/ ( sqrt(var(Y)/12 + var(X)/12) )
t <- t.test(Y,X)$stat  # t-statistic
t.test(X,Y)$p.value  # p-value with t-statistics

Z <- ( mean(Y) - mean(X) ) / sqrt( var(X)/12 + var(Y)/12)
2*( 1-pnorm(Z))  # p-value with CLT

