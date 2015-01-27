---
title: "Tyler1"
output: html_document
---

Simmulate a social network 
1a for i in 1,…,n respondents simulate a position in a 2-dimensional space 
using a draw from a 2d normal distribution

```r
rm(list = ls())
setwd('/Users/yuezhou/Desktop/TYLER')
n <- 100
x <- rnorm(n, mean = 0, sd = 1)
y <- rnorm(n, mean= 0, sd = 1)
plot(x, y)
```

![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-1.png) 
#1b compute the n*n matrix of Euclidean distances between individuals 

```r
A <- cbind(x, y)
distance <- dist(A)
distance <- as.matrix(distance)
```
For each pair of individuals (say individual i and individual j), 
simulate an edge using a draw from a bernoulli distribution with a success probability 
equal to exp(a*d_ij)/(1+exp(a*d_ij)) where d_ij is the distance between i and j computed in 1b 
and a is a constant that controls the overall expected tie frequency. 

```r
a = 0.1
probability <- exp(a * distance)/ (1+exp(a*distance))
rvariable <- rbinom(n*n, 1, as.vector(probability))
rvariable = matrix(rvariable, ncol = n, nrow = n)
plot(x, y)
for(i in 1:n){
  for(j in 1:n){
    if (rvariable[i, j] == 1){
    segments(x[i],y[i],x[j],y[j])
    }
  }
}
```

![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3-1.png) 
Simulate covariates for each individual in the network 
2a simulate each individual’s age and gender 
by taking a draw from a distribution that looks like the us population age distribution 
and flipping a fair coin, respectively. 

```r
age <- runif(n, min = 0, max =100)
gender <- rbinom(n, 1, 1/2)
```
2b Now simulate an outcome for each person. 
Take a draw from a standard normal with mean equal to beta0+beta1*gender+beta2*age.

```r
beta0 <- 10
beta1 <- 10
beta2 <- 10
dontknow <- rnorm(n, mean = beta0 + beta1*gender + beta2*age, sd = 1)
```
Take a sample and compute a confidence interval
3a Randomly choosing a set of individuals 
then running a regression model to estimate beta1 and beta2. 

```r
index <-1:length(dontknow)
m <- 50 
somesample <- sample(index, m, replace=FALSE, prob = NULL)
dontknow[somesample]
```

```
##  [1]  770.19296  852.83328  404.10088   24.56174  373.95926  724.55825
##  [7]  115.85478  389.46874  216.34460  282.52243  619.69382  284.63713
## [13]  439.79948  687.55733  795.42568   24.51116  965.28599  774.51612
## [19]  456.70026  767.68715   37.51872  649.45083   12.76021   89.88596
## [25]  297.99231  238.19120  896.30894   15.54972  957.84196  276.42997
## [31]   76.93026  553.05644  602.56881  281.85349  781.08429  652.24837
## [37]  433.33428  359.74406  602.56678  230.18235 1020.24341  144.12991
## [43]  192.40994  797.91889  532.09980  935.48948  194.33613  346.86761
## [49]   49.66997 1017.67991
```

```r
reg <- lm(dontknow[somesample] ~ age[somesample] + gender[somesample])
reg
```

```
## 
## Call:
## lm(formula = dontknow[somesample] ~ age[somesample] + gender[somesample])
## 
## Coefficients:
##        (Intercept)     age[somesample]  gender[somesample]  
##              9.538              10.007               9.915
```

```r
summary(reg)
```

```
## 
## Call:
## lm(formula = dontknow[somesample] ~ age[somesample] + gender[somesample])
## 
## Residuals:
##      Min       1Q   Median       3Q      Max 
## -1.59298 -0.59554 -0.03784  0.59636  2.49127 
## 
## Coefficients:
##                     Estimate Std. Error t value Pr(>|t|)    
## (Intercept)         9.537996   0.219432   43.47   <2e-16 ***
## age[somesample]    10.007306   0.003762 2660.20   <2e-16 ***
## gender[somesample]  9.915437   0.233462   42.47   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 0.8106 on 47 degrees of freedom
## Multiple R-squared:      1,	Adjusted R-squared:      1 
## F-statistic: 3.588e+06 on 2 and 47 DF,  p-value: < 2.2e-16
```
Check if the confidence interval contains the true values of beta1 and beta2.

```r
confint(reg, "(Intercept)")
```

```
##                2.5 %   97.5 %
## (Intercept) 9.096555 9.979436
```

```r
confint(reg, "age[somesample]")
```

```
##                    2.5 %   97.5 %
## age[somesample] 9.999738 10.01487
```

```r
confint(reg, "gender[somesample]")
```

```
##                       2.5 %  97.5 %
## gender[somesample] 9.445771 10.3851
```
Yes the confidence interval contains the true value of beta1 and beta 2

3b Now randomly sample *edges* from the social network
and evaluate in the same way as 3b.  
After doing this multiple times compute the coverage of the confidence intervals

```r
trials <- 100
countbeta1 <- 0
countbeta2 <- 0
index <-1:length(dontknow)
m <- 50 
for(i in 1:trials){
  somesample <- sample(index, m, replace=FALSE, prob = NULL)
  reg <- lm(dontknow[somesample] ~ age[somesample] + gender[somesample])
  summary(reg)
  interval1 = confint(reg, "age[somesample]")
  interval2 = confint(reg, "gender[somesample]")
  if (beta1 < interval1[2] && beta1 > interval1[1]) {
    countbeta1 <- countbeta1 + 1
  }
  if (beta2 < interval2[2] && beta2 > interval2[1]) {
    countbeta2 <- countbeta2 + 1
  }
}
countbeta1
```

```
## [1] 100
```

```r
countbeta2
```

```
## [1] 98
```
After computing a hundred times, 
the fraction of times the true parameter beta1
is contained in the estimated interval is around 99/100
and the fraction of times the true parameter beta2
is contained in the estimated interval is around 60/100
```
