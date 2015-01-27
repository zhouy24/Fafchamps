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
##  [1] 515.45302 183.95931 628.59850 887.53845 142.10945 863.56907 162.78695
##  [8] 496.01463 410.87574  82.25616 141.06774 230.08028 859.54948 371.02938
## [15] 765.08084  46.80418 495.88993 656.68190 889.34892 390.33008 578.71230
## [22] 290.03223 187.06135 716.07174 850.09596 585.39172 634.53958 679.93663
## [29] 160.42309 747.11880 481.95525 191.24637 145.93851 120.02227 395.47743
## [36] 180.71949 602.11474 245.70372 267.10239 419.52973  79.23692 953.82908
## [43] 968.07923 300.98878 902.62343 808.95057 613.90774 133.89977 125.66490
## [50] 638.84348
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
##              9.740               9.993              10.534
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
## -2.33483 -0.59191 -0.02788  0.65989  2.23012 
## 
## Coefficients:
##                     Estimate Std. Error t value Pr(>|t|)    
## (Intercept)         9.740073   0.358278   27.19   <2e-16 ***
## age[somesample]     9.992654   0.005496 1818.23   <2e-16 ***
## gender[somesample] 10.533908   0.313569   33.59   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 1.086 on 47 degrees of freedom
## Multiple R-squared:      1,	Adjusted R-squared:      1 
## F-statistic: 1.67e+06 on 2 and 47 DF,  p-value: < 2.2e-16
```
Check if the confidence interval contains the true values of beta1 and beta2.

```r
confint(reg, "(Intercept)")
```

```
##               2.5 %   97.5 %
## (Intercept) 9.01931 10.46084
```

```r
confint(reg, "age[somesample]")
```

```
##                    2.5 %   97.5 %
## age[somesample] 9.981598 10.00371
```

```r
confint(reg, "gender[somesample]")
```

```
##                       2.5 %   97.5 %
## gender[somesample] 9.903088 11.16473
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
## [1] 99
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
