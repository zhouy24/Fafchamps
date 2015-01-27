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
##  [1] 917.79386 230.19167 436.54081 248.44031 160.51412 615.84445 664.05953
##  [8] 850.02571 442.79959 213.23401 827.72347 435.07857 371.95230  80.08189
## [15] 648.21904 294.45161 406.84769 292.48212 993.65162 100.60596  99.16315
## [22] 678.24064  27.58611 705.79853 829.99811  97.49354 475.45019 181.10204
## [29] 262.26738 724.59939 379.26265 862.96437 275.60205  55.93078 324.07150
## [36] 533.96649 538.74854  97.50017 182.08662 498.81579 471.26478 823.97078
## [43] 949.63593 473.92164 506.66931 313.05863 926.54682 467.86445  91.97510
## [50] 606.22260
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
##              10.02               10.00               10.00
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
##     Min      1Q  Median      3Q     Max 
## -1.9822 -0.8109 -0.0143  0.5387  2.0355 
## 
## Coefficients:
##                     Estimate Std. Error t value Pr(>|t|)    
## (Intercept)        10.019469   0.312489   32.06   <2e-16 ***
## age[somesample]    10.000844   0.005169 1934.67   <2e-16 ***
## gender[somesample] 10.004978   0.283370   35.31   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 0.9977 on 47 degrees of freedom
## Multiple R-squared:      1,	Adjusted R-squared:      1 
## F-statistic: 1.873e+06 on 2 and 47 DF,  p-value: < 2.2e-16
```
Check if the confidence interval contains the true values of beta1 and beta2.

```r
confint(reg, "(Intercept)")
```

```
##                2.5 %   97.5 %
## (Intercept) 9.390822 10.64812
```

```r
confint(reg, "age[somesample]")
```

```
##                    2.5 %   97.5 %
## age[somesample] 9.990445 10.01124
```

```r
confint(reg, "gender[somesample]")
```

```
##                      2.5 %   97.5 %
## gender[somesample] 9.43491 10.57505
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
## [1] 94
```

```r
countbeta2
```

```
## [1] 100
```
After computing a hundred times, 
the fraction of times the true parameter beta1
is contained in the estimated interval is around 99/100
and the fraction of times the true parameter beta2
is contained in the estimated interval is around 60/100
```
