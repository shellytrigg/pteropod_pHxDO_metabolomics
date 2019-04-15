2016 pteropod survival proportions analysis
================
Shelly Trigg
4/15/2019

Load libraries

    ## ## FSA v0.8.22. See citation('FSA') if used in publication.
    ## ## Run fishR() for related website and fishR('IFAR') for related book.

    ## Loading required package: ggpubr

    ## Loading required package: magrittr

    ## 
    ## Attaching package: 'magrittr'

    ## The following object is masked from 'package:tidyr':
    ## 
    ##     extract

    ## Loading required package: Matrix

    ## 
    ## Attaching package: 'Matrix'

    ## The following object is masked from 'package:tidyr':
    ## 
    ##     expand

    ## 
    ## Attaching package: 'lmerTest'

    ## The following object is masked from 'package:lme4':
    ## 
    ##     lmer

    ## The following object is masked from 'package:stats':
    ## 
    ##     step

Read in data

    ## New names:
    ## * `` -> `..1`
    ## * `` -> `..2`
    ## * `` -> `..15`

Format survival data

remove samples with 2L, M, or LL (2L = 2 animals in one jar (misplacement), M = missing, LL = live but lost..not sure what that actually means)

Add Duration column

Add survival column

merge treatment data

**survival analysis** ![](2016_pteropod_survival_proportions_analysis_files/figure-markdown_github/unnamed-chunk-5-1.png)

**treatment effect on duration to death**

``` r
durTreat <- aov(duration ~ Treatment_abbv,data = dSub[which(dSub$isDead == 1),])
summary(durTreat)
```

    ##                 Df Sum Sq Mean Sq F value Pr(>F)
    ## Treatment_abbv   3   17.2   5.736   1.286  0.281
    ## Residuals      144  642.1   4.459

**MOATS effect on duration to death**

``` r
durMOATS <- aov(duration ~ MOATS,data = dSub[which(dSub$isDead == 1),])
summary(durMOATS)
```

    ##              Df Sum Sq Mean Sq F value Pr(>F)
    ## MOATS        11   72.3   6.571   1.522   0.13
    ## Residuals   136  587.0   4.317

**Test significance of treatment on duration to death with MOATS as a random effect**

``` r
fitME <- lmer(duration~ Treatment_abbv + (1|MOATS), data = dSub[which(dSub$isDead == 1),])
summary(fitME)
```

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: duration ~ Treatment_abbv + (1 | MOATS)
    ##    Data: dSub[which(dSub$isDead == 1), ]
    ## 
    ## REML criterion at convergence: 637.6
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -2.6886 -0.7520  0.2094  0.6895  1.1938 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  MOATS    (Intercept) 0.1924   0.4387  
    ##  Residual             4.3339   2.0818  
    ## Number of obs: 148, groups:  MOATS, 12
    ## 
    ## Fixed effects:
    ##                  Estimate Std. Error      df t value Pr(>|t|)    
    ## (Intercept)        7.4386     0.4145  5.9280  17.946 2.16e-06 ***
    ## Treatment_abbvHL  -0.4796     0.6158  7.0971  -0.779    0.461    
    ## Treatment_abbvLH  -0.8514     0.6211  7.5902  -1.371    0.210    
    ## Treatment_abbvLL  -0.7196     0.5734  5.4898  -1.255    0.260    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) Trt_HL Trt_LH
    ## Trtmnt_bbHL -0.673              
    ## Trtmnt_bbLH -0.667  0.449       
    ## Trtmnt_bbLL -0.723  0.487  0.482

**plot survival proportions for treatments** ![](2016_pteropod_survival_proportions_analysis_files/figure-markdown_github/unnamed-chunk-9-1.png)

**plot survival proportions (fraction) for treatments** ![](2016_pteropod_survival_proportions_analysis_files/figure-markdown_github/unnamed-chunk-10-1.png)

**plot survival proportions for MOATS** ![](2016_pteropod_survival_proportions_analysis_files/figure-markdown_github/unnamed-chunk-11-1.png)

**calculate chi square and p values for proportions**

``` r
kable(surv_stats, caption = "Treatment effect on survival proportions test table")
```

| Comparison | ChiSq            | P.value             |  P.value\_bonferroni|
|:-----------|:-----------------|:--------------------|--------------------:|
| HH\_HL     | 1.52051612180922 | 0.217541405340122   |                    1|
| HH\_LH     | 3.35220932658209 | 0.0671148785326453  |                    1|
| HH\_LL     | 1.24945485257985 | 0.263656623025348   |                    1|
| HL\_LH     | 0.33296118511538 | 0.563920615529862   |                    1|
| HL\_LL     | 5.36592618879012 | 0.0205338811344844  |                    1|
| LH\_LL     | 8.50433362477545 | 0.00354301640506706 |                    1|

EXCLUDING MOATS 2
-----------------

**survival analysis** ![](2016_pteropod_survival_proportions_analysis_files/figure-markdown_github/surv_analysis_noM2-1.png)

**treatment effect on duration to death (without M2)**

``` r
durTreat <- aov(duration ~ Treatment_abbv,data = dSub[which(dSub$isDead == 1),])
summary(durTreat)
```

    ##                 Df Sum Sq Mean Sq F value Pr(>F)
    ## Treatment_abbv   3   24.3   8.099   1.985  0.119
    ## Residuals      130  530.5   4.081

**MOATS effect on duration to death (without M2)**

``` r
durMOATS <- aov(duration ~ MOATS,data = dSub[which(dSub$isDead == 1),])
summary(durMOATS)
```

    ##              Df Sum Sq Mean Sq F value Pr(>F)
    ## MOATS        10   56.7    5.67     1.4  0.188
    ## Residuals   123  498.1    4.05

**Test significance of treatment on duration to death with MOATS as a random effect (without M2)**

``` r
fitME <- lmer(duration~ Treatment_abbv + (1|MOATS), data = dSub[which(dSub$isDead == 1),])
summary(fitME)
```

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: duration ~ Treatment_abbv + (1 | MOATS)
    ##    Data: dSub[which(dSub$isDead == 1), ]
    ## 
    ## REML criterion at convergence: 565.5
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -2.8217 -0.7412  0.2147  0.7093  1.2077 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  MOATS    (Intercept) 0.03648  0.191   
    ##  Residual             4.05705  2.014   
    ## Number of obs: 134, groups:  MOATS, 11
    ## 
    ## Fixed effects:
    ##                  Estimate Std. Error      df t value Pr(>|t|)    
    ## (Intercept)        7.4562     0.3342  5.2626  22.313 2.05e-06 ***
    ## Treatment_abbvHL   0.2004     0.6075  9.3322   0.330    0.749    
    ## Treatment_abbvLH  -0.8697     0.5136  7.5809  -1.693    0.131    
    ## Treatment_abbvLL  -0.7339     0.4586  4.7853  -1.600    0.173    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) Trt_HL Trt_LH
    ## Trtmnt_bbHL -0.550              
    ## Trtmnt_bbLH -0.651  0.358       
    ## Trtmnt_bbLL -0.729  0.401  0.474

**plot survival proportions for treatments (without M2)** ![](2016_pteropod_survival_proportions_analysis_files/figure-markdown_github/unnamed-chunk-17-1.png)

**plot survival proportions (as fractions) for treatments (without M2)** ![](2016_pteropod_survival_proportions_analysis_files/figure-markdown_github/unnamed-chunk-18-1.png)

**plot survival proportions for MOATS (without M2)** ![](2016_pteropod_survival_proportions_analysis_files/figure-markdown_github/unnamed-chunk-19-1.png)

**calculate chi square and p values for proportions (without M2)**

``` r
kable(surv_stats, caption = "Treatment effect on survival proportions test table (without M2)")
```

| Comparison | ChiSq              | P.value             |  P.value\_bonferroni|
|:-----------|:-------------------|:--------------------|--------------------:|
| HH\_HL     | 2.91246517501663   | 0.0878974270396332  |                    1|
| HH\_LH     | 3.35220932658209   | 0.0671148785326453  |                    1|
| HH\_LL     | 1.24945485257985   | 0.263656623025348   |                    1|
| HL\_LH     | 0.0164108278078504 | 0.898066074487066   |                    1|
| HL\_LL     | 7.03450751815149   | 0.00799538572174349 |                    1|
| LH\_LL     | 8.50433362477545   | 0.00354301640506706 |                    1|
