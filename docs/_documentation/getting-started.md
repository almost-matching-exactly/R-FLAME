---
layout: default
title: Getting Started
nav_order: 2
description: ""
permalink: /getting-started
---

# Getting Started
{: .no_toc }

<details open markdown="block">
  <summary>
    Table of contents
  </summary>
  {: .text-delta }
1. TOC
{:toc}
</details>

## Dependencies
This package requires prior installation of
- R (>= 3.3)
- Python (>= 2.7)

## Installation
The FLAME R Package can be installed directly from CRAN:
<div class="code-example" markdown="1">
```r
install.packages('FLAME')
```
</div>

Alternatively, this package can be downloaded from [the author's Github](https://github.com/vittorioorlandi/FLAME):
<div class="code-example" markdown="1">
```r
devtools::install_github('https://github.com/vittorioorlandi/FLAME')
```
</div>

## Quickstart Example
To generate sample data for exploring FLAMEs functionality, use the function `gen_data` as shown below. 
Remember to load the `FLAME` package a shown in line 1 before calling any of the functions discussed 
in this section. This example generates a data frame with n = 250 units and p = 5 covariates:
<div class="code-example" markdown="1">
```r
library('FLAME')

data <- gen_data(n = 250, p = 5)
```
</div>

Note that you can also write the generated dataset to a file by adjusting parameters detailed in the 
next section.

To run the algorithm, use the `FLAME` function as shown in line 3. The required data parameter can 
either be a path to a .csv file or a dataframe. In this example, a .csv file path is used:
<div class="code-example" markdown="1">
```r
library('FLAME')

FLAME_out <- FLAME(data = "data.csv", treated_column_name="treated", outcome_column_name="outcome")
print(FLAME_out$data)
```
</div>
The object FLAME_out is a list of six entries:

| FLAME_out$data:          | a data frame containing the original data with an extra logical column denoting whether a unit was matched and an extra numeric column denoting how many times a unit was matched. The covariates that each unit was not matched on are denoted with asterisks. |
| FLAME_out$MGs:           | a list of every matched group formed by the algorithm.                                                                                                                                                                                                          |
| FLAME_out$CATE:          | a vector containing the conditional average treatment effect (CATE) for every matched group formed.                                                                                                                                                             |
| FLAME_out$matched_on:    | a list corresponding to MGs that gives the covariates, and their values, on which units in each matched group were matched.                                                                                                                                     |
| FLAME_out$matching_covs: | a list containing the covariates that were used for matched on each iteration of the algorithm.                                                                                                                                                                 |
| FLAME_out$dropped:       | a vector of the covariate dropped at each iteration. 

To find the matched groups of particular units after running `FLAME`, use the function `MG` as 
shown below. In this example, the function would return the matched groups of units 1 and 2:

<div class="code-example" markdown="1">
```r
MG(c(1,2),FLAME_out)
```
</div>

To find the CATEs of particular units, use the function `CATE` as shown below. In this example, the 
function would return the matched groups of units 1 and 2:
<div class="code-example" markdown="1">
```r
CATE(c(1,2),FLAME_out)
```
</div>

To find the average treatment effect (ATE) or average treatment effect on the treated (ATT), use 
the functions `ATE` and `ATT`, respectively, as shown below:
<div class="code-example" markdown="1">
```r
ATE(FLAME_out = FLAME_out)
ATT(FLAME_out = FLAME_out)
```
</div>