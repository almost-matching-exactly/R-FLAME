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

## Installation
The FLAME R Package can be downloaded from [the author's Github](https://github.com/vittorioorlandi/FLAME), via:
<div class="code-example" markdown="1">
```r
devtools::install_github('https://github.com/vittorioorlandi/FLAME')
```
</div>

## Quickstart Example
To generate sample data for exploring the functionality of FLAME and DAME, we can use the `gen_data` function included in the package, as shown below. 
Remember to load the `FLAME` package, as shown in line 1, before calling any of the functions discussed in this section. This example generates a data frame with n = 250 units and p = 5 categorical covariates, whose outcome is a linear function of treatment and 
covariate values:
<div class="code-example" markdown="1">
```r
library(FLAME)

data <- gen_data(n = 250, p = 5)
```
</div>

To run either the FLAME or DAME algorithm, use the `FLAME` or `DAME` functions, respectively. Here, we will focus on `FLAME`, but the syntax and output for `DAME` is typically identical. The data to match on can be supplied either as a data frame or
as a .csv file to be read into memory. Here, we'll use the `data` data frame we created above.
<div class="code-example" markdown="1">
```r
FLAME_out <- FLAME(data = data)
```
The names of the columns denoting treatment and outcome have default values of '
treated' and 'outcome', respectively, although you can supply whatever values you
have in your data using the `treated_column_name` and `outcome_column_name` arguments. 
</div>
The object output by `FLAME` is an object of class `ame`. By default, this consists of 
three entries:

| FLAME_out$data:          | a data frame containing the original data with an extra logical column denoting whether a unit was matched and an extra numeric column denoting how many times a unit was matched. The covariates that each unit was not matched on are denoted with asterisks. |
| FLAME_out$MGs:           | a list whose i'th entry contains the main matched group for unit i.                                                                                                                                                                   |
| FLAME_out$CATE:          | a vector whose i'th entry contains the conditional average treatment effect (CATE) estimate for the i'th unit.                                         

To find the units contained in a given matched group, we can simply look at the relevant entry of the `MG` list in the output of `FLAME`. In this example, we retrieve the units in the main matched group of unit 1:

<div class="code-example" markdown="1">
```r
FLAME_out$MG[[1]]
```
</div>

If we are interested in the covariate values that these units matched on, we can 
use the `MG` function as follows:

<div class="code-example" markdown="1">
```r
MG(1, FLAME_out)
```
</div>


To find the CATE estimates of particular units, we can simply look at the relevant entry of the `CATE` vector in the output of `FLAME`. For example, the code below extracts all the CATE estimates of the units in the main matched group of unit 1. Note that these not all be the same in case some subset of those units first matched on more covariates than they share with unit 1.
<div class="code-example" markdown="1">
```r
FLAME_out$CATE[FLAME_out$MGs[[1]]]
```
</div>

To estimate the average treatment effect (ATE) or average treatment effect on the treated (ATT), use 
the functions `ATE` and `ATT`, respectively, as shown below:
<div class="code-example" markdown="1">
```r
ATE(FLAME_out)
ATT(FLAME_out)
```
</div>
