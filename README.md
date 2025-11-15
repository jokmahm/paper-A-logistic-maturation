# Age–Length Structured Population Model (Paper A)

This repository contains the R and C++ (TMB) code used in **Paper A of my PhD thesis**, where an age–length structured model is developed for analysing capelin stock dynamics. The work combines R for data processing and estimation routines with **TMB (Template Model Builder)** C++ templates for fast likelihood evaluation. A small amount of MATLAB code is included for optional figure generation.

All scripts required to reproduce the results in the paper are included here.

---

## How to Run the R + TMB Model

1. Install required R packages:
    ```r
    install.packages(c("TMB", "ggplot2", "dplyr"))
    ```

2. Compile the TMB C++ template:
    ```r
    library(TMB)
    compile("mature.cpp")
    dyn.load(dynlib("mature"))
    ```

3. Run the main estimation script:
    ```r
    source("Main.R")
    ```

4. Run different methods as:
    ```r
    source("Method1.R")
    source("Method2.R")
    ```

These scripts load the input data, construct the model, run optimization via TMB, and produce the results presented in the paper.

---

## Markov properties

Results related to Markov properties can be produced using a peice of MATLAB code. To run these:

```matlab
cd matlab
run("markov.m")
