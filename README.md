# rHCL

## Description
This package implements the Homeostatic-Circadian-Light (HCL) model in R, as described by Skeldon et al. (2023) <https://doi.org/10.1371/journal.pcbi.1011743>. 
It estimates sleep pressure (i.e., mu) and circadian period (i.e., tau) parameters from light exposure and sleep data. Also estimates the daily 
circadian effect of light with respect to the speeding up or slowing down of the circadian pacemaker.

## Usage
The primary function is `rhcl()`, which estimates $\\mu$ and $\\tau$ from observed light data and target sleep duration/midpoint. Predicted changes in the circadian pacemaker to light (i.e., advancement or delay) for each day can be estimated using `circLight()`. See the respective documentation for each function for further details. 
