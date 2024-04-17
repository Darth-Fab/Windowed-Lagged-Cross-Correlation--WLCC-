# Windowed-Lagged Cross-Correlation for Dyadic Studies
This repository hosts an R script for conducting Windowed-Lagged Cross-Correlation analysis in dyadic studies exploring trust and prosocial behavior. The script enables the assessment of synchrony strength between heart rate (HR) and skin conductance (SC) signals at different time lags. Future updates will include support for analyzing facial expressions.

## Features
- Dyadic Interaction Synchrony Estimation with WLCC function: Evaluate HR-SC synchrony in real dyadic interactions
- Temporal Dynamics: Explore synchrony strength at varying time lags by customizable parameters (Window-Size, Window-increase, Lag-size, Lag-increase)

## Stuff that I am fixing
- The script is still hard to customize. I am working on making it more accessible and easy to tailor
- Trial-by-trial analysis is still quite a mess. The script unnecessarily saves several data files for each trial, which is a problem with a large number of trials 
- Still need to implement the section for Facial Synchrony - estimated with OpenFace output 
