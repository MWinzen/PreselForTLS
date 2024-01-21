# AUTOMATING THE TRANSIT LEAST SQUARES WITH A PRESELECTION ALGORITHM

The repository contains a set of python files that can be used to automatically download and pre-scan Kepler Light Curves. Afterwards the Transit Least Squares algorithm of Heller and Hippke is executed with promising settings for each light curve.

# Summary

In 2019 René Heller and Michael Hippke published the so-called Transit Least Squares (TLS), which compared to the Box Fitting Least Squares also considers the stellar limb darkening and planetary ingress and egress. This improvement allowed finding another 18 exoplanet candidates by scanning the Kepler light curve data. Also in 2019 Hippke showed that detrending with Tukey’s biweight filter and using a running window of three times the duration of a central transit yields the highest recovery rates for injected transits of simulated data. 

So far, transit durations derived from manually chosen trial periods of circular orbits were used as input to the biweight filter (as for Kepler 160 TLS analysis [1]). The Preselection Algorithm directly identifies potential transits in the light curve in order to use three times the transit duration as running window for the detrended light curve, which is analyzed by the TLS afterwards. The advantage of firstly applying the biweight filter to enable the TLS to find also planets that usually would be hidden by noise is automized by using durations from transit candidates of Preselection Algorithm as input for the running window.

Transit candidates of the light curve adaptive Preselection Algorithm are also compared with each other by depth and duration in order to derive potential exoplanet periods. The result of the TLS can be compared to the potential exoplanet periods of the Preselection Algorithm. Other than the TLS the Preselection Algorithm is not sensitive to data gaps, which permits to automate the validation of TLS results regarding false positives due to data gaps. 

# Requirements

In order to use this repository, please install first the TLS and wotan of René Heller and Michael Hippke:

https://github.com/hippke/tls

https://github.com/hippke/wotan

The following packages are also required

-progressbar

# How to run it:
Here is an example:

import TransitScanConfirmed as TSC

TSC.main("LLC",1,0.005)

Arguments are as follows:

FileFormat:  of Kepler Light Curves ->  LLC or SLC (Please check Kepler Data Release 7 Notes)

min_kid: The Kepler ID to start with

max_f_1: Please refer to PreselectionAlgorithm.pdf
