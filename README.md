# uCTDread

This repository contains a MATLAB function (uCTDread.m) to import and process  data from underway CTD ASCII files. The function isolates downcast measurements and transforms them in MATLAB format after computing absolute salinity and potential density anomaly. Temperature and conductivity measurements are aligned before calculating salinity by applying a negative time lag to temperature measurements. The effect of the time lag correction can be checked in the plot produced by the function where uncorrected and corrected salinity are shown in the same panel.

Benedetto Barone - Oct 2017
