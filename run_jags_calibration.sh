#!/bin/bash -l
#$ -V
#$ -N Run_model_calibration

module load R/3.5.2
Rscript /usr3/graduate/tmccabe/mccabete/Fire_forecast_509/scripts/tess_geo_fork/fire_area_forecast/03.1_Run_Historical_calibration.R
