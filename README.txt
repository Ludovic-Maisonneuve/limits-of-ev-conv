Programs related to the following preprint article : Maisonneuve, L., Elias, M., Smadi, C. and Llaurens, V., The limits of evolutionary convergence in sympatry: reproductive interference and historical constraints leading to local diversity in warning traits. Evolution, 75: 149-165. https://doi.org/10.1101/2021.01.22.427743
Corresponding author: ludovic.maisonneuve.2015@polytechnique.org

In this study we developped a mathematical model to investigate how reproductive interference impacts the mimetic relationship between sympatric species.

Ludovic Maisonneuve wrote all the scripts.

Content:

- Model 1: contains the scripts to obtain article figure assuming that trait and preference are fixed in a model species
  - results: to strore data from numerical simulation
  - launch_figXX: generates data for figure XX in the manuscript. It saves data in the folder results
  - open_figXX: generates the figure XX from data in results.
  
- Model 2: contains the scripts to obtain article figure assuming that traits and preferences co-evolve in both species
  - results: to strore data from numerical simulation
  - launch_figXX: generates data for figure XX in the manuscript. It saves data in the folder results
  - open_figXX: generates the figure XX from data in results.

Runs with python3.9
Required package:
- numpy
- os
- matplotlib
