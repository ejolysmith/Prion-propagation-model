# Prion-propagation-model
Code for generating numerical simulation results and modelling figures from 'Measuring prion propagation in single bacteria elucidates mechanism of loss', PNAS 2023.  


Here is the code we used to simulated the stochastic model, along with code to make Figs 5 and 
all figures in the SI that have to do with the numerical modelling. Simulations were performed 
on a 2020 M1 MacBook Air, and C code was compiled used gcc Apple clang version 13.0.0 (clang-1300.0.29.30).
All C code was compiled using the following command: 

gcc code.c -o executable -lm

Figures were generated using the simulation outputs with python 3.9.7 .

Each directory corresponds to a different figure and contains the C code and python files used 
to make the respective figures. There is a readme file in each directory to sort things out. 
All files with "tloss" in the name are to generate simulations and compute the time-of-loss curves 
as shown in Fig 5d. All files with "time-rel-to-loss" are to generate simulations and compute 
the distribution of aggregates prior to loss, as shown in Fig. S12. Both "tloss" and "time-rel-
to-loss" executables take in two inputs as follows: 

./P a b 

where P is the executable, a is the name (an integer) to put at the end of the saved filename, 
and b is the seed for the random number generation (an integer).

All files with "divs-rel-to-loss" in file name are to compute the average partitioning errors 
prior to the loss, as shown in Figs 5h,i. The respective executables do not take any inputs. 

The figure outputs from the python files were edited manually using inkscape to change colors, 
linewidths, legends, and fontsizes, as they are then shown in the paper. The numerical results 
from the python figures are the same as those in the paper.  
