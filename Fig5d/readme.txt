prion-model_tloss.c is the C code that simulates the stochastic model in 
Fig. 5a while saving the time-of-loss to generate the data plotted in 
Fig. 5d. The program called P is the compiled code using the following bash 
command:

gcc prion-model_tloss.c -o P -lm

To run P, you use the following bash command: 

./P a b 

where a is the name (an integer) to put at the end of the saved filename, 
and b is the seed for the random number generation (an integer).
The parameters for the model are all set in the C code. gamma and alpha 
are varied manually at line 180 of the code. The code is run 5 times, each 
time with a different gamma, alpha parameter pair, to create 5 output files 
(saved in the Data directory). The plot_tloss.py python code then takes in 
these data files to generate Fig. 5d. 
