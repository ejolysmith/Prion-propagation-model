prion-model_time_relative_to_loss.c is the C code that simulates the stochastic model in 
Fig. 5a while saving the time-of-loss to generate the data plotted in 
Fig. 5d. The program called P is the compiled code using the following bash 
command:

gcc prion-model_time_relative_to_loss.c -o P -lm

To run P, you use the following bash command: 

./P a b 

where a is the name (an integer) to put at the end of the saved filename, 
and b is the seed for the random number generation (an integer).
The parameters for the model are all set in the C code. The code outputs a 
txt file, which is taken by the python code plot_time_relative-to-loss.py 
to generate the figure. 
