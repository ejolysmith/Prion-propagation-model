prion-model_equilibrium.c is the C code that simulates the stochastic model in 
Fig. 5d to generate the curve in Fig S12. 
 The program called P is the compiled code using the following bash 
command:

gcc prion-model_equilibrium.c -o P -lm

To run P, you use the following bash command: 

./P a b 

where a is the name (an integer) to put at the end of the saved filename, 
and b is the seed for the random number generation (an integer).
The parameters for the model are all set in the C code.
The python file plot_equilibrium.py generates Fig S12. 

