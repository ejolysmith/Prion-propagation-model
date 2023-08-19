prion-model_cell-cycle-position.c is the C code that simulates the 
stochastic model with parameters given by orange dot in Fig. 5b, 
1000 times and saves the cell-cycle position at the time of loss 
each time. 
The program called P is the compiled code using the following bash 
command:

gcc prion-model_cell-cycle-position.c -o P -lm

To run P, you use the following bash command: 

./P a b 

where a is the name (an integer) to put at the end of the saved filename, 
and b is the seed for the random number generation (an integer).
The parameters for the model are all set in the C code. gamma and alpha 
are set at line 180. The python file plot_cell-cycle-position-at-time-of-loss.py 
takes in the data from the C code and makes Fig. 5g. 
