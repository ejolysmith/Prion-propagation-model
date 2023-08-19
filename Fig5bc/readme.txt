prion-model_big-sweep.c is the C code that simulates the stochastic model over 
many different gamma and alpha parameter values. See the caption in Fig. S14.  
The program called P is the compiled code using the following bash 
command:

gcc prion-model_big-sweep.c -o P -lm

To run P, you use the following bash command: 

./P a b c d 

where a and b determine the range of alpha and gamma to vary (these are 
alpha_den and gamma_den respectively used in lines 174 to 177 of the code),  
c is the name (an integer) to put at the end of the saved filename
d is the seed for the random number generation (an integer), and e
is the number of parameter to vary randomly in the range given by a and b. 
The program saves the average time-of-loss and absolute partitioning errors 
over the range of parameters simulated and the two python codes generate 
contour plots from this scattered data. To find the order of saved 
measurements in the data file, see line 648 of C code. 
