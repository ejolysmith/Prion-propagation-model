prion-model_divisions_relative_to_loss.c is the C code that simulates 
the stochastic model in Fig. 5a while saving the partitioning errors to 
generate the data plotted in Fig. 5h,i. 
The program called P is the compiled code using the following bash 
command:

gcc prion-mode_divisions_relative_to_loss.c -o P -lm

To run P, you use the following bash command: 

./P

The parameters for the model are all set in the C code. gamma and alpha 
are varied manually at line 180 of the code. The plot_fig5h.py and 
plot_fig5i.py python codes then takes in this data file to generate 
the figures. 
