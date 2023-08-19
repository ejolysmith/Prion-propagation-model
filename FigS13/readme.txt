The C program compiled from prion-model_lambda.c is run for different lambda values, 
which are set manually in the code on line 53. The output is an array of partitioning errors, 
we take the average absolute partitioning errors and save them in the python file 
plot_lambda.py for each lambda along with their error taken as the standard error of the mean. 
We compile prion-model_lambda.c as follows
gcc prion-model_lambda.c -o P -lm

where P is run as ./P 

The python file generates the figure
