prion-model_tloss-initial-seed.c is compiled into a program that runs the stochastic model 
1000 times while saving the time of prion loss. Parameters are set in the code. 
The initial seed size is set by the parameter seed_factor on line 51. Varying 
this parameter changes the initial seed size. 
To compile we do the following bash command

gcc prion-model_tloss-initial-seed.c Pseed -lm

and to run the process we do ./P a b
where a is the name (an integer) to put at the end of the saved filename, 
and b is the seed for the random number generation (an integer). The python file 
plot_tloss-initial-seed.py generates Fig S15a. 

prion-model_tloss-initial-size.c does the same thing but is used to vary the 
initial cell size distribution. This is done manually three times on lines 200 and 433. 
It's compiled similarly to above to Psize. The python code plot_tloss-initial-volume.py 
generates Fig S15b. 
