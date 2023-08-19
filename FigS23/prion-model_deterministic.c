#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

// defining some constants used for the random number generator
#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

// defining the random number generator
float ran1(int *idum)
{
	int j;
	long k;
	static long iy=0;
	static long iv[NTAB];
	float temp;

	if (*idum<=0 || !iy)
	{
        if (-(*idum)<1) *idum=1;
        else *idum=-(*idum);
        for (j=NTAB+7;j>=0;j--)
		{
			k=(*idum)/IQ;
			*idum=IA*(*idum-k*IQ)-IR*k;
			if (*idum<0) *idum+=IM;
			if (j<NTAB) iv[j]=*idum;
		}
        iy=iv[0];
	}

	k=(*idum)/IQ;
	*idum=IA*(*idum-k*IQ)-IR*k;
	if (*idum<0) *idum+=IM;
	j=iy/NDIV;
	iy=iv[j];
	iv[j]=*idum;
	if ((temp=AM*iy)>RNMX) return RNMX;
	else return temp;
}

double beta_x = 0;
double beta_y = 0;
double beta_c = 0;

double step = 0.001;

double cbar = 100;
double ccbar;
double lambda_c; //Change if beta_c not zero

double lambda = 1.75;//10; //Prion production rate
double gammma  ;
double alpha  ;

double alpha_2  ;


long double prion_production_rate(double Volume){
    return lambda*Volume*1/step;
}

long double prion_degradation_rate(double x){
    return beta_x*x/step;
}

long double total_elongation_rate( double x, double Y_tot, double Volume ){
    return gammma*x*Y_tot/Volume/1/step;
}


long double total_fragmentation_rate( double yy , double y0, double Volume, int c){
    return alpha*yy/step;//*(c/Volume)/ccbar ;
    
}

long double total_aggregate_degradation_rate( double Y_tot ){
    return beta_y*Y_tot/step;
}

long double chaperone_production_rate(double Volume){
    return 0*lambda_c*Volume*1/step;
}


//long double chaperone_degradation_rate(int c){
//    return beta_c*c;
//}

double rate_total_function(double x, double Y_tot, double y0){
    return ( prion_production_rate(2) + prion_degradation_rate(x) + total_elongation_rate(x,Y_tot, 1) + total_aggregate_degradation_rate(Y_tot) + chaperone_production_rate(2)  );
}

main(  )
{
	
	clock_t begin = clock();
    srand ( time(NULL) );



    FILE *ff = fopen("prion-simulation-data-model-6-divisions-relative-to-loss.txt", "w");


    int seed = -15;
    //int seed_int = -2;//atoi(argv[3]);
	
    double wait_time = 0;
	int max_aggregate_size = 2000;
	double rate_total;

    int sumsum;
    double x, c;  //Individual prion
	double y[max_aggregate_size]; //The prion aggregate
	double Y_tot; //Total number of aggregate (sum(y))

    int z;
        
    int kk;
	long N_steps;  //Number of steps in the gillespie simulation
	N_steps = 1e9;
	int number_of_components = 7;
    int N_parameters = 1;
	int N_sims = 1   ; //Number of times to run the gillespie simulation
	int sim = 0;
	double y0; //The number of prions not attached to anything (Y1 in slide)
	float ysum = 0;

	double random_number;
	int random_integer;


    double t, t_tot, t_tot_0;
    double tau_c = 50;
    double Volume = 1;
	double Probs[number_of_components];
    double Probs_fragmentation[max_aggregate_size];
	
	double Eta[number_of_components][number_of_components];
	double Means[number_of_components];

    double before, after;
	int i,j,k;
	
	int f = 0;
    int ii;
    int jj;
    int ww;
    int run;
    double xx, yy, YY_tot, cc;
    int nd = 0;

    double average_t_loss;
    double average_number_of_aggregates;
    double average_aggregate_size;
    double average_partitioning_error;
    int total_cell_divisions;
    int total_number_of_sims;
    
    
    double average_concentration_t[100];
    double time_bins[100];
    double number_of_points_in_bins[100];
    double total_fluorescence;
    double predicted_average;
    double t_trunc;
    double t_tot_bin;
    int bin_index;
	

    int is_it_loss;
    
    int cell_div_index;
    
    lambda_c = cbar/( tau_c/log(2)/log(2) );
    ccbar = lambda_c*tau_c/log(2);

    
    for(run = 0; run < N_parameters; run++)
    {
        average_aggregate_size = 0;
        average_number_of_aggregates = 0;
        average_partitioning_error = 0;
        average_t_loss = 0;
        total_cell_divisions = 0;
        total_number_of_sims = 0;
        

        gammma =  1.12200000e-03;//- log(random_number)/2/(pow(10, (double)(gamma_den) ))  ;
        alpha = 1.97000000e-04;//- log(random_number)/(pow(10, (double)(alpha_den) ))  ;
        
        for( i = 0; i<100; i++)
        {
            average_concentration_t[i] = 0;
            time_bins[i] = 10*i;
            number_of_points_in_bins[i] = 0;
        }
        
            
            
            
            
            wait_time = 0;
            sim = 0;
            while(sim < N_sims)
                {
                    

                    x = 0; //prion initial conditon
                    c = (int)(ccbar*Volume);
                    
                    Y_tot = 0;
                    random_number = (ran1(&seed)) ;
                    Volume = 1 + random_number;   //Uniform
                    
                    cell_div_index = 0;
                    is_it_loss = 0;
                    
                    //random_number = 2*log(2)*(1 + (ran1(&seed)));
                    //Volume = - log(random_number/4/log(2))/log(2);    //Exponential
                
                for (i = 0; i<max_aggregate_size; i++) //aggregate initial condition
                    {
                        y[i] = 0;
                        Y_tot = Y_tot + y[i];
                    }
                    y[(int)(126*6*Volume)] = 1;

                    
                    Y_tot = 1;

                    t_tot = 0;
                    t_tot_0 = 0;

                    fprintf(ff, "%lf %d %lf %lf %lf \n", (double)(after - (before - after))/( before )    ,  cell_div_index, t_tot - wait_time, Y_tot/Volume , Y_tot );
                for (i = 0; i < N_steps - 1; i++)
                    {
                        y0 = y[0];
                        rate_total = rate_total_function(x, Y_tot, y0);
                        
                        for (ii = 1; ii < max_aggregate_size; ii++)
                        {
                            rate_total = rate_total + ii*total_fragmentation_rate( y[ii] , y0, 1, c);
                       
                        }

                        random_number = (ran1(&seed)) ;
                        t = 0;
                        t = - log(random_number)/rate_total; //Figure out when a reaction occurs
                        
                        if (t_tot > wait_time){
                            average_number_of_aggregates = average_number_of_aggregates + Y_tot*t;
                        
                            for (ii = 1; ii < max_aggregate_size; ii++)
                            {
                                average_aggregate_size = average_aggregate_size + (ii+1)*y[ii]*t;
                       
                            }
                        }
                        
                        t_tot = t_tot + t;
                        
                        //Cell division:
                        
                        Volume = Volume + exp((t_tot - t_tot_0)*log(2)/tau_c) - exp((t_tot - t_tot_0-t)*log(2)/tau_c);
                        
                        if ( Volume > 2){
                            
                            Volume = Volume - exp((t_tot - t_tot_0)*log(2)/tau_c) - exp((t_tot - t_tot_0-t)*log(2)/tau_c);
                    
                    	    t_tot = tau_c/log(2)*log(2 - Volume -exp((t_tot - t_tot_0 - t)*log(2)/tau_c)) + t_tot_0;
                    
                            t_tot_0 = t_tot;

                            Volume = 1;

                          
                            
                            xx = x;
                            x = x/2;
                            
                            YY_tot = Y_tot;
                            before = xx;
                            after  = x;
                            
                            for(jj = 0; jj < max_aggregate_size; jj++)
                            {
                                yy = y[jj];
                                before = before + y[jj]*(jj+1);
                                y[jj] = y[jj]/2;
                                
                                after = after + y[jj]*(jj+1);
                            }
                            
                            Y_tot = Y_tot/2;
                            
                            if (is_it_loss > 0)
                            {
                                cell_div_index = cell_div_index + 1;
                            }
                            
                            if (t_tot > 0)
                            {
                                fprintf(ff, "%lf %d %lf %lf %lf \n", (double)(after - (before - after))/( before )    ,  cell_div_index, t_tot - wait_time, Y_tot/Volume , Y_tot );
                            }
                            nd = nd + 1;

                    }
                    
                    else
                    {
                        
                        
                        //Figure out which reaction occurs:

                        Probs[0] = prion_production_rate(Volume)/rate_total;
                        Probs[1] = prion_degradation_rate(x)/rate_total;
                        Probs[2] = total_elongation_rate(x, Y_tot, Volume)/rate_total;
                        Probs[3] = total_aggregate_degradation_rate(Y_tot)/rate_total;
                        Probs[4] = chaperone_production_rate(Volume)/rate_total;
                        Probs[5] = 0;//chaperone_degradation_rate(c)/rate_total;
                        Probs[6] = 0;
                        
                        for (ii = 1; ii < max_aggregate_size; ii++)
                        {
                            Probs_fragmentation[ii] = ii*total_fragmentation_rate( y[ii] , y0, Volume, c)/rate_total;
                            Probs[6] = Probs[6] + Probs_fragmentation[ii];
                        }


                        random_number = (ran1(&seed)) ;
                        
                            if (random_number <= Probs[0])
                            {
                                x = x + step;
                                 
                            }
                            
                            else if(random_number <= Probs[0] + Probs[1])
                            {
                                x = x - step;
                            }

                            else if (random_number <= Probs[0] + Probs[1] + Probs[2])
                            {
                                random_number = (ran1(&seed));
                                ysum = ((float)y[0])/((float)Y_tot);
                                
                                for (j = 0; j<max_aggregate_size; j++)
                                {
                                    if (random_number <= ysum )
                                    {
                                        y[j] = y[j] - step;
                                        y[j+1] = y[j+1] + step;
                                        x = x - step;
                                        break;
                                     }
                                         ysum = ysum + ((float)y[j+1])/((float)Y_tot);
                                }
                            }
                            
                            
                            else if (random_number <= Probs[0] + Probs[1] + Probs[2] + Probs[3] )
                            {
                                    
                                random_number = (ran1(&seed));
                                ysum = ((float)y[0])/((float)Y_tot);
                                for (j = 0; j<max_aggregate_size; j++)
                                {
                                    if (random_number <= ysum)
                                    {
                                        y[j] = y[j] - step;
                                        Y_tot = Y_tot - step;
                                        
                                        break;
                                     }
                                    ysum = ysum + ((float)y[j+1])/((float)Y_tot);
                                
                                }
                            }
                        
                            else if (random_number <= Probs[0] + Probs[1] + Probs[2] + Probs[3] + Probs[4] )
                            {
                                    
                                c = c + step;
                            }
                        
                            else if (random_number <= Probs[0] + Probs[1] + Probs[2] + Probs[3] + Probs[4] + Probs[5])
                            {
                                    
                                c = c - step;
                            }
                        
                            else if (random_number < Probs[0] + Probs[1] + Probs[2] + Probs[3] + Probs[4] + Probs[5] + Probs[6])
                            {
                                
                                ysum = Probs_fragmentation[1] + Probs[0] + Probs[1] + Probs[2] + Probs[3] + Probs[4] + Probs[5] ;
                                for (j = 1; j<max_aggregate_size; j++)
                                {
                                    if (random_number <= ysum)
                                    {
                                        random_integer =    rand() % (j) ;
                                        y[j] = y[j] - step;
                                        y[random_integer] = y[random_integer] + step;
                                        y[j - random_integer - 1] = y[j - random_integer - 1] + step;
                                        Y_tot = Y_tot + step;
                                        
                                        break;
                                     }
                                    ysum = ysum + Probs_fragmentation[j+1];
                                }
                            }
                        
                        
                    }
                        
                        
                        fprintf(ff, "%lf %d %lf %lf %lf \n", (double)(after - (before - after))/( before )    ,  cell_div_index, t_tot - wait_time, Y_tot/Volume , Y_tot );
                        
                        
                        
                        if (Y_tot < 0.1 && t_tot < wait_time)
                        {
                            break;
                        }
                        if (Y_tot < 1 && is_it_loss < 1 && t_tot > wait_time)
                        {
                            is_it_loss = 1;
                            average_t_loss = average_t_loss + t_tot;
                            total_number_of_sims = total_number_of_sims + 1;
                        }
                        
                        if (cell_div_index > 9)
                        {
                            printf("%d %lf \n", sim, t_tot);
                            sim = sim + 1;
                            break;
                            
                        }
                        
                        if (t_tot > 100000)
                        {
                            break;
                        }
                         
                     }
                
                    
                }

            
           
            
         
        }
    
    
    
    clock_t end = clock();
    printf("The time spent is: %lf \n\n\n", (double)(end - begin)/ CLOCKS_PER_SEC);

	fclose(ff);
 
        

}
