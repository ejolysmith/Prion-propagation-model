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

double cbar = 1;
double ccbar;
double lambda_c; //Change if beta_c not zero

double lambda = 1.75;//10; //Prion production rate
double gammma  = 0.00112;//0.00278; //0.002 + 0*0.00005;//= 0.44;//0.01; //Elongation/conversion rate constant
double alpha  = 0.0003; //+ 20*(0.00015 - 0.0001);//0.001 + 0*0.00005*20;//= 0.1;//0.05  ;//0.02 ; //Fragmentation rate constant


double Vscale = 1;

long double prion_production_rate(double Volume){
    return lambda*Volume*1;
}

long double prion_degradation_rate(int x){
    return beta_x*x;
}

long double total_elongation_rate( int x, int Y_tot, double Volume ){
    return gammma*x*Y_tot/Volume/1;
}


long double total_fragmentation_rate( int yy , int y0, double Volume, int c){
    return alpha*yy*(c/Volume)/ccbar ;
    
}

long double total_aggregate_degradation_rate( int Y_tot ){
    return beta_y*Y_tot;
}

long double chaperone_production_rate(double Volume){
    return lambda_c*Volume*1;   //Here the rate was changed three times to simulate the three dynamics
}


double rate_total_function(int x, int Y_tot, int y0){
    return ( prion_production_rate(2) + prion_degradation_rate(x) + total_elongation_rate(x,Y_tot, 1) + total_aggregate_degradation_rate(Y_tot) + chaperone_production_rate(2)  );
}

main( int argc, char *argv[] )
{
	clock_t begin = clock();
    srand ( time(NULL) );

    
    int name_id   = atoi(argv[1]);

    char str[50];
	//sprintf(str, "prion-simulation-data_-1.txt", name);
    sprintf(str, "prion-data-model-6-tloss-1-v%d.txt" ,name_id);
    double Vscale = 1;
    FILE *ff = fopen( str, "w");

    
    int seed = atoi(argv[2]);
    
    
    double wait_time = 0;
	int max_aggregate_size = 2000;
	double rate_total;

    int sumsum;
	int x;  //Individual prion
    int c;
	int y[max_aggregate_size]; //The prion aggregate	
	int Y_tot; //Total number of aggregate (sum(y))

    int z;
        
    int kk;
	long N_steps;  //Number of steps in the gillespie simulation
	N_steps = 1e9;
	int number_of_components = 7;
        int N_parameters = 1;
	int N_sims = 1000   ; //Number of times to run the gillespie simulation
	int sim = 0;
	int y0; //The number of prions not attached to anything (Y1 in slide)
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
    int xx, yy, YY_tot, cc;
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
        

        gammma =1.12200000e-03;//0.001902;//0.001805;//0.001707;//0.001512 ; //1.12200000e-03;// - log(random_number)/2/(pow(10, (double)(gamma_den) ))  ;
        
        alpha = 1.97000000e-04;//0.0009097;//0.0008194;//0.000729;//0.0005484; //1.97000000e-04;//- log(random_number)/(pow(10, (double)(alpha_den) ))  ;
        
        for( i = 0; i<100; i++)
        {
            average_concentration_t[i] = 0;
            time_bins[i] = 10*i;
            number_of_points_in_bins[i] = 0;
        }
        
    
	for(sim = 0; sim<N_sims; sim++)
		{
            
            x = 0; //prion initial conditon
            c = (int)(ccbar*Volume);
            
            Y_tot = 0;
            random_number = (ran1(&seed)) ;
            Volume = 1 + random_number;   //Uniform
            
        
		for (i = 0; i<max_aggregate_size; i++) //aggregate initial condition
			{
				y[i] = 0;
				Y_tot = Y_tot + y[i];
			}
            y[ (int)(126*6*Volume*Vscale) ] = 1;

            Y_tot = 1;

		    t_tot = 0;
            t_tot_0 = -log(Volume)*tau_c/log(2);
            t_tot_bin = 0;
            bin_index = 0;
            
		for (i = 0; i < N_steps - 1; i++)
			{	
				y0 = y[0];
				rate_total = rate_total_function(x, Y_tot, y0);
                
                for (ii = 1; ii < max_aggregate_size; ii++)
                {
                    rate_total = rate_total + ii*total_fragmentation_rate( y[ii] , y0, Volume, c);
               
                }

                random_number = (ran1(&seed)) ;
				t = 0;
				t = - log(random_number)/rate_total; //Figure out when a reaction occurs
                
                
                
                
                t_tot = t_tot + t;
                
                //Cell division:
                
                Volume = Volume + exp((t_tot - t_tot_0)*log(2)/tau_c) - exp((t_tot - t_tot_0-t)*log(2)/tau_c);
                
                if ( Volume > 2){
                    Volume = 1;
                    
                    t_tot = t_tot_0 + tau_c;
                    
                    t_tot_0 = t_tot;
                    
                    cc = c;
                    for(j = 0; j < cc; j++){
                        random_number = ran1(&seed);
                        if(random_number <= 0.5){
                        c = c - 1;
                        }
                    }
                    
                    xx = x;
                    for(j = 0; j < xx; j++){
                        random_number = ran1(&seed);
                        if(random_number <= 0.5){
                        x = x - 1;
                        }
                    }
                    
                    YY_tot = Y_tot;
                    before = xx;
                    after  = x;
                    
                    for(jj = 0; jj < max_aggregate_size; jj++)
                    {
                        yy = y[jj];
                        before = before + y[jj]*(jj+1);
                        for(j = 0; j < yy; j++){
                            random_number = ran1(&seed);
                            if(random_number <= 0.5){
                                y[jj] = y[jj] - 1;
                                Y_tot = Y_tot - 1;
                            }
                        }
                        after = after + y[jj]*(jj+1);
                    }
                    
                    

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
                    Probs_fragmentation[ii] = ii*total_fragmentation_rate( y[ii] , y0, 1, c)/rate_total;
                    Probs[6] = Probs[6] + Probs_fragmentation[ii];
                }

                //printf("%f \n", Probs[] );
                
                //Y_tot = 0;
                //for (ii = 0; ii < max_aggregate_size; ii++){
                //    Y_tot = Y_tot + y[ii];
                //}
                //printf("%d \n", Y_tot);

                random_number = (ran1(&seed)) ;
                
                    if (random_number <= Probs[0])
                    {
                        x = x + 1;
                        //printf("%d", 1);
                         
                    }
                    
                    else if(random_number <= Probs[0] + Probs[1])
                    {
                        x = x - 1;
                        //printf("%d", 2);
                    }

                    else if (random_number <= Probs[0] + Probs[1] + Probs[2])
                    {   //printf("  %lf  \n", t_tot);
                        random_number = (ran1(&seed));
                        ysum = ((float)y[0])/((float)Y_tot);
                        
                        for (j = 0; j<max_aggregate_size; j++)
                        {
                            if (random_number <= ysum )
                            {
                                y[j] = y[j] - 1;
                                y[j+1] = y[j+1] + 1;
                                x = x - 1;
                                //printf("%d  %lf %lf  \n", j, t_tot, ysum);
                                break;
                             }
                                 ysum = ysum + ((float)y[j+1])/((float)Y_tot);
                        }
                    }
                    
                    
                    else if (random_number <= Probs[0] + Probs[1] + Probs[2] + Probs[3] )
                    {    //printf("%d", 5);
                            
                        random_number = (ran1(&seed));
                        ysum = ((float)y[0])/((float)Y_tot);
                        for (j = 0; j<max_aggregate_size; j++)
                        {
                            if (random_number <= ysum)
                            {
                                y[j] = y[j] - 1;
                                Y_tot = Y_tot - 1;
                                
                                break;
                             }
                            ysum = ysum + ((float)y[j+1])/((float)Y_tot);
                        
                        }
                    }
                
                    else if (random_number <= Probs[0] + Probs[1] + Probs[2] + Probs[3] + Probs[4] )
                    {    //printf("%d", 5);
                            
                        c = c + 1;
                    }
                
                    else if (random_number <= Probs[0] + Probs[1] + Probs[2] + Probs[3] + Probs[4] + Probs[5])
                    {    //printf("%d", 5);
                            
                        c = c - 1;
                    }
                
                    else if (random_number < Probs[0] + Probs[1] + Probs[2] + Probs[3] + Probs[4] + Probs[5] + Probs[6])
                    {   // printf("%d", 4);
                        //random_number = (ran1(&seed));
                        ysum = Probs_fragmentation[1] + Probs[0] + Probs[1] + Probs[2] + Probs[3] + Probs[4] + Probs[5] ;
                        for (j = 1; j<max_aggregate_size; j++)
                        {
                            if (random_number <= ysum)
                            {
                                random_integer =    rand() % (j) ;//+1  //Use a better seed or random number generator?
                                //printf("%d \n", Y_tot);  (max_number + 1 - minimum_number) + minimum_number
//printf("%d", random_integer);
                                y[j] = y[j] - 1;
                                y[random_integer] = y[random_integer] + 1;
                                y[j - random_integer - 1] = y[j - random_integer - 1] + 1;
                                Y_tot = Y_tot + 1;
                                
                                break;
                             }
                            ysum = ysum + Probs_fragmentation[j+1];
                        }
                    }
                
                total_fluorescence = (float)x;
                for (j = 0; j<max_aggregate_size; j++)
                {
                    total_fluorescence = total_fluorescence + (float)((j+1)*y[j]);
                }
                
                if (  t_tot_bin < t_tot )
                {
                    average_concentration_t[bin_index] = average_concentration_t[bin_index] + total_fluorescence/Volume/Vscale;
                    number_of_points_in_bins[bin_index] = number_of_points_in_bins[bin_index] + 1;
                    t_tot_bin = t_tot_bin + 10;
                    bin_index = bin_index + 1;
                }
                
            }
                
                
                
                //printf("%d \n", x);

                //if (y[max_aggregate_size-1] > 0)
                //{
                //    printf("Max size!");
                //    break;
                //}
                
                //fprintf(ff, "%lf ", t_tot);
                //for (kk = 0; kk<max_aggregate_size; kk++)
                //{
               //     fprintf(ff, "%d ", y[kk]);
                //}
				//fprintf(ff, "\n");
                
                
                //sumsum = 0;
                //for (z = 0; z<max_aggregate_size; z++ ){
                //    sumsum = sumsum + y[z]*(z+1);
                //}
                //printf("%d \n", sumsum );
                
                 
                if (t_tot > 1000)
                {   printf("%d %lf \n", sim, t_tot);
                    break;
                }
                
                
            
		 	}
            
		

            
        }
        
        
        
            for( j = 0; j<100; j++)
            {
                average_concentration_t[j] = average_concentration_t[j]/number_of_points_in_bins[j];
            }
        
            
            predicted_average = lambda*tau_c/log(2)*1;
            t_trunc = 0;
            for( j = 0; j<98; j++)
            {
                if ( fabs(average_concentration_t[j] - predicted_average) < 0.025*predicted_average && fabs(average_concentration_t[j+1] - predicted_average) < 0.025*predicted_average && fabs(average_concentration_t[j+2] - predicted_average) < 0.025*predicted_average )
                { t_trunc = time_bins[j];
                    break;
                }
            }
            //printf("%lf  %lf \n", fabs(average_concentration_t[80] - predicted_average), number_of_points_in_bins[60] );
            //t_loss = t_loss + t_tot;
         
            
            
            //t_loss = t_loss/N_sims;
            //fprintf(ff, "%lf %lf %lf \n", t_loss, gamma, alpha);
            //printf("%lf", t_loss);
            
            
            
            
            wait_time = t_trunc;
            sim = 0;
            while(sim < 1000)
                {
                    

                    x = 0; //prion initial conditon
                    c = (int)(ccbar*Volume);
                    
                    Y_tot = 0;
                    random_number = (ran1(&seed)) ;
                    Volume = 1 + random_number;   //Uniform
                    
                    //random_number = 0.5*log(2)*log(2)*( 1 + (ran1(&seed)));
                    //Volume = - log(random_number/2/log(2)/log(2))/log(2);    //Exponential
                    
                    //Volume = 1; //Identical

                for (i = 0; i<max_aggregate_size; i++) //aggregate initial condition
                    {
                        y[i] = 0;
                        Y_tot = Y_tot + y[i];
                    }
                    y[(int)(126*6*Volume*Vscale)] = 1;
                    //y[100] = 1;
                    //y[300] = 1;
                    //y[200] = 1;

                    Y_tot = 1;

                    t_tot = 0;
                    t_tot_0 = -log(Volume)*tau_c/log(2);
                    //t_tot_bin = 0;
                    //bin_index = 0;
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
                            
                            Volume = 1;//Volume - exp((t_tot - t_tot_0)*log(2)/tau_c) - exp((t_tot - t_tot_0-t)*log(2)/tau_c);
                            
                            t_tot = t_tot_0 + tau_c;//tau_c/log(2)*log(2 - Volume -exp((t_tot - t_tot_0 - t)*log(2)/tau_c)) + t_tot_0;
                            
                            t_tot_0 = t_tot;
                            
                            cc = c;
                            for(j = 0; j < cc; j++){
                                random_number = ran1(&seed);
                                if(random_number <= 0.5){
                                c = c - 1;
                                }
                            }
                            
                            xx = x;
                            for(j = 0; j < xx; j++){
                                random_number = ran1(&seed);
                                if(random_number <= 0.5){
                                x = x - 1;
                                }
                            }
                            
                            YY_tot = Y_tot;
                            before = xx;
                            after  = x;
                            
                            for(jj = 0; jj < max_aggregate_size; jj++)
                            {
                                yy = y[jj];
                                before = before + y[jj]*(jj+1);
                                for(j = 0; j < yy; j++){
                                    random_number = ran1(&seed);
                                    if(random_number <= 0.5){
                                        y[jj] = y[jj] - 1;
                                        Y_tot = Y_tot - 1;
                                    }
                                }
                                after = after + y[jj]*(jj+1);
                            }
                            
                            
                            if (t_tot > wait_time){
                            average_partitioning_error = average_partitioning_error + fabs( (double)(after - (before - after)) )/( before );
                            total_cell_divisions = total_cell_divisions + 1;
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
                                x = x + 1;
                                 
                            }
                            
                            else if(random_number <= Probs[0] + Probs[1])
                            {
                                x = x - 1;
                            }

                            else if (random_number <= Probs[0] + Probs[1] + Probs[2])
                            {
                                random_number = (ran1(&seed));
                                ysum = ((float)y[0])/((float)Y_tot);
                                
                                for (j = 0; j<max_aggregate_size; j++)
                                {
                                    if (random_number <= ysum )
                                    {
                                        y[j] = y[j] - 1;
                                        y[j+1] = y[j+1] + 1;
                                        x = x - 1;
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
                                        y[j] = y[j] - 1;
                                        Y_tot = Y_tot - 1;
                                        
                                        break;
                                     }
                                    ysum = ysum + ((float)y[j+1])/((float)Y_tot);
                                
                                }
                            }
                        
                            else if (random_number <= Probs[0] + Probs[1] + Probs[2] + Probs[3] + Probs[4] )
                            {
                                    
                                c = c + 1;
                            }
                        
                            else if (random_number <= Probs[0] + Probs[1] + Probs[2] + Probs[3] + Probs[4] + Probs[5])
                            {
                                    
                                c = c - 1;
                            }
                        
                            else if (random_number < Probs[0] + Probs[1] + Probs[2] + Probs[3] + Probs[4] + Probs[5] + Probs[6])
                            {
                                ysum = Probs_fragmentation[1] + Probs[0] + Probs[1] + Probs[2] + Probs[3] + Probs[4] + Probs[5] ;
                                for (j = 1; j<max_aggregate_size; j++)
                                {
                                    if (random_number <= ysum)
                                    {
                                        random_integer =    rand() % (j) ;
                                        y[j] = y[j] - 1;
                                        y[random_integer] = y[random_integer] + 1;
                                        y[j - random_integer - 1] = y[j - random_integer - 1] + 1;
                                        Y_tot = Y_tot + 1;
                                        
                                        break;
                                     }
                                    ysum = ysum + Probs_fragmentation[j+1];
                                }
                            }
                        
                        
                    }
                        
                        
                        
                        
                        
                        
                        if (Y_tot < 1 && t_tot < wait_time)
                        {
                            break;
                        }
                        
                        
                        if (Y_tot < 1 && t_tot > wait_time)
                        {
                            fprintf(ff, "%lf \n", t_tot - wait_time);
                            printf("%d %lf \n", sim, t_tot - wait_time);
                            average_t_loss = average_t_loss + t_tot - wait_time;
                            total_number_of_sims = total_number_of_sims + 1;
                            sim = sim + 1;
                            break;
                        }
                        
                        if (t_tot > 3000)
                        {
                            fprintf(ff, "%lf \n", t_tot - wait_time);
                            printf("%d %lf \n", sim, t_tot - wait_time);
                            average_t_loss = average_t_loss + t_tot - wait_time;
                            total_number_of_sims = total_number_of_sims + 1;
                            break;
                        }
                         
                     }
                
                    
                }

            
            average_number_of_aggregates = average_number_of_aggregates/average_t_loss;
            average_aggregate_size = average_aggregate_size/average_t_loss/average_number_of_aggregates;
            
            average_t_loss = average_t_loss/total_number_of_sims;
            
            average_partitioning_error = average_partitioning_error/total_cell_divisions;
            
       
    }
    
 
    clock_t end = clock();
    printf("The time spent is: %lf \n\n\n", (double)(end - begin)/ CLOCKS_PER_SEC);

	fclose(ff);
 
        

}
