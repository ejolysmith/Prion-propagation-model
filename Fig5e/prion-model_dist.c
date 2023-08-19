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

double beta_x = 0;//Individual prion degradation constant (if set at zero then there is no degradation)
double beta_y = 0;//aggregate degradation constant (if set at zero then there is no degradation)

double lambda = 1.75; //Prion production rate
double gammma; //Elongation/conversion rate constant, varied below
double alpha; //Fragmentation rate constant, varied below


long double prion_production_rate(double Volume){
    return lambda*Volume*1;
}

long double prion_degradation_rate(int x){
    return beta_x*x;
}

long double total_elongation_rate( int x, int Y_tot, double Volume ){
    return gammma*x*Y_tot/Volume/1;
}


long double total_fragmentation_rate( int yy , int y0 , double Volume ){
    return alpha*yy;
    
}

long double total_aggregate_degradation_rate( int Y_tot ){
    return beta_y*Y_tot;
}

double rate_total_function(int x, int Y_tot, int y0, double Volume){
    return (prion_production_rate(2) + prion_degradation_rate(x) + total_elongation_rate(x,Y_tot, 1) + total_aggregate_degradation_rate(Y_tot));
}

main( int argc, char *argv[] )
{
    
    clock_t begin = clock();
    srand ( time(NULL) );

    int name_id   = atoi(argv[1]);

    char str[50];
    sprintf(str, "prion-data-model-dist-lam1p75-uni-v%d-v2.txt" ,name_id);
    FILE *ff = fopen( str, "w");

    int seed = atoi(argv[2]);
    
    double wait_time = 0;
    int max_aggregate_size = 2000;
    double rate_total;

    int sumsum;
    int x;  //Individual prion
    int y[max_aggregate_size]; //The prion aggregate
    int Y_tot; //Total number of aggregate (sum(y))

    int z;
    
    int number_of_timesteps = 100;
    float timestep = 10;
    float dist[number_of_timesteps][max_aggregate_size];
    float dist_final[number_of_timesteps][max_aggregate_size];
    float dist_x[number_of_timesteps];
    float dist_numb_sims[number_of_timesteps];
    float dist_x_final[number_of_timesteps];
    float timesteps[number_of_timesteps];
    float dist_time[number_of_timesteps];
    float dist_time_final[number_of_timesteps];

    
    int kk;
    long N_steps;  //Number of steps in the gillespie simulation
    N_steps = 1e9;
    int number_of_components = 5;
        int N_parameters = 1;
    int N_sims = 100000   ; //Number of times to run the gillespie simulation
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
    int xx, yy, YY_tot;
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
    double tt = 0;
    int bin_index;
    int time_of_loss_index = 0;
    

    for (ii=0; ii<number_of_timesteps; ii++)
    {   dist_x[ii] = 0;
        dist_time[ii] = 0;
        for (jj=0; jj<max_aggregate_size; jj++)
        {
            dist[ii][jj] = 0;
        }
    }
    
    for (ii=0; ii<number_of_timesteps; ii++)
    {   dist_x_final[ii] = 0;
        dist_numb_sims[ii] = 0;
        dist_time_final[ii] = 0;
        for (jj=0; jj<max_aggregate_size; jj++)
        {
            dist_final[ii][jj] = 0;
        }
    }

    
    for(run = 0; run < N_parameters; run++)
    {
        average_aggregate_size = 0;
        average_number_of_aggregates = 0;
        average_partitioning_error = 0;
        average_t_loss = 0;
        total_cell_divisions = 0;
        total_number_of_sims = 0;
        

        gammma = 1.12200000e-03;//0.001902;//0.001805;//0.001707;//0.001512 ; //1.12200000e-03;//
        alpha = 1.97000000e-04;//0.0009097;//0.0008194;//0.000729;//0.0005484; //1.97000000e-04;//
        
        for( i = 0; i<100; i++)
        {
            average_concentration_t[i] = 0;
            time_bins[i] = 10*i;
            number_of_points_in_bins[i] = 0;
        }
        
        
        
    for(sim = 0; sim<1000; sim++)
        {
            
            x = 0; //prion initial conditon
        
            Y_tot = 0;
            random_number = (ran1(&seed)) ;
            Volume = 1 + random_number;   //Uniform
        
        for (i = 0; i<max_aggregate_size; i++) //aggregate initial condition
            {
                y[i] = 0;
                Y_tot = Y_tot + y[i];
            }
            y[ (int)(126*6*Volume) ] = 1;
  
            Y_tot = 1;

            t_tot = 0;
            t_tot_0 = -log(Volume)*tau_c/log(2);
            t_tot_bin = 0;
            bin_index = 0;
        for (i = 0; i < N_steps - 1; i++)
            {
                y0 = y[0];
                rate_total = rate_total_function(x, Y_tot, y0, Volume);
                
                for (ii = 1; ii < max_aggregate_size; ii++)
                {
                    rate_total = rate_total + ii*total_fragmentation_rate( y[ii] , y0, Volume);
               
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
                Probs[4] = 0;
                
                for (ii = 1; ii < max_aggregate_size; ii++)
                {
                    Probs_fragmentation[ii] = ii*total_fragmentation_rate( y[ii] , y0, Volume)/rate_total;
                    Probs[4] = Probs[4] + Probs_fragmentation[ii];
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
                
                    else if (random_number < Probs[0] + Probs[1] + Probs[2] + Probs[3] + Probs[4])
                    {
                        ysum = Probs_fragmentation[1] + Probs[0] + Probs[1] + Probs[2] + Probs[3] ;
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
                
                total_fluorescence = (float)x;
                for (j = 0; j<max_aggregate_size; j++)
                {
                    total_fluorescence = total_fluorescence + (float)((j+1)*y[j]);
                }
                
                if (  t_tot_bin < t_tot )
                {
                    average_concentration_t[bin_index] = average_concentration_t[bin_index] + total_fluorescence/Volume;
                    number_of_points_in_bins[bin_index] = number_of_points_in_bins[bin_index] + 1;
                    t_tot_bin = t_tot_bin + 10;
                    bin_index = bin_index + 1;
                }
                
            }
                
                
                
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
                if ( fabs(average_concentration_t[j] - predicted_average) < 0.05*predicted_average && fabs(average_concentration_t[j+1] - predicted_average) < 0.05*predicted_average && fabs(average_concentration_t[j+2] - predicted_average) < 0.05*predicted_average )
                { t_trunc = time_bins[j];
                    break;
                }
            }
            
            
            
            
            wait_time = t_trunc;
            sim = 0;
            while(sim < N_sims)
                {
                    for (ii=0; ii<number_of_timesteps; ii++)
                    {
                        dist_x[ii] = 0;
                        dist_time[ii] = 0;
                        for (jj=0; jj<max_aggregate_size; jj++)
                        {
                            dist[ii][jj] = 0;
                        }
                    }

                    x = 0; //prion initial conditon
                
                    Y_tot = 0;
                    random_number = (ran1(&seed)) ;
                    Volume = 1 + random_number;   //Uniform
                    
                
                for (i = 0; i<max_aggregate_size; i++) //aggregate initial condition
                    {
                        y[i] = 0;
                        Y_tot = Y_tot + y[i];
                    }
                    y[(int)(126*6*Volume)] = 1;

                    Y_tot = 1;

                    t_tot = 0;
                    t_tot_0 = -log(Volume)*tau_c/log(2);
                for (i = 0; i < N_steps - 1; i++)
                    {
                        y0 = y[0];
                        rate_total = rate_total_function(x, Y_tot, y0, Volume);
                        
                        for (ii = 1; ii < max_aggregate_size; ii++)
                        {
                            rate_total = rate_total + ii*total_fragmentation_rate( y[ii] , y0, Volume);
                       
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
                            
                            
                            t_tot = t_tot - t;
                            t = t_tot_0 + tau_c - t_tot;
                            
                            
                            if ( t_tot > wait_time )
                            {
                                for (ii=0; ii<number_of_timesteps; ii++)
                                {
                                    
                                    if (ii*timestep < t_tot - wait_time && t_tot - wait_time < (ii+1)*timestep && t_tot + t - wait_time < (ii+1)*timestep)
                                    {
                                        dist_x[ii] = dist_x[ii] + (float)x * tau_c/log(2)*( exp(-log(2)/tau_c*(t_tot - t_tot_0)) - exp(-log(2)/tau_c*(t_tot - t_tot_0 + t)) ) ;
                                        dist_time[ii] = dist_time[ii] + t;
                                        
                                        for (jj=0; jj<max_aggregate_size; jj++)
                                        {
                                            dist[ii][jj] = dist[ii][jj] + (float)y[jj] * tau_c/log(2)*(exp(-log(2)/tau_c*(t_tot - t_tot_0)) - exp(-log(2)/tau_c*(t_tot - t_tot_0 + t)));
                                            
                                        }
                                    }
                                    else if (ii*timestep < t_tot - wait_time&& t_tot - wait_time < (ii+1)*timestep)
                                    {
                                        tt = ((ii+1)*timestep - t_tot + wait_time);
                                        
                                        dist_x[ii] = dist_x[ii] + (float)x * tau_c/log(2)*(exp(-log(2)/tau_c*(t_tot - t_tot_0)) - exp(-log(2)/tau_c*(t_tot - t_tot_0 + tt))) ;
                                        dist_time[ii] = dist_time[ii] + tt;
                                        
                                        
                                        tt = (t_tot - wait_time + t - (ii+1)*timestep);
                                        if (ii < number_of_timesteps-1)
                                        {
                                            dist_x[ii+1] = dist_x[ii+1] + (float)x * tau_c/log(2)*(exp(-log(2)/tau_c*((ii+1)*timestep + wait_time - t_tot_0)) - exp(-log(2)/tau_c*((ii+1)*timestep + wait_time - t_tot_0 + tt))) ;
                                            dist_time[ii+1] = dist_time[ii+1] + tt;
                                        }
                                        
                                        tt = ((ii+1)*timestep - t_tot + wait_time);
                                        for (jj=0; jj<max_aggregate_size; jj++)
                                        {
                                            dist[ii][jj] = dist[ii][jj] + (float)y[jj] * tau_c/log(2)*(exp(-log(2)/tau_c*(t_tot - t_tot_0)) - exp(-log(2)/tau_c*(t_tot - t_tot_0 + tt)));
                                            tt = (t_tot - wait_time + t - (ii+1)*timestep);
                                            if (ii < number_of_timesteps-1)
                                            {
                                                dist[ii+1][jj] = dist[ii+1][jj] + (float)y[jj] * tau_c/log(2)*(exp(-log(2)/tau_c*((ii+1)*timestep + wait_time - t_tot_0)) - exp(-log(2)/tau_c*((ii+1)*timestep + wait_time - t_tot_0 + tt))) ;
                                            }
                                        }
                                    }
                                }
                            }
                            Volume = 1;
                            
                            t_tot = t_tot_0 + tau_c;
                            t_tot_0 = t_tot;
                            
                            Volume = 1;

                            
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
                        t_tot = t_tot - t;
                        if ( t_tot > wait_time )
                        {
                            for (ii=0; ii<number_of_timesteps; ii++)
                            {
                                
                                if (ii*timestep < t_tot - wait_time && t_tot - wait_time < (ii+1)*timestep && t_tot + t - wait_time < (ii+1)*timestep)
                                {
                                    dist_x[ii] = dist_x[ii] + (float)x * tau_c/log(2)*(exp(-log(2)/tau_c*(t_tot - t_tot_0)) - exp(-log(2)/tau_c*(t_tot - t_tot_0 + t))) ;
                                    dist_time[ii] = dist_time[ii] + t;
                                    
                                    for (jj=0; jj<max_aggregate_size; jj++)
                                    {
                                        dist[ii][jj] = dist[ii][jj] + (float)y[jj] * tau_c/log(2)*(exp(-log(2)/tau_c*(t_tot - t_tot_0)) - exp(-log(2)/tau_c*(t_tot - t_tot_0 + t)));
                                       
                                    }
                                }
                                else if (ii*timestep < t_tot - wait_time&& t_tot - wait_time < (ii+1)*timestep)
                                {
                                    tt = ((ii+1)*timestep - t_tot + wait_time);
                                    
                                    dist_x[ii] = dist_x[ii] + (float)x * tau_c/log(2)*(exp(-log(2)/tau_c*(t_tot - t_tot_0)) - exp(-log(2)/tau_c*(t_tot - t_tot_0 + tt))) ;
                                    dist_time[ii] = dist_time[ii] + tt;
                                    
                                    
                                    tt = (t_tot - wait_time + t - (ii+1)*timestep);
                                    if (ii < number_of_timesteps-1)
                                    {
                                        dist_x[ii+1] = dist_x[ii+1] + (float)x * tau_c/log(2)*(exp(-log(2)/tau_c*((ii+1)*timestep + wait_time - t_tot_0)) - exp(-log(2)/tau_c*((ii+1)*timestep + wait_time - t_tot_0 + tt))) ; //* (t_tot - wait_time + t - (ii+1)*timestep)/Volume;
                                        dist_time[ii+1] = dist_time[ii+1] + tt;
                                    }
                                    
                                    tt = ((ii+1)*timestep - t_tot + wait_time);
                                    for (jj=0; jj<max_aggregate_size; jj++)
                                    {
                                        dist[ii][jj] = dist[ii][jj] + (float)y[jj] * tau_c/log(2)*(exp(-log(2)/tau_c*(t_tot - t_tot_0)) - exp(-log(2)/tau_c*(t_tot - t_tot_0 + tt)));
                                        tt = (t_tot - wait_time + t - (ii+1)*timestep);
                                        if (ii < number_of_timesteps-1)
                                        {
                                            dist[ii+1][jj] = dist[ii+1][jj] + (float)y[jj] * tau_c/log(2)*(exp(-log(2)/tau_c*((ii+1)*timestep + wait_time - t_tot_0)) - exp(-log(2)/tau_c*((ii+1)*timestep + wait_time - t_tot_0 + tt))) ;
                                        }
                                    }
                                }
                            }
                        }
                        t_tot = t_tot + t;
                        
                        //Figure out which reaction occurs:

                        Probs[0] = prion_production_rate(Volume)/rate_total;
                        Probs[1] = prion_degradation_rate(x)/rate_total;
                        Probs[2] = total_elongation_rate(x, Y_tot, Volume)/rate_total;
                        Probs[3] = total_aggregate_degradation_rate(Y_tot)/rate_total;
                        Probs[4] = 0;
                        
                        for (ii = 1; ii < max_aggregate_size; ii++)
                        {
                            Probs_fragmentation[ii] = ii*total_fragmentation_rate( y[ii] , y0, Volume)/rate_total;
                            Probs[4] = Probs[4] + Probs_fragmentation[ii];
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
                        
                            else if ( Probs[0] + Probs[1] + Probs[2] + Probs[3] + Probs[4])
                            {
                                ysum = Probs_fragmentation[1] + Probs[0] + Probs[1] + Probs[2] + Probs[3] ;
                                for (j = 1; j<max_aggregate_size; j++)
                                {
                                    if (random_number <= ysum)
                                    {
                                        random_integer =    rand() % (j) ;
                                        y[j] = y[j] - 1;
                                        y[random_integer] = y[random_integer] + 1;
                                        y[j - random_integer - 1] = y[j - random_integer - 1] + 1;
                                        
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
                             for (ii=0; ii<number_of_timesteps; ii++)
                             {
                                 if (ii*timestep < t_tot - wait_time  && t_tot - wait_time  < (ii+1)*timestep)
                                 {
                                     time_of_loss_index = ii;
                                     break;
                                 }
                             }
                             for (ii=0; ii<time_of_loss_index+1; ii++)
                             {   dist_x_final[number_of_timesteps - 1 - ii] = dist_x_final[number_of_timesteps - 1 - ii] + dist_x[time_of_loss_index - ii];
                                 dist_numb_sims[number_of_timesteps - 1 - ii] = dist_numb_sims[number_of_timesteps - 1 - ii] + 1;
                                 dist_time_final[number_of_timesteps - 1 - ii] = dist_time_final[number_of_timesteps - 1 - ii] + dist_time[time_of_loss_index - ii];
                                 for (jj=0; jj<max_aggregate_size; jj++)
                                 {
                                     dist_final[number_of_timesteps -1 - ii][jj] = dist_final[number_of_timesteps -1 - ii][jj] + dist[time_of_loss_index - ii][jj];
                                 }
                             }
                             printf("%d %lf \n", sim, t_tot - wait_time);
                             sim = sim + 1;
                             break;
                         }
                        
                        if (t_tot + wait_time > 3000)
                        {
                            
                            break;
                        }
                         
                     }
                
                    
                }

            
           
         
        
        
        for (ww = 0; ww<number_of_timesteps; ww++)
        {
            fprintf(ff, "%f ", ww*(timestep+0.5) );
            fprintf(ff, "%f ", dist_time_final[ww]);
            fprintf(ff, "%f ", dist_x_final[ww]);
            for (kk = 0; kk<max_aggregate_size; kk++)
            {
                fprintf(ff, "%f ", dist_final[ww][kk]);
                
            }
            
            fprintf(ff, "\n");
            
        }
    }
    
    
    clock_t end = clock();
    printf("The time spent is: %lf \n\n\n", (double)(end - begin)/ CLOCKS_PER_SEC);

    fclose(ff);
 
        

}



