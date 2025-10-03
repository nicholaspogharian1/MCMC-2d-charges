#include<stdio.h>
#include<math.h>

#define PI 3.1415926536

// generate gaussian random_numbers with mean 0 and variance 1
double gaussran()
{                                                                                                                                                                           
    int g;
    static double ran1, ran2;
    static int phase=0;
    double res;
    if (phase == 0) {
    	ran1 = (rand()+1.)/(RAND_MAX+2.);
    	ran2 = (double)rand()/(RAND_MAX+1.);
    	res = sqrt(-2 * log(ran1)) * sin(2 * PI * ran2);
    }
    else {
        res = sqrt(-2 * log(ran1)) * cos(2 * PI * ran2);
    }
    phase = 1 - phase;
    
    
    return res;
}
