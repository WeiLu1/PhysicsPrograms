/*
A GAUSSIAN SAMPLE 3 WAYS
VERSION 3.1
WEI LU
14/12/2018

This program utilises a number of different methods to generate random numbers (RNGs).
To utilise specific methods require the user to input a command line argument.
(RNG = " -a", Approximate method = " -b", Box-Muller method = " -c" and Rejection Method = " -d".) The Rejection Method also utilises Inverse Transforms.
To shift the Gaussian distribution needed for the second part of the work, new constants were added: "MU" is the mean and "VAR" is the variance.
The program will write the random number results to the a text file (same file is used independent of user choice for method), whereupon gnuplot is used
to plot the histogram. The script "run_plot.plt" was written to automate the plot of histograms in gnuplot. It is ran in gnuplot with the command:
' load "run_plot.plt" '.
Care must be taken in remembering which method the user has used to generate the random numbers, as the text file is over written each time.
On the command line, information will be presented on the data set, including mean, standard deviation etc. An execution timer will also be shown that is
used for analysis on the speeds of each methods.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#define RANDOM_MAX 0x7FFFFFFF
#define LOWER 3.0
#define UPPER 7.0
#define MU 5.0
#define VAR 2.0
#define N 1e4   //Number of Samples
#define APPROX_METHOD_LIMIT 12
#define PI M_PI
#define MAXEVAL 1000

static void error_message();
static void * xmalloc(size_t bytes);
static double RNG(void);    /* Variates from probability density f(X) */
static double PDF(double x); /* Probability density f(X) as in RNG() above */
static void analysis(double *result);
static double approx_method(void);
static double box_mull_method(void);
static double rejection_method(void);

typedef enum {Accept, Reject} condition;
static int box_mull_counter = 0;

int main(int argc, char *argv[])
{
    double (*func_ptr)(void);
    int i;
    double time_took, PDF_out;
    double *randnum_results = xmalloc(N * sizeof(double));
    FILE *rand_num_list = fopen("RandomNumbers.txt", "w");
    FILE *pdf_list = fopen("PDF.txt", "w");

    if (argc == 2){
        if (strcmp(argv[1], "-a") == 0){
            func_ptr = &RNG;
        }
        else if (strcmp(argv[1], "-b") == 0){
            func_ptr = &approx_method;
        }
        else if (strcmp(argv[1], "-c") == 0){
            func_ptr = &box_mull_method;
        }
        else if (strcmp(argv[1], "-d") == 0){
            func_ptr = &rejection_method;
        }
        else{
            error_message();
            return 0;
        }
    }
    else{
        error_message();
        return 0;
    }
    clock_t t;
    t = clock();
    srand(time(NULL));

    for (i = 0; i < N; i++){
        randnum_results[i] = func_ptr();
        PDF_out = PDF(randnum_results[i]);
        fprintf(rand_num_list, "%f\n", randnum_results[i]);
        fprintf(pdf_list, "%f \t %f\n", randnum_results[i], PDF_out);
    }
    analysis(randnum_results);

    t = clock() - t;
    time_took = ((double)t)/CLOCKS_PER_SEC;
    printf("Execution time: %f seconds. \n", time_took);
    free(randnum_results);
    fclose(rand_num_list);
    fclose(pdf_list);
    return 0;
}

static void error_message(){
    printf("Incorrect argument. \n");
    printf("To choose RNG, enter: ' -a' as an argument. \nTo choose Approximate Method, enter: ' -b' as an argument. \n");
    printf("To choose Box-Muller Method, enter: ' -c' as an argument. \nTo choose Rejection Method, enter: ' -d' as an argument. \n");
}

static void * xmalloc(size_t bytes){
    void *retVal = malloc(bytes);
    if (retVal) {
        return retVal;
    }
    printf("ERROR: Memory not allocated.");
    exit(99);
}

static double RNG(void){
    return((rand() / (double)RANDOM_MAX) * (UPPER - LOWER) + LOWER);
}

static double approx_method(void){
    int num;
    double total, u_rand_num, z_rand_num;

    total = 0.0;
    for (num = 0; num < APPROX_METHOD_LIMIT; num++){
        u_rand_num = rand() / (double)RANDOM_MAX ;
        total += u_rand_num;
    }
    z_rand_num = total - 6;
    return MU + (z_rand_num * sqrt(VAR));
}

static double box_mull_method(void){
    double u_rand_num1 = rand() / (double)RANDOM_MAX;
    double u_rand_num2 = rand() / (double)RANDOM_MAX;
    double z_rand_num;

    if (box_mull_counter % 2 == 0){
        z_rand_num = sqrt(-2 * log(u_rand_num1)) * cos(2 * PI * u_rand_num2);
    }
    else{
        z_rand_num = sqrt(-2 * log(u_rand_num1)) * sin(2 * PI * u_rand_num2);
    }
    box_mull_counter++;
    return MU + (z_rand_num * sqrt(VAR));
}

static double rejection_method(void){
    double p_x, a_x, u_rand_num1, u_rand_num2, Y;
    int count = 0;

    do{
        u_rand_num1 = rand() / (double)RANDOM_MAX;
        u_rand_num2 = rand() / (double)RANDOM_MAX;
        Y = -tan((PI / 2) - (PI * u_rand_num1)); // normalised inverse CDF equation
        p_x = 1 / sqrt(2 * PI) * exp(-0.5 * pow(Y, 2));
        a_x = 0.5 / (1 + pow(Y, 2));
        if (u_rand_num2 <= p_x / a_x){
            return MU + (sqrt(VAR) * Y);
        }
        if (count > MAXEVAL){
            printf("Error: conditions not met.");
            exit(99);
        }
        count++;
    }while(Reject == 1);
}

static double PDF(double x){
    if ((LOWER < x) && (x < UPPER)){
        return 1 / (UPPER - LOWER);
    }
    else if (x == LOWER || x == UPPER){
        return 0.5 / (UPPER - LOWER);
    }
    return 0.0; // Only occurs if x is NaN
}

static void analysis(double *result){
    int counter;
    double maximum, minimum, mean, variance, stan_dev, SEM;
    double *holder = xmalloc(N * sizeof(double));

    maximum = result[0];
    minimum = result[0];

    for (counter = 0; counter < N; counter++){
        if (result[counter] > maximum){
            maximum = result[counter];
        }
        if (result[counter] < minimum){
            minimum = result[counter];
        }
        mean += (result[counter] / N);
    }

    for (counter = 0; counter < N; counter++){
        holder[counter] = pow((result[counter] - mean), 2);
        variance += holder[counter] / N;
    }
    stan_dev = sqrt(variance);
    SEM = stan_dev / sqrt(N);

    free(holder);
    printf("Number of Samples: %.0f \n", N);
    printf("Maximum: %f \n", maximum);
    printf("Minimum: %f \n", minimum);
    printf("Mean: %f \n", mean);
    printf("Variance: %f \n", variance);
    printf("Standard deviation: %f\n", stan_dev);
    printf("Standard error of the mean: %f\n", SEM);
}
