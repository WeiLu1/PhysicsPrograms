/*
BESSEL FUNCTION ROOT FUNDING PROGRAM
VERSION 1.1
WEI LU
02/11/2018

This program outputs the roots of the Bessel function of the first kind, from a range between 0 and 15.
The initial estimates of the roots are found by using a for loop to iterate through the range. From the two Bessel function results given from the
for loop outputs, there is an inequality check. If the product of the two Bessel functions are negative, there is a root present in that interval.
To initiate finding of the root, the user must pass a command line argument to choose secant or bisection method for this. Failure
to provide an argument will result in the program not running at all.
The total tolerance calculation of the root utilises the absolute and relative tolerance, which are constants in this program.
*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#define ABSTOL 1e-15
#define RELTOL 1e-12


int solver(double bessel(double x), double x1, double x2, int maxeval, int choice);
void bisect(double input1, double input2, double *mid);
void secant(double in1, double in2, double res1, double res2, double *mid);
double tolerance_calc(double x_est, double x_true, double *uncertain);


int main(int argc, char *argv[])
{
    double j0(double x);
    double guess1, guess2, j_out1, j_out2;
    int option;

    if (argc == 2){
        if (strcmp(argv[1], "-b") == 0 ){
            option = 1;
        }
        else if (strcmp(argv[1], "-s") == 0 ){
            option = 2;
        }
        else{
            printf("Incorrect argument. \n");
            printf("To choose secant, enter: ' -s' as an argument. \nTo choose bisection, enter: ' -b' as an argument.");
            return 0;
        }
        for (guess1 = 0.0; guess1 < 15.0; guess1++){
            guess2 = guess1 + 1.0;
            j_out1 = j0(guess1);
            j_out2 = j0(guess2);
                if (j_out1 * j_out2 < 0){
                   solver(j0, guess1, guess2, 100, option);
                }
        }
    }
    else{
        printf("Incorrect argument. \n");
        printf("To choose secant, enter: ' -s' as an argument. \nTo choose bisection, enter: ' -b' as an argument.");
    }
    return 0;
}


int solver(double bessel(double x), double x1, double x2, int maxeval, int choice){

    double result1, result2, x_out, x_root, tolerance, uncertainty;
    int iteration = 0;

    do{
        result1 = bessel(x1);
        result2 = bessel(x2);

        if (choice == 1){   // BISECTION
            bisect(x1, x2, &x_out);
            if (bessel(x_out) * result1 > 0.0){
                x1 = x_out;
            }
            else{
                x2 = x_out;
            }
            bisect(x1, x2, &x_root);
        }

        else if (choice == 2){  // SECANT
            secant(x1, x2, result1, result2, &x_out);
            x1 = x2;
            x2 = x_out;
            secant(x1, x2, bessel(x1), bessel(x2), &x_root);
            }

        tolerance = tolerance_calc(x_out, x_root, &uncertainty);
        if (fabs(x_root - x_out) < tolerance){
            printf("Root: %.20f \tUncertainty: %.5le \tIteration number: %d\n", x_root, uncertainty, iteration);
            return 0;
            }
        iteration++;
    }
    while(iteration < maxeval);
    printf("Maximum evaluation limit reached. \n");
    return 0;
}


void bisect(double in1, double in2, double *mid){
    *mid = (in1 + in2) / 2.0;
}


void secant(double in1, double in2, double res1, double res2, double *mid){
    *mid = in2 - ((res2 * (in2 - in1)) / (res2 - res1));
}


double tolerance_calc(double x_est, double x_true, double *uncertain){
    double tol_tot;

    *uncertain = fabs((x_true - x_est) / x_est);
    tol_tot = ABSTOL + RELTOL * fabs(x_est);
    return tol_tot;
}
