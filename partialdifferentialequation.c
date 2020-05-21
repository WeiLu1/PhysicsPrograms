/*
HW5: SOLVING A PARTIAL DIFFERENTIAL EQUATION
VERSION 1.6
WEI LU
18/01/2019
For use of GSL libraries, run with: gcc main.c -Wall -I/usr/local/include -L/usr/local/lib -lgsl -lgslcblas -lm -o test.exe

This program calculates the temperature change over time from the radioactive decay around a nuclear waste rod. This is done by solving the partial
differential equation that is provided by the works of Dr Louise Olsen-Kettle:
(page 43 of http://espace.library.uq.edu.au/eserv/UQ:239427/Lectures_Book.pdf).

Solving this heat equation utilises the backward Euler Method subject to Neumann and Dirichlet boundary conditions (differential of temperature with time
zero radius is zero and temperature at the end of the rod is 300K, respectively).

This program uses the GSL libraries for solving the linear algebra, LU decomposition to find the new time forward temperature vector (essentially Ax = b).
As well as this, it is used for matrix and vectors.

Constants can be altered for the initial set up of the problem, the user can change the length of the rod, the time step between each temperature value,
the number of points used to make up the different sections of the rod etc.

Information is written to a text file under three columns: time, radius and temperature. Then GNUplot is used to plot the resulting graph, where only the
second and third columns are considered.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

#define R_C 100.0
#define N 99.0
#define M 1000.0
#define TIME_LENGTH_SOL 100.0
#define KAPPA 2.0e7
#define T_ROD 1.0
#define T_BOUND 300.0
#define TAU_0 100.0
#define a 25.0

static void setup_matrix_vec(gsl_matrix* A_0, gsl_vector* T_k, gsl_vector* S_k, gsl_vector* r_0, gsl_vector* b_0, double s, double dr);
static void solver(double time_step, gsl_matrix* A_0, gsl_permutation* p_0, gsl_vector* T_k, gsl_vector* S_k, gsl_vector* b_0);

int main()
{
    double delta_r, delta_time, gain_parameter;
    double time_stepper = 0.0;
    int count1, count2, s;
    const int oneyear = 9;
    const int tenyear = 99;
    const int fiftyyear = 499;
    const int hundyear = 999;

    FILE * f = fopen("decayresultsdata.txt", "w");
    FILE * Init = fopen("Initialconditions.txt", "w");

    gsl_vector * r = gsl_vector_alloc(N);
    gsl_vector * b = gsl_vector_alloc(N);
    gsl_vector * Sk = gsl_vector_alloc(N);
    gsl_vector * Tk = gsl_vector_alloc(N);
    gsl_matrix * A = gsl_matrix_alloc(N, N);

    delta_r = R_C / (N + 1);
    delta_time = TIME_LENGTH_SOL / M;
    gain_parameter = (KAPPA * delta_time) / pow(delta_r, 2);

    setup_matrix_vec(A, Tk, Sk, r, b, gain_parameter, delta_r);

    for (count1 = 0; count1 < N; count1++){
        fprintf(Init, "%f \t %f \t %f \n", time_stepper, gsl_vector_get(r, count1), gsl_vector_get(Tk, count1));
    }

    gsl_permutation * p = gsl_permutation_alloc(N);
    gsl_linalg_LU_decomp(A, p, &s);     //sets up system to be solved in solver

    for (count1 = 0; count1 < M; count1++){
        time_stepper = time_stepper + delta_time;
        solver(time_stepper, A, p, Tk, Sk, b);
        if (count1 == oneyear || count1 == tenyear || count1 == fiftyyear || count1 == hundyear){
            for (count2 = 0; count2 < N; count2++){
                fprintf(f, "%f \t %f \t %f \n", time_stepper, gsl_vector_get(r, count2), gsl_vector_get(Tk, count2));
            }
            fprintf(f, "\n");
        }
    }
    gsl_permutation_free(p);
    gsl_matrix_free(A);
    gsl_vector_free(Tk);
    gsl_vector_free(Sk);
    gsl_vector_free(b);
    gsl_vector_free(r);
    fclose(f);
    return 0;
}

// Initialises the matrix and vectors used in main to pass to solver function.
static void setup_matrix_vec(gsl_matrix* A_0, gsl_vector* T_k, gsl_vector* S_k, gsl_vector* r_0, gsl_vector* b_0, double s, double dr){
    int count1, count2;
    double radius = 0.0;

    for (count1 = 0; count1 < N; count1++){
        gsl_vector_set(r_0, count1, radius);
        gsl_vector_set(T_k, count1, T_BOUND);
        radius = radius + dr;
        if (count1 == N - 1){
            gsl_vector_set(b_0, count1, -(-s - s / (2 * N)) * T_BOUND);
        }
        else{
            gsl_vector_set(b_0, count1, 0);
        }
        if (count1 < a){
            gsl_vector_set(S_k, count1, 1.0);
        }
        else{
            gsl_vector_set(S_k, count1, 0);
        }
        for (count2 = 0; count2 < N; count2++){
            if (count1 == count2){
                gsl_matrix_set(A_0, count1, count2, 1 + 2*s);
            }
            else if (count1 == count2 - 1){
                gsl_matrix_set(A_0, count1, count2, -s - s / (2 * (count1 + 1)));
            }
            else if (count1 == count2 + 1){
                gsl_matrix_set(A_0, count1, count2, -s + s / (2 * (count1 + 1)));
            }
            else{
                gsl_matrix_set(A_0, count1, count2, 0);
            }
        }
    }
    gsl_matrix_set(A_0, 0, 0, 1 + s + (s / 2));
}

// Solves for the new temperature vector at the later time. Replaces old temperature vector with new values.
// Creates one vector of values in "hold" to be solved along with the matrix in LU decomposition to obtain the new temperature vector.
static void solver(double time_step, gsl_matrix* A_0, gsl_permutation* p_0, gsl_vector* T_k, gsl_vector* S_k, gsl_vector* b_0){
    int counter;
    double  value_holder = 0.0;
    double source_value = 0.0;
    gsl_vector * hold = gsl_vector_alloc(N);
    gsl_vector * T_k_plusone = gsl_vector_alloc(N);

    source_value = (T_ROD * exp(-(time_step / TAU_0))) / pow(a, 2);
    for (counter = 0; counter < N; counter++){
        value_holder = 0.0;
        value_holder = gsl_vector_get(T_k, counter) + (KAPPA * (TIME_LENGTH_SOL/M) * source_value * gsl_vector_get(S_k, counter)) + gsl_vector_get(b_0, counter);
        gsl_vector_set(hold, counter, value_holder);
    }
    gsl_linalg_LU_solve(A_0, p_0, hold, T_k_plusone);

    for(counter = 0; counter < N; counter++){
        gsl_vector_set(T_k, counter, gsl_vector_get(T_k_plusone, counter));
    }
    gsl_vector_free(T_k_plusone);
    gsl_vector_free(hold);
}
