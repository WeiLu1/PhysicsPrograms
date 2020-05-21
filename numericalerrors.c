/*

This program outputs the value calculated from two functions, that are mathematically
equal, over a range of values spanning between minus 2pi and 2pi. The aim is to show that the outputs from the functions
are equal to an extent that is dominated by the precision of the input data type (double).
Estimated values of these outputs are also found using inputs that include the relative error
of the C double type. The forward error is found by calculating the difference between the true
output using a more precise long double type and the estimated output.
The backward error is found by dividing the relative forward error, for a specific input value, by the
Condition number.
Condition number here is the derivative of the function multiplied by the input and
divided by the function evaluated at the input. The first derivative of the function produces a
division by 0 error.
The derivative of the condition number with respect to the input value is also found for further analysis,
but need not be displayed unless intended.

*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define REL_ERR 0.5e-16
#define TWO_PI 2 * M_PI

static double func_cos(double x);
static double func_sin(double x);
static long double cosl_f_err(long double x);
static long double sinl_f_err(long double x);
static long double cosl_b_err(long double x, long double *ptr, int num);
static long double cosl_condition_deriv(long double x);
static long double sinl_b_err(long double x, long double *ptr, int num);


int main()
{
    double interval = TWO_PI / 12;
    double counter = -TWO_PI;
    int size = (2 * TWO_PI) / interval;
    double *cos_ans_array = malloc(size * sizeof(double)), *sin_ans_array= malloc(size * sizeof(double)), *difference = malloc((size * sizeof(double)));
    int i;
    long double input = 0.0;
    long double *cos_f_err = malloc(11 * sizeof(long double)), *sin_f_err = malloc(11 * sizeof(long double));
    long double cos_b_err, sin_b_err, cos_condition_derivative;

    if(cos_ans_array == NULL || sin_ans_array == NULL || difference == NULL || cos_f_err == NULL || sin_f_err == NULL)
    {
        printf("Error: memory not allocated.");
        exit(99);
    }

    printf("Input: \t\t f(x) output: \t g(x) output: \t difference: \n");
    for (i = 0; i <= size; i++){	// "<=" means you are returning one more element to a predefined array size that is 1 less. You are overrunning array boundaries. 
        cos_ans_array[i] = func_cos(counter);
        sin_ans_array[i] = func_sin(counter);
        difference[i] = fabs(sin_ans_array[i] - cos_ans_array[i]);
        printf("%f \t %f \t %f \t %e \n", counter, cos_ans_array[i], sin_ans_array[i], difference[i]);
        counter += interval;
    }

    printf("\n\n");
    printf("Input:\t\t f(x) forward error: \tg(x) forward error: \n");
    for (i = 0; i <= 10; i++){
        cos_f_err[i] = cosl_f_err(input);
        sin_f_err[i] = sinl_f_err(input);
        printf("%.1Le \t %Le\t\t\t%Le \n", input, cos_f_err[i], sin_f_err[i]);
        input += 1e-8;
    }

    input = 0.0;
    printf("\n\n");
    printf("Input: \t\t f(x) backward error: \tg(x) backward error: \tf(x) condition derivative:\n");
    for (i = 0; i <= 10; i++){
        cos_b_err = cosl_b_err(input, cos_f_err, i);
        sin_b_err = sinl_b_err(input, sin_f_err, i);
        cos_condition_derivative = cosl_condition_deriv(input);
        printf("%.1Le \t %Le \t\t\t%Le \t\t\t%Le\n", input, cos_b_err, sin_b_err, cos_condition_derivative);
        input += 1e-8;
    }

    free(cos_ans_array);
    free(sin_ans_array);
    free(difference);
    free(cos_f_err);
    free(sin_f_err);

    return 0;
}


static double func_cos(double x){
    double result;
    result = (1 - cos(x)) / pow(x, 2);
    return result;
}


static double func_sin(double x){
    double result;
    result = 2 * pow((sin(0.5 * x) / x), 2);
    return result;
}


static long double cosl_f_err(long double x){
    long double result_true, f_err;
    double result_estimate, x_hat;
    x_hat = (double) x * (1 + REL_ERR);
    result_estimate = func_cos(x_hat);
    result_true = (1 - cosl(x)) / powl(x, 2);
    f_err = fabsl(result_true - result_estimate);
    return f_err;
}


static long double sinl_f_err(long double x){
    long double result_true, f_err;
    double result_estimate, x_hat;
    x_hat = (double) x * (1 + REL_ERR);
    result_estimate = func_sin(x_hat);
    result_true = 2 * powl((sinl(0.5 * x) / x), 2);
    f_err = fabsl(result_true - result_estimate);
    return f_err;
}


static long double cosl_b_err(long double x, long double *ptr, int num){
    long double cosl_derivative, cosl_condition, result_true, f_err, b_err;
    f_err = *(ptr + num);
    result_true = (1 - cosl(x)) / powl(x, 2);
    cosl_derivative = (x * sinl(x) + 2 * cosl(x) - 2) / powl(x, 3);
    cosl_condition = fabsl((cosl_derivative * x) / result_true);
    b_err = (f_err / fabsl(result_true)) / fabsl(cosl_condition);
    return b_err;
}

static long double cosl_condition_deriv(long double x){
    long double cosl_condition_derivative;
    cosl_condition_derivative = (-powl(x,2)*powl(sin(x),2)+sin(x)-x*powl(cos(x),2)+(x-sin(x))*cos(x)) / (1-powl(cos(x),2));
	return cosl_condition_derivative;
}


static long double sinl_b_err(long double x, long double *ptr, int num){
    long double sinl_derivative, sinl_condition, result_true, f_err, b_err;
    f_err = *(ptr + num);
    result_true = fabsl(2 * powl((sinl(0.5 * x) / x), 2));
    sinl_derivative = (sinl(0.5*x) * ((2 * x * cosl(0.5*x)) - (4 * sinl(0.5*x)))) / powl(x, 3);
    sinl_condition = fabsl((sinl_derivative * x) / result_true);
    b_err = (f_err / result_true) / sinl_condition;
    return b_err;
}


