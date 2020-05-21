/*
PR2 - FINAL PROJECT - PROJECT OPTION C: ISING MODEL
VERSION 7.2
WEI LU
08/03/2019

This program simulates the 2D Ising Model using the Metropolis-Hastings Algorithm, where the time of the operating system is used to seed the random number generators. 
The outputs of the program are quantities of interest: magnetisation, heat capacity and magnetic susceptibility. 
These quantities are written to text file "Ising_data,txt", in the order: temperature, average magnetisation per spin, heat capacity, susceptibility and average energy per spin. 
While identical data is output onto the command line. 

To change the size of the 2D square system, one must change the variable "D", since the multi-dimensional array is a DxD sized.
The amount of iterations that the program is executed for depends on the macros "EQUIL" (the number of iterations needed to let the system reach equilibrium) and "MCS" (Monte-Carlo
Steps, the number of iterations of the Metropolis-Hastings algorithm after the system has reached equilibrium). 
The system calculates quantities of interest for a minimum temperature up to a maximum temperature, in intervals of a step defined by a macro, however special care must be taken in 
considering the temperature intervals as it affects the size of the memory allocated arrays to store the quantities of interest. Small steps such as 0.25 are preferred to avoid memory
leaks and undefined behaviour, the start and end temperature must also be considered in this way to avoid over running array boundaries. 

The units of temperature for this program are J/k_B, where J is the interaction strength and k_B is the Boltzmann factor. All subsequent calculations are made with these units in 
mind.

Periodic boundary conditions have been introduced in the energy calculations. Therefore, when considering the spins in the lattice, when the program encounters an edge or any first or
last row or column, it will take the spin at the opposite edge as its nearest neighbour. 

To further investigate the quantities of interest, the "investigator" function writes to a file "Investigation_data.txt", the difference between the temperature and the critical 
temperature found by Onsager to the power of the exponents specifically for the three quantiites of interest, and the corresponding values of magnetisation, specific heat and 
susceptibility. It also writes the values for the natural logarithm of the heat capacity as well as the logarithm for the size value D (in steps up to D). 
To ensure that the number of log values of D and C_v are kept the same, there is an appropriate stepping variable used to calculate the neccessary values of ln(D) from 0 to D. 

GNUplot will be used to manually plot the graphs for each relation from the their data files.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#define RANDOM_MAX 0x7FFFFFFF
#define S(x,y) spin[(D + (x)) % D][(D + (y)) % D]
#define MAX_TEMP 5.0
#define START_TEMP 0.5
#define TEMP_STEP 0.25
#define UP 1.0
#define DOWN -1.0
#define D 50				// Lattice size
#define N pow(D, 2.0)		
#define EQUIL (6 * D * D)	// Number of Monte-Carlo sweeps until the system reaches equilibrium
#define MCS (8 * D * D)		// Number of Monte-Carlo steps after the system has reached equilibrium

// Critical exponents for phase transtion analysis
#define ALPHA 0.0
#define BETA (1.0/8.0)
#define GAMMA (-7.0/4.0)
#define TC 2.269			// Analytical solution for the 2D critical temperature

static void * xmalloc(size_t bytes);
static double rng(void);
static int rand_pos(void);
static void initialise_grid(double spin[][D], double * Energy);
static void metropolis_move(double spin[][D], double * Energy, double temp);
static double delta_E(double spin[][D], int x, int y);
static double magnetisation(double spin[][D]);
static double heat_cap(double E_av, double E_squared_av, double temp);
static double susceptibility(double M_av, double M_squared_av, double temp);
static void output(double T, double M, double C, double X);
static void investigator(double * Cv, double * Magn, double * X, int stepweight);

int main()
{	
	double spin[D][D];
	double Temp;
	int sweep;
	double av_E_per_site = 0.0;
	int array_weight = (MAX_TEMP / TEMP_STEP);				// This variable determines the size of the array in which values for the qauntities of interest are stored, ensure the MACRO defintions of temperature variables are set with this in mind
	double * C_v = xmalloc(array_weight * sizeof(double));	// Considerations have been made to ensure there is enough memory allocated
	double * chi = xmalloc(array_weight * sizeof(double));
	double * av_M_lattice = xmalloc(array_weight * sizeof(double));
	int counter = 0;
	double normalize =  1.0 / ((double)MCS);				// Averaging factor
	
	if (TEMP_STEP > 0.50){
		printf("Change TEMP_STEP to be between 0 and 0.5.");// This ensures that there is enough memory allocated in the arrays for the quantities of interest
		exit(99);
	}
	
	FILE * results = fopen("Ising_data.txt", "w");
	if (results == NULL){
		printf("Failed to open file.");
		exit(99);
	}
	
    srand(time(NULL));
	clock_t t;
    t = clock();
	
	for (Temp = START_TEMP; Temp <= MAX_TEMP; Temp += TEMP_STEP){
		double E_out = 0.0, E = 0.0, E_sqrd = 0.0, av_E = 0.0, av_E_sqrd = 0.0;
		double M_out = 0.0, M = 0.0, M_sqrd = 0.0, av_M = 0.0, av_M_sqrd = 0.0;
		
		initialise_grid(spin, &E);
		for (sweep = 0; sweep < EQUIL; sweep++){
			metropolis_move(spin, &E, Temp);
		}
		
		for (sweep = 0; sweep < MCS; sweep++){
			metropolis_move(spin, &E, Temp);
			M = magnetisation(spin);
			
			E_out += E;
			E_sqrd += pow(E, 2.0);
			M_out += M;
			M_sqrd += pow(M, 2.0);
		}
		av_M_lattice[counter] = (fabs(M) / (double)N);		// Average magnetisation of the lattice per spin at the end of the sweeps
		av_E = E_out * normalize;							// Average energy
		av_E_per_site = (av_E * normalize) / (double)N;		// Average energy per site of the lattice
		av_E_sqrd = E_sqrd * normalize;						// Average of the energy squared
		av_M = M_out * normalize;							// Average magnetisation
		av_M_sqrd = M_sqrd * normalize;						// Average of the magnetisation squared
		
		C_v[counter] = heat_cap(av_E, av_E_sqrd, Temp) / (double)N;
		chi[counter] = susceptibility(av_M, av_M_sqrd, Temp) / (double)N;
		
		fprintf(results, "%f \t %f \t %f \t %f \t %f \n", Temp, av_M_lattice[counter], C_v[counter], chi[counter], av_E_per_site);
		output(Temp, av_M_lattice[counter], C_v[counter], chi[counter]);
		counter += 1;
	}
	t = clock() - t;
    double time_took = ((double)t)/CLOCKS_PER_SEC;
    printf("Execution time: %f seconds. \n", time_took);
	
	investigator(C_v, av_M_lattice, chi, counter);
	
	fclose(results);
	free(av_M_lattice);
	free(chi);
	free(C_v);
    return 0;
}

static void * xmalloc(size_t bytes){
    void *retVal = malloc(bytes);
    if (retVal) {
        return retVal;
    }
    printf("ERROR: Memory not allocated.");
    exit(99);
}

// Random number generator between 0 and 1
static double rng(void){
    return((rand() / (double)RANDOM_MAX));
}

// Used in finding random positions in the 2D grid
static int rand_pos(void){
	return(rand() % (int)D);
}

// Sets up the grid with random spins and calculates initial energy of this configuration subject to periodic boundary conditions
static void initialise_grid(double spin[][D], double * Energy){
    int row, col, south, east;
    double rand_num = 0.0;
	double sum = 0.0;
	
    for (row = 0; row < D; row++){
        for (col = 0; col < D; col++){
            rand_num = rng();
            if (rand_num < 0.5){
                S(row, col) = UP;
            }
            else{
                S(row, col) = DOWN;
            }
        }
    }
	for (row = 0; row < D; row++){
		if (row == D - 1){
			south = 0;
		}
		else{
			south = row + 1;
		}
		
		for (col = 0; col < D; col++){
			if (col == D - 1){
				east = 0;
			}
			else{
				east = col + 1;
			}
			sum = (S(south, col) + S(row, east));
			*Energy -= (S(row, col) * sum);
		}
	}
}

// The Metropolis-Hastings Algorithm
// Called iteratively from main from equilibrium and MCS loops, to randomly flip spins
static void metropolis_move(double spin[][D], double * Energy, double temp){
	int spin_step, a, b;
	double dE, rand_num, w;
	
	for (spin_step = 0; spin_step < N; spin_step++){
		a = rand_pos();
		b = rand_pos();
		dE = delta_E(spin, a, b);
		rand_num = rng();
		w = exp(-dE / temp);
			
		if (dE <= 0.0){
			S(a, b) = -S(a, b);
			*Energy += -dE;
		}
		else if (rand_num <= w){
			S(a, b) = -S(a, b);
			*Energy += -dE;
		}
		else{
			S(a, b) = S(a, b);
		}
	}
}

// Calculation of the change in energy if the spin were to be flipped, subject to periodic boundary conditions
// The return value of this function dictates whether or not the spin is flipped in the "metropolis_move" function
static double delta_E(double spin[][D], int row, int col){
	int north = 0;
	int east = 0;
	int south = 0;
	int west = 0;
	double dE = 0.0;
	
	if (row == 0){
		north = S(D - 1, col);
	}
	else if (row == D - 1){
		south = S(0, col);
	}
	else{
		north = S(row - 1, col);
		south = S(row + 1, col);
	}
	
	if (col == 0){
		west = S(row, D - 1);
	}
	else if ( col == D - 1){
		east = S(row, 0);
	}
	else{
		west = S(row, col - 1);
		east = S(row, col + 1);
	}
	dE = S(row, col) * (2.0 * (north + east + south + west));
	
	return dE;
}

// Calculation of the total spin of the system
static double magnetisation(double spin[][D]){
	int row, col; 
	double spin_sum = 0.0;
	double total_mag = 0.0;
	
	for (row = 0; row < D; row++){
		for (col = 0; col < D; col++){
			spin_sum += S(row, col);
		}
	}
	total_mag = spin_sum;
	return total_mag;
}

// Calculation of the heat capacity
static double heat_cap(double E_av, double E_squared_av, double temp){
	double C_v = 0.0;
	
	C_v = (E_squared_av - pow(E_av, 2.0)) / pow(temp, 2.0);
	return C_v;
}

// Calculation of the magnetic susceptibility 
static double susceptibility(double M_av, double M_squared_av, double temp){
	double chi = 0.0;
	
	chi = (M_squared_av - pow(M_av, 2.0)) / temp;
	return chi;
}

static void output(double T, double M, double C, double X){
	printf("Temperature: %f \n", T);
	printf("Magnetisation: %.15f \n", M);
	printf("Heat Capacity: %f \n", C);
	printf("Susceptibility: %f \n", X);
	printf("\n");
}

// Further analysis of the quanities of interest, values to be written to file and plotted in GNUplot 
static void investigator(double * Cv, double * Magn, double * X, int stepweight){
	int counter = 0;
	double temp;
	double Cv_temp_check = 0.0;
	double M_temp_check = 0.0;
	double X_temp_check = 0.0;
	double D_step = ((double)D / stepweight);	// Used to step through lattize size variable D equally
	double D_stepper = 0.0;
	double log_Cv = 0.0;
	double log_D = 0.0;
	
	FILE * investigate = fopen("Investigation_data.txt", "w");
	if (investigate == NULL){
		printf("Failed to open file.");
		exit(99);
	}
	
	for (temp = START_TEMP; temp <= MAX_TEMP; temp += TEMP_STEP){
		if (D_stepper == 0.0){
			D_stepper += D_step;
			continue;							// This will avoid taking the logarithm of 0
		}
		else{
			log_Cv = log(Cv[counter]);
			log_D = log(D_stepper);
		}
		Cv_temp_check = pow(fabs(temp - TC), ALPHA);
		M_temp_check = 	pow(fabs(temp - TC), BETA);
		X_temp_check = pow(fabs(temp - TC), GAMMA);
		
		fprintf(investigate, "%f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \n", Cv_temp_check, Cv[counter], M_temp_check, Magn[counter], X_temp_check, X[counter], log_D, log_Cv);
		
		counter += 1;
		D_stepper += D_step;
	}
	
	fclose(investigate);
}

