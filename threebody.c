/*
PROJECT 1: N-BODY PROBLEM
VERSION 2.1
WEI LU
30/11/2018

This program calculates the orbits of a 3-body system. The initial conditions used for this program must be read from a text file.
The values provided here are for the solutions for the Earth, Sun and Moon. As well as this, the program solves a system of three identical particles.
Command line arguments are taken to choose between the two stable orbit data outputs. To choose the Earth, Sun and Moon, "-a" is needed and for the figure
eight orbit, "-b" is needed.

This work includes the contribution from Dr Williams' "ReadOrbits.c", which is used to read initial conditions of the bodies from a file, where upon the
information is used for further calculations. (lines 100 - 120)

The algorithm used to solve the 3-body problem is the Verlet method, which is a numerical solution to integrating the differential equation of motion.
It is called iteratively in a loop until a maximum number of evaluations is reached.
This loop provides the logic for the constant updating of position, acceleration, velocity and force vectors required for the simulation.
The convention with the 1D arrays for this work is that the elements 0, 1, 2 correspond to X, Y and Z directions.
Energy considerations are also met here and used to confirm the validity of the program. The total energy is calculated and printed. The potential energy
is the gravitational potential energy. The kinetic energy is for the translational velocity of the body in 3 dimensions.

X and Y components of position will be written into a text file for each of the three bodies provided. GNUplot is then used to plot the graphs manually.
The logic and calculation loop of this program will work for n-bodies as long as there are relevant initial conditions provided from a file.
One must change the "MAX_BODIES" constant to accommodate as well as change the number of text files needed to write out to.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MAX_FILE_LINE_SIZE 250
#define ITEMS_PER_LINE 8
#define MAX_BODIES 3
#define MAX_NAME_SIZE 32
#define DIMENSIONS 3    // Spatial Dimensions


typedef enum coords {X, Y, Z, N_COORDS} Coords;

typedef struct body{
    char name[MAX_NAME_SIZE];
    double mass;
    double r[N_COORDS];
    double v[N_COORDS];
    double a[N_COORDS];

} Body;

static void error_message();
static void do_energy(int BodyN, Body *bodies, double G);
static double kinetic_energy_calc(Body body);
static double potential_energy_calc(Body body1, Body body2, double G);
static void * xmalloc(size_t bytes);
static double * v_acc_calc(Body body1, Body body2, double G);
static double * v_position(Body body, double DT);
static double * v_velocity(Body body, double DT);


int main(int argc, char *argv[])
{
    double G, DT, MAX_ITR;
    double *new_positions, *new_acceleration, *new_velocity;
    double time = 0.0;
    int body_num, count;
    char line[MAX_FILE_LINE_SIZE];
    char nameBuf[MAX_FILE_LINE_SIZE];
    FILE *input;
    FILE *Body1_info, *Body2_info, *Body3_info;

    if (argc == 2){
        if (strcmp(argv[1], "-a") == 0){
            G = 6.67e-11;
            DT = 60;
            MAX_ITR = 1051200;
            Body1_info = fopen("EarthData.txt", "w");
            Body2_info = fopen("SunData.txt", "w");
            Body3_info = fopen("MoonData.txt", "w");
            input = fopen("selfOrbitData.txt", "r");
        }
        else if (strcmp(argv[1], "-b") == 0){
            G = 1;
            DT = 0.0005;
            MAX_ITR = 25000;
            Body1_info = fopen("mass1data.txt", "w");
            Body2_info = fopen("mass2data.txt", "w");
            Body3_info = fopen("mass3data.txt", "w");
            input = fopen("figureeightdata.txt", "r");
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

    if (!input) {
        fprintf(stderr, "Error: Could not open orbit data file.\n");
        exit(1);
    }
    Body *bodies = malloc(MAX_BODIES * sizeof(Body));
    int bodyN = 0;
    while ( bodyN < MAX_BODIES && fgets(line, MAX_FILE_LINE_SIZE, input) ) {
        if (line[0] != '#') {
            int nFound = sscanf(line,"%s %lg %lg %lg %lg %lg %lg %lg",
                        nameBuf, &bodies[bodyN].mass,
                        &bodies[bodyN].r[X], &bodies[bodyN].r[Y], &bodies[bodyN].r[Z],
                        &bodies[bodyN].v[X], &bodies[bodyN].v[Y], &bodies[bodyN].v[Z] );
            if (nFound == ITEMS_PER_LINE) {
                strncpy(bodies[bodyN].name,nameBuf,MAX_NAME_SIZE);
                bodyN++;
            }
            else {
                fprintf(stderr, "Unknown format: %s\n",line);
            }
        }
    }
    printf("Initial energy of system (Joules): ");
    do_energy(bodyN, bodies, G);

    do{
        for (body_num = 0; body_num < bodyN; body_num++){
            bodies[body_num].a[X] = bodies[body_num].a[Y] = bodies[body_num].a[Z] = 0.0;
            for (count = 0; count < bodyN; count++){
                if (body_num == count){
                    continue;
                }
                else{
                        new_acceleration = v_acc_calc(bodies[body_num], bodies[count], G);
                        bodies[body_num].a[X] -= new_acceleration[0];
                        bodies[body_num].a[Y] -= new_acceleration[1];
                        bodies[body_num].a[Z] -= new_acceleration[2];
                }
            }
        }
        free(new_acceleration);

        for (body_num = 0; body_num < bodyN; body_num++){
            new_positions = v_position(bodies[body_num], DT);
            bodies[body_num].r[X] = new_positions[0];
            bodies[body_num].r[Y] = new_positions[1];
            bodies[body_num].r[Z] = new_positions[2];
        }

        for (body_num = 0; body_num < bodyN; body_num++){
            bodies[body_num].a[X] = bodies[body_num].a[Y] = bodies[body_num].a[Z] = 0.0;
            for (count = 0; count < bodyN; count++){
                if (body_num == count){
                    continue;
                }
                else{
                    new_acceleration = v_acc_calc(bodies[body_num], bodies[count], G);
                    bodies[body_num].a[X] -= new_acceleration[0];
                    bodies[body_num].a[Y] -= new_acceleration[1];
                    bodies[body_num].a[Z] -= new_acceleration[2];
                }
            }
        }
        for (body_num = 0; body_num < bodyN; body_num++){
            new_velocity = v_velocity(bodies[body_num], DT);
            bodies[body_num].v[X] = new_velocity[0];
            bodies[body_num].v[Y] = new_velocity[1];
            bodies[body_num].v[Z] = new_velocity[2];
        }

        fprintf(Body1_info, "%f \t %f \n", bodies[0].r[X], bodies[0].r[Y]);
        fprintf(Body2_info, "%f \t %f \n", bodies[1].r[X], bodies[1].r[Y]);
        fprintf(Body3_info, "%f \t %f \n", bodies[2].r[X], bodies[2].r[Y]);

        free(new_velocity);
        free(new_acceleration);
        free(new_positions);
        time += DT;
    } while(time < (DT * MAX_ITR));
    printf("Final energy of system (Joules): ");
    do_energy(bodyN, bodies, G);

    free(bodies);

    fclose(Body1_info);
    fclose(Body2_info);
    fclose(Body3_info);

    fclose(input);
    return 0;
}

static void error_message(){
printf("Incorrect argument. \n");
printf("To choose Earth, Sun and Moon orbit, enter: ' -a' as an argument. \nTo choose figure eight orbit, enter: ' -b' as an argument. \n");
}

/*Energy calculations for the total energy, which calls kinetic_energy_calc and potential_energy_calc.*/
static void do_energy(int bodyN, Body *bodies, double G){
    int body_num, count;
    double k_energy, p_energy = 0.0;

    for (body_num = 0; body_num < bodyN; body_num++){
		k_energy += kinetic_energy_calc(bodies[body_num]);
        for (count = 0; count < bodyN; count++){
            if (body_num == count){
                continue;
            }
            else{
                p_energy += potential_energy_calc(bodies[body_num], bodies[count], G);
            }
        }
    }
    printf("%le \n", k_energy - p_energy);
}

static double kinetic_energy_calc(Body body){
    double energy;
    double velocity_mag;

    velocity_mag = sqrt(pow(body.v[X], 2) + pow(body.v[Y], 2) + pow(body.v[Z], 2));
    energy = 0.5 * body.mass * pow(velocity_mag, 2);
    return energy;
}

static double potential_energy_calc(Body body1, Body body2, double G){
    double x, y, z, dist, energy;

    x = body1.r[X] - body2.r[X];
    y = body1.r[Y] - body2.r[Y];
    z = body1.r[Z] - body2.r[Z];
    dist = sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));
    energy = -(G * body1.mass * body2.mass) / dist;
    return energy;
}

static void * xmalloc(size_t bytes){
    void *retVal = malloc(bytes);
    if (retVal) {
        return retVal;
    }
    printf("ERROR: Memory not allocated.");
    exit(99);
}

/* VERLET ALGORTIHM FUNCTIONS */
static double * v_acc_calc(Body body1, Body body2, double G){	//calculates the acceleration at the given point provided by values in the struct.
    double x, y, z, dist_cube, f_scalar;
    double *acceleration = xmalloc(DIMENSIONS * sizeof(double));

    x = body1.r[X] - body2.r[X];
    y = body1.r[Y] - body2.r[Y];
    z = body1.r[Z] - body2.r[Z];
    dist_cube = pow(sqrt(pow(x,2) + pow(y,2) + pow(z,2)), 3);
    f_scalar = (G * body1.mass * body2.mass) / dist_cube;
    acceleration[0] = (f_scalar * x) / body1.mass;
    acceleration[1] = (f_scalar * y) / body1.mass;
    acceleration[2] = (f_scalar * z) / body1.mass;
    return acceleration;
}

static double * v_position(Body body, double DT){	//calculates new position of the body.
    double *pos = xmalloc(DIMENSIONS * sizeof(double));

    pos[0] = body.r[X] + (body.v[X] * DT) + (body.a[X] * pow(DT, 2));
    pos[1] = body.r[Y] + (body.v[Y] * DT) + (body.a[Y] * pow(DT, 2));
    pos[2] = body.r[Z] + (body.v[Z] * DT) + (body.a[Z] * pow(DT, 2));
    return pos;
}

static double * v_velocity(Body body, double DT){	//calculates new velocity of the body.
    double *vel = xmalloc(DIMENSIONS * sizeof(double));

    vel[0] = body.v[X] + (body.a[0] * DT);
    vel[1] = body.v[Y] + (body.a[1] * DT);
    vel[2] = body.v[Z] + (body.a[2] * DT);
    return vel;
}
