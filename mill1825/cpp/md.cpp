#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

using namespace std;

int n = 4;
int N = 4*pow(n, 3); 
double rho = 1.0;   
double T = 5.0;    
double N_tau = 1;
double nu = 100.0;

double rc2 = 1e20;
double ecut;
double L;

void initialize();
double Temperature();
void BerendsenThermostat(double dt);
void AndersenThermostat(double dt);

double **r; 
double **v; 
double **a; 

gsl_rng * rr = gsl_rng_alloc(gsl_rng_mt19937);
unsigned long int Seed = 23410981;

void initialize() {
    
    gsl_rng_set(rr,Seed);    

    double rr3 = 1.0/(rc2*rc2*rc2);
    ecut = 4*(rr3*rr3*rr3*rr3-rr3*rr3);

    cout << "ecut = " << ecut << endl;

    r = new double* [N];
    v = new double* [N];
    a = new double* [N];
    for (int i = 0; i < N; i++) {
        r[i] = new double [3];
        v[i] = new double [3];
        a[i] = new double [3];
    }

    L = pow(N / rho, 1.0/3);
    int M = 1;
    while (4 * M * M * M < N)
        ++M;
    double a = L / M; 

    cout << "Particles = " << N << endl;
    cout << "Boxes =  " << M << endl;
    
    double xCell[4] = {0.25,0.75, 0.75, 0.25};
    double yCell[4] = {0.25,0.75, 0.25, 0.75};
    double zCell[4] = {0.25,0.25, 0.75, 0.75};

    int n = 0;
    for (int x = 0; x < M; x++)
        for (int y = 0; y < M; y++)
            for (int z = 0; z < M; z++)
                for (int k = 0; k < 4; k++)
                    if (n < N) {
                        r[n][0] = (x + xCell[k]) * a;
                        r[n][1] = (y + yCell[k]) * a;
                        r[n][2] = (z + zCell[k]) * a;
                        ++n;
                    }

    for (int n = 0; n < N; n++)
        for (int i = 0; i < 3; i++)
            v[n][i] = gsl_ran_gaussian(rr,1.0);

    double vCM[3] = {0, 0, 0};
    for (int n = 0; n < N; n++)
        for (int i = 0; i < 3; i++)
            vCM[i] += v[n][i];
    for (int i = 0; i < 3; i++)
        vCM[i] /= N;
    for (int n = 0; n < N; n++)
        for (int i = 0; i < 3; i++)
            v[n][i] -= vCM[i];

    double vSqdSum = 0;
    for (int n = 0; n < N; n++)
        for (int i = 0; i < 3; i++)
            vSqdSum += v[n][i] * v[n][i];
    double lambda = sqrt( 3 * (N-1) * T / vSqdSum );
    for (int n = 0; n < N; n++)
        for (int i = 0; i < 3; i++)
            v[n][i] *= lambda;
}

void computeAccelerations() {

    for (int i = 0; i < N; i++) 
        for (int k = 0; k < 3; k++) 
            a[i][k] = 0; 

    for (int i = 0; i < N-1; i++) 
        for (int j = i+1; j < N; j++) { 
            double rij[3]; 
            double rSqd = 0; 
            for (int k = 0; k < 3; k++) {
                rij[k] = r[i][k] - r[j][k];
                if (abs(rij[k]) > 0.5 * L) {
                    if (rij[k] > 0)
                        rij[k] -= L;
                    else
                        rij[k] += L;
                }
                rSqd += rij[k] * rij[k];
            }
            double f = 24 * (2 * pow(rSqd, -7) - pow(rSqd, -4));
            for (int k = 0; k < 3; k++) {
                a[i][k] += rij[k] * f;
                a[j][k] -= rij[k] * f;
            }
        }
}

void BerendsenThermostat(double dt) {
    double tau = N_tau*dt;
    double T_i = Temperature();
    double lambda = sqrt( 1 + dt/tau*(T_i/T-1) );
    for (int n = 0; n < N; n++)
        for (int i = 0; i < 3; i++)
            v[n][i] *= lambda;
}

void AndersenThermostat(double dt) {
    for (int i = 0; i < N; i++) { 
        if ( gsl_rng_uniform(rr) <= 1-exp(-nu*dt)) {
            for (int k = 0; k < 3; k++)
                v[i][k] = sqrt(T) * gsl_ran_gaussian(rr,1.0);
        }
    }

}

void velocityVerlet(double dt) {
    computeAccelerations();
    for (int i = 0; i < N; i++)
        for (int k = 0; k < 3; k++) {
            r[i][k] += v[i][k] * dt + 0.5 * a[i][k] * dt * dt;
            if (r[i][k] < 0)
                r[i][k] += L;
            if (r[i][k] >= L)
                r[i][k] -= L;
            v[i][k] += 0.5 * a[i][k] * dt;
        }
    computeAccelerations();
    for (int i = 0; i < N; i++)
        for (int k = 0; k < 3; k++)
            v[i][k] += 0.5 * a[i][k] * dt;
}

double Temperature() {
    double sum = 0;
    for (int i = 0; i < N; i++)
        for (int k = 0; k < 3; k++)
            sum += v[i][k] * v[i][k];
    return sum / (3 * (N - 1));
}

double KineticEnergy() {
    double KE = 0;
    for (int i = 0; i < N; i++)
        for (int k = 0; k < 3; k++)
            KE += v[i][k] * v[i][k]; 
    return 0.5*KE;
}

double PotentialEnergy() {
    double PE = 0, P = 0;
    for (int i = 0; i < N; i++) {
        for (int j = i+1; j < N; j++) { 
            double rij[3]; 
            double rSqd = 0; 
            for (int k = 0; k < 3; k++) {
                rij[k] = r[i][k] - r[j][k];
                if (abs(rij[k]) > 0.5 * L) {
                    if (rij[k] > 0)
                        rij[k] -= L;
                    else
                        rij[k] += L;
                }
                rSqd += rij[k] * rij[k];
            }
            if (rSqd < rc2) 
                PE += 4 * (pow(rSqd, -6) - pow(rSqd, -3)) - ecut;
                P += 24 * (2 * pow(rSqd, -7) - pow(rSqd, -4));
        }
        
    }
    return PE;
}

int main() {
    initialize();
    double dt = 1e-4;
    int Nsteps = 20000;
    int Nrescale = Nsteps/200;
    int AndersenEnd = Nsteps/2;

    ofstream fileE("./results/E.dat"), fileT("./results/T_nu100.dat");
    for (int i = 0; i < Nsteps; i++) {
        velocityVerlet(dt);
        fileE << i << " " << KineticEnergy() << " " << PotentialEnergy() << endl;
        fileT << i << " " << Temperature() << endl;
        if (i <= AndersenEnd) { AndersenThermostat(dt); }
    }
    fileE.close();
    fileT.close();
}

