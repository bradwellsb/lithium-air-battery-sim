// main.cpp
// This program simulates the electrochemical behavior of a lithium-air battery system.
// It uses numerical methods to solve differential equations describing the battery's state variables.

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <cmath>
#include <iostream>
#include <fstream>
using namespace std;

// Constants for the battery model
double const L = 0.075;       // Length of the battery electrode in meters
int const N = 50;             // Number of spatial discretization points
double h = L / N;             // Spatial step size
double tau = 86400;           // Time step in seconds (1 day)

// Physical constants
double const F = 96485;       // Faraday's constant (C/mol)
double const I = 0.1e-3;      // Applied current (A)
double const A = 1;           // Cross-sectional area of the electrode (m^2)
double const eps0 = 0.75;     // Initial porosity
double const reactionBeta = 0.5; // Transfer coefficient (changed name from beta to avoid ambiguity with std)
double const r0 = 2e-6;       // Initial radius of lithium peroxide particles (m)
double const rho_li2o2 = 2.14; // Density of Li2O2 (g/cm^3)
double const rho_c = 2.26;    // Density of carbon (g/cm^3)
double const M_li2o2 = 45.88; // Molar mass of Li2O2 (g/mol)
double const c0 = 3.26e-6;    // Initial oxygen concentration (mol/m^3)
double const k = 1.76e-10;    // Rate constant for oxygen reduction (m/s)
double const I0 = 1.3e-8;     // Exchange current density (A/m^2)
double const sigma = 1.2;     // Conductivity of the electrode (S/m)
double const kappa = 0.0124;  // Ionic conductivity of the electrolyte (S/m)
double const D_li = 1e-5;     // Diffusion coefficient of Li+ in the electrolyte (m^2/s)
double const tplus = 0.25;    // Transference number of Li+
double const D_o2 = 7e-6;     // Diffusion coefficient of O2 in the electrolyte (m^2/s)
double const Vt = 0.025;      // Thermal voltage (V)
double const kappa_d = Vt * kappa * (1 - tplus);

// Function to solve a system of linear equations A * X = rhs using the Gauss elimination technique
// A, rhs, and n are input variables
// n is the number of equations in the system
// X is the output vector
void Solve(double **A, double *rhs, double *X, int n)
{
    double **w = new double*[n];
    for (int i = 0; i < n; i++) w[i] = new double[n];
    int **index = new int*[n];
    for (int i = 0; i < n; i++) index[i] = new int[3];

    int nv  = 1;
    int irow = 0;
    int icol = 0;
    double pivot;
    
    for (int i = 0; i < n; i++){
        w[i][0] = rhs[i];
        index[i][2] = 0;
    }

    double determ = 1;
    for (int i = 0; i < n; i++){
        double big = 0;
        for (int j = 0; j < n; j++){
            if (index[j][2] != 1) 
                for (int k = 0; k < n; k++)
                {
                    if (index[k][2] > 1)
                    {
                        printf("Error at Gauss-Jordan: Matrix is singular\n");
                        exit(1);
                    }
                    if (index[k][2] < 1)
                        if (fabs(A[j][k]) > big)
                        {
                            irow = j;
                            icol = k;
                            big = fabs(A[j][k]);
                        }
                }
        }

        index[icol][2] += 1;
        index[icol][0] = irow;
        index[icol][1] = icol;

        if (irow != icol)
        {
            determ = -determ;
            for (int l = 0; l < n; l++)
            {
                double temp = A[irow][l];
                A[irow][l] = A[icol][l];
                A[icol][l] = temp;
            }
            if (nv > 0)
                for (int l = 0; l < nv; l++)
                {
                    double temp = w[irow][l];
                    w[irow][l] = w[icol][l];
                    w[icol][l] = temp;
                }
        }

        pivot = A[icol][icol];
        determ *= pivot;
        A[icol][icol] = 1;

        for (int l = 0; l < n; l++)
            A[icol][l] /= pivot;
        if (nv != 0)
            for (int l = 0; l < nv; l++)
                w[icol][l] /= pivot;
        for (int l1 = 0; l1 < n; l1++){
            if (l1 != icol)
            {
                double t = A[l1][icol];
                A[l1][icol] = 0;
                for (int l = 0; l < n; l++)
                    A[l1][l] -= A[icol][l]*t;
                if (nv != 0)
                    for (int l = 0; l < nv; l++)
                        w[l1][l] -= w[icol][l]*t;
            }
        }
    }


    for (int i = 0; i < n; i++){
        int l = n-i-1;
        if (index[l][0] != index[l][1]){
            irow = index[l][0];
            icol = index [l][1];
            for (int k = 0; k < n; k++)
            {
                double temp = A[k][irow];
                A[k][irow] = A[k][icol];
                A[k][icol] = temp;
            }
        }
    }

    for (int k = 0; k < n; k++)
        if (index[k][2] != 1){
            printf("Error at Gauss-Jordan: Matrix is singular\n");
            exit(1);
        }

    for (int i = 0; i < n; i++)
        X[i] = w[i][0];

    for (int i = 0; i < n; i++) delete [](w[i]);
    for (int i = 0; i < n; i++) delete [](index[i]);
    delete []w;
    delete []index;
}

// Functions to compute the reaction rate and its derivatives for the oxygen reduction
double RC(double phi, double phi_li, double c_o2, double eps)
{
    return (k * (2 * c_o2 * sqrt(eps * eps0)) / r0) * (exp(reactionBeta * (phi_li - phi) / Vt) - exp(-reactionBeta * (phi_li - phi) / Vt));
}

double dRC_dphi(double phi, double phi_li, double c_o2, double eps)
{
    return (k * (2 * c_o2 * sqrt(eps * eps0)) / r0) * (reactionBeta / Vt) * (-exp(reactionBeta * (phi_li - phi) / Vt) - (exp(-reactionBeta * (phi_li - phi) / Vt)));
}

double dRC_dphi_li(double phi, double phi_li, double c_o2, double eps)
{
    return (k * (2 * c_o2 * sqrt(eps * eps0)) / r0) * (reactionBeta / Vt) * ((exp(reactionBeta * (phi_li - phi) / Vt)) + (exp(-reactionBeta * (phi_li - phi) / Vt)));
}

double dRC_dc_o2(double phi, double phi_li, double c_o2, double eps)
{
    return (k * (2 * sqrt(eps * eps0)) / r0) * (exp(reactionBeta * (phi_li - phi) / Vt) - exp(-reactionBeta * (phi_li - phi) / Vt));
}

double dRC_deps(double phi, double phi_li, double c_o2, double eps)
{
    return (k * (2 * c_o2 * eps0) / (2 * r0 * sqrt(eps * eps0))) * (exp(reactionBeta * (phi_li - phi) / Vt) - exp(-reactionBeta * (phi_li - phi) / Vt));
}

double RA(double phi_li)
{
    return I0 * (exp(reactionBeta * ((-phi_li + 2.95) / Vt)) - exp(-reactionBeta * ((-phi_li + 2.95) / Vt)));
}

double dRA_dphi_li(double phi_li)
{
    return I0 * (reactionBeta / Vt) * (-exp(reactionBeta * ((-phi_li + 2.95) / Vt)) - exp(-reactionBeta * ((-phi_li + 2.95) / Vt)));
}

// Compute the Right Hand Side (RHS) of the differential equations
void ComputeRHS(double *rhs, double *x, double *x_old)
{
	//electron conduction
	rhs[0]=(x[1]-x[0])/h;
	for(int i=1;i<N;i++)
		rhs[i]=sigma*(x[i+1]+x[i-1]-2*x[i])/pow(h,2.0) + F*RC(x[i],x[i+(N+1)],x[i+3*(N+1)],x[i+4*(N+1)]);
	rhs[N]=-sigma*A*(x[N]-x[N-1])/h -I;

	//lithium conductivity
	rhs[0+(N+1)]=-kappa*(x[1+(N+1)]-x[0+(N+1)])/h + kappa_d*(log(x[1+2*(N+1)])-log(x[0+2*(N+1)]))/h - F*RA(x[0+(N+1)]);
	for(int i=1;i<N;i++)
		rhs[i+(N+1)]=kappa*(x[i+1+(N+1)]+x[i-1+(N+1)]-2*x[i+(N+1)])/pow(h,2.0) - kappa_d*(log(x[i+1+2*(N+1)]) + log(x[i-1+2*(N+1)]) - 2*log(x[i+2*(N+1)]))/pow(h,2.0) - F*RC(x[i],x[i+(N+1)],x[i+3*(N+1)],x[i+4*(N+1)]);
	rhs[N+(N+1)]=-kappa*(x[N+(N+1)] - x[N-1+(N+1)])/h + kappa_d*(log(x[N+2*(N+1)]) - log(x[N-1+2*(N+1)]))/h;

	//lithium diffusion
	rhs[0+2*(N+1)] = x[0+2*(N+1)]-1e-3;
	for(int i=1; i<N; i++)
		rhs[i+2*(N+1)] = (x[i+4*(N+1)]*x[i+2*(N+1)]-x_old[i+4*(N+1)]*x_old[i+2*(N+1)])/tau - D_li*(x[i+1+2*(N+1)]+x[i-1+2*(N+1)]-2*x[i+2*(N+1)])/pow(h,2.0) + (1-tplus)*RC(x[i],x[i+(N+1)],x[i+3*(N+1)],x[i+4*(N+1)]);
	rhs[N+2*(N+1)]=(x[N+2*(N+1)]-x[N-1+2*(N+1)])/h;

	//oxygen diffusion
	rhs[0+3*(N+1)]=(x[1+3*(N+1)]-x[0+3*(N+1)])/h;
	for(int i=1;i<N;i++)
		rhs[i+3*(N+1)]=(x[i+4*(N+1)]*x[i+3*(N+1)]- x_old[i+4*(N+1)]*x_old[i+3*(N+1)])/tau - D_o2*(x[i+1+3*(N+1)]+x[i-1+3*(N+1)]-2*x[i+3*(N+1)])/pow(h,2.0) + (RC(x[i],x[i+(N+1)],x[i+3*(N+1)],x[i+4*(N+1)]))/2;
	rhs[N+3*(N+1)]=x[N+3*(N+1)]-c0;

	//porosity change
	for(int i=0;i<=N;i++)
		rhs[i+4*(N+1)]=(x[i+4*(N+1)]-x_old[i+4*(N+1)])/tau + (M_li2o2/(2*rho_li2o2))*RC(x[i],x[i+(N+1)],x[i+3*(N+1)],x[i+4*(N+1)]);
}

// Compute the Jacobian matrix for the system of equations
void ComputeJ(double **J, double *x, double *x_old)
{
	for (int i=0; i<5*(N+1);i++)
	{
		for(int j=0;j<5*(N+1);j++)
			J[i][j]=0;
	}

	//electron conduction
	J[0][0]=-1/h;
	J[0][1]=1/h; //0th equation in terms of x[1]
	for(int i=1;i<N;i++)
	{
		//														phi, phi_li,	c_o2,		 eps
		//rhs[i]=sigma*(x[i+1]+x[i-1]-2*x[i])/pow(h,2.0) + F*RC(x[i],x[i+(N+1)],x[i+3*(N+1)],x[i+4*(N+1)]);
		J[i][i-1]=sigma/pow(h,2.0);
		J[i][i]=-2*sigma/pow(h,2.0) + F*dRC_dphi(x[i],x[i+(N+1)],x[i+3*(N+1)],x[i+4*(N+1)]);
		J[i][i+1]=sigma/pow(h,2.0);
		J[i][i+(N+1)]=F*dRC_dphi_li(x[i],x[i+(N+1)],x[i+3*(N+1)],x[i+4*(N+1)]);
		J[i][i+3*(N+1)]=F*dRC_dc_o2(x[i],x[i+(N+1)],x[i+3*(N+1)],x[i+4*(N+1)]);
		J[i][i+4*(N+1)]=F*dRC_deps(x[i],x[i+(N+1)],x[i+3*(N+1)],x[i+4*(N+1)]);
	}
	J[N][N]=-sigma*A/h;
	J[N][N-1]=sigma*A/h;

	//lithium conductivity
	//																										 phi_li
	//rhs[0+(N+1)]=-kappa*(x[1+(N+1)]-x[0+(N+1)])/h + kappa_d*(log(x[1+2*(N+1)])-log(x[0+2*(N+1)]))/h - F*RA(x[0+(N+1)]);
	J[0+(N+1)][1+(N+1)]=-kappa/h;
	J[0+(N+1)][0+(N+1)]=kappa/h - F*dRA_dphi_li(x[0+(N+1)]);
	J[0+(N+1)][1+2*(N+1)]=(kappa_d/h)*(1/x[1+2*(N+1)]);
	J[0+(N+1)][0+2*(N+1)]=(-kappa_d/h)*(1/x[0+2*(N+1)]);
	for(int i=1; i<N; i++)
	{
		//																																										phi, phi_li,	c_o2,		eps
		//rhs[i+(N+1)]=kappa*(x[i+1+(N+1)]+x[i-1+(N+1)]-2*x[i+(N+1)])/pow(h,2.0) - kappa_d*(log(x[i+1+2*(N+1)]) + log(x[i-1+2*(N+1)]) - 2*log(x[i+2*(N+1)]))/pow(h,2.0) - F*RC(x[i],x[i+(N+1)],x[i+3*(N+1)],x[i+4*(N+1)]);
		J[i+(N+1)][i+1+(N+1)]=kappa/pow(h,2.0);
		J[i+(N+1)][i-1+(N+1)]=kappa/pow(h,2.0);
		J[i+(N+1)][i+(N+1)]=-2*kappa/pow(h,2.0) -F*dRC_dphi_li(x[i],x[i+(N+1)],x[i+3*(N+1)],x[i+4*(N+1)]);
		J[i+(N+1)][i+1+2*(N+1)]=(-kappa_d/pow(h,2.0))*(1/x[i+1+2*(N+1)]);
		J[i+(N+1)][i-1+2*(N+1)]=(-kappa_d/pow(h,2.0))*(1/x[i-1+2*(N+1)]);
		J[i+(N+1)][i+2*(N+1)]=(2*kappa_d/pow(h,2.0))*(1/x[i+2*(N+1)]);
		J[i+(N+1)][i]=-F*dRC_dphi(x[i],x[i+(N+1)],x[i+3*(N+1)],x[i+4*(N+1)]);
		J[i+(N+1)][i+3*(N+1)]=-F*dRC_dc_o2(x[i],x[i+(N+1)],x[i+3*(N+1)],x[i+4*(N+1)]);
		J[i+(N+1)][i+4*(N+1)]=-F*dRC_deps(x[i],x[i+(N+1)],x[i+3*(N+1)],x[i+4*(N+1)]);
	}
	//rhs[N+(N+1)]=-kappa*(x[N+(N+1)] - x[N-1+(N+1)])/h + kappa_d*(log(x[N+2*(N+1)]) - log(x[N-1+2*(N+1)]))/h;
	J[N+(N+1)][N+(N+1)]=-kappa/h;
	J[N+(N+1)][N-1+(N+1)]=kappa/h;
	J[N+(N+1)][N+2*(N+1)]=(kappa_d/h)*(1/x[N+2*(N+1)]);
	J[N+(N+1)][N-1+2*(N+1)]=(-kappa_d/h)*(1/x[N-1+2*(N+1)]);

	//lithium diffusion
	J[0+2*(N+1)][0+2*(N+1)]=1;
	for(int i=1; i<N; i++)
	{
		//																																									 phi, phi_li,	 c_o2,		  eps
		//rhs[i+2*(N+1)] = (x[i+4*(N+1)]*x[i+2*(N+1)]-x_old[i+4*(N+1)]*x_old[i+2*(N+1)])/tau - D_li*(x[i+1+2*(N+1)]+x[i-1+2*(N+1)]-2*x[i+2*(N+1)])/pow(h,2.0) + (1-tplus)*RC(x[i],x[i+(N+1)],x[i+3*(N+1)],x[i+4*(N+1)]);
		J[i+2*(N+1)][i+4*(N+1)]=x[i+2*(N+1)]/tau + (1-tplus)*dRC_deps(x[i],x[i+(N+1)],x[i+3*(N+1)],x[i+4*(N+1)]);
		J[i+2*(N+1)][i+2*(N+1)]=x[i+4*(N+1)]/tau + 2*D_li/pow(h,2.0);
		J[i+2*(N+1)][i+1+2*(N+1)]=-D_li/pow(h,2.0);
		J[i+2*(N+1)][i-1+2*(N+1)]=-D_li/pow(h,2.0);
		J[i+2*(N+1)][i]=(1-tplus)*dRC_dphi(x[i],x[i+(N+1)],x[i+3*(N+1)],x[i+4*(N+1)]);
		J[i+2*(N+1)][i+(N+1)]=(1-tplus)*dRC_dphi_li(x[i],x[i+(N+1)],x[i+3*(N+1)],x[i+4*(N+1)]);
		J[i+2*(N+1)][i+3*(N+1)]=(1-tplus)*dRC_dc_o2(x[i],x[i+(N+1)],x[i+3*(N+1)],x[i+4*(N+1)]);
	}
	J[N+2*(N+1)][N+2*(N+1)]=1/h;
	J[N+2*(N+1)][N-1+2*(N+1)]=-1/h;

	//oxygen diffusion
	J[0+3*(N+1)][1+3*(N+1)]=1/h;
	J[0+3*(N+1)][0+3*(N+1)]=-1/h;
	for(int i=1; i<N; i++)
	{
		//																																						   phi, phi_li,	   c_o2,		eps
		//rhs[i+3*(N+1)]=(x[i+4*(N+1)]*x[i+3*(N+1)]- x_old[i+4*(N+1)]*x_old[i+3*(N+1)])/tau - D_o2*(x[i+1+3*(N+1)]+x[i-1+3*(N+1)]-2*x[i+3*(N+1)])/pow(h,2.0) + (RC(x[i],x[i+(N+1)],x[i+3*(N+1)],x[i+4*(N+1)]))/2;
		J[i+3*(N+1)][i+4*(N+1)]=x[i+3*(N+1)]/tau + dRC_deps(x[i],x[i+(N+1)],x[i+3*(N+1)],x[i+4*(N+1)])/2;
		J[i+3*(N+1)][i+3*(N+1)]=x[i+4*(N+1)]/tau + 2*D_o2/pow(h,2.0) + dRC_dc_o2(x[i],x[i+(N+1)],x[i+3*(N+1)],x[i+4*(N+1)])/2;
		J[i+3*(N+1)][i+1+3*(N+1)]=-D_o2/pow(h,2.0);
		J[i+3*(N+1)][i-1+3*(N+1)]=-D_o2/pow(h,2.0);
		J[i+3*(N+1)][i]=dRC_dphi(x[i],x[i+(N+1)],x[i+3*(N+1)],x[i+4*(N+1)])/2;
		J[i+3*(N+1)][i+(N+1)]=dRC_dphi_li(x[i],x[i+(N+1)],x[i+3*(N+1)],x[i+4*(N+1)])/2;
	}
	J[N+3*(N+1)][N+3*(N+1)]=1;

	//porosity change
	for(int i=0; i<=N; i++)
	{
		//																				 phi, phi_li,	  c_o2,		   eps
		//rhs[i+4*(N+1)]=(x[i+4*(N+1)]-x_old[i+4*(N+1)])/tau + (M_li2o2/(2*rho_li2o2))*RC(x[i],x[i+(N+1)],x[i+3*(N+1)],x[i+4*(N+1)]);
		J[i+4*(N+1)][i+4*(N+1)]=1/tau + (M_li2o2/(2*rho_li2o2))*dRC_deps(x[i],x[i+(N+1)],x[i+3*(N+1)],x[i+4*(N+1)]);
		J[i+4*(N+1)][i]=(M_li2o2/(2*rho_li2o2))*dRC_dphi(x[i],x[i+(N+1)],x[i+3*(N+1)],x[i+4*(N+1)]);
		J[i+4*(N+1)][i+(N+1)]=(M_li2o2/(2*rho_li2o2))*dRC_dphi_li(x[i],x[i+(N+1)],x[i+3*(N+1)],x[i+4*(N+1)]);
		J[i+4*(N+1)][i+3*(N+1)]=(M_li2o2/(2*rho_li2o2))*dRC_dc_o2(x[i],x[i+(N+1)],x[i+3*(N+1)],x[i+4*(N+1)]);
	}
}

// Calculate the Euclidean norm of a vector
double Norm(double *x, int n)
{
    double res = 0;
    for(int i = 0; i < n; i++)
        res += x[i] * x[i];
    return sqrt(res);
}

// Check if porosity is still positive
bool porosity(double *x) 
{
    for(int i = 4 * (N + 1); i < 5 * (N + 1); i++)
    {
        if (x[i] <= 0.0)
            return false;
    }
    return true;
}

int main()
{
    // Allocate memory for solution vectors and Jacobian matrix
    double *x = new double[5 * (N + 1)];
    double *x_old = new double[5 * (N + 1)];
    double *y = new double[5 * (N + 1)];
    double *rhs = new double[5 * (N + 1)];

    double **J = new double*[5 * (N + 1)];
    for(int i = 0; i < 5 * (N + 1); i++) 
        J[i] = new double[5 * (N + 1)];

    // Initialize J matrix to zero
    for (int i = 0; i < 5 * (N + 1); i++)
    {
        for(int j = 0; j < 5 * (N + 1); j++)
            J[i][j] = 0;
    }

    int count;
    double t = 0;
    double err = 1e-6;

    // Open file for writing battery voltage data
    ofstream myfile;
    myfile.open("output_voltage.txt");

    // Set initial conditions
    for(int i = 0; i < N + 1; i++)
        x[i] = x_old[i] = 2.95;  // Initial electric potential
    for(int i = N + 1; i < 2 * (N + 1); i++)
        x[i] = x_old[i] = 2.95;  // Initial lithium potential
    for(int i = 2 * (N + 1); i < 3 * (N + 1); i++)
        x[i] = x_old[i] = 1.e-3; // Initial lithium concentration
    for(int i = 3 * (N + 1); i < 4 * (N + 1); i++)
        x[i] = x_old[i] = c0;    // Initial oxygen concentration
    for(int i = 4 * (N + 1); i < 5 * (N + 1); i++)
        x[i] = x_old[i] = 0.75;  // Initial porosity

    // Main time-stepping loop
    do{  
        t += tau;
        count = 0;
        do{
            ComputeRHS(rhs, x, x_old);
            ComputeJ(J, x, x_old);
            Solve(J, rhs, y, 5 * (N + 1));
            
            for(int i = 0; i < 5 * (N + 1); i++)
                x[i] -= y[i];

            count++;

            if(!porosity(x)) break;  // Check for battery degradation

        } while(Norm(y, 5 * (N + 1)) > err && count < 100);
        
        for(int i = 0; i < 5 * (N + 1); i++)
            x_old[i] = x[i];

        // Write voltage data to file
        myfile << t / 86400 << "\t" << x[N] << endl;

    } while(porosity(x));

    myfile.close();

    // Write final state variables to separate files
    myfile.open("output_phi.txt");
    for (int i = 0; i < (N + 1); i++)
        myfile << i << "\t" << x[i] << endl;
    myfile.close();

    myfile.open("output_phi_li.txt");
    for (int i = 0; i < (N + 1); i++)
        myfile << i << "\t" << x[i + (N + 1)] << endl;
    myfile.close();

    myfile.open("output_c_li.txt");
    for (int i = 0; i < (N + 1); i++)
        myfile << i << "\t" << x[i + 2 * (N + 1)] << endl;
    myfile.close();

    myfile.open("output_c_o2.txt");
    for (int i = 0; i < (N + 1); i++)
        myfile << i << "\t" << x[i + 3 * (N + 1)] << endl;
    myfile.close();

    myfile.open("output_eps.txt");
    for (int i = 0; i < (N + 1); i++)
        myfile << i << "\t" << x[i + 4 * (N + 1)] << endl;
    myfile.close();

    printf("Simulation Complete");

    getchar();
}
