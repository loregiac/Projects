/***********************************************************************
 *                                                                     *
 * Code for project of "Fundamentals of CFD Methods" course            *
 *                                                                     *
 * Author:   Lorenzo Giacomel                                          *
 * Date:    13/12/2017                                                 *
 *                                                                     *
 * This code can be freely used for non-commercial purposes as long    *
 * as this header is left intact.                                      *
 * Copyright © 2017 Lorenzo Giacomel. All rights reserved.             *
 ***********************************************************************/

#include <iostream>
#include <string>
#include <fstream>
#include <math.h>

#define pi 3.14159265358979323846

/**
    Compute the source term of the fluid phase at point x for the Method of Manufactured Solutions with manufactured solutions Tf = cos(2*pi/L) and Ts = cos(4*pi/L)
 
    @param x The x coordinate of point at which the source term is computed
    @param alpha The diffusion coefficient for the fluid phase
    @param u The advection speed for the fluid phase
    @param dx The cell width
    @param h The coupling coefficient for the fluid phase
 
    @return the source term of the fluid phase at point x
 */
double s_Tf(double x, double alpha, double u, double dx, double h, int N){
    double k1 = 3 * pi * 1/(N*dx);
    double k2 = 2 * pi * 1/(N*dx);
    
    //double k1 = 2 * pi * 1 / (N * dx);
    //double k2 = 2 * pi * 2 / (N * dx);
    
    return (alpha * k1 + h / k1) * (sin(k1 * (x + 0.5 * dx)) - sin(k1 * (x - 0.5 * dx))) + u * (cos(k1 * (x + 0.5 * dx))
            - cos(k1 * (x - 0.5 * dx))) - h / k2 * (sin(k2 * (x + 0.5 * dx)) - sin(k2 * (x - 0.5 * dx)));
}

/**
    Perform one step of the FTBS method during the charging periods for the fluid phase

    @param Tf_old The fluid phase temperature distribution at time n
    @param Tf_new The fluid phase temperature distribution at time n+1
    @param Ncells The number of cells for the Finite Volumes discretization
    @param dt The time step
    @param dx The cell width
    @param alpha_f The diffusion coefficient for the fluid phase
    @param u_f The advectin speed for the fluid phase
    @param TB The temperature at the left boundary
    @param x_cells the x coordinates of the cell centroids
    @param h The coupling coefficient for the fluid phase
    @param MMS The flag which indicates wheter the simulation should include the MMS source term
 */
void step_TfCharge(double* Tf_old, double* Tf_new, int Ncells, double dt, double dx,
                   double alpha_f, double u_f, double TB, double* x_cells, double h, bool MMS){
    double diff = alpha_f * dt / (dx*dx);
    double adv = u_f * dt / dx;
    //double k = 2 * pi * 1/(Ncells*dx);
    
    Tf_new[0] = Tf_old[0] - adv * (Tf_old[0] - TB) + diff * (Tf_old[1] - Tf_old[0]) +
    (dt / dx) * s_Tf(x_cells[0], alpha_f , u_f, dx, h, Ncells) * MMS;
    
    Tf_new[Ncells-1] = Tf_old[Ncells-1] - adv * (Tf_old[Ncells - 1] - Tf_old[Ncells - 2]) +
    diff * (Tf_old[Ncells - 2] - Tf_old[Ncells - 1]) + (dt/dx) * s_Tf(x_cells[Ncells-1], alpha_f, u_f, dx, h, Ncells) * MMS;
    
    for (int i = 1; i < Ncells - 1; i++) {
        Tf_new[i] = Tf_old[i]  - adv * (Tf_old[i] - Tf_old[i-1]) + diff * (Tf_old[i+1] - 2*Tf_old[i] + Tf_old[i-1]) +
        (dt / dx) * s_Tf(x_cells[i], alpha_f, u_f, dx, h, Ncells) * MMS;
    }
}

/**
    Perform one step of the FTFS method during the discharging periods for the fluid phase
 
    @param Tf_old The fluid phase temperature distribution at time n
    @param Tf_new The fluid phase temperature distribution at time n+1
    @param Ncells The number of cells for the Finite Volumes discretization
    @param dt The time step
    @param dx The cell width
    @param alpha_f The diffusion coefficient for the fluid phase
    @param u_f The advectin speed for the fluid phase
    @param TB The temperature at the right boundary
    @param x_cells the x coordinates of the cell centroids
    @param h The coupling coefficient for the fluid phase
    @param MMS The flag which indicates wheter the simulation should include the MMS source term
 */
void step_TfDischarge(double* Tf_old, double* Tf_new, int Ncells, double dt, double dx,
                      double alpha_f, double u_f, double TB, double* x_cells, double h, bool MMS){
    double diff = alpha_f * dt / (dx * dx);
    double adv = u_f * dt / dx;
    double k = 3 * pi * 1/(Ncells * dx);
    
    Tf_new[0] = Tf_old[0] - adv * (Tf_old[1] - Tf_old[0]) + diff * (Tf_old[1] - Tf_old[0])
    +  (dt / dx)* s_Tf(x_cells[0], alpha_f, u_f, dx, h, Ncells) * MMS;
    
    Tf_new[Ncells - 1] = Tf_old[Ncells - 1] - adv * (TB - Tf_old[Ncells - 1])
    + diff * (Tf_old[Ncells - 2] - Tf_old[Ncells - 1])
    + (dt / dx) * s_Tf(x_cells[Ncells - 1], alpha_f, u_f, dx, h, Ncells) * MMS;
    
    for (int i = 1; i < Ncells - 1; i++) {
        Tf_new[i] = Tf_old[i]  - adv * (Tf_old[i + 1] - Tf_old[i]) + diff * (Tf_old[i + 1] - 2*Tf_old[i] + Tf_old[i - 1])
        + (dt / dx) * s_Tf(x_cells[i], alpha_f, u_f, dx, h, Ncells) * MMS;
    }
}

/**
    Compute the source term of the solid phase for the Method of Manufactured Solutions with manufactured solutions Tf = cos(2*pi/L) and Ts = cos(4*pi/L)
 
    @param x The x coordinate of point at which the source term is computed
    @param alpha The diffusion coefficient for the solid phase
    @param dx The cell width
    @param h The coupling coeeficient for the solid phase
 
    @return the source term of the solid phase at point x
 */
double s_Ts(double x, double alpha, double dx, double h, int N){
    double k1 = 3 * pi * 1 / (N * dx);
    double k2 = 2 * pi * 1 / (N * dx);
    return (alpha * k2 + h / k2) * (sin(k2 * (x + 0.5 * dx)) - sin(k2 * (x - 0.5 * dx)))-
    h / k1 * (sin(k1 * (x + 0.5 * dx)) - sin(k1 * (x - 0.5 * dx)));
}

/**
    Perform one step of the FTCS method for the solid phase
 
    @param Ts_old The solid phase temperature distribution at time n
    @param Ts_new The solid phase temperature distribution at time n+1
    @param Ncells The number of cells for the Finite Volumes discretization
    @param dt The time step
    @param dx The cell width
    @param alpha_s The diffusion coefficient for the solid phase
    @param x_cells the x coordinates of the cell centroids
    @param h The coupling coefficient for the solid phase
    @param MMS The flag which indicates wheter the simulation should include the MMS source term
 */
void step_Ts(double* Ts_old, double* Ts_new, int Ncells, double dt, double dx, double alpha_s,
             double* x_cells, double h, int MMS){
    double diff = alpha_s * dt / (dx*dx);
    
    Ts_new[0] = Ts_old[0] + diff * (Ts_old[1] - Ts_old[0]) + (dt / dx) * s_Ts(x_cells[0], alpha_s, dx, h, Ncells) * MMS;
    Ts_new[Ncells - 1] = Ts_old[Ncells - 1] + diff * (Ts_old[Ncells - 2] - Ts_old[Ncells - 1]) +
    (dt / dx) * s_Ts(x_cells[Ncells-1], alpha_s, dx, h, Ncells) * MMS;
    
    for (int i = 1; i < Ncells - 1; i++) {
        Ts_new[i] = Ts_old[i]  + diff * (Ts_old[i + 1] - 2 * Ts_old[i] + Ts_old[i - 1]) +
        (dt / dx) * s_Ts(x_cells[i], alpha_s, dx, h, Ncells) * MMS;
    }
}


/**
    Compute the temperature distributions of the solid and fluid phases after one time step
    according to the coupled Schumann equations
 
    @param Tf_new The fluid phase temperature distribution at time n+1
    @param Ts_new The solid phase temperature distribution at time n+1
    @param Tf_old The fluid phase temperature distribution at time n
    @param Ts_old The solid phase temperature distribution at time n
    @param h_vf The coupling coefficient for the fluid phase
    @param h_vs The coupling coefficient for the solid phase
    @param Ncells The number of cells for the Finite Volumes discretization
    @param dt The time step
    @param dx The cell width
    @param alpha_f The diffusion coefficient for the fluid phase
    @param alpha_s The diffusion coefficient for the solid phase
    @param TB The boundary temperature, i.e. the temperature of the inflowing fluid
    @param x_cells the x coordinates of the cell centroids
    @param MMS The flag which indicates wheter the simulation should include the MMS source term
 */
void coupledSolver(double* Tf_new, double* Ts_new, double* Tf_old, double* Ts_old,
                   double h_vf, double h_vs, int Ncells, double dt, double dx, double u_f,
                   double alpha_f, double alpha_s, double TB, double* xcells, int MMS){
    
    double Tf_star[Ncells];
    double Ts_star[Ncells];
    
    if (u_f > 0) {
        step_TfCharge(Tf_old, Tf_star, Ncells, dt, dx, alpha_f, u_f, TB, xcells, h_vf, MMS);
    } else {
        step_TfDischarge(Tf_old, Tf_star, Ncells, dt, dx, alpha_f, u_f, TB, xcells, h_vf, MMS);
    }
    
    step_Ts(Ts_old, Ts_star, Ncells, dt, dx, alpha_s, xcells, h_vs, MMS);
    
    
    
    double alpha = 1 + h_vf * dt;
    double beta = - h_vf * dt;
    double gamma = - h_vs * dt;
    double delta = 1 + h_vs * dt;
    
    
    for (int i = 0; i < Ncells; i++) {
        Tf_new[i] = (beta * Ts_star[i] - delta * Tf_star[i]) / (gamma * beta - alpha * delta);
        Ts_new[i] = (alpha * Ts_star[i] - gamma * Tf_star[i]) / (alpha * delta - gamma * beta);
    }
    
}

