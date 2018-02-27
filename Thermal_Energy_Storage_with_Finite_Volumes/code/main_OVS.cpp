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
#include <assert.h>
#include "init.hpp"
#include "write.hpp"
#include "step.hpp"
#include "errors.hpp"
#include "reconstruct.hpp"
#include "create_states.hpp"
#include <utility>

#define pi 3.14159265358979323846

int main(int argc, char* argv[]){
    
    unsigned int N_cells;                               //Number of cells for Finite Volumes discretization
    std::cout<<"Please enter the number of cells:"<<std::endl;
    std::cin>>N_cells;
    double d = 1;
    double height = 1;

    double dx = height/N_cells;                                          //Cells size
    double alpha_f;                                     //Diffusion coefficient for the fluid phase
    double alpha_s;                                     //Diffusion coefficient for the solid phase
    double u_f;                                         //Advection speed of the fluid phase
    double eps = 0.4;
    double d_s = 0.03;
    double rho_s = 2600.0;
    double rho_f = 1835.6;
    double C_s = 900.0;
    double C_pf = 1511.0;
    double k_s = 2.0;
    double k_f = 0.52;
    double mu_f = 2.63;
    double mdot_f = 0.1;
    double h_vf, h_vs;
    double TB = 873.0;
    double LinfErr_Tf, LinfErr_Ts;
    int N_points = N_cells+1;
    double Tf_old[N_cells];                             //array for temperature of the fluid phase at time n
    double Ts_old[N_cells];                             //array for temperature of the solid phase at time n
    double Tf_new[N_cells];                             //array for temperature of the fluid phase at time n
    double Ts_new[N_cells];                             //array for temperature of the solid phase at time n+1
    double Tf_points[N_cells+1];                        //array for temperature of the fluid phase after the reconstruction
    double Ts_points[N_cells+1];                        //array for temperature of the solid phase after the reconstruction
    double x_cells[N_cells];                            //array for the x coordinate of the cell centroids
    double x_points[N_cells+1];                         //array for the x coordinate of the cell interfaces

    //boundary and initial condition
    for(int i = 0; i<N_cells; i++){
        Tf_old[i] = 288.15;
        Ts_old[i] = 288.15;
        x_cells[i] = 0.5*dx + i*dx;
    }

    for(int i = 0; i<N_cells+1; i++){
        x_points[i] = i*dx;
    }

    write_temperature(Tf_old, Ts_old, x_cells, N_cells, "temperatures/temperature"+std::to_string(0)+".dat");

    int n=0;
    double diff = 1;
    u_f = (mdot_f) / ((rho_f) * pi * (d) * (d) / 4 * (eps));                  //advection speed
    double Re = (eps) * (rho_f) * (u_f) * (d_s) / mu_f;                      //Reynolds number
    double Pr = mu_f * (C_pf) / k_f;                                       //Prandtl number
    double Nu_fs = 0.255 / (eps) * pow(Pr,1.0/3) * pow(Re,2.0/3);          //Nusselt number
    double h_fs = Nu_fs * k_f / (d_s);                                      //fluid-solid heat transfer
    double h = 1.0 / (1.0 / h_fs + (d_s) / (10.0 * k_s));                         //overall heat transfer
    double h_v = 6 * (1 - (eps)) * h / (d_s);                                  //volumetric heat transfer
    h_vf = h_v / ((eps) * (rho_f) * (C_pf));
    h_vs = h_v / ((1 - (eps)) * (rho_s) * (C_s));
    alpha_f = k_f / ((eps) * (rho_f) * (C_pf));
    alpha_s = k_f / ((1 - (eps)) * (rho_s) * (C_s));
    //u_f = alpha_f;
    double dt = 0.01;
    double s = u_f * dt/dx;
    double di = alpha_f * dt/(dx*dx);
    assert(s*s <= s + 2*di && s + 2*di <=1);
    int nim = 0;
    while(diff > 1e-10 && n*dt<5000){
        n++;
        //perform one step of the BTCS
        coupledSolver(Tf_new, Ts_new, Tf_old, Ts_old, h_vf, h_vs, N_cells, dt, dx, u_f, alpha_f, alpha_s, TB, x_cells, 0);
        diff = std::max(infNorm(Tf_new,Tf_old,N_cells),infNorm(Ts_new,Ts_old,N_cells));
        for(int j = 0; j < N_cells; j++){
            Tf_old[j] = Tf_new[j];
            Ts_old[j] = Ts_new[j];
        }
    }

    //reconstruction
    for(int i = 1; i < N_points-1; i++){
        Ts_points[i] = 0.5*(Ts_new[i-1]+Ts_new[i]);
        Tf_points[i] = 0.5*(Tf_new[i-1]+Tf_new[i]);
    }
    Ts_points[0] = Ts_new[0];
    Ts_points[N_points-1] = Ts_new[N_cells-1];
    Tf_points[0] = TB;
    Tf_points[N_points-1] = Tf_new[N_cells-1];

    write_temperature(Tf_points, Ts_points, x_points, N_points, "comparison/comparison"+std::to_string(N_cells)+".dat");
    LinfErr_Tf = error_spaceTf(Tf_points, N_points, x_points, dx);
    LinfErr_Ts = error_spaceTs(Ts_points, N_points, x_points, dx);
    double L2Err_Tf = errorL2_spaceTf(Tf_points, N_points, x_points, dx);
    double L2Err_Ts = errorL2_spaceTs(Ts_points, N_points, x_points, dx);
    std::cout<<"L_inf err Tf = "<< LinfErr_Tf <<std::endl;
    std::cout<<"L_inf err Ts = "<< LinfErr_Ts <<std::endl;
    std::cout<<"L_2 err Tf = "<< L2Err_Tf <<std::endl;
    std::cout<<"L_2 err Ts = "<< L2Err_Ts <<std::endl;

}



