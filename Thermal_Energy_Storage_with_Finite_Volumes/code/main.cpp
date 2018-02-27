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
#include "create_states.hpp"
#include "reconstruct.hpp"
#include "postprocessing.hpp"


#define pi 3.14159265358979323846

int main(int argc, char* argv[]){
    
    double h;                                           //Height of the storage
    double d;                                           //Diameter of the storage
    unsigned int N_cells;                               //Number of cells for Finite Volumes discretization
    double dx;                                          //Cells size
    double D1, D2, D3, D4;                              //Duration of charging, idle charging -> discharging,
                                                            //discharging, idle discharging -> charging
    unsigned int N_cycles;                               //Number of cycles
    unsigned int Nt;                                    //Number of time steps per cycle
    double alpha_f;                                     //Diffusion coefficient for the fluid phase
    double alpha_s;                                     //Diffusion coefficient for the solid phase
    double u_f;                                         //Advection speed of the fluid phase
    double h_vf;                                        //Coupling term for the fluid phase
    double h_vs;                                        //Coupling term for the solid phase
    double Tc;                                          //Temperature of the inflowing fluid during charge period
    double Td;                                          //Temperature of the inflowing fluid during discharge period
    double mdot_f;                                      //Mass flow rate for the fluid phase
    double C_pf;                                        //Specific heat at constant pressure of the fluid phase
    double C_s;                                         //Specific heat of the solid phase
    double eps;                                         //Porosity
    double rho_f;                                       //Density of the fluid phase
    double rho_s;                                       //Density of the solid phase

    if (argc <= 1) {
        std::cerr << "Please insert either file or terminal option"<<std::endl;
        abort();
    }
    
    std::string flag = argv[1];
    
    //Initializations from file
    init(&h, &d, &N_cells, &Tc, &Td, &dx, &D1, &D2, &D3, &D4, &N_cycles, &Nt, &alpha_f,
         &alpha_s, &u_f, &h_vf, &h_vs, &mdot_f, &C_pf, &C_s, &eps, &rho_f, &rho_s,flag);
    
    unsigned int N_points = N_cells+1;                  //Number of points for Finite Volumes discretization
    double Tf_old[N_cells];                             //Array for cells temperature of the fluid phase at time n
    double Ts_old[N_cells];                             //Array for cells temperature of the solid phase at time n
    double Tf_new[N_cells];                             //Array for cells temperature of the fluid phase at time n + 1
    double Ts_new[N_cells];                             //Array for cells temperature of the solid phase at time n + 1
    double Tf_points[N_cells + 1];                      //Array for temperature of the fluid phase after the reconstruction
    double Ts_points[N_cells + 1];                      //Array for temperature of the solid phase after the reconstruction
    double x_cells[N_cells];                            //Array for the x coordinate of the cell centroids
    double x_points[N_cells + 1];                       //Array for the x coordinate of the cell interfaces

    //Initialization of boundary and initial conditions and of x_cells
    for(int i = 0; i < N_cells; i++){
        Tf_old[i] = Td;
        Ts_old[i] = Td;
        x_cells[i] = 0.5 * dx + i * dx;
    }

    //Initialization of x_points
    for(int i = 0; i<N_cells + 1; i++){
        x_points[i] = i * dx;
    }

    //create state array
    int Nt_tot = Nt * N_cycles;                          //Total number of time steps
    double dt = (D1 + D2 + D3 + D4) * N_cycles / (Nt_tot - 1);   //Time step
    int* states = new int[Nt_tot];                      //Vector storing the storage state at each time instant
    
    create_states(D1, D2, D3, D4, N_cycles, Nt, states, dt);     //Initialize states vector
    write_temperature(Tf_old, Ts_old, x_cells, N_cells,
                      "temperatures/temperature"+std::to_string(0)+".dat");     //Write initial temperature to external file

    //Check if the chosen dt and dx are within the stability region
    double s = u_f * dt / dx;
    double di = alpha_f * dt / (dx * dx);
    assert(s * s <= s + 2 * di);
    assert(s + 2 * di <= 1);
    
    int nim = 0;                                        //Auxiliary variable for counting the number of saved images
    double u_fn;                                        //Advection speed at n time step
    double TB;
    double Tf_max_outflow = 0;                          //Max temperature of the fluid at the outflow
    double Tf_charge_out[Nt];                     //Temperatures of the fluid at the outflow at each time instant of a charge period
    double Tf_charge_in[Nt];                      //Temperatures of the fluid at the inflow at each time instant of a charge period
    double Tf_lastcharge[N_points];               //Temperature distribution of the fluid phase after the last charge period
    double Ts_lastdischarge[N_points];            //Temperature distribution of the solid phase after the last discharge period
    double Tf_discharge_out[Nt];                  //Temperatures of the fluid at the outflow at each time instant of a discharge
                                                    //period
    double Tf_discharge_in[Nt];                   //Temperatures of the fluid at the outflow at each time instant of a discharge
                                                    //period
    int u=0;
    int v=0;
    for (int n=1; n <= Nt_tot; n++){
        //Set the right advection speed
        u_fn = states[n] * u_f;
        //Set the right boundary temperature
        TB = Tc - (Tc - Td) * (states[n] == -1);
        //Perform one step of the BTCS
        coupledSolver(Tf_new, Ts_new, Tf_old, Ts_old, h_vf, h_vs, N_cells, dt, dx, u_fn, alpha_f, alpha_s, TB, x_cells, 0);

        //Update old solution
        for (int j = 0; j < N_cells; j++) {
            Tf_old[j] = Tf_new[j];
            Ts_old[j] = Ts_new[j];
        }
        
        //Carry out the reconstruction
        reconstruct(Tf_points, Ts_points, Ts_new, Tf_new, TB, N_points, N_cells, states, n);
        
        //Last charging phase
        if(n >= Nt_tot - Nt && n < Nt_tot - 3.0 / 4.0 * Nt) {
            Tf_max_outflow = std::max(Tf_max_outflow, Tf_new[N_cells - 1]);
            Tf_charge_out[u] = Tf_points[0];
            Tf_charge_in[u] = Tf_points[N_points - 1];
            u++;
        }
        
        //Last discharging phase
        if(n >= Nt_tot - Nt / 2.0 && n < Nt_tot - Nt / 4.0){
            Tf_max_outflow = std::max(Tf_max_outflow, Tf_new[N_cells - 1]);
            Tf_discharge_out[v] = Tf_points[N_points - 1];
            Tf_discharge_in[v] = Tf_points[0];
            v++;
        }
        
        //Before last charging phase
        if(n == Nt_tot - 3.0 / 4.0 * Nt){
            for (int i = 0; i < N_points; i++) {
                Tf_lastcharge[i] = Tf_points[i];
            }
        }
        
        //Before last discharging phase
        if (n == Nt_tot - Nt / 4.0) {
            for (int i = 0; i < N_points; i++) {
                Ts_lastdischarge[i] = Ts_points[i];
            }
        }
        
        //Write temperatures to file for animation purposes
        if(n % ((int) (Nt_tot / 1080.0)) == 0){
            nim++;
            write_temperature(Tf_points, Ts_points, x_points, N_points, "temperatures/temperature"+std::to_string(n)+".dat");
        }
    }
    
    //Terminal outputs
    std::cout << "================ results for storage design ================" << std::endl;
    std::cout << "diameter: "<<d<<std::endl;
    std::cout << "max temperature difference at outflow: " << Tf_max_outflow-Td<<std::endl;
    double eta = compute_cycle_exergy_efficiency(mdot_f, C_pf, dx, N_points, Tf_charge_in, Tf_charge_out, Tf_discharge_in, Tf_discharge_out);
    std::cout << "eta: "<<eta<<std::endl;
    double cap = compute_capacity_factor(Tf_lastcharge, Ts_lastdischarge, d, h, eps, rho_f, C_pf, rho_s, C_s, Td, Tc, N_points, dx);
    std::cout << "capacity factor: "<<cap<<std::endl;
    
    std::cout << "================= parameters for animation =================" << std::endl;
    std::cout << "Number of saved images: " << nim << std::endl;
    std::cout << "Time interval between images: " << (int) (Nt_tot / 1080.0) << std::endl;
    std::cout << "dt: " << dt << std::endl;
    
    delete[] states;
}



