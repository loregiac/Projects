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
    Write the temperature distributions for both the fluid and the solid phases to an external file
 
    @param Tf The temperature distribution of the fluid phase
    @param Ts The temperature distribution of the solid phase
    @param x The x coordinates of the points at which the temperatures are given
    @param N The number of elements in Tf, Ts and x arrays
    @param name The name of the file to which the temperatures have to be written
 
 */
void write_temperature(double* Tf, double* Ts, double* x, unsigned int N, std::string name)
{
    std::ofstream out_file(name, std::ios::out);
    
    for (int i = 0; i < N; i++) {
        out_file << x[i] << '\t' << Tf[i] << '\t' << Ts[i] << "\n";
    }
    
    out_file.close();
}

/**
    Write the state of the storage at every time instant to an external file
 
    @param t The time instants
    @param states The state of each time instant
    @param Nt_tot The number of time instants 
 */
void write_states(double* t, double* states, int Nt_tot)
{
    std::ofstream out_file("states.dat", std::ios::out);
    out_file << "t" << '\t' << "state"<< "\n";
    
    for (int i = 0; i < Nt_tot; i++) {
        out_file << t[i] << '\t' << states[i] << "\n";
    }
    
    out_file.close();
}
