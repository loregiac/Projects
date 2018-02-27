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
    Sets the parameters for the simulations.
 
    @param D1 The duration of a charge period
    @param D2 The duration of a idle period between charging and discharging
    @param D3 The duration of a discharge period
    @param D4 The duration of a idle period between discharging and charging
    @param Ncycles The number of charging-idle-discharging-idle cycles
    @param Nt The number of time steps per cycle
    @param states the vector storing the state of the storage at each time instant
    @param dt The time step
 */
void create_states(double D1, double D2, double D3, double D4, int Ncycles, int Nt, int* states, double dt){
    int N1 = D1/dt;
    int N2 = N1 + D2/dt;
    int N3 = N2 + D3/dt;
    int N4 = Nt;
    
    for (int i = 0; i < Ncycles; i++) {
        for (int j = 0; j < N1; j++) {
            states[j + i * Nt] = 1;
        }
    }
    for (int i = 0; i < Ncycles; i++) {
        for (int j = N1; j < N2; j++) {
            states[j + i * Nt] = 0;
        }
    }
    for (int i = 0; i < Ncycles; i++) {
        for (int j = N2; j < N3; j++) {
            states[j + i * Nt] = -1;
        }
    }
    for (int i = 0; i < Ncycles; i++) {
        for (int j = N3; j < N4; j++) {
            states[j + i * Nt] = 0;
        }
    }
}
