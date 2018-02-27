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
    Carry out the reconstruction of the solution for the solid phase through piecewise linear functions.
 
    @param Tf_points The temperature of the fluid phase at the cell interfaces
    @param Ts_points The temperature of the solid phase at the cell interfaces
    @param Tf_cells The temperature of the fluid phase at the cell centroids
    @param Ts_cells The temperature of the solid phase at the cell centroids
    @param TB Temperature at the left/right boundary
    @param Npoints The Number of points fo the Finite Volumes discretization
    @param Ncells The Number of points fo the Finite Volumes discretization
    @param states The vector storing the state of the storage at each time instant
    @param n The current time instant

 */
void reconstruct(double* Tf_points, double* Ts_points, double* Ts_cells, double* Tf_cells, double TB, unsigned int Npoints, unsigned int Ncells, int* states, int n)
{
    for (int i = 1; i < Npoints-1; i++) {
        Ts_points[i] = 0.5 * (Ts_cells[i - 1] + Ts_cells[i]);
        Tf_points[i] = 0.5 * (Tf_cells[i - 1] + Tf_cells[i]);
    }
    
    if (states[n] == 1) {
        Ts_points[0] = Ts_cells[0];
        Ts_points[Npoints - 1] = Ts_cells[Ncells - 1];
        Tf_points[0] = TB;
        Tf_points[Npoints - 1] = Tf_cells[Ncells - 1];
    } else if (states[n] == -1) {
        Ts_points[0] = Ts_cells[0];
        Ts_points[Npoints - 1] = Ts_cells[Ncells - 1];
        Tf_points[0] = Tf_cells[0];
        Tf_points[Npoints - 1] = TB;
    } else {
        Ts_points[0] = Ts_cells[0];
        Ts_points[Npoints - 1] = Ts_cells[Ncells - 1];
        Tf_points[0] = Tf_cells[0];
        Tf_points[Npoints - 1] = Tf_cells[Ncells - 1];
    }
}
