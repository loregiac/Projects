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
    Computes the L^inf error for the fluid phase with respect to the manufactured solution cos(3*pi*x/L)
 
    @param T The temperature of the fluid phase
    @param N The number of points
    @param x The x coordinates of the points
    @param dx The cell width
 
    @return The L^inf error for the fluid phase with respect to the manufactured solution cos(3*pi*x/L)
*/
double error_spaceTf(double* T, unsigned int N, double* x, double dx){
    double err = 0;
    double k = 3 * pi * 1 / ((N - 1) * dx);
    for (int i = 0; i < N; i++) {
        err = std::max(err, fabs(T[i] - cos(k * x[i])));
    }
    return err;
}

/**
    Computes the L^inf error for the solid phase with respect to the manufactured solution cos(2*pi*x/L)
 
    @param T The temperature of the solid phase
    @param N The number of points
    @param x The x coordinates of the points
    @param dx The cell width
 
    @return The L^inf error for the fluid phase with respect to the manufactured solution cos(3*pi*x/L)
 */
double error_spaceTs(double* T, unsigned int N, double* x, double dx){
    double err = 0;
    double k = 2 * pi * 1 / ((N-1) * dx);
    for (int i = 0; i < N; i++) {
        err = std::max(err, fabs(T[i] - cos(k * x[i])));
    }
    return err;
}

/**
    Computes the normalized L^2 error for the fluid phase with respect to the manufactured solution cos(3*pi*x/L)
 
    @param T The temperature of the fluid phase
    @param N The number of points
    @param x The x coordinates of the points
    @param dx The cell width
 
    @return The L^2 error for the fluid phase with respect to the manufactured solution cos(3*pi*x/L)
 */
double errorL2_spaceTf(double* T, unsigned int N, double* x, double dx){
    double err = 0;
    double k = 3 * pi * 1 / ((N - 1) * dx);
    
    for (int i = 0; i < N; i++) {
        err += (T[i] - cos(k * x[i])) * (T[i] - cos(k * x[i]));
    }
    
    return sqrt(err / N);
}

/**
    Computes the normalized L^2 error for the solid phase with respect to the manufactured solution cos(2*pi*x/L)
 
    @param T The temperature of the fluid phase
    @param N The number of points
    @param x The x coordinates of the points
    @param dx The cell width
 
    @return The L^2 error for the solid phase with respect to the manufactured solution cos(2*pi*x/L)
 */
double errorL2_spaceTs(double* T, unsigned int N, double* x, double dx){
    double err = 0;
    double k = 2 * pi * 1 / ((N-1) * dx);
    
    for (int i = 0; i < N; i++) {
        err += (T[i] - cos(k * x[i])) * (T[i] - cos(k * x[i]));
    }
    
    return sqrt(err / N);
}

/**
    Computes the L^inf norm of the difference between two vectors
 
    @param vec1 The first vector to be compared
    @param vec2 The second vector to be compared
    @param N The number of elements in the two vectors
 
    @return the L^inf norm of the difference between two vectors
 */
double infNorm(double* vec1, double* vec2, unsigned int N){
    double res = 0;
    for (int i = 0; i < N; i++) {
        res = std::max(res, fabs(vec1[i] - vec2[i]));
    }
    return res;
}
