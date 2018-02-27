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
    Computes the exergy flux for the fluid phase
 
    @param T The temperature of the fluid phase
    @param Npoints The number of elements in the temperature vector
    @param dx The cell width
    @param mdot_f The mass flow rate of the fluid phase
    @param C_pf The specific heat at constant pressure of the flui phase
 
    @return The exergy flux for the fluid phase
 */
double compute_exergy_flux(double* T, int Npoints, double dx, double m_dot, double C_pf){
    double result = 0;
    double integ[Npoints];
    
    for (int i = 0; i < Npoints; i++) {
        integ[i] = m_dot * C_pf * (T[i] - 288.15 - 288.15 * log(T[i] / 288.15));
    }
    
    for (int i = 1; i <= Npoints - 2; i++) {
        result += 2 * integ[i];
    }
    
    result += integ[0] + integ[Npoints - 1];
    return dx / 2.0 * result;
    
}

/**
    Computes the exergy efficiency for a cycle using the trapezoidal rule for the integrals
 
    @param mdot_f The mass flow rate of the fluid phase
    @param C_pf The specific heat at constant pressure of the flui phase
    @param dx The cell width
    @param Npoints The number of elements in the temperature vectors
    @param T_charge_in The inflow temperature of the fluid phase during a charge period
    @param T_charge_out The outflow temperature of the fluid phase during a charge period
    @param T_discharge_in The inflow temperature of the fluid phase during a discharge period
    @param T_discharge_out The outflow temperature of the fluid phase during a discharge period
 
    @return the exergy efficiency for a cycle
 */
double compute_cycle_exergy_efficiency(double m_dot, double C_pf, double dx, int Npoints, double* T_charge_in, double* T_charge_out, double* T_discharge_in, double* T_discharge_out){
    double ef_dout, ef_din, ef_cout, ef_cin;
    ef_dout = compute_exergy_flux(T_discharge_out, Npoints, dx, m_dot, C_pf);
    ef_din = compute_exergy_flux(T_discharge_in, Npoints, dx, m_dot, C_pf);
    ef_cin = compute_exergy_flux(T_charge_in, Npoints, dx, m_dot, C_pf);
    ef_cout = compute_exergy_flux(T_charge_out, Npoints, dx, m_dot, C_pf);
    return (ef_dout - ef_din) / (ef_cin - ef_cout);
}

/**
    Computes the capacity factor of the storage using the trapezoidal rule for the integrals
 
    @param Tf The temperature of the fluid phase
    @param Ts The temperature of the solid phase
    @param D The diameter of the storage
    @param H The height of the storage
    @param eps The Porosity of the solid phase
    @param rho_f The density of the fluid phase
    @param C_pf The specific heat at constant pressure of the fluid phase
    @param rho_s The density of the fluid phase
    @param C_s The specific heat of the solid phase
    @param Tc The temperature of the inflowing fluid during the charge periods
    @param Td The temperature of the inflowing fluid during the discharge periods
    @param Npoints The number of points for the Finite Volumes discretization
    @param dx The cell width
 
    @return the capacity factor of the storage
 */
double compute_capacity_factor(double* Tf, double* Ts, double D, double H, double eps, double rho_f, double C_pf, double rho_s, double C_s, double Td, double Tc, int Npoints, double dx){
    
    double result1 = 0;
    double result2 = 0;
    double integ1[Npoints], integ2[Npoints];
    
    for(int i = 0; i < Npoints; i++){
        integ1[i] = Tf[i] - Td;
        integ2[i] = Ts[i] - Td;
    }
    
    for(int i = 1; i <= Npoints - 2; i++){
        result1 += 2 * integ1[i];
        result2 += 2 * integ2[i];
    }
    
    result1 += integ1[0] + integ1[Npoints - 1];
    result1 = dx / 2.0 * result1;
    
    result2 += integ2[0] + integ2[Npoints - 1];
    result2 = dx / 2.0 * result2;

    double Q_t = pi / 4 * D * D * (eps * rho_f * C_pf * result1 + (1 - eps) * rho_s * C_s * result2);
    double Q_max_t = (eps * rho_f * C_pf + (1 - eps) * rho_s * C_s) * pi / 4 * D * D * H * (Tc - Td);

    return Q_t / Q_max_t;
}
