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
 
    @param height The height of the storage.
    @param d The diameter of the storage
    @param Ncells The number of cells for the Finite Volumes discretization
    @param Tc The temperature of the inflowing fluid during the charge periods
    @param Td The temperature of the inflowing fluid during the discharge periods
    @param dx The size of each cell
    @param D1 The duration of a charge period
    @param D2 The duration of a idle period between charging and discharging
    @param D3 The duration of a discharge period
    @param D4 The duration of a idle period between discharging and charging
    @param Ncycles The number of charging-idle-discharging-idle cycles
    @param Nt The number of time steps per cycle
    @param alpha_f The diffusion coefficient of the fluid phase
    @param alpha_s The diffusion coefficient of the solid phase
    @param u_f The advection speed of the fluid phase
    @param h_vf The coupling coefficient of the fluid phase
    @param h_vs The coupling coefficient of the solid phase
    @param mdot_f The mass flow rate of the fluid phase
    @param C_pf The specific heat at constant pressure of the fluid phase
    @param C_s The specific heat of the solid phase
    @param eps The Porosity of the solid phase
    @param rho_f The density of the fluid phase
    @param rho_s The density of the solid phase
    @param file The flag that indicates whether to read the parameters from the file "params.dat" or to ask them in the terminal
 */
void init(double* height, double* d, unsigned int* Ncells, double* Tc, double* Td, double* dx, double* D1, double* D2, double* D3, double* D4, unsigned int* Ncycles, unsigned int* Nt, double* alpha_f, double* alpha_s, double* u_f, double *h_vf, double *h_vs, double* mdot_f, double* C_pf, double* C_s, double* eps, double* rho_f, double* rho_s, std::string file){
    
    double k_s;
    double k_f;
    double mu_f;
    double d_s;
    double V;
    
    if(file.compare("file") && file.compare("terminal")){
        std::cerr << "Please insert either file or terminal option"<<std::endl;
        abort();
    }
    
    if(!file.compare("terminal")){                      //if chosen option is to read from terminal
        std::cout << "Please enter the volume of the storage:" << std::endl;
        std::cin >> V;
        *height = 4 * V / (pi * (*d) * (*d));
        
        std::cout << "Please enter the diameter of the storage:" << std::endl;
        std::cin >> *d;
        
        std::cout << "Please enter the number of cells:" << std::endl;
        std::cin >> *Ncells;
    
        *dx = (*height) / (*Ncells);
     
        std::cout << "Please enter the duration of the charge periods:" << std::endl;
        std::cin >> *D1;
     
        std::cout << "Please enter the duration of idle periods between charging and discharging:" << std::endl;
        std::cin >> *D2;
     
        std::cout << "Please enter the duration of the discharge periods:" << std::endl;
        std::cin >> *D3;
     
        std::cout << "Please enter the duration of idle state between discharging and charging:" << std::endl;
        std::cin >> *D4;
     
        std::cout << "Please enter the number of cycles:" << std::endl;
        std::cin >> *Ncycles;
     
        std::cout << "Please enter the number of time steps per cycle:" << std::endl;
        std::cin >> *Nt;
     
        std::cout << "Please enter alpha_f:" << std::endl;
        std::cin >> *alpha_f;
     
        std::cout << "Please enter alpha_s:" << std::endl;
        std::cin >> *alpha_s;
        
        std::cout << "Please enter k_f:" << std::endl;
        std::cin >> k_f;
        
        std::cout << "Please enter k_s:" << std::endl;
        std::cin >> k_s;
        
        std::cout << "Please enter mu_f:" << std::endl;
        std::cin >> mu_f;
        
        std::cout << "Please enter d_s:" << std::endl;
        std::cin >> d_s;
        
        std::cout << "Please enter Tc:" << std::endl;
        std::cin >> *Tc;
        
        std::cout << "Please enter Td:" << std::endl;
        std::cin >> *Td;
        
        std::cout << "Please enter eps:" << std::endl;
        std::cin >> *eps;
        
        std::cout << "Please enter rho_f:" << std::endl;
        std::cin >> *rho_f;
        
        std::cout << "Please enter rho_s:" << std::endl;
        std::cin >> *rho_s;
    }else{                                //if chosen option is to read from file
        
        //open the file "params.dat"
        std::string line;
        std::ifstream myfile ("params.dat");

    
        if(myfile.is_open()){                           //if the file has been opened correctly
            getline(myfile,line);
            V = std::stod(line);
            getline(myfile,line);
            *d = std::stod(line);
            *height = 4 * V / (pi * (*d) * (*d));
            *Ncells = 2000;
            *dx = *height / (*Ncells);
            getline(myfile,line);
            *Tc = std::stod(line);
            getline(myfile,line);
            *Td = std::stod(line);
            getline(myfile,line);
            *D1 = std::stod(line);
            getline(myfile,line);
            *D2 = std::stod(line);
            getline(myfile,line);
            *D3 = std::stod(line);
            getline(myfile,line);
            *D4 = std::stod(line);
            getline(myfile,line);
            *Ncycles = std::stoi(line);
            getline(myfile,line);
            *Nt = std::stoi(line);
            getline(myfile,line);
            k_f = std::stod(line);
            getline(myfile,line);
            k_s = std::stod(line);
            getline(myfile,line);
            *eps = std::stod(line);
            getline(myfile,line);
            d_s = std::stod(line);
            getline(myfile,line);
            *rho_s = std::stod(line);
            getline(myfile,line);
            *rho_f = std::stod(line);
            getline(myfile,line);
            *C_s = std::stod(line);
            getline(myfile,line);
            *C_pf = std::stod(line);
            getline(myfile,line);
            mu_f = std::stod(line);
            getline(myfile,line);
            *mdot_f = std::stod(line);
        }else{                                  //if the program has encountered problems in opening the file
            std::cout << "Unable to open file";
        }
    }
    
    //compute the remaining parameters
    *u_f = (*mdot_f) / ((*rho_f) * pi * (*d) * (*d) / 4 * (*eps));                  //advection speed
    double Re = (*eps) * (*rho_f) * (*u_f) * (d_s) / mu_f;                      //Reynolds number
    double Pr = mu_f * (*C_pf) / k_f;                                       //Prandtl number
    double Nu_fs = 0.255 / (*eps) * pow(Pr,1.0/3) * pow(Re,2.0/3);          //Nusselt number
    double h_fs = Nu_fs * k_f / (d_s);                                      //fluid-solid heat transfer
    double h = 1.0 / (1.0 / h_fs + (d_s) / (10.0 * k_s));                         //overall heat transfer
    double h_v = 6 * (1 - (*eps)) * h / (d_s);                                  //volumetric heat transfer
    *h_vf = h_v / ((*eps) * (*rho_f) * (*C_pf));
    *h_vs = h_v / ((1 - (*eps)) * (*rho_s) * (*C_s));
    *alpha_f = k_f / ((*eps) * (*rho_f) * (*C_pf));
    *alpha_s = k_f / ((1 - (*eps)) * (*rho_s) * (*C_s));
}
