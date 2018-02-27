
#include <iostream>
#include <algorithm>
#include <string>
#include <fstream>
#include <cassert>
#include <array>
#include <cmath>
#include <mkl_vsl.h>
#include <omp.h>

#define PI 3.141592653589793

typedef int size_type;

class Diffusion2D {
public:
    	Diffusion2D(const double D,
                const double L,
                const size_type M,
                const double dt,
                const double tot)
    	: D(D), L_(L), M(M),  dt_(dt), m(M-1), tot(tot)
    	{
     		dh = L_ / m;
        	l = D*dt/(dh*dh);
        	rho = new double[M*M];
        	k= new int[M*M];
        	knew = new int[M*M];
        	vslNewStream(&stream, VSL_BRNG_MT19937, 777);
        	initialize_density();
    	}
	    
    	void advance()
    	{

		size_type i;
        	for(i = 0; i<(m+1)*(m+1); i=i++){
                	knew[i] = 0.;
        	}


        	size_type j;
        	for(i = 1; i<m; i++){
            		for(j = 1; j<m-m%4; j=j+4){
                
                		n =k[i*M+j];
                
				if (n == 0 ) continue;

				int moves[4];
				do{
			        	viRngBinomial(VSL_RNG_METHOD_BINOMIAL_BTPE,stream,4,moves,n,l);               
	        	        	stay = n - moves[0]-moves[1]-moves[2]-moves[3];
				}while(stay<0);
	                
	        	        knew[i*M + j] = knew[i*M + j] + stay;
	        	        knew[(i-1)*M + j] = knew[(i-1)*M + j] + moves[0];
	        	        knew[(i+1)*M + j] = knew[(i+1)*M + j] + moves[1];
	        	        knew[i*M + (j-1)] = knew[i*M + (j-1)] + moves[2];
	        	        knew[i*M + (j+1)] = knew[i*M + (j+1)] + moves[3];
	                
               
	        	        n =k[i*M+j+1];
	        	        if (n == 0 ) continue;
	        	        do{
			        	viRngBinomial(VSL_RNG_METHOD_BINOMIAL_BTPE,stream,4,moves,n,l);               
	        	        	stay = n - moves[0]-moves[1]-moves[2]-moves[3];
				}while(stay<0);
	        	        
	        	        knew[i*M + j+1] = knew[i*M + j+1] + stay;
	        	        knew[(i-1)*M + j+1] = knew[(i-1)*M + j+1] + moves[0];
	             	        knew[(i+1)*M + j+1] = knew[(i+1)*M + j+1] + moves[1];
	                	knew[i*M + (j+1-1)] = knew[i*M + (j+1-1)] + moves[2];
	                	knew[i*M + (j+1+1)] = knew[i*M + (j+1+1)] + moves[3];
                
	                	n =k[i*M+j+2];
	                	if (n == 0 ) continue;
	
		                do{
			        	viRngBinomial(VSL_RNG_METHOD_BINOMIAL_BTPE,stream,4,moves,n,l);               
		                	stay = n - moves[0]-moves[1]-moves[2]-moves[3];
				}while(stay<0);
                
		                knew[i*M + j+2] = knew[i*M + j+2] + stay;
		                knew[(i-1)*M + j+2] = knew[(i-1)*M + j+2] + moves[0];
		                knew[(i+1)*M + j+2] = knew[(i+1)*M + j+2] + moves[1];
		                knew[i*M + (j+2-1)] = knew[i*M + (j+2-1)] + moves[2];
		                knew[i*M + (j+2+1)] = knew[i*M + (j+2+1)] + moves[3];

	        	        n =k[i*M+j+3];
	        	        if (n == 0 ) continue;
	
	        	        do{
			        	viRngBinomial(VSL_RNG_METHOD_BINOMIAL_BTPE,stream,4,moves,n,l);               
	        	        	stay = n - moves[0]-moves[1]-moves[2]-moves[3];
				}while(stay<0);
	                
	        	        knew[i*M + j+3] = knew[i*M + j+3] + stay;
	        	        knew[(i-1)*M + j+3] = knew[(i-1)*M + j+3] + moves[0];
	        	        knew[(i+1)*M + j+3] = knew[(i+1)*M + j+3] + moves[1];
	        	        knew[i*M + (j+3-1)] = knew[i*M + (j+3-1)] + moves[2];
	        	        knew[i*M + (j+3+1)] = knew[i*M + (j+3+1)] + moves[3];

	        	}
	        	for(; j<m;j++){
	                
	        	        n =k[i*M+j];
				if (n == 0 ) continue;
				int moves[4];
	        	        do{
			        	viRngBinomial(VSL_RNG_METHOD_BINOMIAL_BTPE,stream,4,moves,n,l);               
	        	        	stay = n - moves[0]-moves[1]-moves[2]-moves[3];
				}while(stay<0);
                

                
			        knew[i*M + j] = knew[i*M + j] + stay;
	        	    	knew[(i-1)*M + j] = knew[(i-1)*M + j] + moves[0];
	        	    	knew[(i+1)*M + j] = knew[(i+1)*M + j] + moves[1];
	        	    	knew[i*M + (j-1)] = knew[i*M + (j-1)] + moves[2];
	        	    	knew[i*M + (j+1)] = knew[i*M + (j+1)] + moves[3];
            		}
        	}
	 		for(size_type i = 0; i<m+1; i++){
	        	knew[0*M + i] = 0;
	        	knew[M*m + i] = 0;
	        	knew[i*M + 0] = 0;
			knew[i*M + m] = 0;
		}
        
        	std::swap(k,knew);
	}

	void compute_density()
	{
	        for(size_type i = 1; i<m; i++){
			for(size_type j = 1; j<m; j++){
		                rho[i*M + j] = k[i*M + j]*integ/(tot*dh*dh);
			}
	        }
	}


     
	void write_density(std::string const& filename) const
	{
		std::ofstream out_file(filename, std::ios::out);
	        
	        for(size_type i = 0; i < M; ++i) {
	            for(size_type j = 0; j < M; ++j)
	                out_file << (i*dh - L_/2.) << '\t' << (j*dh - L_/2.) << '\t' << rho[i*M + j] << "\n";
	            out_file << "\n";
	        }
	        out_file.close();
	}
	    
	void write_density_diff(std::string const& filename, double current_time) const
	{
	        std::ofstream out_file(filename, std::ios::out);
	        
	        for (size_type i = 0; i < M; ++i) {
	            for (size_type j = 0; j < M; ++j) {
	   	             out_file << (i*dh - L_ / 2.) << '\t' << (j*dh - L_ / 2.) << '\t' << rho[i*M + j] - sin(PI*i*dh) * sin(PI*j*dh) * std::exp(-2 * D*PI*PI*current_time) << "\n";
          		     out_file << "\n";
            		}
        	}	
       		out_file.close();
    	}
	    
    	double variance(double current_time) 
	{
        	double error = 0;
                for (size_type i = 1; i < M-1; ++i){
            		for(size_type j = 1; j < M-1; ++j){
                		double numerical =rho[i*M+j];
                		double analytic  = sin(PI*i*dh) * sin(PI*j*dh) * std::exp(-2*D*PI*PI*current_time);
                		error += k[i*M + j]*(numerical - analytic)*(numerical-analytic);
                        }
        	}
        	return error;
    	}
    
private:
    
    	void initialize_density()
    	{
    		for(size_type i = 0;i<m+1; i++){
            		for(size_type j = 0;j<m+1; j++){
                		k[i*M+j] = floor(sin(PI*i*dh)*sin(PI*j*dh)*tot*dh*dh/integ);
            		}
        	}
    	}
    
    
    	double D, L_, tot,l;
    	size_type M, m;
    	double dh, dt_, k_;
    	double *rho;
    	int *k,*knew;
    	int stay, n;        	


        	
    	const double integ = (2/PI)*(2/PI);

	
    	VSLStreamStatePtr stream;

    
};


int main(int argc, char* argv[])
{
	if (argc !=4) {
     		std::cerr << "Usage: " << argv[0] << " m (grid intervals) Nt (timesteps) Np (number of particles)" << std::endl;
     		return 1;
    	}
     
    	const double D  = 1;
    	const double L  = 1;
    	const size_type m = std::stod(argv[1]);
    	const size_type Nt = std::stod(argv[2]);
    	const size_type  M  = m+1;
    	const double l = 0.15;
    	const double dh = L/m;
    	const double dt = l*dh*dh/D;
    	std::cout<<"dt= "<<dt<<std::endl;
    	const double tmax = Nt*dt;
    	const double tot = std::stod(argv[3]); 

    
	Diffusion2D system(D, L, M, dt,tot);
    
    	double time = 0;
    	double var =0;
    	double norm_stdev;
    	double t1,t2;
    	t1 = omp_get_wtime();
    	while (time < tmax) {
    	    system.advance();
    	    time += dt;
        }
    	t2 = omp_get_wtime();
    	system.compute_density();
    
    	var = system.variance(time);
    	
    	norm_stdev = sqrt(var/tot)/(m*m);
    	std::cout <<"Error: "<< norm_stdev << std::endl;
    
    
    	std::cout << "Timing : "<< t2-t1 << std::endl;
    
    
    	return 0;
}












