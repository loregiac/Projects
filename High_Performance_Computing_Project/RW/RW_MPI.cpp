#include <iostream>
#include <algorithm>
#include <string>
#include <fstream>
#include <cassert>
#include <array>
#include <cmath>
#include <mkl_vsl.h>
#include <omp.h>
#include <mpi.h>
#include <unistd.h>


#define PI 3.141592653589793
#define TAG 0

typedef int size_type;


struct position{
	double x,y;
};

class Diffusion2D {
public:
    	Diffusion2D(const double D,
                	const double L,
                	const size_type M,
                	const double dt,
                	const double tot, 
			const int size,
			const int rank)
    	: D(D), L_(L), M(M),  dt_(dt), m(M-1), tot(tot), size(size), rank(rank)
    	{	
		M_glob = M;
 		Mx_loc = M_glob;
        	My_loc = M_glob / size;
        	dh = L_ / m;
       		int res; 
        	l = D*dt/(dh*dh);
        	rho = new double[Mx_loc*(My_loc)];
        	k= new int[Mx_loc*(My_loc+2)];
        	knew = new int[Mx_loc*(My_loc+2)];
		
        	vslNewStream( &stream, VSL_BRNG_MT2203+rank, 777 );

		

		xmin_loc = 0;
		xmax_loc = xmin_loc + (Mx_loc) * dh;
		ymin_loc = rank * (My_loc) * dh;
		ymax_loc = ymin_loc + (My_loc) * dh;

		MPI_Type_contiguous(Mx_loc-2,MPI_INT,&top_boundary);
		MPI_Type_commit(&top_boundary);

        	MPI_Type_contiguous(Mx_loc-2, MPI_INT, &bottom_boundary);
		MPI_Type_commit(&bottom_boundary);
	        initialize_density();
	
}    

	void advance()
    	{	
		MPI_Request request[4];
                MPI_Status status[4];

        	for(size_type i = 0; i<My_loc; i++){
            		for(size_type j = 0; j<Mx_loc; j++){
                		knew[i*Mx_loc + j] = 0;
            		}
        	}
        	int n;
                
        	for(size_type i = 1 + (rank==0); i<My_loc+1 - (rank==(size-1)); i++){
            		for(size_type j = 1; j<Mx_loc-1; j++){
                		n =k[i*Mx_loc+j];
				int moves[5];
				if (n == 0 ){
					moves[0]= 0;
					moves[1]= 0;
					moves[2]= 0;
					moves[3]= 0;
					moves[4]= 0;
					continue;
				}
		
		
	        		do {
					viRngBinomial(VSL_RNG_METHOD_BINOMIAL_BTPE,stream,4,moves,n,l);
	                		moves[4] = n - moves[0]-moves[1]-moves[2]-moves[3];
				}while (moves[4]< 0)	;
                                knew[i*Mx_loc + j] = knew[i*Mx_loc + j] + moves[4];
                		knew[(i-1)*Mx_loc + j] = knew[(i-1)*Mx_loc + j] + moves[0];
                		knew[(i+1)*Mx_loc + j] = knew[(i+1)*Mx_loc + j] + moves[1];
               		 	knew[i*Mx_loc + (j-1)] = knew[i*Mx_loc + (j-1)] + moves[2];
                		knew[i*Mx_loc + (j+1)] = knew[i*Mx_loc + (j+1)] + moves[3];
            		}
        	}


		int add_bottom[Mx_loc-2],add_top[Mx_loc-2];
		
        	
		//Commmunication
		int next_rank = rank+1;
		int prev_rank = rank-1;

		if (prev_rank >= 0) {
			MPI_Irecv(&add_top[0], 1, top_boundary, prev_rank, 100, MPI_COMM_WORLD, &request[0]);
			MPI_Isend(&knew[1], 1, top_boundary, prev_rank, 100, MPI_COMM_WORLD, &request[1]);
		} else {
			request[0] = MPI_REQUEST_NULL;
			request[1] = MPI_REQUEST_NULL;
		}
	
		if (next_rank < size) {
			MPI_Irecv(&add_bottom[0], 1, bottom_boundary, next_rank, 100, MPI_COMM_WORLD, &request[2]);
			MPI_Isend(&knew[(My_loc+1)*Mx_loc+1], 1, bottom_boundary, next_rank, 100, MPI_COMM_WORLD, &request[3]);
		}else{
			request[2] = MPI_REQUEST_NULL;
			request[3] = MPI_REQUEST_NULL;
		}
		MPI_Waitall(4,request,status);

	

		if(prev_rank>=0){
		        for(size_type j = 1; j<Mx_loc-1; j++){
            			knew[Mx_loc + j] += add_top[j-1];
        		}
		}	

		if(next_rank<size){
		        for(size_type j = 1; j<Mx_loc-1; j++){
            			knew[Mx_loc*(My_loc) + j] += add_bottom[j-1];
        		}	
		}

		for(size_type j = 0; j<Mx_loc; j++){
			knew[(0+rank==0)*Mx_loc+j] = 0;
			knew[(My_loc + 1 -(rank==(size-1)))*Mx_loc+j] = 0;
		}
		
        	for(size_type i = 0; i<My_loc+1; i++){
            		knew[i*Mx_loc + 0] = 0;
            		knew[i*Mx_loc + (Mx_loc-1)] = 0;
        	}
		
		std::swap(k,knew);


    	}

	void compute_density()
	{
        	for(size_type i = 1; i<My_loc+1; i++){
            		for(size_type j = 0; j<Mx_loc; j++){
                		rho[(i-1)*Mx_loc + j] = k[i*Mx_loc + j]*integ/(tot*dh*dh);
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
		double sum_error,x,y;
        	for (size_type i = 0; i < My_loc; ++i){
            		for(size_type j = 0; j < Mx_loc; ++j){
				x = j*dh;
				y = (i+rank*My_loc);
                		double numerical =rho[i*Mx_loc+j];
                		double analytic  = sin(PI*x) * sin(PI*y) * std::exp(-2*D*PI*PI*current_time);
                		error += k[(i+1)*Mx_loc + j]*(numerical - analytic)*(numerical-analytic);
               		}
        	}
		MPI_Allreduce(&error,&sum_error,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        	return sum_error;
    	}
    
private:
    
	void initialize_density()
    	{
        	//initial condition in the internal part
		position current;
        	for(size_type i = 1;i<My_loc+1; i++){
            		for(size_type j = 1;j<Mx_loc-1; j++){
			current = get_position(i,j);
                	k[i*Mx_loc+j] = floor(sin(PI*current.x)*sin(PI*current.y)*tot*dh*dh/integ);
            		}
        	}
		//initialize the ghost cells with zeros

		
    	}
    
    	inline position get_position(size_type i, size_type j) const
	{
	    	position p;
	    	p.x = xmin_loc + j*dh;
	    	p.y = (rank*My_loc + (i-1))*dh;
		return p;
	}
    	double D, L_, tot,l;
    	size_type M, m, Mx_loc, My_loc, M_glob;
	int size, rank;
    	double dh, dt_, k_;
    	double *rho;
    	int *k,*knew;



    	const double integ = (2/PI)*(2/PI);
	uint32_t curr_a = 3454325;
	double xmin_loc,ymin_loc;
    	double xmax_loc, ymax_loc;	

  	VSLStreamStatePtr stream;

	MPI_Comm cart_comm;
    	MPI_Datatype top_boundary, bottom_boundary;
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
    	
    	const double tmax = Nt*dt;
    	const double tot = std::stod(argv[3]);
    
	MPI_Init(&argc,&argv);
	
    	int rank, size;
    	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    	MPI_Comm_size(MPI_COMM_WORLD, &size);




	const size_type M_world = (M-2) % size == 0 ? M : M + (size - (M-2) % size);

    	Diffusion2D system(D, L, M_world, dt,tot,size,rank);
    

    
    	double time = 0;
    	double var =0;
        double norm_stdev;
    	double t1,t2;

    	t1 = MPI_Wtime();
    	while (time < tmax) {
        	system.advance();
        	time += dt;
    	}
	t2 = MPI_Wtime();

    	var = system.variance(time);

        norm_stdev = sqrt(var/tot)/(m*m);


	MPI_Finalize();
	std::cout <<"Error: "<< norm_stdev << std::endl;
    	double timing = t2-t1;
    	
    	std::cout << "Timing : "<< timing << std::endl;        
    	return 0;
}













