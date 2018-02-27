# 2 "ADI_omp_simple.cpp"
#include <iostream>
#include <algorithm>
#include <string>
#include <fstream>
#include <cassert>
#include <vector>
#include <cmath>
#include <malloc.h>
#include <immintrin.h>
#include "timer.hpp"
#include <omp.h>
#define PI 3.141592653589793

typedef int size_type;

class Diffusion2D {
    public:
        Diffusion2D(const double D,
                const double L,
                const size_type M,
                const double dt)
            : D_(D), L_(L), M_(M), Mtot(M_*M_), dt_(dt), m(M_ - 1)
        {

            dh_ = L_ / m;
            M_padded = M_;

            std::cout<<M_padded<<std::endl;

            k_ = dt_ * D_ / (dh_ * dh_);


            k_2 = k_/2;

            kb = 1 + k_;
            kbd = 1 / kb;

            ka = -k_ / 2;
            kc = -k_ / 2;
            nip = M_ - 2;


            k_ = 1 - k_;



            rho_ = (double*) _mm_malloc(M_*M_padded* sizeof(double), 32);
            rho_mid = (double*) _mm_malloc(M_*M_padded * sizeof(double), 32);
            initialize_density();
        }
        ~Diffusion2D()
        {
            _mm_free(rho_);
            _mm_free(rho_mid);
        }



        void advance(double *b, double *c)
        {
            size_type n = nip - 1;



#pragma omp for schedule (dynamic)
            for (size_t j=1; j < m; j++) {

                rho_mid[(0)*M_ + j ] = rho_[0*M_ + (j)];
                for (size_t i = 1; i < m; i++) {
                    b[i - 1] = k_2 * rho_[i*M_ + (j - 1)] + (k_)*rho_[i*M_ + j] + k_2 * rho_[i*M_ + (j + 1)];
                    rho_mid[(i)*M_ + j ] = (b[i-1] - ka*rho_mid[(i-1)*M_ + j ])*denom[i-1];
                }

                b[nip - 1] = b[nip - 1] + k_2 * rho_[m*M_ + j ];
                rho_mid[(n + 1)*M_ + j ] = (b[n] - ka*rho_mid[n*M_ + j ]) *denom[n];

                for (size_t i = n; i-- > 0;) {
                    rho_mid[(i + 1)*M_ + j ] = rho_mid[(i + 1)*M_ + j] -cc[i] * rho_mid[(i + 2)*M_ + j ];
                }


            }






#pragma omp for schedule (dynamic)
            for (size_t i = 1; i < m; ++i) {

                rho_[(i)*M_ +0 ] = rho_mid[i*M_ + 0];
                for (size_t j = 1; j < m; j ++) {
                    c[j-1 ] = k_2 * rho_mid[(i - 1)*M_ + j] + (k_)*rho_mid[i*M_ + j] + k_2 * rho_mid[(i + 1)*M_ + j];
                    rho_[i*M_ + (j)] = (c[j-1] - ka*rho_[i*M_ + j-1 ]) *denom[j-1];
                }

                c[m - 2] = c[m - 2] + k_2 * rho_mid[i*M_ + m];
                rho_[i*M_ + n + 1] = (c[n] - ka*rho_[i*M_ + n]) *denom[n];


                for (size_t ii = n; ii-- > 0; ) {
                    rho_[i*M_ + (ii + 1)] -= cc[ii] * rho_[i*M_ + ii + 2];
                }
            }

        }





        void write_density(std::string const& filename) const
        {
            std::ofstream out_file(filename, std::ios::out);

            for (size_type i = 0; i < M_; ++i) {
                for (size_type j = 0; j < M_; ++j)
                    out_file << (i*dh_) << '\t' << (j*dh_) << '\t' << rho_[i*M_padded + j] << "\n";
                out_file << "\n";
            }
            out_file.close();
        }

        double linf_error(double current_time) {
            double max_error = 0;
            double current_error = 0;

            for (size_type i = 1; i < M_ - 1; ++i) {
                for (size_type j = 1; j < M_ - 1; ++j) {
                    double numerical = rho_[i*M_padded + j];
                    double analytic = sin(PI*i*dh_) * sin(PI*j*dh_) * std::exp(-2 * D_*PI*PI*current_time);
                    current_error = fabs(numerical - analytic);
                    max_error = std::max(max_error, current_error);
                }
            }
            return max_error;
        }
        void thomas_LU() {
            cc = new double[nip];
            denom = new double[nip];

            cc[0] = kc *kbd;
            denom[0] = 1 / (kb);
            for (int i = 1; i < nip; i++) {
                cc[i] = kc / (kb - ka*cc[i - 1]);
                denom[i] = 1 / (kb - ka*cc[i - 1]);
            }
            cc[nip-1] = kc;
        }

    private:




        void initialize_density()
        {


            for (size_type i = 1; i < m; i++) {
                for (size_type j = 1; j < m; j++) {
                    rho_[i*M_padded + j] = sin(PI*i*dh_)*sin(PI*j*dh_);
                }
            }
        }


        double D_, L_;
        size_type M_, Mtot, m, nip;
        size_type M_padded;
        double dh_, dt_, k_, k_2;

        double ka, kb, kc, kbd;

        double *cc, *denom;

        double *rho_, *rho_mid;



};


int main(int argc, char* argv[])
{
    const double D = 1;
    const double L = 1;
    const size_type m = std::stod(argv[1]);
    const size_type M = m + 1;
    const double dt = std::stod(argv[2]);



    Diffusion2D system(D, L, M, dt);



    double time = 0;
    double tmax = 0.1;
    double max_err = 0;
    double curr_err;
    timer t;
    system.thomas_LU();
#pragma omp parallel
    {
        double *b = (double*) _mm_malloc(M*sizeof(double), 64);
        double *c = (double*) _mm_malloc(M*sizeof(double), 64);
#pragma omp single
        {t.start();}

        while (time < tmax) {
            system.advance(b,c);
#pragma omp single
            time += dt;
        }

#pragma omp single
        {
            t.stop();
        }
        _mm_free(b);
        _mm_free(c);
    }
    std::cout << "Timing : " << t.get_timing() << std::endl;
    curr_err = system.linf_error(time);
    max_err = std::max(max_err, curr_err);

    std::cout << "err: " << max_err << std::endl;


    return 0;
}
