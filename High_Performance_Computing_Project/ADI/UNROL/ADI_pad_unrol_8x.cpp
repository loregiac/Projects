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
            if(M_padded%4) M_padded = M_padded + 4 - (M_padded%4);
            std::cout<<M_padded<<std::endl;

            k_ = dt_ * D_ / (dh_ * dh_);


            k_2 = k_/2;

            kb = 1 + k_;
            kbd = 1 / kb;

            ka = -k_ / 2;
            kc = -k_ / 2;
            nip = M_ - 2;


            k_ = 1 - k_;

            b = new double[8*nip];
            c = new double[nip];

            rho_ = (double*) _mm_malloc(M_*M_padded* sizeof(double), 32);
            rho_mid = (double*) _mm_malloc(M_*M_padded * sizeof(double), 32);
            initialize_density();
        }
        ~Diffusion2D()
        {
            delete(b);
            delete(c);
            _mm_free(rho_);
            _mm_free(rho_mid);
        }



        void advance()
        {
            size_type n = nip - 1;




            for (j=8; j < m; j=j+8) {

                for ( i = 1; i < m; i++) {
                    b[8*i - 8] = k_2 *rho_[i*M_padded + (j - 8)] + ( k_)*rho_[i*M_padded + j-7] + k_2 * rho_[i*M_padded + (j - 6)];
                    b[8*i - 7] = k_2 * rho_[i*M_padded + (j - 7)] + ( k_)*rho_[i*M_padded + j-6] + k_2 * rho_[i*M_padded + (j - 5)];
                    b[8*i - 6] = k_2 * rho_[i*M_padded + (j - 6)] + ( k_)*rho_[i*M_padded + j-5] + k_2 * rho_[i*M_padded + (j - 4)];
                    b[8*i - 5] = k_2 * rho_[i*M_padded + (j - 5)] + ( k_)*rho_[i*M_padded + j-4] + k_2 * rho_[i*M_padded + (j-3)];
                    b[8*i - 4] = k_2 *rho_[i*M_padded + (j - 4)] + ( k_)*rho_[i*M_padded + j-3] + k_2 * rho_[i*M_padded + (j - 2)];
                    b[8*i - 3] = k_2 * rho_[i*M_padded + (j - 3)] + ( k_)*rho_[i*M_padded + j-2] + k_2 * rho_[i*M_padded + (j - 1)];
                    b[8*i - 2] = k_2 * rho_[i*M_padded + (j - 2)] + ( k_)*rho_[i*M_padded + j-1] + k_2 * rho_[i*M_padded + (j )];
                    b[8*i - 1] = k_2 * rho_[i*M_padded + (j - 1)] + ( k_)*rho_[i*M_padded + j] + k_2 * rho_[i*M_padded + (j+1)];

                    rho_mid[(i )*M_padded + j-7 ] = (b[8*i-8] - ka*rho_mid[(i-1)*M_padded + j-7 ])*denom[i-1];
                    rho_mid[(i )*M_padded + j-6 ] = (b[8*i-7] - ka*rho_mid[(i-1)*M_padded + j-6 ])*denom[i-1];
                    rho_mid[(i)*M_padded + j-5 ] = (b[8*i-6] - ka*rho_mid[(i-1)*M_padded + j-5 ])*denom[i-1];
                    rho_mid[(i )*M_padded + j-4 ]=(b[8*i-5] - ka*rho_mid[(i-1)*M_padded + j-4 ])*denom[i-1];
                    rho_mid[(i )*M_padded + j-3 ] = (b[8*i-4] - ka*rho_mid[(i-1)*M_padded + j-3 ])*denom[i-1];
                    rho_mid[(i )*M_padded + j-2 ] = (b[8*i-3] - ka*rho_mid[(i-1)*M_padded + j-2 ])*denom[i-1];
                    rho_mid[(i)*M_padded + j-1 ] = (b[8*i-2] - ka*rho_mid[(i-1)*M_padded + j-1 ])*denom[i-1];
                    rho_mid[(i )*M_padded + j ] = (b[8*i-1] - ka*rho_mid[(i-1)*M_padded + j ])*denom[i-1];

                }


                b[8*nip - 8] = b[8*nip - 8] + k_2 * rho_[m*M_padded + j-7];
                b[8*nip - 7] = b[8*nip - 7] + k_2 * rho_[m*M_padded + j-6];
                b[8*nip - 6] = b[8*nip - 6] + k_2 * rho_[m*M_padded + j-5];
                b[8*nip - 5] = b[8*nip - 5] + k_2 * rho_[m*M_padded + j-4];
                b[8*nip - 4] = b[8*nip - 4] + k_2 * rho_[m*M_padded + j-3];
                b[8*nip - 3] = b[8*nip - 3] + k_2 * rho_[m*M_padded + j-2];
                b[8*nip - 2] = b[8*nip - 2] + k_2 * rho_[m*M_padded + j-1];
                b[8*nip - 1] = b[8*nip - 1] + k_2 * rho_[m*M_padded + j ];



                rho_mid[(n + 1)*M_padded + j - 7 ] = (b[8*nip-8] - ka*rho_mid[n*M_padded + j - 7]) *denom[n];
                rho_mid[(n + 1)*M_padded + j - 6 ] = (b[8*nip-7] - ka*rho_mid[n*M_padded + j - 6]) *denom[n];
                rho_mid[(n + 1)*M_padded + j - 5 ] = (b[8*nip-6] - ka*rho_mid[n*M_padded + j - 5]) *denom[n];
                rho_mid[(n + 1)*M_padded + j - 4 ] = (b[8*nip-5] - ka*rho_mid[n*M_padded + j - 4]) *denom[n];
                rho_mid[(n + 1)*M_padded + j - 3 ] = (b[8*nip-4] - ka*rho_mid[n*M_padded + j - 3]) *denom[n];
                rho_mid[(n + 1)*M_padded + j - 2 ] = (b[8*nip-3] - ka*rho_mid[n*M_padded + j - 2]) *denom[n];
                rho_mid[(n + 1)*M_padded + j - 1 ] = (b[8*nip-2] - ka*rho_mid[n*M_padded + j - 1]) *denom[n];
                rho_mid[(n + 1)*M_padded + j ] = (b[8*nip-1] - ka*rho_mid[n*M_padded + j - 0]) *denom[n];

                for ( i = n; i-- > 0;) {
                    rho_mid[(i + 1)*M_padded + j - 7] = rho_mid[(i + 1)*M_padded + j - 7] -cc[i] * rho_mid[(i + 2)*M_padded + j - 7];
                    rho_mid[(i + 1)*M_padded + j - 6] = rho_mid[(i + 1)*M_padded + j - 6] -cc[i] * rho_mid[(i + 2)*M_padded + j - 6];
                    rho_mid[(i + 1)*M_padded + j - 5] = rho_mid[(i + 1)*M_padded + j - 5] -cc[i] * rho_mid[(i + 2)*M_padded + j - 5];
                    rho_mid[(i + 1)*M_padded + j - 4] = rho_mid[(i + 1)*M_padded + j - 4] -cc[i] * rho_mid[(i + 2)*M_padded + j - 4];
                    rho_mid[(i + 1)*M_padded + j - 3] = rho_mid[(i + 1)*M_padded + j - 3] -cc[i] * rho_mid[(i + 2)*M_padded + j - 3];
                    rho_mid[(i + 1)*M_padded + j - 2] = rho_mid[(i + 1)*M_padded + j - 2] -cc[i] * rho_mid[(i + 2)*M_padded + j - 2];
                    rho_mid[(i + 1)*M_padded + j - 1] = rho_mid[(i + 1)*M_padded + j - 1] -cc[i] * rho_mid[(i + 2)*M_padded + j - 1];
                    rho_mid[(i + 1)*M_padded + j ] = rho_mid[(i + 1)*M_padded + j ] -cc[i] * rho_mid[(i + 2)*M_padded + j ];
                }
            }
            j = j - 7;
            for(; j < m; j++){
                for ( i = 1; i < m; i++) {
                    b[i - 1] = k_2 * rho_[i*M_padded+ (j - 1)] + (k_)*rho_[i*M_padded + j] + k_2 * rho_[i*M_padded + (j + 1)];
                    rho_mid[(i)*M_padded + j ] = (b[i-1] - ka*rho_mid[(i-1)*M_padded + j ])*denom[i-1];
                }
                b[nip-1]=b[nip-1] + k_2 * rho_[0*M_padded + j];

                rho_mid[(n + 1)*M_padded + j ] = (b[n] - ka*rho_mid[n*M_padded + j ]) *denom[n];
                for ( i = n; i-- > 0;) {
                    rho_mid[(i + 1)*M_padded + j ] = rho_mid[(i + 1)*M_padded + j] -cc[i] * rho_mid[(i + 2)*M_padded + j ];
                }
            }





            for (i = 1; i < m; ++i) {

                for (j = 1; j < m; j ++) {
                    c[j - 1] = k_2 * rho_mid[(i - 1)*M_padded + j] + (k_)*rho_mid[i*M_padded + j] + k_2 * rho_mid[(i + 1)*M_padded + j];
                    rho_[i*M_padded + (j)] = (c[j-1] - ka*rho_[i*M_padded + j-1]) *denom[j-1];
                }




                c[m - 2] = c[m - 2] + k_2 * rho_mid[i*M_padded + m];
                rho_[i*M_padded + n + 1] = (c[n] - ka*rho_[i*M_padded + n]) *denom[n];

                for (ii = n; ii -- > 0; ) {
                    rho_[i*M_padded + (ii + 1)] -= cc[ii] * rho_[i*M_padded + ii + 2];
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
        size_type M_, Mtot, m, nip, j, i, ii;
        size_type M_padded;
        double dh_, dt_, k_, k_2;

        double ka, kb, kc, kbd;

        double *cc, *denom;

        double *rho_, *rho_mid;

        double *b, *c;


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
    t.start();
    system.thomas_LU();
    while (time < tmax) {
        system.advance();
        time += dt;
    }


    t.stop();

    std::cout << "Timing : " << t.get_timing() << std::endl;
    curr_err = system.linf_error(time);
    max_err = std::max(max_err, curr_err);

    std::cout << "err: " << max_err << std::endl;


    return 0;
}
