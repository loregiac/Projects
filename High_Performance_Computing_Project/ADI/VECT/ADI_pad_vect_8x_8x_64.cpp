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

            b = (double*) _mm_malloc(8*nip*sizeof(double), 64);
            c = new double[8*nip];

            rho_ = (double*) _mm_malloc((M_*M_padded+3)* sizeof(double), 64);
            rho_ = rho_ + 3;
            rho_mid = (double*) _mm_malloc((M_*M_padded+3) * sizeof(double), 64);
            rho_mid = rho_mid + 3;
            initialize_density();
        }
        ~Diffusion2D()
        {
            _mm_free(b);
            delete(c);
            rho_ = rho_-3;
            rho_mid = rho_mid - 3;
            _mm_free(rho_);
            _mm_free(rho_mid);
        }



        void advance()
        {
            size_type n = nip - 1;



            __m256d k2p = _mm256_broadcast_sd(&k_2);
            __m256d kp = _mm256_broadcast_sd(&k_);
            __m256d kap = _mm256_broadcast_sd(&ka);

            for (j=8; j < m; j=j+8) {

                for ( i = 1; i < m; i++) {
                    __m256d rho1 = _mm256_loadu_pd(rho_ + i*M_padded + (j-8));
                    __m256d rho2 = _mm256_load_pd(rho_ + i*M_padded + (j-7));
                    __m256d rho3 = _mm256_loadu_pd(rho_ + i*M_padded + (j-6));
                    __m256d b1 = _mm256_mul_pd(k2p,rho1);
                    __m256d b2 = _mm256_fmadd_pd(kp,rho2,b1);
                    __m256d b3 = _mm256_fmadd_pd(k2p,rho3,b2);
                    _mm256_store_pd(b + 8*i-8,b3);

                    __m256d rho4 = _mm256_loadu_pd(rho_ + i*M_padded + (j-4));
                    __m256d rho5 = _mm256_load_pd(rho_ + i*M_padded + (j-3));
                    __m256d rho6 = _mm256_loadu_pd(rho_ + i*M_padded + (j-2));
                    __m256d b4 = _mm256_mul_pd(k2p,rho4);
                    __m256d b5 = _mm256_fmadd_pd(kp,rho5,b4);
                    __m256d b6 = _mm256_fmadd_pd(k2p,rho6,b5);
                    _mm256_store_pd(b + 8*i-4,b6);



                    __m256d rho7 = _mm256_load_pd(rho_mid + (i-1)*M_padded + (j-7));
                    __m256d rho8 = _mm256_load_pd(rho_mid + (i-1)*M_padded + (j-3));
                    __m256d denomp = _mm256_broadcast_sd(denom+i-1);

                    __m256d temp = _mm256_fnmadd_pd(kap,rho7, b3);
                    __m256d res = _mm256_mul_pd(temp, denomp);

                    _mm256_store_pd(rho_mid + (i)*M_padded + j-7,res);



                    __m256d temp2 = _mm256_fnmadd_pd(kap,rho8, b6);
                    __m256d res2 = _mm256_mul_pd(temp2, denomp);
                    _mm256_store_pd(rho_mid + (i )*M_padded + j-3,res2);



                }





                __m256d rho1 = _mm256_load_pd(rho_ + m*M_padded + j-7);
                __m256d b1 = _mm256_load_pd(b + 8*nip -8 );
                __m256d res1 = _mm256_fmadd_pd(k2p,rho1,b1);
                _mm256_store_pd(b + 8*nip -8,res1);

                __m256d rho2 = _mm256_load_pd(rho_ + m*M_padded + j-3);
                __m256d b2 = _mm256_load_pd(b + 8*nip -4 );
                __m256d res2 = _mm256_fmadd_pd(k2p,rho2,b2);
                _mm256_store_pd(b + 8*nip -4,res2);




                __m256d rho3 = _mm256_load_pd(rho_mid + n*M_padded + (j-7));
                __m256d denomp = _mm256_broadcast_sd(denom+n);


                __m256d temp1 = _mm256_fnmadd_pd(kap,rho3, res1);
                __m256d res3 = _mm256_mul_pd(temp1, denomp);
                _mm256_store_pd(rho_mid + (n + 1)*M_padded + j-7,res3);


                __m256d rho4 = _mm256_load_pd(rho_mid + n*M_padded + (j-3));

                __m256d temp2 = _mm256_fnmadd_pd(kap,rho4, res2);
                __m256d res4 = _mm256_mul_pd(temp2, denomp);
                _mm256_store_pd(rho_mid + (n + 1)*M_padded + j-3,res4);






                for ( i = n; i-- > 0;) {
                    __m256d rho2 = _mm256_load_pd(rho_mid + (i+1)*M_padded + (j-7));
                    __m256d rho3 = _mm256_loadu_pd(rho_mid + (i+2)*M_padded + (j-7));
                    __m256d ccp = _mm256_broadcast_sd(cc + i);

                    __m256d res = _mm256_fnmadd_pd(ccp,rho3, rho2);
                    _mm256_store_pd(rho_mid + (i + 1)*M_padded + j-7,res);


                    rho2 = _mm256_load_pd(rho_mid + (i+1)*M_padded + (j-3));
                    rho3 = _mm256_loadu_pd(rho_mid + (i+2)*M_padded + (j-3));
                    ccp = _mm256_broadcast_sd(cc + i);

                    res = _mm256_fnmadd_pd(ccp,rho3, rho2);
                    _mm256_store_pd(rho_mid + (i + 1)*M_padded + j-3,res);




                }
            }

            for(j=m-m%8-7; j < m; j++){
                for ( i = 1; i < m; i++) {
                    b[i - 1] = k_2 * rho_[i*M_padded+ (j - 1)] + (k_)*rho_[i*M_padded + j] + k_2 * rho_[i*M_padded + (j + 1)];
                    rho_mid[(i)*M_padded + j ] = (b[i-1] - ka*rho_mid[(i-1)*M_padded + j ])*denom[i-1];
                }
                b[nip-1]=b[nip-1] + k_2 * rho_[m*M_padded + j];

                rho_mid[(n + 1)*M_padded + j ] = (b[n] - ka*rho_mid[n*M_padded + j ]) *denom[n];
                for ( i = n; i-- > 0;) {
                    rho_mid[(i + 1)*M_padded + j ] = rho_mid[(i + 1)*M_padded + j] -cc[i] * rho_mid[(i + 2)*M_padded + j ];
                }
            }





            for (i = 8; i < m; i+=8) {

                for (size_type j = 1; j < m; ++j) {
                    c[8*j-8] = k_2*rho_mid[(i - 8)*M_padded + j] + (k_)*rho_mid[(i-7)*M_padded + j] + k_2 * rho_mid[(i - 6)*M_padded + j];
                    c[8*j-7] = k_2*rho_mid[(i - 7)*M_padded + j] + (k_)*rho_mid[(i-6)*M_padded + j] + k_2 * rho_mid[(i - 5)*M_padded + j];
                    c[8*j-6] = k_2 * rho_mid[(i - 6)*M_padded + j] + (k_)*rho_mid[(i-5)*M_padded + j] + k_2 * rho_mid[(i-4 )*M_padded + j];
                    c[8*j-5] = k_2 * rho_mid[(i - 5)*M_padded + j] + (k_)*rho_mid[(i-4)*M_padded + j] + k_2 * rho_mid[(i - 3)*M_padded + j];

                    c[8*j-4] = k_2*rho_mid[(i - 4)*M_padded + j] + (k_)*rho_mid[(i-3)*M_padded + j] + k_2 * rho_mid[(i - 2)*M_padded + j];
                    c[8*j-3] = k_2*rho_mid[(i - 3)*M_padded + j] + (k_)*rho_mid[(i-2)*M_padded + j] + k_2 * rho_mid[(i - 1)*M_padded + j];
                    c[8*j-2] = k_2 * rho_mid[(i - 2)*M_padded + j] + (k_)*rho_mid[(i-1)*M_padded + j] + k_2 * rho_mid[(i )*M_padded + j];
                    c[8*j-1] = k_2 * rho_mid[(i - 1)*M_padded + j] + (k_)*rho_mid[i*M_padded + j] + k_2 * rho_mid[(i + 1)*M_padded + j];

                    rho_[(i-7)*M_padded+j]= (c[8*j-8] - ka*rho_[(i-7)*M_padded+j-1]) *denom[j-1];
                    rho_[(i-6)*M_padded+j]= (c[8*j-7] - ka*rho_[(i-6)*M_padded+j-1]) *denom[j-1];
                    rho_[(i-5)*M_padded+j]= (c[8*j-6] - ka*rho_[(i-5)*M_padded+j-1]) *denom[j-1];
                    rho_[(i-4)*M_padded+j]= (c[8*j-5] - ka*rho_[(i-4)*M_padded+j-1]) *denom[j-1];
                    rho_[(i-3)*M_padded+j]= (c[8*j-4] - ka*rho_[(i-3)*M_padded+j-1]) *denom[j-1];
                    rho_[(i-2)*M_padded+j]= (c[8*j-3] - ka*rho_[(i-2)*M_padded+j-1]) *denom[j-1];
                    rho_[(i-1)*M_padded+j]= (c[8*j-2] - ka*rho_[(i-1)*M_padded+j-1]) *denom[j-1];
                    rho_[(i-0)*M_padded+j]= (c[8*j-1] - ka*rho_[(i-0)*M_padded+j-1]) *denom[j-1];
                }

                c[8*nip-8] += k_2 * rho_mid[(i-7)*M_padded + m];
                c[8*nip-7] += k_2 * rho_mid[(i-6)*M_padded + m];
                c[8*nip-6] += k_2 * rho_mid[(i-5)*M_padded + m];
                c[8*nip-5] += k_2 * rho_mid[(i-4)*M_padded + m];

                c[8*nip-4] += k_2 * rho_mid[(i-3)*M_padded + m];
                c[8*nip-3] += k_2 * rho_mid[(i-2)*M_padded + m];
                c[8*nip-2] += k_2 * rho_mid[(i-1)*M_padded + m];
                c[8*nip-1] += k_2 * rho_mid[i*M_padded + m];


                rho_[(i-7)*M_padded+nip] = (c[8*nip-8] - ka*rho_[(i-7)*M_padded+nip - 1]) *denom[n];
                rho_[(i-6)*M_padded+nip] = (c[8*nip-7] - ka*rho_[(i-6)*M_padded+ nip- 1]) *denom[n];
                rho_[(i-5)*M_padded+nip] = (c[8*nip-6] - ka*rho_[(i-5)*M_padded+nip-1]) *denom[n];
                rho_[(i-4)*M_padded+nip] = (c[8*nip-5] - ka*rho_[(i-4)*M_padded+nip - 1]) *denom[n];

                rho_[(i-3)*M_padded+nip] = (c[8*nip-4] - ka*rho_[(i-3)*M_padded+nip - 1]) *denom[n];
                rho_[(i-2)*M_padded+nip] = (c[8*nip-3] - ka*rho_[(i-2)*M_padded+ nip- 1]) *denom[n];
                rho_[(i-1)*M_padded+nip] = (c[8*nip-2] - ka*rho_[(i-1)*M_padded+nip-1]) *denom[n];
                rho_[(i)*M_padded+nip] = (c[8*nip-1] - ka*rho_[i*M_padded+nip - 1]) *denom[n];


                for (int ii = n; ii-- > 0;) {
                    rho_[(i-7)*M_padded+ii+1] = rho_[(i-7)*M_padded+ii+1]- cc[ii] * rho_[(i-7)*M_padded +ii+ 2];
                    rho_[(i-6)*M_padded+ii+1] = rho_[(i-6)*M_padded+ii+1] - cc[ii] * rho_[(i-6)*M_padded+ii + 2];
                    rho_[(i-5)*M_padded+ii+1] = rho_[(i-5)*M_padded+ ii+1]- cc[ii] * rho_[(i-5)*M_padded+ii + 2];
                    rho_[(i-4)*M_padded+ii+1] = rho_[(i-4)*M_padded+ ii+1] - cc[ii] * rho_[(i-4)*M_padded+ii + 2];

                    rho_[(i-3)*M_padded+ii+1] = rho_[(i-3)*M_padded+ii+1]- cc[ii] * rho_[(i-3)*M_padded +ii+ 2];
                    rho_[(i-2)*M_padded+ii+1] = rho_[(i-2)*M_padded+ii+1] - cc[ii] * rho_[(i-2)*M_padded+ii + 2];
                    rho_[(i-1)*M_padded+ii+1] = rho_[(i-1)*M_padded+ ii+1]- cc[ii] * rho_[(i-1)*M_padded+ii + 2];
                    rho_[i*M_padded+ii+1] = rho_[i*M_padded+ ii+1] - cc[ii] * rho_[i*M_padded+ii + 2];
                }

            }
            i = i - 7;
            for (; i < m; ++i) {

                for (j = 1; j < m; j ++) {
                    c[j - 1] = k_2 * rho_mid[(i - 1)*M_padded + j] + (k_)*rho_mid[i*M_padded + j] + k_2 * rho_mid[(i + 1)*M_padded + j];
                    rho_[i*M_padded + (j)] = (c[j-1] - ka*rho_[i*M_padded + j-1]) *denom[j-1];
                }

                c[m - 2] = c[m - 2] + k_2 * rho_mid[i*M_padded + m];




                rho_[i*M_padded + n + 1] = (c[n] - ka*rho_[i*M_padded + n]) *denom[n];

                for (ii = n; ii-- > 0; ) {
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
    system.thomas_LU();
    t.start();
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
