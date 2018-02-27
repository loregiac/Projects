#include <iostream>
#include <algorithm>
#include <string>
#include <fstream>
#include <cassert>
#include <vector>
#include <cmath>
#include <malloc.h>
#include <immintrin.h>
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

            b = (double*) _mm_malloc(8*nip*sizeof(double), 32);
            c = new double[nip];

            rho_ = (double*) _mm_malloc((M_*M_padded+3)* sizeof(double), 32);
            rho_ = rho_ + 3;
            rho_mid = (double*) _mm_malloc((M_*M_padded+3) * sizeof(double), 32);
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




            for (j=8; j < m; j=j+8) {

                __m256d k2p = _mm256_broadcast_sd(&k_2);
                __m256d kp = _mm256_broadcast_sd(&k_);
                for ( i = 1; i < m; i++) {
                    __m256d rho1 = _mm256_loadu_pd(rho_ + i*M_padded + (j-8));
                    __m256d rho2 = _mm256_load_pd(rho_ + i*M_padded + (j-7));
                    __m256d rho3 = _mm256_loadu_pd(rho_ + i*M_padded + (j-6));
                    __m256d b1 = _mm256_mul_pd(k2p,rho1);
                    __m256d b2 = _mm256_fmadd_pd(kp,rho2,b1);
                    __m256d b3 = _mm256_fmadd_pd(k2p,rho3,b2);
                    _mm256_store_pd(b + 8*i-8,b3);

                    rho1 = _mm256_loadu_pd(rho_ + i*M_padded + (j-4));
                    rho2 = _mm256_load_pd(rho_ + i*M_padded + (j-3));
                    rho3 = _mm256_loadu_pd(rho_ + i*M_padded + (j-2));
                    b1 = _mm256_mul_pd(k2p,rho1);
                    b2 = _mm256_fmadd_pd(kp,rho2,b1);
                    b3 = _mm256_fmadd_pd(k2p,rho3,b2);
                    _mm256_store_pd(b + 8*i-4,b3);





                }

                __m256d rho1 = _mm256_load_pd(rho_ + (j-7));
                __m256d b1 = _mm256_load_pd(b);
                __m256d res1 = _mm256_fmadd_pd(k2p,rho1,b1);
                _mm256_store_pd(b ,res1);

                rho1 = _mm256_load_pd(rho_ + (j-3));
                b1 = _mm256_load_pd(b);
                res1 = _mm256_fmadd_pd(k2p,rho1,b1);
                _mm256_store_pd(b ,res1);






                __m256d rho2 = _mm256_load_pd(rho_ + m*M_padded + j-7);
                __m256d b2 = _mm256_load_pd(b + 8*nip -8 );
                __m256d res2 = _mm256_fmadd_pd(k2p,rho2,b2);
                _mm256_store_pd(b + 8*nip -8,res2);

                rho2 = _mm256_load_pd(rho_ + m*M_padded + j-3);
                b2 = _mm256_load_pd(b + 8*nip -4 );
                res2 = _mm256_fmadd_pd(k2p,rho2,b2);
                _mm256_store_pd(b + 8*nip -4,res2);





                rho_mid[M_padded + j-7] = b[0] * kbd;
                rho_mid[M_padded + j-6] = b[1] * kbd;
                rho_mid[M_padded + j-5] = b[2] * kbd;
                rho_mid[M_padded + j-4] = b[3] * kbd;

                rho_mid[M_padded + j-3] = b[4] * kbd;
                rho_mid[M_padded + j-2] = b[5] * kbd;
                rho_mid[M_padded + j-1] = b[6] * kbd;
                rho_mid[M_padded + j ] = b[7] * kbd;
                __m256d kap = _mm256_broadcast_sd(&ka);
                for ( i = 1; i < n; i++) {
                    __m256d rho2 = _mm256_load_pd(rho_mid + i*M_padded + (j-7));
                    __m256d denomp = _mm256_broadcast_sd(denom+i);
                    __m256d b1 = _mm256_load_pd(b+8*i);

                    __m256d temp = _mm256_fnmadd_pd(kap,rho2, b1);
                    __m256d res = _mm256_mul_pd(temp, denomp);
                    _mm256_store_pd(rho_mid + (i + 1)*M_padded + j-7,res);

                    rho2 = _mm256_load_pd(rho_mid + i*M_padded + (j-3));
                    denomp = _mm256_broadcast_sd(denom+i);
                    b1 = _mm256_load_pd(b+8*i+4);

                    temp = _mm256_fnmadd_pd(kap,rho2, b1);
                    res = _mm256_mul_pd(temp, denomp);
                    _mm256_store_pd(rho_mid + (i + 1)*M_padded + j-3,res);





                }
                rho2 = _mm256_load_pd(rho_mid + n*M_padded + (j-7));
                __m256d denomp = _mm256_broadcast_sd(denom+n);
                b1 = _mm256_load_pd(b+8*nip-8);

                __m256d temp = _mm256_fnmadd_pd(kap,rho2, b1);
                __m256d res = _mm256_mul_pd(temp, denomp);
                _mm256_store_pd(rho_mid + (n + 1)*M_padded + j-7,res);


                rho2 = _mm256_load_pd(rho_mid + n*M_padded + (j-3));
                denomp = _mm256_broadcast_sd(denom+n);
                b1 = _mm256_load_pd(b+8*nip-4);

                temp = _mm256_fnmadd_pd(kap,rho2, b1);
                res = _mm256_mul_pd(temp, denomp);
                _mm256_store_pd(rho_mid + (n + 1)*M_padded + j-3,res);






                cc[n] = kc;
                for ( i = n; i-- > 0;) {
                    rho2 = _mm256_load_pd(rho_mid + (i+1)*M_padded + (j-7));
                    __m256d rho3 = _mm256_loadu_pd(rho_mid + (i+2)*M_padded + (j-7));
                    __m256d ccp = _mm256_broadcast_sd(cc + i);

                    res = _mm256_fnmadd_pd(ccp,rho3, rho2);
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
                }
                b[0] = b[0] + k_2 * rho_[0 * M_padded + j ];
                b[nip-1]=b[nip-1] + k_2 * rho_[0*M_padded + j];

                rho_mid[M_padded + j ] = b[0] * kbd;
                for ( i = 1; i < n; i++) {
                    rho_mid[(i + 1)*M_padded + j ] = (b[i] - ka*rho_mid[(i)*M_padded + j ])*denom[i];
                }
                rho_mid[(n + 1)*M_padded + j ] = (b[n] - ka*rho_mid[n*M_padded + j ]) *denom[n];
                cc[n] = kc;
                for ( i = n; i-- > 0;) {
                    rho_mid[(i + 1)*M_padded + j ] = rho_mid[(i + 1)*M_padded + j] -cc[i] * rho_mid[(i + 2)*M_padded + j ];
                }
            }





            for (i = 1; i < m; ++i) {

                for (j = 1; j < m; j ++) {
                    c[j - 1] = k_2 * rho_mid[(i - 1)*M_padded + j] + (k_)*rho_mid[i*M_padded + j] + k_2 * rho_mid[(i + 1)*M_padded + j];
                }
                c[0] = c[0] + k_2 * rho_mid[i*M_padded + 0];
                c[m - 2] = c[m - 2] + k_2 * rho_mid[i*M_padded + m];


                size_type n = nip - 1;
                rho_[i*M_padded + 1] = c[0] * kbd;

                for (ii = 1; ii < n; ii++) {
                    rho_[i*M_padded + (ii + 1)] = (c[ii] - ka*rho_[i*M_padded + ii]) *denom[ii];
                }
                rho_[i*M_padded + n + 1] = (c[n] - ka*rho_[i*M_padded + n]) *denom[n];
                cc[n] = kc;

                for (ii = n; ii >= 0; ii = ii-1) {
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

            for (int i = 1; i < nip; i++) {
                cc[i] = kc / (kb - ka*cc[i - 1]);
                denom[i] = 1 / (kb - ka*cc[i - 1]);
            }
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
    double t1, t2;

    system.thomas_LU();
    t1=omp_get_wtime();
    while (time < tmax) {
        system.advance();
        time += dt;
    }


    t2 = omp_get_wtime();

    std::cout << "Timing : " << t2-t1 << std::endl;
    curr_err = system.linf_error(time);
    max_err = std::max(max_err, curr_err);

    std::cout << "err: " << max_err << std::endl;


    return 0;
}
