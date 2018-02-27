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
            if (M_padded % 4) M_padded = M_padded + 4 - (M_padded % 4);
            std::cout << M_padded << std::endl;

            k_ = dt_ * D_ / (dh_ * dh_);


            k_2 = k_ / 2;

            kb = 1 + k_;
            kbd = 1 / kb;

            ka = -k_ / 2;
            kc = -k_ / 2;
            nip = M_ - 2;


            k_ = 1 - k_;


            rho_ = (double*)_mm_malloc((M_*M_padded + 3) * sizeof(double), 32);
            rho_mid = (double*)_mm_malloc((M_*M_padded + 3) * sizeof(double), 32);
            rho_trans = (double*)_mm_malloc((M_*M_padded + 3) * sizeof(double), 32);
            rho_ = rho_ + 3;
            rho_mid = rho_mid + 3;
            rho_trans = rho_trans + 3;
            initialize_density();
        }
        ~Diffusion2D()
        {
            rho_ = rho_ - 3;
            rho_mid = rho_mid - 3;
            _mm_free(rho_);
            _mm_free(rho_mid);
        }



        double run(double dt, double tmax)
        {
            double time;
            double *bb = (double*)_mm_malloc((8 * nip) * sizeof(double), 32);
            double *b;
            thomas_LU();
            cc[nip - 1] = kc;
            size_type n = nip - 1;
            double t1 = omp_get_wtime();

            for (time = 0; time < tmax; time += dt) {




                b = bb ;

                for (int j = 8; j < m; j = j + 8) {

                    __m256d zero = _mm256_set1_pd(0);

                    _mm256_store_pd(rho_mid + (j-7),zero);
                    _mm256_store_pd(rho_mid + (j-3),zero);
                    __m256d k2p = _mm256_broadcast_sd(&k_2);
                    __m256d kp = _mm256_broadcast_sd(&k_);
                    __m256d kap = _mm256_broadcast_sd(&ka);
                    for (int i = 1; i < m; i++) {
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



                    __m256d rho3 = _mm256_load_pd(rho_mid + n*M_padded + (j-7));
                    __m256d denomp = _mm256_broadcast_sd(denom+n);
                    __m256d b3 = _mm256_load_pd(b+8*nip-8);

                    __m256d temp1 = _mm256_fnmadd_pd(kap,rho3, b3);
                    __m256d res3 = _mm256_mul_pd(temp1, denomp);



                    __m256d rho4 = _mm256_load_pd(rho_mid + n*M_padded + (j-3));
                    __m256d b4 = _mm256_load_pd(b+8*nip-4);

                    __m256d temp2 = _mm256_fnmadd_pd(kap,rho4, b4);
                    __m256d res4 = _mm256_mul_pd(temp2, denomp);




                    __m128d tmp0 = _mm256_extractf128_pd(res3, 0);
                    __m128d tmp1 = _mm256_extractf128_pd(res3, 1);

                    _mm_store_sd(rho_trans + (n + 1) + (j - 7)*M_padded, tmp0);
                    _mm_storeh_pd(rho_trans + (n + 1) + (j - 6)*M_padded, tmp0);
                    _mm_store_sd(rho_trans + (n + 1) + (j - 5)*M_padded, tmp1);
                    _mm_storeh_pd(rho_trans + (n + 1) + (j - 4)*M_padded, tmp1);


                    tmp0 = _mm256_extractf128_pd(res4, 0);
                    tmp1 = _mm256_extractf128_pd(res4, 1);

                    _mm_store_sd(rho_trans + (n + 1) + (j - 3)*M_padded, tmp0);
                    _mm_storeh_pd(rho_trans + (n + 1) + (j - 2)*M_padded, tmp0);
                    _mm_store_sd(rho_trans + (n + 1) + (j - 1)*M_padded, tmp1);
                    _mm_storeh_pd(rho_trans + (n + 1) + (j)*M_padded, tmp1);

                    for (int i = n; i-- > 0;) {
                        rho_trans[(i + 1) + (j - 7)*M_padded ] = rho_mid[(i + 1)*M_padded + j - 7] -cc[i] * rho_trans[(i + 2) + (j - 7)*M_padded ];
                        rho_trans[(i + 1) + (j - 6)*M_padded ] = rho_mid[(i + 1)*M_padded + j - 6] -cc[i] * rho_trans[(i + 2) + (j - 6)*M_padded ];
                        rho_trans[(i + 1) + (j - 5)*M_padded ] = rho_mid[(i + 1)*M_padded + j - 5] -cc[i] * rho_trans[(i + 2) + (j - 5)*M_padded ];
                        rho_trans[(i + 1) + (j - 4)*M_padded ] = rho_mid[(i + 1)*M_padded + j - 4] -cc[i] * rho_trans[(i + 2) + (j - 4)*M_padded ];

                        rho_trans[(i + 1) + (j - 3)*M_padded ] = rho_mid[(i + 1)*M_padded + j - 3] -cc[i] * rho_trans[(i + 2) + (j - 3)*M_padded ];
                        rho_trans[(i + 1) + (j - 2)*M_padded ] = rho_mid[(i + 1)*M_padded + j - 2] -cc[i] * rho_trans[(i + 2) + (j - 2)*M_padded ];
                        rho_trans[(i + 1) + (j - 1)*M_padded ] = rho_mid[(i + 1)*M_padded + j - 1] -cc[i] * rho_trans[(i + 2) + (j - 1)*M_padded ];
                        rho_trans[(i + 1) + j *M_padded ] = rho_mid[(i + 1)*M_padded + j ] -cc[i] * rho_trans[(i + 2) + (j - 0)*M_padded ];

                    }

                }


                for (int j = m - 7 + m % 8; j < m; j++) {
                    for ( size_type i = 1; i < m; i++) {
                        b[i - 1] = k_2 * rho_[i*M_padded+ (j - 1)] + (k_)*rho_[i*M_padded + j] + k_2 * rho_[i*M_padded + (j + 1)];
                        rho_mid[(i)*M_padded + j ] = (b[i-1] - ka*rho_mid[(i-1)*M_padded + j ])*denom[i-1];
                    }
                    b[nip-1]=b[nip-1] + k_2 * rho_[0*M_padded + j];

                    rho_trans[(n + 1)+ j*M_padded ] = (b[n] - ka*rho_mid[n*M_padded + j ]) *denom[n];
                    for (size_type i = n; i-- > 0;) {
                        rho_trans[(i + 1)+ j*M_padded] = rho_mid[(i + 1)*M_padded + j] -cc[i] * rho_trans[(i + 2) + j*M_padded ];
                    }


                }

                __m256d zero = _mm256_set1_pd(0);
                _mm256_store_pd(rho_mid + (1),zero);
                _mm256_store_pd(rho_mid + (5),zero);

                for (int j = 8; j < m; j = j + 8) {



                    __m256d k2p = _mm256_broadcast_sd(&k_2);
                    __m256d kp = _mm256_broadcast_sd(&k_);
                    __m256d kap = _mm256_broadcast_sd(&ka);
                    for (int i = 1; i < m; i++) {
                        __m256d rho1 = _mm256_loadu_pd(rho_trans+ i*M_padded + (j-8));
                        __m256d rho2 = _mm256_load_pd(rho_trans+ i*M_padded + (j-7));
                        __m256d rho3 = _mm256_loadu_pd(rho_trans+ i*M_padded + (j-6));
                        __m256d b1 = _mm256_mul_pd(k2p,rho1);
                        __m256d b2 = _mm256_fmadd_pd(kp,rho2,b1);
                        __m256d b3 = _mm256_fmadd_pd(k2p,rho3,b2);
                        _mm256_store_pd(b + 8*i-8,b3);

                        __m256d rho4 = _mm256_loadu_pd(rho_trans+ i*M_padded + (j-4));
                        __m256d rho5 = _mm256_load_pd(rho_trans+ i*M_padded + (j-3));
                        __m256d rho6 = _mm256_loadu_pd(rho_trans+ i*M_padded + (j-2));
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



                    __m256d rho3 = _mm256_load_pd(rho_mid + n*M_padded + (j-7));
                    __m256d denomp = _mm256_broadcast_sd(denom+n);
                    __m256d b3 = _mm256_load_pd(b+8*nip-8);

                    __m256d temp1 = _mm256_fnmadd_pd(kap,rho3, b3);
                    __m256d res3 = _mm256_mul_pd(temp1, denomp);



                    __m256d rho4 = _mm256_load_pd(rho_mid + n*M_padded + (j-3));
                    __m256d b4 = _mm256_load_pd(b+8*nip-4);

                    __m256d temp2 = _mm256_fnmadd_pd(kap,rho4, b4);
                    __m256d res4 = _mm256_mul_pd(temp2, denomp);




                    __m128d tmp0 = _mm256_extractf128_pd(res3, 0);
                    __m128d tmp1 = _mm256_extractf128_pd(res3, 1);

                    _mm_store_sd(rho_ + (n + 1) + (j - 7)*M_padded, tmp0);
                    _mm_storeh_pd(rho_ + (n + 1) + (j - 6)*M_padded, tmp0);
                    _mm_store_sd(rho_ + (n + 1) + (j - 5)*M_padded, tmp1);
                    _mm_storeh_pd(rho_ + (n + 1) + (j - 4)*M_padded, tmp1);


                    __m128d tmp2 = _mm256_extractf128_pd(res4, 0);
                    __m128d tmp3 = _mm256_extractf128_pd(res4, 1);

                    _mm_store_sd(rho_ + (n + 1) + (j - 3)*M_padded, tmp2);
                    _mm_storeh_pd(rho_ + (n + 1) + (j - 2)*M_padded, tmp2);
                    _mm_store_sd(rho_ + (n + 1) + (j - 1)*M_padded, tmp3);
                    _mm_storeh_pd(rho_ + (n + 1) + (j)*M_padded, tmp3);

                    for (int i = n; i-- > 0;) {
                        rho_[(i + 1) + (j - 7)*M_padded ] = rho_mid[(i + 1)*M_padded + j - 7] -cc[i] * rho_[(i + 2) + (j - 7)*M_padded ];
                        rho_[(i + 1) + (j - 6)*M_padded ] = rho_mid[(i + 1)*M_padded + j - 6] -cc[i] * rho_[(i + 2) + (j - 6)*M_padded ];
                        rho_[(i + 1) + (j - 5)*M_padded ] = rho_mid[(i + 1)*M_padded + j - 5] -cc[i] * rho_[(i + 2) + (j - 5)*M_padded ];
                        rho_[(i + 1) + (j - 4)*M_padded ] = rho_mid[(i + 1)*M_padded + j - 4] -cc[i] * rho_[(i + 2) + (j - 4)*M_padded ];

                        rho_[(i + 1) + (j - 3)*M_padded ] = rho_mid[(i + 1)*M_padded + j - 3] -cc[i] * rho_[(i + 2) + (j - 3)*M_padded ];
                        rho_[(i + 1) + (j - 2)*M_padded ] = rho_mid[(i + 1)*M_padded + j - 2] -cc[i] * rho_[(i + 2) + (j - 2)*M_padded ];
                        rho_[(i + 1) + (j - 1)*M_padded ] = rho_mid[(i + 1)*M_padded + j - 1] -cc[i] * rho_[(i + 2) + (j - 1)*M_padded ];
                        rho_[(i + 1) + j *M_padded ] = rho_mid[(i + 1)*M_padded + j ] -cc[i] * rho_[(i + 2) + (j - 0)*M_padded ];

                    }

                }


                for (int j = m - 7 + m % 8; j < m; j++) {
                    rho_mid[(0)*M_padded + j] = 0;
                    for ( size_type i = 1; i < m; i++) {
                        b[i - 1] = k_2 * rho_trans[i*M_padded+ (j - 1)] + (k_)*rho_trans[i*M_padded + j] + k_2 * rho_trans[i*M_padded + (j + 1)];
                        rho_mid[(i)*M_padded + j ] = (b[i-1] - ka*rho_mid[(i-1)*M_padded + j ])*denom[i-1];
                    }
                    b[nip-1]=b[nip-1] + k_2 * rho_[0*M_padded + j];

                    rho_[(n + 1)+ j*M_padded ] = (b[n] - ka*rho_mid[n*M_padded + j ]) *denom[n];
                    for (size_type i = n; i-- > 0;) {
                        rho_[(i + 1)+ j*M_padded] = rho_mid[(i + 1)*M_padded + j] -cc[i] * rho_[(i + 2) + j*M_padded ];
                    }


                }
            }
            double t2 = omp_get_wtime();
            std::cout << "time without LU: " << t2-t1 << std::endl;
            return time;
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


        double *rho_, *rho_mid, *rho_trans, *cc, *denom;



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
    double t1 = omp_get_wtime();
    time = system.run(dt, tmax);
    double t2 = omp_get_wtime();

    std::cout << "Timing : " << t2 - t1 << std::endl;
    curr_err = system.linf_error(time);
    max_err = std::max(max_err, curr_err);

    std::cout << "err: " << max_err << std::endl;


    return 0;
}
