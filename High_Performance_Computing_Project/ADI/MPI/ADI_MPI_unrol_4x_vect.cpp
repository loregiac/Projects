# 3 "ADI_MPI_unrol_4x_vect.cpp"
#include <iostream>
#include <algorithm>
#include <string>
#include <fstream>
#include <cassert>
#include <vector>
#include <cmath>
#include <malloc.h>
#include <immintrin.h>
#include <mpi.h>

#define PI 3.141592653589793
#define TAG 0
typedef int size_type;

struct world_info
{

    int size;
    int dims_x;
    int dims_y;


    int left_proc;
    int right_proc;
    int top_proc;
    int bottom_proc;

    int rank;
    int cart_rank;
    int coord_x;
    int coord_y;


} world;

struct position
{
    double x,y;
};

class Diffusion2D {
    public:
        Diffusion2D(const double D,
                const double L,
                const size_type M,
                const double dt )
            : D_(D), L_(L), M_(M), Mtot(M_*M_), dt_(dt), m(M_ - 1)
        {


            M_glob = M_;


            Mx_loc = M_glob / world.dims_x;
            My_loc = M_glob / world.dims_y;

            M_block = (M_glob-2) / world.size ;


            dh_ = L_ / m;
            M_padded = Mx_loc;
            int periods[2] = {false, false};
            int dims[2] = {world.dims_x, world.dims_y};
            MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, true, &cart_comm);

            MPI_Comm_rank(cart_comm,&world.cart_rank);

            MPI_Cart_shift(cart_comm, 1, 1, &world.bottom_proc, &world.top_proc);


            int coords[2];
            MPI_Cart_coords(cart_comm, world.cart_rank, 2, coords);

            world.coord_x = coords[0];
            world.coord_y = coords[1];


            MPI_Type_contiguous(M_block, MPI_DOUBLE, &block_row);
            MPI_Type_commit(&block_row);


            MPI_Type_hvector(M_block, 1, (world.size*(M_block)+2)*sizeof(double), block_row, &block_send);
            MPI_Type_commit(&block_send);
            MPI_Type_create_resized(block_send, 0, (M_block)*sizeof(double),
                    &block_resized_send);
            MPI_Type_commit(&block_resized_send);


            MPI_Type_create_resized(block_row, 0, (M_block+2)*sizeof(double), &block_recv);
            MPI_Type_commit(&block_recv);
            MPI_Type_contiguous(M_block,block_recv, &block_resized_recv);

            MPI_Type_commit(&block_resized_recv);





            MPI_Type_contiguous(world.size*M_block, MPI_DOUBLE, &bottom_boundary);
            MPI_Type_commit(&bottom_boundary);

            MPI_Type_contiguous(world.size*M_block, MPI_DOUBLE, &top_boundary);
            MPI_Type_commit(&top_boundary);

            MPI_Type_vector(world.size*M_block, 1, M_block + 2, MPI_DOUBLE, &left_boundary);
            MPI_Type_commit(&left_boundary);

            MPI_Type_vector(world.size*M_block, 1, M_block + 2, MPI_DOUBLE, &right_boundary);
            MPI_Type_commit(&right_boundary);

            MPI_Cart_shift(cart_comm, 1, 1, &world.top_proc, &world.bottom_proc);



            xmin_loc = world.coord_x * Mx_loc * dh_;
            xmax_loc = xmin_loc + (Mx_loc - 1) * dh_;
            ymin_loc = world.coord_y * (M_block) * dh_;
            if( world.coord_y == 0) ymin_loc = 0;
            ymax_loc = ymin_loc + (My_loc-1) * dh_;
            M_padded = Mx_loc;

            k_ = dt_ * D_ / (dh_ * dh_);


            k_2 = k_/2;

            kb = 1 + k_;
            kbd = 1 / kb;

            ka = -k_ / 2;
            kc = -k_ / 2;
            nip = M_ - 2;


            k_ = 1 - k_;

            b = new double[4*nip+1];
            c = new double[8*nip+1];


            rho_ = (double*) _mm_malloc((Mx_loc*(M_block+2)+100)*sizeof(double), 32);
            rho_mid = (double*) _mm_malloc((Mx_loc*(My_loc+2)+100)*sizeof(double), 32);
            rho_recv = (double*) _mm_malloc((Mx_loc*(My_loc+2)+100)*sizeof(double), 32);
        }
        ~Diffusion2D()
        {
            delete(b);
            delete(c);
            _mm_free(rho_);
            _mm_free(rho_mid);
            _mm_free(rho_recv);
            MPI_Type_free(&left_boundary);
            MPI_Type_free(&right_boundary);
            MPI_Type_free(&bottom_boundary);
            MPI_Type_free(&top_boundary);





            MPI_Type_free(&block_send);
            MPI_Type_free(&block_resized_recv);
            MPI_Type_free(&block_resized_send);
            MPI_Type_free(&block_recv);

            MPI_Type_free(&block_row);

            MPI_Comm_free(&cart_comm);
            MPI_Finalize();
        }



        void advance()
        {
            size_type n = nip - 1;




            MPI_Request request[4];
            MPI_Status status[4];
# 198 "ADI_MPI_unrol_4x_vect.cpp"
            if (world.coord_y % 2 == 0) {

                MPI_Isend(&rho_[(Mx_loc)*(M_block)+1],1, bottom_boundary, world.bottom_proc, TAG, cart_comm, &request[0]);
                MPI_Irecv(&rho_[(M_block+1)*(Mx_loc)+1], 1, bottom_boundary, world.bottom_proc, TAG, cart_comm, &request[1]);
                MPI_Isend(&rho_[(Mx_loc)+1], 1, top_boundary, world.top_proc, TAG, cart_comm, &request[2]);
                MPI_Irecv(&rho_[1], 1, top_boundary, world.top_proc, TAG, cart_comm, &request[3]);


            }
            else {

                MPI_Irecv(&rho_[1], 1, top_boundary, world.top_proc, TAG, cart_comm, &request[0]);
                MPI_Isend(&rho_[(Mx_loc)+1], 1, top_boundary, world.top_proc, TAG, cart_comm, &request[1]);
                MPI_Irecv(&rho_[(Mx_loc)*(M_block+1)+1], 1, bottom_boundary, world.bottom_proc, TAG, cart_comm, &request[2]);
                MPI_Isend(&rho_[(Mx_loc)*(M_block)+1], 1, bottom_boundary, world.bottom_proc, TAG, cart_comm, &request[3]);

            }
# 243 "ADI_MPI_unrol_4x_vect.cpp"
            for (size_t i = 9; i < M_block; i+=8) {

                for (size_type j = 1; j < m; ++j) {
                    c[j-1] = k_2*rho_[(i - 8)*M_padded + j] + (k_)*rho_[(i-7)*M_padded + j] + k_2 * rho_[(i - 6)*M_padded + j];
                    c[nip+j-1] = k_2*rho_[(i - 7)*M_padded + j] + (k_)*rho_[(i-6)*M_padded + j] + k_2 * rho_[(i - 5)*M_padded + j];
                    c[2* nip+ j-1] = k_2 * rho_[(i - 6)*M_padded + j] + (k_)*rho_[(i-5)*M_padded + j] + k_2 * rho_[(i-4 )*M_padded + j];
                    c[3* nip+j-1] = k_2 * rho_[(i - 5)*M_padded + j] + (k_)*rho_[(i-4)*M_padded + j] + k_2 * rho_[(i - 3)*M_padded + j];

                    c[4*nip+j-1] = k_2*rho_[(i - 4)*M_padded + j] + (k_)*rho_[(i-3)*M_padded + j] + k_2 * rho_[(i - 2)*M_padded + j];
                    c[5*nip+j-1] = k_2*rho_[(i - 3)*M_padded + j] + (k_)*rho_[(i-2)*M_padded + j] + k_2 * rho_[(i - 1)*M_padded + j];
                    c[6*nip+j-1] = k_2 * rho_[(i - 2)*M_padded + j] + (k_)*rho_[(i-1)*M_padded + j] + k_2 * rho_[(i )*M_padded + j];
                    c[7*nip+j-1] = k_2 * rho_[(i - 1)*M_padded + j] + (k_)*rho_[i*M_padded + j] + k_2 * rho_[(i + 1)*M_padded + j];
                }


                size_type n = nip - 1;
                rho_mid[(i-7)*M_padded+1] = c[0]*kbd;
                rho_mid[(i-6)*M_padded+1] = c[nip] * kbd;
                rho_mid[(i-5)*M_padded+1] = c[2*nip] * kbd;
                rho_mid[(i-4)*M_padded+1] = c[3*nip] * kbd;
                rho_mid[(i-3)*M_padded+1] = c[4*nip]*kbd;
                rho_mid[(i-2)*M_padded+1] = c[5*nip] * kbd;
                rho_mid[(i-1)*M_padded+1] = c[6*nip] * kbd;
                rho_mid[(i-0)*M_padded+1] = c[7*nip] * kbd;

                for (size_t ii = 1; ii < n; ii++) {
                    rho_mid[(i-7)*M_padded+ii+1]= (c[ii] - ka*rho_mid[(i-7)*M_padded+ii]) *denom[ii];
                    rho_mid[(i-6)*M_padded+ii+1]= (c[nip+ii] - ka*rho_mid[(i-6)*M_padded+ii]) *denom[ii];
                    rho_mid[(i-5)*M_padded+ii+1]= (c[2*nip+ii] - ka*rho_mid[(i-5)*M_padded+ii]) *denom[ii];
                    rho_mid[(i-4)*M_padded+ii+1]= (c[3*nip+ii] - ka*rho_mid[(i-4)*M_padded+ii]) *denom[ii];
                    rho_mid[(i-3)*M_padded+ii+1]= (c[4*nip+ii] - ka*rho_mid[(i-3)*M_padded+ii]) *denom[ii];
                    rho_mid[(i-2)*M_padded+ii+1]= (c[5*nip+ii] - ka*rho_mid[(i-2)*M_padded+ii]) *denom[ii];
                    rho_mid[(i-1)*M_padded+ii+1]= (c[6*nip+ii] - ka*rho_mid[(i-1)*M_padded+ii]) *denom[ii];
                    rho_mid[(i-0)*M_padded+ii+1]= (c[7*nip+ii] - ka*rho_mid[(i-0)*M_padded+ii]) *denom[ii];
                }
                rho_mid[(i-7)*M_padded+nip] = (c[n] - ka*rho_mid[(i-7)*M_padded+nip - 1]) *denom[n];
                rho_mid[(i-6)*M_padded+nip] = (c[nip+n] - ka*rho_mid[(i-6)*M_padded+ nip- 1]) *denom[n];
                rho_mid[(i-5)*M_padded+nip] = (c[2*nip+n] - ka*rho_mid[(i-5)*M_padded+nip-1]) *denom[n];
                rho_mid[(i-4)*M_padded+nip] = (c[3*nip+n] - ka*rho_mid[(i-4)*M_padded+nip - 1]) *denom[n];

                rho_mid[(i-3)*M_padded+nip] = (c[4*nip+n] - ka*rho_mid[(i-3)*M_padded+nip - 1]) *denom[n];
                rho_mid[(i-2)*M_padded+nip] = (c[5*nip+n] - ka*rho_mid[(i-2)*M_padded+ nip- 1]) *denom[n];
                rho_mid[(i-1)*M_padded+nip] = (c[6*nip+n] - ka*rho_mid[(i-1)*M_padded+nip-1]) *denom[n];
                rho_mid[(i)*M_padded+nip] = (c[7*nip+n] - ka*rho_mid[i*M_padded+nip - 1]) *denom[n];


                for (size_t ii = n; ii-- > 0;) {
                    rho_mid[(i-7)*M_padded+ii+1] = rho_mid[(i-7)*M_padded+ii+1]- cc[ii] * rho_mid[(i-7)*M_padded +ii+ 2];
                    rho_mid[(i-6)*M_padded+ii+1] = rho_mid[(i-6)*M_padded+ii+1] - cc[ii] * rho_mid[(i-6)*M_padded+ii + 2];
                    rho_mid[(i-5)*M_padded+ii+1] = rho_mid[(i-5)*M_padded+ ii+1]- cc[ii] * rho_mid[(i-5)*M_padded+ii + 2];
                    rho_mid[(i-4)*M_padded+ii+1] = rho_mid[(i-4)*M_padded+ ii+1] - cc[ii] * rho_mid[(i-4)*M_padded+ii + 2];

                    rho_mid[(i-3)*M_padded+ii+1] = rho_mid[(i-3)*M_padded+ii+1]- cc[ii] * rho_mid[(i-3)*M_padded +ii+ 2];
                    rho_mid[(i-2)*M_padded+ii+1] = rho_mid[(i-2)*M_padded+ii+1] - cc[ii] * rho_mid[(i-2)*M_padded+ii + 2];
                    rho_mid[(i-1)*M_padded+ii+1] = rho_mid[(i-1)*M_padded+ ii+1]- cc[ii] * rho_mid[(i-1)*M_padded+ii + 2];
                    rho_mid[i*M_padded+ii+1] = rho_mid[i*M_padded+ ii+1] - cc[ii] * rho_mid[i*M_padded+ii + 2];
                }




            }

            for (size_t i = M_block-(M_block -9)%8; i < M_block; ++i) {

                for (size_t j = 1; j < m; j ++) {
                    c[j - 1] = k_2 * rho_[(i - 1)*M_padded + j] + (k_)*rho_[i*M_padded + j] + k_2 * rho_[(i + 1)*M_padded + j];
                }


                size_type n = nip - 1;
                rho_mid[i*M_padded + 1] = c[0] * kbd;

                for (size_t ii = 1; ii < n; ii++) {
                    rho_mid[i*M_padded + (ii + 1)] = (c[ii] - ka*rho_mid[i*M_padded + ii]) *denom[ii];
                }
                rho_mid[i*M_padded + n + 1] = (c[n] - ka*rho_mid[i*M_padded + n]) *denom[n];
                cc[n] = kc;

                for (size_t ii = n; ii-- > 0; ) {
                    rho_mid[i*M_padded + (ii + 1)] -= cc[ii] * rho_mid[i*M_padded + ii + 2];
                }
            }


            MPI_Waitall(4,request,status);

            {size_t i = 1;
                for (size_t j = 1; j < Mx_loc-1; j ++) {
                    c[j - 1] = k_2 * rho_[(i - 1)*M_padded + j] + (k_)*rho_[i*M_padded + j] + k_2 * rho_[(i + 1)*M_padded + j];
                }



                size_type n = nip - 1;
                rho_mid[i*M_padded + 1] = c[0] * kbd;

                for (size_t j = 1; j < n; j++) {
                    rho_mid[i*M_padded + (j + 1)] = (c[j] - ka*rho_mid[i*M_padded + j]) *denom[j];
                }
                rho_mid[i*M_padded + n + 1] = (c[n] - ka*rho_mid[i*M_padded + n]) *denom[n];
                cc[n] = kc;

                for (size_t j = n; j-->0;) {
                    rho_mid[i*M_padded + (j + 1)] -= cc[j] * rho_mid[i*M_padded + j + 2];
                }
            }

            {size_t i = M_block;
                for (size_t j = 1; j < Mx_loc-1; j ++) {
                    c[j - 1] = k_2 * rho_[(i - 1)*M_padded + j] + (k_)*rho_[i*M_padded + j] + k_2 * rho_[(i + 1)*M_padded + j];
                }



                size_type n = nip - 1;
                rho_mid[i*M_padded + 1] = c[0] * kbd;

                for (size_t j = 1; j < n; j++) {
                    rho_mid[i*M_padded + (j + 1)] = (c[j] - ka*rho_mid[i*M_padded + j]) *denom[j];
                }
                rho_mid[i*M_padded + n + 1] = (c[n] - ka*rho_mid[i*M_padded + n]) *denom[n];
                cc[n] = kc;

                for (size_t j = n; j-->0;) {
                    rho_mid[i*M_padded + (j + 1)] -= cc[j] * rho_mid[i*M_padded + j + 2];
                }
            }


            MPI_Alltoall(&rho_mid[(M_block)*world.size+3],1 , block_resized_send, &rho_recv[M_block+2+1],1, block_resized_recv, MPI_COMM_WORLD);



            MPI_Request requests[4];
            if (world.coord_y % 2 == 0) {

                MPI_Isend(&rho_recv[2*(M_block+2)-2],1, right_boundary, world.bottom_proc, TAG, cart_comm, &requests[0]);
                MPI_Irecv(&rho_recv[2*(M_block+2)-1], 1, right_boundary, world.bottom_proc, TAG, cart_comm, &requests[1]);
                MPI_Isend(&rho_recv[(M_block+2)+1], 1, left_boundary, world.top_proc, TAG, cart_comm, &requests[2]);
                MPI_Irecv(&rho_recv[(M_block+2)], 1, left_boundary, world.top_proc, TAG, cart_comm, &requests[3]);


            }
            else {

                MPI_Irecv(&rho_recv[(M_block+2)], 1, left_boundary, world.top_proc, TAG, cart_comm, &requests[0]);
                MPI_Isend(&rho_recv[(M_block+2)+1], 1, left_boundary, world.top_proc, TAG, cart_comm, &requests[1]);
                MPI_Irecv(&rho_recv[2*(M_block+2)-1], 1, right_boundary, world.bottom_proc, TAG, cart_comm, &requests[2]);
                MPI_Isend(&rho_recv[2*(M_block+2)-2], 1, right_boundary, world.bottom_proc, TAG, cart_comm, &requests[3]);

            }

            M_padded = M_block + 2;

            for (size_t j=5; j < M_block; j=j+4) {

                __m256d k2p = _mm256_broadcast_sd(&k_2);
                __m256d kp = _mm256_broadcast_sd(&k_);
                for (size_t i = 1; i < M_block *world.size +1; i++) {
                    __m256d rho1 = _mm256_loadu_pd(rho_recv + i*M_padded + (j-4));
                    __m256d rho2 = _mm256_loadu_pd(rho_recv + i*M_padded + (j-3));
                    __m256d rho3 = _mm256_loadu_pd(rho_recv + i*M_padded + (j-2));
                    __m256d b1 = _mm256_mul_pd(k2p,rho1);
                    __m256d b2 = _mm256_fmadd_pd(kp,rho2,b1);
                    __m256d b3 = _mm256_fmadd_pd(k2p,rho3,b2);
                    _mm256_storeu_pd(b + 4*i-4,b3);

                }


                rho_mid[M_padded + j-3] = b[0] * kbd;
                rho_mid[M_padded + j-2] = b[1] * kbd;
                rho_mid[M_padded + j-1] = b[2] * kbd;
                rho_mid[M_padded + j ] = b[3] * kbd;
                __m256d kap = _mm256_broadcast_sd(&ka);
                for (size_t i = 1; i < n; i++) {
                    __m256d rho2 = _mm256_loadu_pd(rho_mid + i*M_padded + (j-3));
                    __m256d denomp = _mm256_broadcast_sd(denom+i);
                    __m256d b1 = _mm256_loadu_pd(b+4*i);

                    __m256d temp = _mm256_fnmadd_pd(kap,rho2, b1);
                    __m256d res = _mm256_mul_pd(temp, denomp);
                    _mm256_storeu_pd(rho_mid + (i + 1)*M_padded + j-3,res);

                }
                __m256d rho2 = _mm256_loadu_pd(rho_mid + n*M_padded + (j-3));
                __m256d denomp = _mm256_broadcast_sd(denom+n);
                __m256d b1 = _mm256_loadu_pd(b+4*nip-4);

                __m256d temp = _mm256_fnmadd_pd(kap,rho2, b1);
                __m256d res = _mm256_mul_pd(temp, denomp);
                _mm256_storeu_pd(rho_mid + (n + 1)*M_padded + j-3,res);

                for (size_t i = n; i-- > 0;) {
                    rho2 = _mm256_loadu_pd(rho_mid + (i+1)*M_padded + (j-3));
                    __m256d rho3 = _mm256_loadu_pd(rho_mid + (i+2)*M_padded + (j-3));
                    __m256d ccp = _mm256_broadcast_sd(cc + i);

                    res = _mm256_fnmadd_pd(ccp,rho3, rho2);
                    _mm256_storeu_pd(rho_mid + (i + 1)*M_padded + j-3,res);

                }
            }





            for(size_t j = M_block-(M_block -5)%4; j < M_block; j++){
                for ( size_t i = 1; i < m; i++) {
                    b[i - 1] = k_2 * rho_recv[i*M_padded+ (j - 1)] + (k_)*rho_recv[i*M_padded + j] + k_2 * rho_recv[i*M_padded + (j + 1)];
                }


                rho_mid[M_padded + j ] = b[0] * kbd;
                for (size_t i = 1; i < n; i++) {
                    rho_mid[(i + 1)*M_padded + j ] = (b[i] - ka*rho_mid[(i)*M_padded + j ])*denom[i];
                }
                rho_mid[(n + 1)*M_padded + j ] = (b[n] - ka*rho_mid[n*M_padded + j ]) *denom[n];
                cc[n] = kc;
                for (size_t i = n; i-- > 0;) {
                    rho_mid[(i + 1)*M_padded + j ] = rho_mid[(i + 1)*M_padded + j] -cc[i] * rho_mid[(i + 2)*M_padded + j ];
                }
            }


            MPI_Waitall(4,requests,status);

            {size_t j=1;

                for (size_t i = 1; i < M_block *world.size +1; i++) {
                    b[i - 1] = k_2 * rho_recv[i*(M_block+2) + (j - 1)] + (k_)*rho_recv[i*(M_block+2) + j] + k_2 * rho_recv[i*(M_block+2) + (j + 1)];
                }

                size_type n = nip - 1;
                rho_mid[(M_block+2) + j ] = b[0] * kbd;

                for (size_t i = 1; i < n; i++) {
                    rho_mid[(i + 1)*(M_block+2) + j ] = (b[i] - ka*rho_mid[(i)*(M_block+2) + j ])*denom[i];
                }

                rho_mid[(n + 1)*(M_block+2) + j ] = (b[n] - ka*rho_mid[n*(M_block+2) + j ]) *denom[n];
                cc[n] = kc;
                for (size_t i = n; i-- > 0;) {
                    rho_mid[(i + 1)*(M_block+2) + j ] = rho_mid[(i + 1)*(M_block+2) + j] -cc[i] * rho_mid[(i + 2)*(M_block+2) + j ];

                }


            }

            {size_t j=M_block;

                for (size_t i = 1; i < M_block *world.size +1; i++) {
                    b[i - 1] = k_2 * rho_recv[i*(M_block+2) + (j - 1)] + (k_)*rho_recv[i*(M_block+2) + j] + k_2 * rho_recv[i*(M_block+2) + (j + 1)];
                }

                size_type n = nip - 1;
                rho_mid[(M_block+2) + j ] = b[0] * kbd;

                for (size_t i = 1; i < n; i++) {
                    rho_mid[(i + 1)*(M_block+2) + j ] = (b[i] - ka*rho_mid[(i)*(M_block+2) + j ])*denom[i];
                }

                rho_mid[(n + 1)*(M_block+2) + j ] = (b[n] - ka*rho_mid[n*(M_block+2) + j ]) *denom[n];
                cc[n] = kc;
                for (size_t i = n; i-- > 0;) {
                    rho_mid[(i + 1)*(M_block+2) + j ] = rho_mid[(i + 1)*(M_block+2) + j] -cc[i] * rho_mid[(i + 2)*(M_block+2) + j ];

                }


            }



            MPI_Alltoall(&rho_mid[(M_block+2)+1],1 ,block_resized_recv , &rho_[(M_block)*world.size+3],1,block_resized_send , MPI_COMM_WORLD);

            M_padded = Mx_loc;
            write_density();
        }





        void write_density() const
        {
            position current;
            std::string filename = "density" + std::to_string(world.rank) + ".dat";
            std::ofstream out_file(filename, std::ios::out);

            for (size_type i = 1; i < M_block+1; ++i) {
                for (size_type j = 1; j < Mx_loc-1; ++j){
                    current = get_position((i),j);
                    out_file << (current.x) << '\t' << (current.y) << '\t' << rho_[i*M_padded + j] << "\n";
                }
                out_file << "\n";
            }
            out_file.close();
        }

        void write_density_int() const
        {
            position current;
            std::string filename = "test.dat";
            if (world.rank == 0){filename = "density_int_0.dat";}

            if (world.rank ==1){filename = "density_int_1.dat";}
            std::ofstream out_file(filename, std::ios::out);

            for (size_type i = 1; i < Mx_loc-1; ++i) {
                for (size_type j = 1; j < M_block+1; ++j){
                    current = get_position((i),j);
                    out_file << i*dh_ << '\t' << j*dh_ << '\t' << rho_[i*(M_block+2) + j] << "\n";
                }
                out_file << "\n";
            }
            out_file.close();
        }

        double linf_error(double current_time) {
            double max_error = 0;
            double global_error = 0;
            double current_error = 0;
            int i_max, j_max;
            position current;
            for (size_type i = 1; i < M_block+1; ++i) {
                for (size_type j = 1; j < Mx_loc-1; ++j) {
                    current = get_position((i),j);
                    double numerical = rho_[i*M_padded + j];
                    double analytic = sin(PI*current.x) * sin(PI*current.y) * std::exp(-2 * D_*PI*PI*current_time);
                    current_error = fabs(numerical - analytic);
                    if (max_error < current_error ) {
                        i_max = i;
                        j_max = j;
                    }
                    max_error = std::max(max_error, current_error);
                }
            }
            MPI_Allreduce(&max_error, &global_error, 1, MPI_DOUBLE, MPI_MAX,
                    MPI_COMM_WORLD);

            return global_error;
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
        void initialize_density()
        {

            position current;
            for (size_type i = 0; i < M_block+2; i++) {
                for (size_type j = 0; j < Mx_loc; j++) {
                    current = get_position((i),j);
                    rho_[i*Mx_loc + j] = sin(PI*current.x)*sin(PI*current.y);
                    rho_mid[i*Mx_loc + j] =0;
                    rho_recv[i*Mx_loc + j] = 0;
                }
            }
        }

        inline position get_position(size_type i, size_type j) const
        {
            position p;
            p.x = xmin_loc + j*dh_;
            p.y = ymin_loc + i*dh_;
            return p;
        }


    private:







        double D_, L_;
        size_type M_, Mtot, m, nip, Mx_loc,My_loc, M_glob, M_block;
        size_type M_padded;
        double dh_, dt_, k_, k_2;
        int rank, prcos;
        double ka, kb, kc, kbd;
        double xmin_loc,ymin_loc;
        double *cc, *denom;

        double *rho_, *rho_mid, *rho_recv;
        MPI_Datatype block_send, block_resized_send,block_resized_recv, block_recv, block_row;
        double xmax_loc, ymax_loc;

        double *b, *c;

        MPI_Datatype left_boundary, right_boundary, bottom_boundary, top_boundary;
        MPI_Comm cart_comm;

};


int main(int argc, char* argv[])
{
    const double D = 1;
    const double L = 1;
    const size_type m = std::stod(argv[1]);
    const size_type M = m + 1;
    const double dt = std::stod(argv[2]);

    MPI_Init(&argc,&argv);

    int rank, procs;
    MPI_Comm_rank(MPI_COMM_WORLD, &world.rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world.size);



    int dims[2] = {1,world.size};
    MPI_Dims_create(world.size, 2, dims);
    world.dims_x = dims[0];
    world.dims_y = dims[1];

    const size_type M_world = (M-2) % world.size == 0 ? M : M + (world.size - (M-2) % world.size);


    Diffusion2D system(D, L, M_world, dt);
    system.initialize_density();



    double time = 0;
    double tmax = 0.1;
    double max_err = 0;
    double curr_err;
    double t1, t2;
    system.thomas_LU();
    t1 = MPI_Wtime();
    while (time < tmax) {
        system.advance();
        time += dt;
    }


    t2 = MPI_Wtime();

    std::cout << "Timing : " << t2-t1 << std::endl;
    curr_err = system.linf_error(time);
    max_err = std::max(max_err, curr_err);

    if (world.rank == 0) {
        std::cout << "err: " << max_err << std::endl;
    }


}
