#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <nd_malloc.h>
#include <imageio.h>


#define BS      15
#define PSF_HS  11
#define SGMC    0.6
#define DELTA   0.001
#define NX      320
#define NY      240
#define NZ      3
#define HANTEI  0.001
#define CENTER_X 159
#define CENTER_Y 119
#define NDM     0


int main(void)
{
  FILE *fp;
  IMAGE in_img;
  char flnm[256];
  double **p, **dpdd, **dpdk, tmp1;
  double k, d, sgmc, RMSE;
  double **pd, **pk, **pdk;
  double **p0d, **p0k, **p1d, **p1k, **p2d, **p2k;
  double **kekkad, **kekkak;
  double **f, **fd, **fk;
  double **psfd, **psfk, **psff, **sabun;
  double t1 = 0.0, t2 = 0.0, t3 = 0.0, t4 = 0.0, t5 = 0.0, t6 = 0.0, t7 = 0.0, t8 = 0.0, t9 = 0.0, t10, t11, t12;
  double tf0p1d, tf1p0d, tf1p2d, tf2p1d, tf0p2d, tf2p0d;
  double tf0p1dk, tf1p0dk, tf1p2dk, tf2p1dk, tf0p2dk, tf2p0dk;
  double tf0p1k, tf1p0k, tf1p2k, tf2p1k, tf0p2k, tf2p0k;
  double te1f0p1d, te1f1p0d, te2f1p2d, te2f2p1d, te3f0p2d, te3f2p0d;
  double te1f0p1k, te1f1p0k, te2f1p2k, te2f2p1k, te3f0p2k, te3f2p0k;
  double a00, a01, a10, a11, b0, b1;
  double **f0p1, **f1p0, **f1p2, **f2p1, **f0p2, **f2p0;
  double **f0pd1, **f1pd0, **f1pd2, **f2pd1, **f0pd2, **f2pd0;
  double **f0pk1, **f1pk0, **f1pk2, **f2pk1, **f0pk2, **f2pk0;
  double ***fn;
  double **p0, **p1, **p2;
  double data[NX * NY], **f0p1delta;
  double **f0p1delta_d, **f1p0delta_d, **f1p2delta_d, **f2p1delta_d, **f0p2delta_d, **f2p0delta_d;
  double **f0p1delta_k, **f1p0delta_k, **f1p2delta_k, **f2p1delta_k, **f0p2delta_k, **f2p0delta_k;
  double **f01delta_d, **f10delta_d, **f12delta_d, **f21delta_d, **f02delta_d, **f20delta_d;
  double **f01delta_k, **f10delta_k, **f12delta_k, **f21delta_k, **f02delta_k, **f20delta_k;
  double ***p_delta;
  //double fp;
  double e[NZ], dedd[NZ], dedk[NZ], dekk[NZ], Ddd, Ddk, Dkk, Did, Dik;
  int z, nr, nw, zn;
  int i, j, n, m;
  double k0, d0, deltad, deltak, dprev, kprev, eprev;
  int count;
  void psf(double k, int z, double d, double sgmc, int psf_hs, double **p);
  void dpsfdd(double k, int z, double d, double sgmc, int psf_hs, double **p);
  void dpsfdk(double k, int z, double d, double sgmc, int psf_hs, double **p);
  void make_fnpn(int bs, int center_x, int center_y, double **fn, double **pn, double **fnpn);
  void make_fnpdn(int bs, int center_x, int center_y, double **fn, double **fnpdn, double k, int z, double d, double sgmc);
  void make_fnpkn(int bs, int center_x, int center_y, double **fn, double **fnpkn, double k, int z, double d, double sgmc);  
  void num_diff_method(int bs, double ** fnpndelta, double **fnpn, double **f);
  void num_diff_method_p(double ** fnpndelta, double **fnpn, double **f);
  
  k = 0.5;
  //z = 0;
  d = 1.0;
  sgmc = SGMC;
  
  //-----------------------------------------------------------------------------
  // malloc
  p    = malloc_double_2d(2 * PSF_HS + 1, PSF_HS,  2 * PSF_HS + 1, PSF_HS);
  p0   = malloc_double_2d(2 * PSF_HS + 1, PSF_HS,  2 * PSF_HS + 1, PSF_HS);
  p1   = malloc_double_2d(2 * PSF_HS + 1, PSF_HS,  2 * PSF_HS + 1, PSF_HS);
  p2   = malloc_double_2d(2 * PSF_HS + 1, PSF_HS,  2 * PSF_HS + 1, PSF_HS);  
  dpdd = malloc_double_2d(2 * PSF_HS + 1, PSF_HS,  2 * PSF_HS + 1, PSF_HS);
  dpdk = malloc_double_2d(2 * PSF_HS + 1, PSF_HS,  2 * PSF_HS + 1, PSF_HS);
  pd = malloc_double_2d(2 * PSF_HS + 1, PSF_HS,  2 * PSF_HS + 1, PSF_HS);
  pk = malloc_double_2d(2 * PSF_HS + 1, PSF_HS,  2 * PSF_HS + 1, PSF_HS);
  pdk = malloc_double_2d(2 * PSF_HS + 1, PSF_HS,  2 * PSF_HS + 1, PSF_HS);
  kekkad = malloc_double_2d(2 * PSF_HS + 1, PSF_HS,  2 * PSF_HS + 1, PSF_HS);
  kekkak = malloc_double_2d(2 * PSF_HS + 1, PSF_HS,  2 * PSF_HS + 1, PSF_HS);
  f = malloc_double_2d(2 * PSF_HS + BS, PSF_HS,  2 * PSF_HS + BS, PSF_HS);
  p_delta = malloc_double_3d(6, 0, 2 * PSF_HS + BS, PSF_HS,  2 * PSF_HS + BS, PSF_HS);
  
  fd = malloc_double_2d(BS, 0, BS, 0);
  fk = malloc_double_2d(BS, 0, BS, 0);
  psfd = malloc_double_2d(BS, 0, BS, 0);
  psfk = malloc_double_2d(BS, 0, BS, 0);
  psff = malloc_double_2d(BS, 0, BS, 0);
  sabun = malloc_double_2d(BS, 0, BS, 0);

  f0p1 = malloc_double_2d(BS, 0, BS, 0);
  f1p0 = malloc_double_2d(BS, 0, BS, 0);
  f1p2 = malloc_double_2d(BS, 0, BS, 0);
  f2p1 = malloc_double_2d(BS, 0, BS, 0);
  f0p2 = malloc_double_2d(BS, 0, BS, 0);
  f2p0 = malloc_double_2d(BS, 0, BS, 0);

  f0pd1 = malloc_double_2d(BS, 0, BS, 0);
  f1pd0 = malloc_double_2d(BS, 0, BS, 0);
  f1pd2 = malloc_double_2d(BS, 0, BS, 0);
  f2pd1 = malloc_double_2d(BS, 0, BS, 0);
  f0pd2 = malloc_double_2d(BS, 0, BS, 0);
  f2pd0 = malloc_double_2d(BS, 0, BS, 0);

  f0pk1 = malloc_double_2d(BS, 0, BS, 0);
  f1pk0 = malloc_double_2d(BS, 0, BS, 0);
  f1pk2 = malloc_double_2d(BS, 0, BS, 0);
  f2pk1 = malloc_double_2d(BS, 0, BS, 0);
  f0pk2 = malloc_double_2d(BS, 0, BS, 0);
  f2pk0 = malloc_double_2d(BS, 0, BS, 0); 
  
  fn = malloc_double_3d(NZ, 0, 2 * PSF_HS + NX, PSF_HS, 2 * PSF_HS + NY, PSF_HS);
  f0p1delta_d = malloc_double_2d(BS, 0, BS, 0);
  f1p0delta_d = malloc_double_2d(BS, 0, BS, 0);
  f1p2delta_d = malloc_double_2d(BS, 0, BS, 0);
  f2p1delta_d = malloc_double_2d(BS, 0, BS, 0);
  f0p2delta_d = malloc_double_2d(BS, 0, BS, 0);
  f2p0delta_d = malloc_double_2d(BS, 0, BS, 0);
  
  f0p1delta_k = malloc_double_2d(BS, 0, BS, 0);
  f1p0delta_k = malloc_double_2d(BS, 0, BS, 0);
  f1p2delta_k = malloc_double_2d(BS, 0, BS, 0);
  f2p1delta_k = malloc_double_2d(BS, 0, BS, 0);
  f0p2delta_k = malloc_double_2d(BS, 0, BS, 0);
  f2p0delta_k = malloc_double_2d(BS, 0, BS, 0);

  p0d = malloc_double_2d(2 * PSF_HS + 1, PSF_HS,  2 * PSF_HS + 1, PSF_HS);
  p0k = malloc_double_2d(2 * PSF_HS + 1, PSF_HS,  2 * PSF_HS + 1, PSF_HS);
  p1d = malloc_double_2d(2 * PSF_HS + 1, PSF_HS,  2 * PSF_HS + 1, PSF_HS);
  p1k = malloc_double_2d(2 * PSF_HS + 1, PSF_HS,  2 * PSF_HS + 1, PSF_HS);
  p2d = malloc_double_2d(2 * PSF_HS + 1, PSF_HS,  2 * PSF_HS + 1, PSF_HS);
  p2k = malloc_double_2d(2 * PSF_HS + 1, PSF_HS,  2 * PSF_HS + 1, PSF_HS);

  f01delta_d = malloc_double_2d(BS, 0, BS, 0);
  f10delta_d = malloc_double_2d(BS, 0, BS, 0);
  f12delta_d = malloc_double_2d(BS, 0, BS, 0);
  f21delta_d = malloc_double_2d(BS, 0, BS, 0);
  f02delta_d = malloc_double_2d(BS, 0, BS, 0);
  f20delta_d = malloc_double_2d(BS, 0, BS, 0);

  f01delta_k = malloc_double_2d(BS, 0, BS, 0);
  f10delta_k = malloc_double_2d(BS, 0, BS, 0);
  f12delta_k = malloc_double_2d(BS, 0, BS, 0);
  f21delta_k = malloc_double_2d(BS, 0, BS, 0);
  f02delta_k = malloc_double_2d(BS, 0, BS, 0);
  f20delta_k = malloc_double_2d(BS, 0, BS, 0);
  //-----------------------------------------------------------------------------
  // image input
  for(zn = 0 ; zn < NZ ; zn++){
    sprintf(flnm, "./k05/d0515/d0515_ra_0%d.bin", zn);
    fp = fopen(flnm, "r");
    nr = fread(data, sizeof(double), NX*NY, fp);
    //printf("nr = %d\n", nr);
    nw = NX*NY;
    //printf("nw = %d\n", nw);
    if(nr != nw){
      fprintf(stderr, "\nError : fread\n\n");exit(1);
    }
    for(m = 0 ; m < NY ; m++){
      for(n = 0 ; n < NX ; n++){
	fn[zn][n][m] = data[m * NX + n];
      }
    }
    fclose(fp);
    
    // fn SAVE
    cIMAGE(NX, NY, &in_img, MONO);
    sprintf(flnm, "d0515_%d.dat", zn);
    fp = fopen(flnm, "w");
    for(m = 0 ; m < NY ; m++){
      for(n = 0 ; n < NX ; n++){
	fprintf(fp, "%4d %4d  %f\n", n, m, fn[zn][n][m]);
	D(in_img, n, m) = (unsigned char)(fn[zn][n][m] + 0.5);
      }
    }
    fclose(fp);
    sprintf(flnm, "d0515_%d.ppm", zn);
    wIMAGE(in_img, flnm);
  }

  // orikaeshi
  

  //-----------------------------------------------------------------------------
  /*
  psf(k, z, d, sgmc, PSF_HS, p);
  psf(k, z, d+DELTA, sgmc, PSF_HS, pd);
  psf(k+DELTA, z, d, sgmc, PSF_HS, pk);
  psf(k+DELTA, z, d+DELTA, sgmc, PSF_HS, pdk);

  dpsfdd(k, z, d, sgmc, PSF_HS, dpdd);
  dpsfdk(k, z, d, sgmc, PSF_HS, dpdk);
  */
  //-----------------------------------------------------------------------------
  /*
  for(m = -PSF_HS ; m < PSF_HS + BS ; m++) {
    for(n = -PSF_HS ; n < PSF_HS + BS; n++) {
      f[n][m] = (double)n  + 128.0;
    }
  }
  
  for(j = 0 ; j < BS ; j++) {
    for(i = 0 ; i < BS ; i++) {
	  fd[i][j] = 0.0;
	  psff[i][j] = 0.0;
	  psfd[i][j] = 0.0;
    }
  }
   
  for(j = 0 ; j < BS ; j++) {
    for(i = 0 ; i < BS ; i++) {
      for(m = -PSF_HS ; m <= PSF_HS ; m++) {
	for(n = -PSF_HS ; n <= PSF_HS ; n++) {
	  tmp1 += f[i-n][j-m] * dpdd[n][m];
	  tmp2 += f[i-n][j-m] * p[n][m];
	  tmp3 += f[i-n][j-m] * pd[n][m];
	  //printf("%f ", tmp1);
	}
      }
      fd[i][j] = tmp1;
      psff[i][j] = tmp2;
      psfd[i][j] = tmp3;
      tmp1 = tmp2 = tmp3 = 0.0;
      //printf("%f ", fd[i][j]);
    }
   }

  for(j = 0 ; j < BS ; j++) {
    for(i = 0 ; i < BS ; i++) {
      sabun[i][j] = fabs(psff[i][j] - psfd[i][j]);
    }
   }

  for(j = -PSF_HS ; j <= PSF_HS ; j++) {
    for(i = -PSF_HS ; i <= PSF_HS ; i++) {
      kekkad[i][j] = (pd[i][j] - p[i][j]) / DELTA;
      kekkak[i][j] = (pk[i][j] - p[i][j]) / DELTA;
    }
  }
  */
  //-----------------------------------------------------------------------------
  /*
  make_fnpn(BS, (BS-1)/2, (BS-1)/2, fn[0], p, f0p1);
  make_fnpn(BS, (BS-1)/2, (BS-1)/2, fn[0], pk, f0p1delta);
  make_fnpkn(BS, (BS-1)/2, (BS-1)/2, fn[0], f0pk1, k, z, d, SGMC);
  
  for(j = 0 ; j < BS ; j++) {
    for(i = 0 ; i < BS ; i++) {
      sabun[i][j] = (f0p1delta[i][j] - f0p1[i][j])/DELTA;
    }
  }
  tmp1 = RMSE = 0.0;
  for(j = 0 ; j < BS ; j++) {
    for(i = 0 ; i < BS ; i++) {
      tmp1 = sabun[i][j] - f0pk1[i][j];
      RMSE += tmp1 * tmp1;
    }
  }
  RMSE = sqrt(RMSE / (double)(BS*BS));
  printf("RMSE = %f\n", RMSE);
  */
  //-----------------------------------------------------------------------------

  k0 = 0.6, d0 = 1.2;  

  count = 0;
  
  do{
    psf(k0, 0, d0, sgmc, PSF_HS, p0);
    psf(k0, 1, d0, sgmc, PSF_HS, p1);
    psf(k0, 2, d0, sgmc, PSF_HS, p2);
    
    make_fnpn(BS, CENTER_X, CENTER_Y, fn[0], p1, f0p1);
    make_fnpn(BS, CENTER_X, CENTER_Y, fn[1], p0, f1p0);
    make_fnpn(BS, CENTER_X, CENTER_Y, fn[1], p2, f1p2);
    make_fnpn(BS, CENTER_X, CENTER_Y, fn[2], p1, f2p1);
    make_fnpn(BS, CENTER_X, CENTER_Y, fn[0], p2, f0p2);
    make_fnpn(BS, CENTER_X, CENTER_Y, fn[2], p0, f2p0);    
    
    make_fnpdn(BS, CENTER_X, CENTER_Y, fn[0], f0pd1, k0, 1, d0, SGMC);
    make_fnpdn(BS, CENTER_X, CENTER_Y, fn[1], f1pd0, k0, 0, d0, SGMC);
    make_fnpdn(BS, CENTER_X, CENTER_Y, fn[1], f1pd2, k0, 2, d0, SGMC);
    make_fnpdn(BS, CENTER_X, CENTER_Y, fn[2], f2pd1, k0, 1, d0, SGMC);
    make_fnpdn(BS, CENTER_X, CENTER_Y, fn[0], f0pd2, k0, 2, d0, SGMC);
    make_fnpdn(BS, CENTER_X, CENTER_Y, fn[2], f2pd0, k0, 0, d0, SGMC);
        
    make_fnpkn(BS, CENTER_X, CENTER_Y, fn[0], f0pk1, k0, 1, d0, SGMC);
    make_fnpkn(BS, CENTER_X, CENTER_Y, fn[1], f1pk0, k0, 0, d0, SGMC);
    make_fnpkn(BS, CENTER_X, CENTER_Y, fn[1], f1pk2, k0, 2, d0, SGMC);
    make_fnpkn(BS, CENTER_X, CENTER_Y, fn[2], f2pk1, k0, 1, d0, SGMC);
    make_fnpkn(BS, CENTER_X, CENTER_Y, fn[0], f0pk2, k0, 2, d0, SGMC);
    make_fnpkn(BS, CENTER_X, CENTER_Y, fn[2], f2pk0, k0, 0, d0, SGMC);
    
    // NDM
 #if NDM == 1
    psf(k0, 0, d0 + DELTA, sgmc, PSF_HS, p0d);
    psf(k0 + DELTA, 0, d0, sgmc, PSF_HS, p0k);
    psf(k0, 1, d0 + DELTA, sgmc, PSF_HS, p1d);
    psf(k0 + DELTA, 1, d0, sgmc, PSF_HS, p1k);
    psf(k0, 2, d0 + DELTA, sgmc, PSF_HS, p2d);
    psf(k0 + DELTA, 2, d0, sgmc, PSF_HS, p2k);
    
    /*
    num_diff_method_p(p0d, p0, p_delta[0]);
    num_diff_method_p(p1d, p1, p_delta[1]);
    num_diff_method_p(p2d, p2, p_delta[2]);
    num_diff_method_p(p0k, p0, p_delta[3]);
    num_diff_method_p(p1k, p1, p_delta[4]);
    num_diff_method_p(p2k, p2, p_delta[5]);
    
    make_fnpn(BS, CENTER_X, CENTER_Y, fn[0], p_delta[1], f01delta_d);
    make_fnpn(BS, CENTER_X, CENTER_Y, fn[1], p_delta[0], f10delta_d);
    make_fnpn(BS, CENTER_X, CENTER_Y, fn[1], p_delta[2], f12delta_d);
    make_fnpn(BS, CENTER_X, CENTER_Y, fn[2], p_delta[1], f21delta_d);
    make_fnpn(BS, CENTER_X, CENTER_Y, fn[0], p_delta[2], f02delta_d);
    make_fnpn(BS, CENTER_X, CENTER_Y, fn[2], p_delta[0], f20delta_d);

    make_fnpn(BS, CENTER_X, CENTER_Y, fn[0], p_delta[4], f01delta_k);
    make_fnpn(BS, CENTER_X, CENTER_Y, fn[1], p_delta[3], f10delta_k);
    make_fnpn(BS, CENTER_X, CENTER_Y, fn[1], p_delta[5], f12delta_k);
    make_fnpn(BS, CENTER_X, CENTER_Y, fn[2], p_delta[4], f21delta_k);
    make_fnpn(BS, CENTER_X, CENTER_Y, fn[0], p_delta[5], f02delta_k);
    make_fnpn(BS, CENTER_X, CENTER_Y, fn[2], p_delta[3], f20delta_k); 
    */
    
    make_fnpn(BS, CENTER_X, CENTER_Y, fn[0], p1d, f0p1delta_d);
    make_fnpn(BS, CENTER_X, CENTER_Y, fn[1], p0d, f1p0delta_d);
    make_fnpn(BS, CENTER_X, CENTER_Y, fn[1], p2d, f1p2delta_d);
    make_fnpn(BS, CENTER_X, CENTER_Y, fn[2], p1d, f2p1delta_d);
    make_fnpn(BS, CENTER_X, CENTER_Y, fn[0], p2d, f0p2delta_d);
    make_fnpn(BS, CENTER_X, CENTER_Y, fn[2], p0d, f2p0delta_d);

    make_fnpn(BS, CENTER_X, CENTER_Y, fn[0], p1k, f0p1delta_k);
    make_fnpn(BS, CENTER_X, CENTER_Y, fn[1], p0k, f1p0delta_k);
    make_fnpn(BS, CENTER_X, CENTER_Y, fn[1], p2k, f1p2delta_k);
    make_fnpn(BS, CENTER_X, CENTER_Y, fn[2], p1k, f2p1delta_k);
    make_fnpn(BS, CENTER_X, CENTER_Y, fn[0], p2k, f0p2delta_k);
    make_fnpn(BS, CENTER_X, CENTER_Y, fn[2], p0k, f2p0delta_k);
    
    num_diff_method(BS, f0p1delta_d, f0p1, f01delta_d);
    num_diff_method(BS, f1p0delta_d, f1p0, f10delta_d);
    num_diff_method(BS, f1p2delta_d, f1p2, f12delta_d);
    num_diff_method(BS, f2p1delta_d, f2p1, f21delta_d);
    num_diff_method(BS, f0p2delta_d, f0p2, f02delta_d);
    num_diff_method(BS, f2p0delta_d, f2p0, f20delta_d);

    num_diff_method(BS, f0p1delta_k, f0p1, f01delta_k);
    num_diff_method(BS, f1p0delta_k, f1p0, f10delta_k);
    num_diff_method(BS, f1p2delta_k, f1p2, f12delta_k);
    num_diff_method(BS, f2p1delta_k, f2p1, f21delta_k);
    num_diff_method(BS, f0p2delta_k, f0p2, f02delta_k);
    num_diff_method(BS, f2p0delta_k, f2p0, f20delta_k);
    
    tmp1 = RMSE = 0.0;
    for(j = 0 ; j < BS ; j++) {
      for(i = 0 ; i < BS ; i++) {
	tmp1 = f10delta_k[i][j] * f10delta_d[i][j] - f1pk0[i][j] * f1pd0[i][j];
	RMSE += tmp1 * tmp1;
      }
    }
    RMSE = sqrt(RMSE / (double)(BS*BS));
    printf("                  RMSE = %f\n", RMSE);
 #endif
    for(zn = 0 ; zn > NZ ; zn++){
      dedd[zn] = 0.0;
      dedk[zn] = 0.0;
      e[zn] = 0.0;
    }

    Did = Dik = 0.0;
    t1 = 0.0, t2 = 0.0, t3 = 0.0, t4 = 0.0, t5 = 0.0, t6 = 0.0, t7 = 0.0, t8 = 0.0, t9 = 0.0;
    t10 = t11 = t12 = 0.0;
    tf0p1d = tf1p0d = tf1p2d = tf2p1d = tf0p2d = tf2p0d = 0.0;
    tf0p1dk = tf1p0dk = tf1p2dk = tf2p1dk = tf0p2dk = tf2p0dk = 0.0;
    tf0p1k = tf1p0k = tf1p2k = tf2p1k = tf0p2k = tf2p0k = 0.0;
    te1f0p1d = te1f1p0d = te2f1p2d = te2f2p1d = te3f0p2d = te3f2p0d = 0.0;
    te1f0p1k = te1f1p0k = te2f1p2k = te2f2p1k = te3f0p2k = te3f2p0k = 0.0;
    for(j = 0 ; j < BS ; j++) {
      for(i = 0 ; i < BS ; i++) {
	t1 = (f0p1[i][j] - f1p0[i][j]);
	t2 = (f1p2[i][j] - f2p1[i][j]);
	t3 = (f0p2[i][j] - f2p0[i][j]);
	
	t4 += (f0pd1[i][j] - f1pd0[i][j]) * (f0pd1[i][j] - f1pd0[i][j]);
	t5 += (f1pd2[i][j] - f2pd1[i][j]) * (f1pd2[i][j] - f2pd1[i][j]);  
	t6 += (f0pd2[i][j] - f2pd0[i][j]) * (f0pd2[i][j] - f2pd0[i][j]);
	
	t7 += (f0pk1[i][j] - f1pk0[i][j]) * (f0pk1[i][j] - f1pk0[i][j]);
	t8 += (f1pk2[i][j] - f2pk1[i][j]) * (f1pk2[i][j] - f2pk1[i][j]);  
	t9 += (f0pk2[i][j] - f2pk0[i][j]) * (f0pk2[i][j] - f2pk0[i][j]);

	t10 += (f0pd1[i][j] - f1pd0[i][j]) * (f0pk1[i][j] - f1pk0[i][j]);
	t11 += (f1pd2[i][j] - f2pd1[i][j]) * (f1pk2[i][j] - f2pk1[i][j]);  
	t12 += (f0pd2[i][j] - f2pd0[i][j]) * (f0pk2[i][j] - f2pk0[i][j]);

	Did += t1 * (f0pd1[i][j] - f1pd0[i][j]) + t2 * (f1pd2[i][j] - f2pd1[i][j]) + t3 * (f0pd2[i][j] - f2pd0[i][j]);
	Dik += t1 * (f0pk1[i][j] - f1pk0[i][j]) + t2 * (f1pk2[i][j] - f2pk1[i][j]) + t3 * (f0pk2[i][j] - f2pk0[i][j]);
	/*
	// a00
	tf0p1d += f0pd1[i][j] * f0pd1[i][j];
	tf1p0d += f1pd0[i][j] * f1pd0[i][j];
	tf1p2d += f1pd2[i][j] * f1pd2[i][j];
	tf2p1d += f2pd1[i][j] * f2pd1[i][j];
	tf0p2d += f0pd2[i][j] * f0pd2[i][j];
	tf2p0d += f2pd0[i][j] * f2pd0[i][j];

	// a01 a10
	tf0p1dk += f0pk1[i][j] * f0pd1[i][j];
	tf1p0dk += f1pk0[i][j] * f1pd0[i][j];
	tf1p2dk += f1pk2[i][j] * f1pd2[i][j];
	tf2p1dk += f2pk1[i][j] * f2pd1[i][j];
	tf0p2dk += f0pk2[i][j] * f0pd2[i][j];
	tf2p0dk += f2pk0[i][j] * f2pd0[i][j];

	// a11
	tf0p1k += f0pk1[i][j] * f0pk1[i][j];
	tf1p0k += f1pk0[i][j] * f1pk0[i][j];
	tf1p2k += f1pk2[i][j] * f1pk2[i][j];
	tf2p1k += f2pk1[i][j] * f2pk1[i][j];
	tf0p2k += f0pk2[i][j] * f0pk2[i][j];
	tf2p0k += f2pk0[i][j] * f2pk0[i][j];

	// b0
	te1f0p1d += f0p1[i][j] * f0pd1[i][j];
	te1f1p0d += f1p0[i][j] * f1pd0[i][j];
	te2f1p2d += f1p2[i][j] * f1pd2[i][j];
	te2f2p1d += f2p1[i][j] * f2pd1[i][j];
	te3f0p2d += f0p2[i][j] * f0pd2[i][j];
	te3f2p0d += f2p0[i][j] * f2pd0[i][j];

	// b1
	te1f0p1k += f0p1[i][j] * f0pk1[i][j];
	te1f1p0k += f1p0[i][j] * f1pk0[i][j];
	te2f1p2k += f1p2[i][j] * f1pk2[i][j];
	te2f2p1k += f2p1[i][j] * f2pk1[i][j];
	te3f0p2k += f0p2[i][j] * f0pk2[i][j];
	te3f2p0k += f2p0[i][j] * f2pk0[i][j];
	*/
	/*
	// NDM
	//a00
	tf0p1d += f01delta_d[i][j] * f01delta_d[i][j];
	tf1p0d += f10delta_d[i][j] * f10delta_d[i][j];
	tf1p2d += f12delta_d[i][j] * f12delta_d[i][j];
	tf2p1d += f21delta_d[i][j] * f21delta_d[i][j];
	tf0p2d += f02delta_d[i][j] * f02delta_d[i][j];
	tf2p0d += f20delta_d[i][j] * f20delta_d[i][j];

	// a01, a10
	tf0p1dk += f01delta_k[i][j] * f01delta_d[i][j];
	tf1p0dk += f10delta_k[i][j] * f10delta_d[i][j];
	tf1p2dk += f12delta_k[i][j] * f12delta_d[i][j];
	tf2p1dk += f21delta_k[i][j] * f21delta_d[i][j];
	tf0p2dk += f02delta_k[i][j] * f02delta_d[i][j];
	tf2p0dk += f20delta_k[i][j] * f20delta_d[i][j];

	// a11
	tf0p1k += f01delta_k[i][j] * f01delta_k[i][j];
	tf1p0k += f10delta_k[i][j] * f10delta_k[i][j];
	tf1p2k += f12delta_k[i][j] * f12delta_k[i][j];
	tf2p1k += f21delta_k[i][j] * f21delta_k[i][j];
	tf0p2k += f02delta_k[i][j] * f02delta_k[i][j];
	tf2p0k += f20delta_k[i][j] * f20delta_k[i][j];

	// b0 
	te1f0p1d += f0p1[i][j] * f01delta_d[i][j];
	te1f1p0d += f1p0[i][j] * f10delta_d[i][j];
	te2f1p2d += f1p2[i][j] * f12delta_d[i][j];
	te2f2p1d += f2p1[i][j] * f21delta_d[i][j];
	te3f0p2d += f0p2[i][j] * f02delta_d[i][j];
	te3f2p0d += f2p0[i][j] * f20delta_d[i][j];

	// b1
	te1f0p1k += f0p1[i][j] * f01delta_k[i][j];
	te1f1p0k += f1p0[i][j] * f10delta_k[i][j];
	te2f1p2k += f1p2[i][j] * f12delta_k[i][j];
	te2f2p1k += f2p1[i][j] * f21delta_k[i][j];
	te3f0p2k += f0p2[i][j] * f02delta_k[i][j];
	te3f2p0k += f2p0[i][j] * f20delta_k[i][j];
	*/
      }
    }
    
    e[0] = t1;
    e[1] = t2;
    e[2] = t3;
    a00 = a01 = a10 = a11 = b0 = b1 = 0.0;
    a00 = (tf0p1d - tf1p0d) + (tf1p2d - tf2p1d) + (tf0p2d - tf2p0d);
    a01 = (tf0p1dk - tf1p0dk) + (tf1p2dk - tf2p1dk) + (tf0p2dk - tf2p0dk);
    a10 = a01;
    a11 = (tf0p1k - tf1p0k) + (tf1p2k - tf2p1k) + (tf0p2k - tf2p0k);
    b0 = (te1f0p1d - te1f1p0d) + (te2f1p2d - te2f2p1d) + (te3f0p2d - te3f2p0d);
    b1 = (te1f0p1k - te1f1p0k) + (te2f1p2k - te2f2p1k) + (te3f0p2k - te3f2p0k);

    dedd[0] = t4;//d01
    dedd[1] = t5;//d12
    dedd[2] = t6;//d02
    dedk[0] = t10;
    dedk[1] = t11;
    dedk[2] = t12;
    dekk[0] = t7;
    dekk[1] = t8;
    dekk[2] = t9;
    
    printf("E[%d] = %f\n", count, e[0] + e[1] + e[2]);
    printf("e0: %f\n", e[0]);
    printf("e1: %f\n", e[1]);
    printf("e2: %f\n", e[2]);

    eprev = 0.0;

    //printf("a00 :%f, a10 :%f, a11 :%f\n", a00, a10, a11);
    
    Ddd =  dedd[0] + dedd[1] + dedd[2];//a
    Ddk =  dedk[0] + dedk[1] + dedk[2];//b,c
    Dkk =  dekk[0] + dekk[1] + dekk[2];//d
    
    //Did = -(e[0] * dedd[0] + e[1] * dedd[1] + e[2] * dedd[2]);//e
    //Dik = -(e[0] * dedk[0] + e[1] * dedk[1] + e[2] * dedk[2]);//f
    Did = (-1)*Did;
    Dik = (-1)*Dik;
   
    
    deltad = (Dkk * Did - Ddk * Dik) / (Ddd * Dkk - Ddk * Ddk);
    
    //deltad = 0.0;
    //deltad = (a01 * b1 - a11 * b0) / (a00 * a11 - a01 * a10);
    dprev = d0;
    d0 = dprev + deltad;
    
    printf("deltad: %f\n", deltad);
    printf("d0: %f\n", d0);
    
    deltak = 0.0;
    deltak = (Ddd * Dik - Ddk * Did) / (Ddd * Dkk - Ddk * Ddk);
    //deltak = (a10 * b0 - a00 * b1) / (a00 * a11 - a01 * a10);
    kprev = k0;
    k0 = kprev + deltak;
    
    printf("deltak: %f\n", deltak);
    printf("k0: %f\n", k0);
    /*
    //----------------------------------------------------
    psf(k0, 0, d0, sgmc, PSF_HS, p0);
    psf(k0, 1, d0, sgmc, PSF_HS, p1);
    psf(k0, 2, d0, sgmc, PSF_HS, p2);	
    psf(k0 + DELTA, 0, d0, sgmc, PSF_HS, p0k);
    psf(k0 + DELTA, 1, d0, sgmc, PSF_HS, p1k);
    psf(k0 + DELTA, 2, d0, sgmc, PSF_HS, p2k);

    make_fnpn(BS, CENTER_X, CENTER_Y, fn[0], p1, f0p1);
    make_fnpn(BS, CENTER_X, CENTER_Y, fn[1], p0, f1p0);
    make_fnpn(BS, CENTER_X, CENTER_Y, fn[1], p2, f1p2);
    make_fnpn(BS, CENTER_X, CENTER_Y, fn[2], p1, f2p1);
    make_fnpn(BS, CENTER_X, CENTER_Y, fn[0], p2, f0p2);
    make_fnpn(BS, CENTER_X, CENTER_Y, fn[2], p0, f2p0);     
    
    make_fnpn(BS, CENTER_X, CENTER_Y, fn[0], p1k, f0p1delta_k);
    make_fnpn(BS, CENTER_X, CENTER_Y, fn[1], p0k, f1p0delta_k);
    make_fnpn(BS, CENTER_X, CENTER_Y, fn[1], p2k, f1p2delta_k);
    make_fnpn(BS, CENTER_X, CENTER_Y, fn[2], p1k, f2p1delta_k);
    make_fnpn(BS, CENTER_X, CENTER_Y, fn[0], p2k, f0p2delta_k);
    make_fnpn(BS, CENTER_X, CENTER_Y, fn[2], p0k, f2p0delta_k);   

    num_diff_method(BS, f0p1delta_k, f0p1, f01delta_k);
    num_diff_method(BS, f1p0delta_k, f1p0, f10delta_k);
    num_diff_method(BS, f1p2delta_k, f1p2, f12delta_k);
    num_diff_method(BS, f2p1delta_k, f2p1, f21delta_k);
    num_diff_method(BS, f0p2delta_k, f0p2, f02delta_k);
    num_diff_method(BS, f2p0delta_k, f2p0, f20delta_k);


    double ea, eb;

    for(j = 0 ; j < BS ; j++) {
      for(i = 0 ; i < BS ; i++) {   
	ea = (f0p1[i][j] - f1p0[i][j]) * (f0p1[i][j] - f1p0[i][j]) +
	     (f1p2[i][j] - f2p1[i][j]) * (f1p2[i][j] - f2p1[i][j]) +
	     (f0p2[i][j] - f2p0[i][j]) * (f0p2[i][j] - f2p0[i][j]);
    
	eb = (f0p1delta_k[i][j] - f1p0delta_k[i][j]) * (f0p1delta_k[i][j] - f1p0delta_k[i][j]) +
	     (f1p2delta_k[i][j] - f2p1delta_k[i][j]) * (f1p2delta_k[i][j] - f2p1delta_k[i][j]) +
	     (f0p2delta_k[i][j] - f2p0delta_k[i][j]) * (f0p2delta_k[i][j] - f2p0delta_k[i][j]);

      }
    }
    //printf("E_bibun(k) = %f\n", ((eb - ea) / DELTA));
	       
    //----------------------------------------------------

    for(d0 = 0.8 ; d0 <= 1.2 ; d0 += 0.01){

    psf(k0, 0, d0, sgmc, PSF_HS, p0);
    psf(k0, 1, d0, sgmc, PSF_HS, p1);
    psf(k0, 2, d0, sgmc, PSF_HS, p2);	
    psf(k0, 0, d0 + DELTA, sgmc, PSF_HS, p0k);
    psf(k0, 1, d0 + DELTA, sgmc, PSF_HS, p1k);
    psf(k0, 2, d0 + DELTA, sgmc, PSF_HS, p2k);

    make_fnpn(BS, CENTER_X, CENTER_Y, fn[0], p1, f0p1);
    make_fnpn(BS, CENTER_X, CENTER_Y, fn[1], p0, f1p0);
    make_fnpn(BS, CENTER_X, CENTER_Y, fn[1], p2, f1p2);
    make_fnpn(BS, CENTER_X, CENTER_Y, fn[2], p1, f2p1);
    make_fnpn(BS, CENTER_X, CENTER_Y, fn[0], p2, f0p2);
    make_fnpn(BS, CENTER_X, CENTER_Y, fn[2], p0, f2p0);     
    
    make_fnpn(BS, CENTER_X, CENTER_Y, fn[0], p1k, f0p1delta_k);
    make_fnpn(BS, CENTER_X, CENTER_Y, fn[1], p0k, f1p0delta_k);
    make_fnpn(BS, CENTER_X, CENTER_Y, fn[1], p2k, f1p2delta_k);
    make_fnpn(BS, CENTER_X, CENTER_Y, fn[2], p1k, f2p1delta_k);
    make_fnpn(BS, CENTER_X, CENTER_Y, fn[0], p2k, f0p2delta_k);
    make_fnpn(BS, CENTER_X, CENTER_Y, fn[2], p0k, f2p0delta_k);   

    num_diff_method(BS, f0p1delta_k, f0p1, f01delta_k);
    num_diff_method(BS, f1p0delta_k, f1p0, f10delta_k);
    num_diff_method(BS, f1p2delta_k, f1p2, f12delta_k);
    num_diff_method(BS, f2p1delta_k, f2p1, f21delta_k);
    num_diff_method(BS, f0p2delta_k, f0p2, f02delta_k);
    num_diff_method(BS, f2p0delta_k, f2p0, f20delta_k);


    for(j = 0 ; j < BS ; j++) {
      for(i = 0 ; i < BS ; i++) {   
	ea = (f0p1[i][j] - f1p0[i][j]) * (f0p1[i][j] - f1p0[i][j]) +
	     (f1p2[i][j] - f2p1[i][j]) * (f1p2[i][j] - f2p1[i][j]) +
	     (f0p2[i][j] - f2p0[i][j]) * (f0p2[i][j] - f2p0[i][j]);
    
	eb = (f0p1delta_k[i][j] - f1p0delta_k[i][j]) * (f0p1delta_k[i][j] - f1p0delta_k[i][j]) +
	     (f1p2delta_k[i][j] - f2p1delta_k[i][j]) * (f1p2delta_k[i][j] - f2p1delta_k[i][j]) +
	     (f0p2delta_k[i][j] - f2p0delta_k[i][j]) * (f0p2delta_k[i][j] - f2p0delta_k[i][j]);

      }
    }
    //printf("E_bibun(k) = %f\n", ((eb - ea) / DELTA));
    printf("d = %f, ea = %f, eb = %f\n", d0, ea, eb);
    //printf("d = %f, ea = %f\n",d0, ea);
    
    }*/

    /*----------------------------------------------------*/

    //   if(count == 1)break;   
    count++;
    //} while(fabs(e[0]+e[1]+e[2]) > 1.0E-5);   
    //} while(fabs(deltak) > 1.0E-5);
    printf("\n");
    } while( fabs( (dprev - d) - (d0 - d) ) / fabs(dprev - d) > 0.1);
    //}while(fabs(deltad) > 1.0E-8);
    //}while(fabs());
    //} while( fabs( (kprev - k) - (k0 - k) ) / fabs(kprev - k) > HANTEI);
  
  printf("--------------------------------------------\n");
  printf("loop = %d\n", count);
  printf("d = %f, k = %f\n", d0, k0);
  printf("dprev = %f\n", dprev);
  printf("deltad = %f, deltak = %f\n", deltad, deltak);
  printf("gosa = %f\n", d0 - d);

  
  //-----------------------------------------------------------------------------
  
  fp = fopen("kekkad.dat", "w");
  for(j = -PSF_HS ; j <= PSF_HS ; j++) {
    for(i = -PSF_HS ; i <= PSF_HS ; i++) {
      fprintf(fp, "%2d %2d  %f\n", i, j, kekkad[i][j]);
    }
    fprintf(fp, "\n");
  }
  fclose(fp);

  fp = fopen("kekkak.dat", "w");
  for(j = -PSF_HS ; j <= PSF_HS ; j++) {
    for(i = -PSF_HS ; i <= PSF_HS ; i++) {
      fprintf(fp, "%2d %2d  %f\n", i, j, kekkak[i][j]);
    }
    fprintf(fp, "\n");
  }
  fclose(fp);

    fp = fopen("dpdd.dat", "w");
  for(j = -PSF_HS ; j <= PSF_HS ; j++) {
    for(i = -PSF_HS ; i <= PSF_HS ; i++) {
      fprintf(fp, "%2d %2d  %f\n", i, j, dpdd[i][j]);
    }
    fprintf(fp, "\n");
  }
  fclose(fp);

      fp = fopen("dpdk.dat", "w");
  for(j = -PSF_HS ; j <= PSF_HS ; j++) {
    for(i = -PSF_HS ; i <= PSF_HS ; i++) {
      fprintf(fp, "%2d %2d  %f\n", i, j, dpdk[i][j]);
    }
    fprintf(fp, "\n");
  }
  fclose(fp);
  
  fp = fopen("psf.dat", "w");
  for(j = -PSF_HS ; j <= PSF_HS ; j++) {
    for(i = -PSF_HS ; i <= PSF_HS ; i++) {
      fprintf(fp, "%2d %2d  %f\n", i, j, p[i][j]);
    }
    fprintf(fp, "\n");
  }
  fclose(fp);

  fp = fopen("f.dat", "w");
  for(j = -PSF_HS ; j < PSF_HS + BS ; j++) {
    for(i = -PSF_HS ; i < PSF_HS + BS ; i++) {
      fprintf(fp, "%2d %2d  %f\n", i, j, f[i][j]);
    }
    fprintf(fp, "\n");
  }
  fclose(fp);

  fp = fopen("sabun.dat", "w");
  for(j = 0 ; j < BS ; j++) {
    for(i = 0 ; i < BS ; i++) {
      fprintf(fp, "%2d %2d  %f\n", i, j, sabun[i][j]);
    }
    fprintf(fp, "\n");
  }
  fclose(fp);

  fp = fopen("fd.dat", "w");
  for(j = 0 ; j < BS ; j++) {
    for(i = 0 ; i < BS ; i++) {
      fprintf(fp, "%2d %2d  %f\n", i, j, fd[i][j]);
    }
    fprintf(fp, "\n");
  }
  fclose(fp);

  fp = fopen("psff.dat", "w");
  for(j = 0 ; j < BS ; j++) {
    for(i = 0 ; i < BS ; i++) {
      fprintf(fp, "%2d %2d  %f\n", i, j, psff[i][j]);
    }
    fprintf(fp, "\n");
  }
  fclose(fp);

    fp = fopen("psfd.dat", "w");
  for(j = 0 ; j < BS ; j++) {
    for(i = 0 ; i < BS ; i++) {
      fprintf(fp, "%2d %2d  %f\n", i, j, psfd[i][j]);
    }
    fprintf(fp, "\n");
  }
  fclose(fp);

  fp = fopen("f0p1.dat", "w");
  for(j = 0 ; j < BS ; j++) {
    for(i = 0 ; i < BS ; i++) {
      fprintf(fp, "%2d %2d  %f\n", i, j, f0p1[i][j]);
    }
    fprintf(fp, "\n");
  }
  fclose(fp);
  /*
  fp = fopen("f0p1delta.dat", "w");
  for(j = 0 ; j < BS ; j++) {
    for(i = 0 ; i < BS ; i++) {
      fprintf(fp, "%2d %2d  %f\n", i, j, f0p1delta[i][j]);
    }
    fprintf(fp, "\n");
  }
  fclose(fp);
  *
  fp = fopen("f0pk1.dat", "w");
  for(j = 0 ; j < BS ; j++) {
    for(i = 0 ; i < BS ; i++) {
      fprintf(fp, "%2d %2d  %f\n", i, j, f0pk1[i][j]);
    }
    fprintf(fp, "\n");
  }
  fclose(fp);
  */
  free_double_2d(p, 2 * PSF_HS + 1, PSF_HS,  2 * PSF_HS + 1, PSF_HS);
  
  return 0;
}
  
/*---------------------------------------------------------------------*/
void psf(double k, int z, double d, double sgmc, int psf_hs, double **p)
{
  int i, j;
  double c, K, K2;
 
  K = k * ((double)z - d);
  K2 = K * K + sgmc * sgmc;

  c = 0.0;
  for(j = -PSF_HS ; j <= PSF_HS ; j++) {
    for(i = -PSF_HS ; i <= PSF_HS ; i++) {
      p[i][j] = exp( -(double)(i * i + j * j) / 2.0 / K2);
      c += p[i][j];
    }
  }

  c = 2.0 * M_PI * (K * K + sgmc * sgmc);
  for(j = -PSF_HS ; j <= PSF_HS ; j++) {
    for(i = -PSF_HS ; i <= PSF_HS ; i++) {
      p[i][j] /= c;
    }
  }
}

/*---------------------------------------------------------------------*/
void dpsfdd(double k, int z, double d, double sgmc, int psf_hs, double **dpdd)
{
  int i, j;
  double K, K2;
  double P0, P1, P2, P3, i2j2;


  K = k * ((double)z - d);
  K2 = K * K + sgmc * sgmc;

  P0 = 1.0 / (M_PI * K2 * K2);
  P3 = -k * k * (z - d);
  
  for(j = -PSF_HS ; j <= PSF_HS ; j++) {
    for(i = -PSF_HS ; i <= PSF_HS ; i++) {
      i2j2 = (double)(i * i + j * j);
      P1 = exp(-i2j2 / 2.0 / K2);
      P2 = i2j2 / 2.0 / K2 - 1.0;
      dpdd[i][j] = P0 * P1 * P2 * P3;
    }
  }
}

/*---------------------------------------------------------------------*/
void dpsfdk(double k, int z, double d, double sgmc, int psf_hs, double **dpdk)
{
  int i, j;
  double K, K2;
  double P0, P1, P2, P3, i2j2;


  K = k * ((double)z - d);
  K2 = K * K + sgmc * sgmc;

  P0 = 1.0 / (M_PI * K2 * K2); 
  P3 = k * (z - d) * (z - d);
  
  for(j = -PSF_HS ; j <= PSF_HS ; j++) {
    for(i = -PSF_HS ; i <= PSF_HS ; i++) {
      i2j2 = (double)(i * i + j * j);
      P1 = exp(-i2j2 / 2.0 / K2);
      P2 = i2j2 / 2.0 / K2 - 1.0;
      dpdk[i][j] = P0 * P1 * P2 * P3;
    }
  }
}

/*---------------------------------------------------------------------*/
void make_fnpn(int bs, int center_x, int center_y, double **fn, double **pn, double **fnpn)
{
  int i, j, n, m;
  int hbs = (bs - 1) / 2; // 7
  
  for(j = 0 ; j < bs ; j++) {
    for(i = 0 ; i < bs ; i++) {
      fnpn[i][j] = 0.0;
      for(m = -PSF_HS ; m <= PSF_HS ; m++) {
	for(n = -PSF_HS ; n <= PSF_HS ; n++) {
	  fnpn[i][j] += fn[center_x - hbs + i - n][center_y - hbs + j - m] * pn[n][m];
	}
      }
    }
  }
  
}

/*---------------------------------------------------------------------*/
void make_fnpdn(int bs, int center_x, int center_y, double **fn, double **fnpdn, double k, int z, double d, double sgmc)
{
  int i, j, n, m;
  double **dpdd;
  int hbs = (bs - 1) / 2;
  
  dpdd = malloc_double_2d(2 * PSF_HS + 1, PSF_HS,  2 * PSF_HS + 1, PSF_HS);
  dpsfdd(k, z, d, sgmc, PSF_HS, dpdd);
   
  for(j = 0 ; j < bs ; j++) {
    for(i = 0 ; i < bs ; i++) {
      fnpdn[i][j] = 0.0;
      for(m = -PSF_HS ; m <= PSF_HS ; m++) {
	for(n = -PSF_HS ; n <= PSF_HS ; n++) {
	  fnpdn[i][j] += fn[center_x - hbs + i - n][center_y - hbs + j - m] * dpdd[n][m];
	}
      }
    }
  }
  
}

/*---------------------------------------------------------------------*/
void make_fnpkn(int bs, int center_x, int center_y, double **fn, double **fnpkn, double k, int z, double d, double sgmc)
{
  int i, j, n, m;
  double **dpdk;
  int hbs = (bs - 1) / 2;
  
  dpdk = malloc_double_2d(2 * PSF_HS + 1, PSF_HS,  2 * PSF_HS + 1, PSF_HS);
  dpsfdk(k, z, d, sgmc, PSF_HS, dpdk);
  
  for(j = 0 ; j < bs ; j++) {
    for(i = 0 ; i < bs ; i++) {
      fnpkn[i][j] = 0.0;
      for(m = -PSF_HS ; m <= PSF_HS ; m++) {
	for(n = -PSF_HS ; n <= PSF_HS ; n++) {	  
	  fnpkn[i][j] += fn[center_x  - hbs + i - n][center_y - hbs + j - m] * dpdk[n][m];	  
	}
      }
    }
  }
}

/*---------------------------------------------------------------------*/
void num_diff_method(int bs, double ** fnpndelta, double **fnpn, double **f){
  int i, j;

  for(j = 0 ; j < bs ; j++){
    for(i = 0 ; i < bs ; i++){
      f[i][j] = (fnpndelta[i][j] - fnpn[i][j]) / DELTA;
    }
  }
}

void num_diff_method_p(double ** fnpndelta, double **fnpn, double **f){
  int i, j;

  for(j = -PSF_HS ; j <= PSF_HS ; j++){
    for(i = -PSF_HS ; i <= PSF_HS ; i++){
      f[i][j] = (fnpndelta[i][j] - fnpn[i][j]) / DELTA;
    }
  }
}
