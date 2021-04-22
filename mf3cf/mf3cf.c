/* Dey, Subhadip & Bhattacharya, Avik & Frery, Alejandro & López-Martínez, Carlos. (2020).
A Model-free Four Component Scattering Power Decomposition for Polarimetric SAR Data.

To compile and run, e.g.:
   gcc mf3cf.c -o mf3cf.exe -lm
   ./mf3cf.exe T3

C impl. 20210422 by Ash Richardson, Senior Data Scientist, BC Wildfire Service */
#include<math.h>
#include<stdio.h>
#include<float.h>
#include<stdlib.h>
#include<memory.h>
#include<string.h>
#define for0(i,n) for(i = 0; i < n; i++) /* for loop shorthand */

#define N_IN 9 /* number of input files */
#define T11 0 /* T3 matrix input data indexing */
#define T12_re 1
#define T12_im 2
#define T13_re 3
#define T13_im 4
#define T22 5
#define T23_re 6
#define T23_im 7
#define T33 8

const char* T_fn[] = {"T11.bin", /* T3 matrix input filenames */
                      "T12_real.bin",
                      "T12_imag.bin",
                      "T13_real.bin",
                      "T13_imag.bin",
                      "T22.bin",
                      "T23_real.bin",
                      "T23_imag.bin",
                      "T33.bin"};

float ** T; /* buffers for T3 matrix elements */
float ** T_f; /* buffers for filtered T3 matrix elements */
float ** out_d; /* output data buffers */

#define N_OUT 4
#define _theta_f 0 /* output data indexing */
#define _pd_f 1
#define _ps_f 2
#define _pv_f 3

const char * out_fn[] = {"Theta_FP.bin", /* output filenames */
                         "Pd_FP.bin",
                         "Ps_FP.bin",
                         "Pv_FP.bin"};

char sep(){
  #ifdef _WIN32
    return '\\'; /* windows path separator */
  #else
    return '/'; /* mac/linux/unix path sep */
  #endif
}

void err(const char * m){
  printf("Error: %s\n", m);
  exit(1); /* print message and bail */
}

#define MAX_ARRAYS 1024  /* large number */
int n_arrays = 0; /* count arrays initialized */
void ** arrays; /* track mallocs, free at end */

void * alloc(size_t n){
  void * d = malloc(n); /* create array */
  if(!d) err("failed to allocate memory");
  memset(d, '\0', n); /* must touch memory on windows */
  arrays[n_arrays ++] = d; /* save pointer to free later */
  return d;
}

float * falloc(size_t n){
  return (float *)alloc(n * sizeof(float)); /* float32 array */
}

#define READ 1
#define WRITE 0

FILE * open(const char * fn, int mode){
  printf("+%s %s\n", mode?"r":"w", fn); /* open a file for read or write */
  FILE * f = fopen(fn, mode?"rb":"wb");
  if(!f) err("file access failed");
  return f;
}

void read_config(char * file_name, int * nrow, int * ncol){
  size_t x;
  char tmp[4096]; /* based on PolSARPro by Eric POTTIER and Laurent FERRO-FAMIL */
  FILE * f = open(file_name, READ);
  x = fscanf(f, "%s\n", tmp);
  x = fscanf(f, "%s\n", tmp); // number of rows
  *nrow = atoi(tmp);
  x = fscanf(f, "%s\n", tmp);
  x = fscanf(f, "%s\n", tmp);
  x = fscanf(f, "%s\n", tmp); // number of cols
  *ncol = atoi(tmp);
  fclose(f);
  printf("nrow %d ncol %d\n", *nrow, *ncol);
}

float * read(const char * file_name, size_t n_float){
  FILE * f = open(file_name, READ);
  float * d = falloc(n_float);
  size_t nr = fread(d, sizeof(float), n_float, f);
  if(nr != n_float){
    printf("Expected number of floats: %zu\n", n_float);
    printf("Number of floats read: %zu\n", nr);
    err("unexpected float read count");
  }
  fclose(f);
  return d; /* return array of floats we read in */
}

#define STR_MAX 4096 /* string variable length */
void hwrite(char * bfn, size_t nrow, size_t ncol, size_t nband){
  size_t i;
  char hfn[STR_MAX];
  size_t L = strlen(bfn);

  strcpy(hfn, bfn); /* change ext from bin to hdr */
  hfn[L - 3] = 'h';
  hfn[L - 2] = 'd';
  hfn[L - 1] = 'r';

  FILE * f = open(hfn, WRITE);
  fprintf(f, "ENVI\n");
  fprintf(f, "samples = %zu\n", ncol);
  fprintf(f, "lines = %zu\n", nrow);
  fprintf(f,"bands = %zu\n", nband);
  fprintf(f, "header offset = 0\n");
  fprintf(f, "file type = ENVI Standard\n");
  fprintf(f, "data type = 4\n");
  fprintf(f, "interleave = bsq\n");
  fprintf(f, "byte order = 0\n");
  fprintf(f, "band names = {band 1");
  for0(i, nband - 1) fprintf(f, ",\nband %zu", i + 2);
  fprintf(f, "}\n");
  fclose(f);
}

float nan_to_num(float x){
  return isinf(x) ? 0 : (isnan(x)? FLT_MAX : x); /* replace NAN with zero, and infinity with large number */
}


int main(int argc, char ** argv){

  if(argc < 2) err("M4FC.exe [input T3 directory]");
  char * path = argv[1]; /* T3 matrix data path */
  int i, j, k, np, nrow, ncol, di, dj, ii, jj, x, ix, jx, nw;

  char fn[STR_MAX];
  strcpy(fn, path);
  fn[strlen(path)] = sep();
  strcpy(fn + strlen(path) + 1, "config.txt"); /* path to config.txt */
  read_config(fn, &nrow, &ncol); /* read image dimensions */
  np = nrow * ncol; /* number of px */

  arrays = (void *) malloc(sizeof(void *) * MAX_ARRAYS); /* array of pointers to free later */
  memset(arrays, '\0', sizeof(void *) * MAX_ARRAYS); /* always touch memory on win OS */
  n_arrays = 0; /* start from the beginning */

  T = (float **) alloc(sizeof(float *) * N_IN); /* input file buffers */
  for0(k, N_IN){
    strcpy(fn, path);
    fn[strlen(path)] = sep();
    strcpy(fn + strlen(path) + 1, T_fn[k]); /* [path][sep][filename] e.g. T3/T11.bin */
    T[k] = read(fn, np); /* read each input data band */
  }

  out_d = (float **) alloc(sizeof(float *) * N_OUT); /* output bands buffers */
  for0(i, N_OUT) out_d[i] = falloc(np); /* allocate output space */

  double h, g, span, val, theta_f, pd_f, ps_f, pv_f; /* intermediary variables */
  double t11, t12_r, t12_i, t13_r, t13_i, t22, t23_r, t23_i, t33;
  float * out_d_theta_f, * out_d_pd_f, * out_d_ps_f, * out_d_pv_f;
  double det_im, det_re, trace, trace3, m1_re, m1_im, r, theta, m1;

  float * t11_p = T[T11];
  float * t12_r_p = T[T12_re];
  float * t12_i_p = T[T12_im];
  float * t13_r_p = T[T13_re];
  float * t13_i_p = T[T13_im];
  float * t22_p = T[T22];
  float * t23_r_p = T[T23_re];
  float * t23_i_p = T[T23_im];
  float * t33_p = T[T33];

  out_d_theta_f = out_d[_theta_f];
  out_d_pd_f = out_d[_pd_f];
  out_d_ps_f = out_d[_ps_f];
  out_d_pv_f = out_d[_pv_f];

  for0(i, np){
    t11   = (double)t11_p[i];
    t12_r = (double)t12_r_p[i];
    t12_i = (double)t12_i_p[i];
    t13_r = (double)t13_r_p[i];
    t13_i = (double)t13_i_p[i];
    t22   = (double)t22_p[i];
    t23_r = (double)t23_r_p[i];
    t23_i = (double)t23_i_p[i];
    t33   = (double)t33_p[i];

    det_re = - t33 * t12_i *t12_i  - 2. * t12_i * t13_r * t23_i + t22 * t13_i *t13_i
             - t33 * t12_r * t12_r + 2. * t12_r * t13_r * t23_r - t22 * t13_r * t13_r
             - t11 * t23_i * t23_i - t11 * t23_r * t23_r + t11 * t22 * t33;

    det_im =  - 2. * (t12_i * t13_i * t23_i - t13_i * t12_r * t23_r  + t22 * t13_i * t13_r);

    trace = t11 + t22 + t33;
    trace3 = trace * trace * trace;
    m1_re = 1. - 27. * det_re / trace3;
    m1_im = 0. - 27. * det_im / trace3;

    r = sqrt(m1_re * m1_re + m1_im * m1_im);
    theta = (float)atan2((double)m1_im, (double)m1_re); /* convert to polar */
    m1 = sqrt(r) * cos(theta / 2); /* take square root and real part */

    g = t22 + t33;
    h = t11 - t22 - t33;
    span = t11 + t22 + t33;

    val = (m1 * span * h) / ( t11 * g + m1 * m1 * span * span);
    theta = atan(val);
    theta_f = (180. / M_PI) * theta;

    ps_f = nan_to_num(m1 * span * (1 + sin(2. * theta)) / 2.);
    pd_f = nan_to_num(m1 * span * (1 - sin(2. * theta)) / 2.);
    pv_f = nan_to_num(span * (1 - m1));

    out_d_theta_f[i] = (float)theta_f;
    out_d_pd_f[i] = (float)pd_f;
    out_d_ps_f[i] = (float)ps_f;
    out_d_pv_f[i] = (float)pv_f;
  }

  FILE * out_f[N_OUT];
  for0(i, N_OUT){
    strcpy(fn, path);
    fn[strlen(path)] = sep();
    strcpy(fn + strlen(path) + 1, out_fn[i]);
    out_f[i] = open(fn, WRITE);
    hwrite(fn, nrow, ncol, 1); /* write envi header */
    nw = fwrite(out_d[i], sizeof(float), np, out_f[i]);
    if(nw != np) err("failed to write expected number of floats");
    fclose(out_f[i]);
  }

  for0(k, n_arrays) free(arrays[n_arrays]); /* free anything we malloc'ed */
  free(arrays);
  return 0;
}
