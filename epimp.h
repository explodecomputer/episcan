#include <math.h>
#include "id4117.h"
#include "readplink.h"
#include <omp.h>

void ftest8df(char *X, float *Y, int n, float *F, int *df1, int *df2);
void ftest4df(char *X, float *Y, int n, float *F, float *Fi, int *df1, int *df2);
int standard_scan(char *geno, int nid, int nsnp, int nchr, map *genmap, chromosome *chrstat, ped *dat, char *filename);
int squareomp(char *geno, ped *dat, int nid, int chr1, int chr2, chromosome *chrstat, map *genmap, char *filename);
int squareompi(char *geno, ped *dat, int nid, int chr1, int chr2, chromosome *chrstat, map *genmap, char *filename);
void permute(int nid, ped *dat);
