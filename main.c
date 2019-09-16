#include "epimp.h"

#define ARGS printf("%s",instructions)
float UTHRESHOLD1, UTHRESHOLD2;
int UPERM;
char UTEST;

// usage: 
// epigpu -D [data stuff...]
// epigpu -A [analysis stuff...]

// data
// -s
// -q
// -p - change to r
// -g
// -c
// -i - change to m
// -r - change to n
// -u

// analysis

char *instructions =
	"\n\nepiMP functions in two modes, data management (D) and analysis (A)\n\n" \
	"DATA MANAGEMENT MODE:\n\n" \
	"<epiMP> -D[ arguments ] [ filenames ... ]\n\n\n" \
	"Arguments:\n\n" \
	"r\tRead PLINK data\n" \
	"c\tClean data, remove low call rate SNPs and individuals\n" \
	"m\tImpute missing genotypes based on allele frequency\n" \
	"n\tReplace phenotype with random normally distributed sample\n" \
	"q\tSimulate AxA interactions at specific SNP pairs\n" \
	"s\tSimulate entire dataset\n" \
	"e\tExtract binary data to PLINK format\n\n" \
	"Example:\n" \
	"<epiMP> -Drm [ .ped file ] [ .map file ] [ epigpu file ]\n" \
	"<epiMP> -De [ .ped file ] [ .map file ] [ epigpu file ]\n" \
	"<epiMP> -Dsq [ epigpu file ]\n\n\n\n" \
	"ANALYSIS MODE:\n\n" \
	"<epiMP> -A [ epigpu file] [ output file] [options ... ]\n\n" \
	"Options:\n\n" \
	"p [n]\tPermutation seed\n" \
	"t [f/i]\tFull test or full + interaction test (default = f)\n" \
	"F [n]\tF value threshold for full test, default [n] = 6.5\n" \
	"I [n]\tF value threshold for interaction test, default [n] = 10.5\n" \
	"1 [n]\tSpecify which chr x chr scan to perform (chromosome 1)\n" \
	"2 [n]\tChromosome 2\n\n" \
	"f [s]\tPath to fam file which has a different phenotype to the one in the epigpu file\n\n" \
	"Example:\n" \
	"<epiMP> -A [ epigpu file ] [ output file ]\n" \
	"<epiMP> -A [ epigpu file ] [ output file ] -1 5 -2 23\n" \
	"<epiMP> -A [ epigpu file ] [ output file ] -t i -F 7 -I 11\n\n\n\n" \

	"For additional help please see README.pdf for epiGPU software\n\n\n";


int main(int argc, char **argv)
{
	double tres;
	char *timetaken;
	time_t oval,nval;


	char *binfile, *outfile, *famfile;

	int nsnp, nid, npack, remain, nchr, tothits, UA, UB;
	ped *dat;
	map *genmap;
	int *genop;
	char *geno;
	chromosome *chrstat;
	int i,j,k;
	
	//plinkbin
	int flag;
	char *pedfile;
	char *mapfile;
	int *gpack;

	timetaken = malloc(sizeof(char) * 50);

	if(argc < 3)
	{
		ARGS;
		exit(1);
	}
	if(argv[1][0] != '-')
	{
		printf("Arguments\n");
		exit(1);
	}
	if(argv[1][1] == 'A') // analysis mode
	{
		flag = 2;
	} else if(argv[1][1] == 'D') // data mode
	{
		flag = 3;
	} else {
		ARGS;
		exit(1);
	}


	// DATA MODE
	if(flag == 3)
	{
		flag = datamode(argc, argv, &nid, &nsnp, &npack, &remain, &nchr, &dat, &genmap, &chrstat, &binfile, &pedfile, &mapfile, &genop, &gpack, &geno);
		return 0;
	}

	// ANALYSIS MODE
	// must have at least input and output files
	if(argc < 4 || (((int)argc & 1) != 0))
	{
		printf("Incorrect arguments\n");
		ARGS;
		exit(1);
	}

	// begin timing
	oval = time(NULL);
	binfile = argv[2];
	outfile = argv[3];

	silent = 0;
	UTHRESHOLD1 = 6.5;
	UTHRESHOLD2 = 10.5;
	UPERM = 0;
	UTEST = 'f';
	UA = 0;
	UB = 0;
	if(argc >= 4)
	{

		for(i = 4; i < argc; i+=2)
		{
			if(argv[i][0] != '-')
			{
				printf("Incorrect optional arguments %c\n",argv[i][0]);
				exit(1);
			}
			switch(argv[i][1])
			{
				case 'f' :
					famfile = argv[i+1];
					printf("Using %s as fam file\n", famfile);
					break;
				case 'p' :
					j = 0; k = 0;
					while(argv[i+1][j])
					{
						k += isdigit(argv[i+1][j++]) ? 0 : 1;
					}
					if(k != 0)
					{
						printf("Please enter a numeric value for permutations\nDefault = 0 (no permutation)\n\n"); exit(1);
					}
					UPERM = atoi(argv[i+1]);
					break;
				case 't' :
					if(argv[i+1][0] != 'f' && argv[i+1][0] != 'i')
					{
						printf("Test selection:\nf\tFull 8df test\ni\tFull 8df test + 4df interaction test\n\n");
						exit(1);
					}
					UTEST = argv[i+1][0];
					break;
				case 'F' :
					UTHRESHOLD1 = (float)atof(argv[i+1]);
					if(UTHRESHOLD1 == 0)
					{
						printf("Please enter a numeric F-test threshold value.\nRecommended values 5.5 - 10\nDefault = 7\n\n");
						exit(1);
					}
					break;
				case 'I' :
					UTHRESHOLD2 = (float)atof(argv[i+1]);
					if(UTHRESHOLD2 == 0)
					{
						printf("Please enter a numeric F-test threshold value.\nRecommended values 5.5 - 10\nDefault = 7\n\n");
						exit(1);
					}
					if(UTEST == 'f')
					{
						printf("\n\nWarning: Interaction test threshold has been set, but interaction test isn't being performed\n\n\n");
					}
					break;
				case '1' :
					j = 0; k = 0;
					while(argv[i+1][j])
					{
						k += isdigit(argv[i+1][j++]) ? 0 : 1;
					}
					if(k != 0)
					{
						printf("Please enter a chromosome number\nDefault = 0 (perform all chromosomes)\n\n"); exit(1);
					}
					UA = atoi(argv[i+1]);
					break;
				case '2' :
					j = 0; k = 0;
					while(argv[i+1][j])
					{
						k += isdigit(argv[i+1][j++]) ? 0 : 1;
					}
					if(k != 0)
					{
						printf("Please enter a chromosome number\nDefault = 0 (perform all chromosomes)\n\n"); exit(1);
					}
					UB = atoi(argv[i+1]);
					break;
				default :
					printf("Incorrect optional arguments %c\n",argv[i][0]);
					exit(1);
			}
		}
	}

	// Read in epiGPU format binary file
	readpackedbinary(&nid, &nsnp, &npack, &remain, &nchr, &genmap, &dat, &genop, &chrstat, binfile);
	// Unpack genotypes to char - faster on CPU than packing
	unpack(nid,nsnp,npack,remain,genop,&geno);
	// Begin scan

	ped *dat2;
	if(famfile)
	{
		FILE *FAM;

		dat2 = (ped *)malloc(sizeof(ped) * (nid));

		FAM = fopen(famfile,"r");
		printf("Reading in %s.",famfile);fflush(stdout);
		for(i = 0; i < nid; i++)
		{
			if(i % ((int)nid/10) == 1) {printf(".");fflush(stdout);}
			(void)fscanf(FAM,"%s",dat2[i].family);
			(void)fscanf(FAM,"%s",dat2[i].id);
			(void)fscanf(FAM,"%s",dat2[i].paternal);
			(void)fscanf(FAM,"%s",dat2[i].maternal);
			(void)fscanf(FAM,"%d",&j);
			(void)fscanf(FAM,"%f%*[ /t]",&dat2[i].phen);
			dat2[i].sex = j;
		}
		fclose(FAM);
		free(dat);
		dat = dat2;
		printf("done!\n\n");fflush(stdout);
	}

	if(UPERM > 0)
	{
		printf("Permutation %d\n",UPERM);
		srand((size_t)UPERM);
		permute(nid, dat);
		for(i = 0; i < 10; i++)
		{
			printf("%f\n",dat[i].phen);
		}
	}

	if(UA)
	{
		if(UB)
		{
			printf("Performing scan on %s vs %s...\n\n",chrstat[UA-1].chrname, chrstat[UB-1].chrname);

			if(UTEST == 'f')
			{
				printf("Performing 8df test\n");
				tothits = squareomp(geno, dat, nid, UA-1, UB-1, chrstat, genmap, outfile);
			} else{
				printf("Performing 4df test\n");
				tothits = squareompi(geno, dat, nid, UA-1, UB-1, chrstat, genmap, outfile);
			}

		} else {
			printf("Only one chromosome specified!\n");
			ARGS;
			exit(1);
		}
	} else {
		tothits = standard_scan(geno, nid, nsnp, nchr, genmap, chrstat, dat, outfile);
	}
	nval = time(NULL);
	timediff(&tres,nval,oval);
	readtime(timetaken,tres);
	printf("\nScan completed in %s\n%d epistatic interactions discovered\n",timetaken,tothits);

	return 0;
}



// Function to perform F test 
// Input is preformatted
// - missing values from genotype and phenotype are removed
// - joint genotype coding is stored as X (groups 0-8)
void ftest8df(char *X, float *Y, int n, float *F, int *df1, int *df2)
{
	int i, nfac, factor_count[9];
	float MSB, MSW, SSB = 0, SSW = 0, SSD, mY = 0;
	float mY_fac[9];

	for(i = 0; i < 9; i++)
	{
		factor_count[i] = 0;
		mY_fac[i] = 0;
	}
	
	// calculate class means
	
	for(i = 0; i < n; i++)
	{
		mY += Y[i];
		factor_count[(int)X[i]]++;
		mY_fac[(int)X[i]] += Y[i];
	}

	mY = mY / n;
	nfac = 0;
	// calculate SSB
	for(i = 0; i < 9; i++)
	{
		if (factor_count[i]>0){
			nfac++;
			mY_fac[i] /= factor_count[i];
			SSB += factor_count[i] * pow(mY_fac[i] - mY, 2);
		}
	}
	
	// calculate SSW
	for (i = 0; i < n; i++)
	{
		SSW += pow(Y[i] - mY, 2);
	}

	*df1 = nfac - 1;
	*df2 = n - nfac;
	SSD = SSW - SSB;

	MSB = SSB / *df1;
	MSW = SSD / *df2;
	*F = MSB / MSW;
}

void ftest4df(char *X, float *Y, int n, float *F, float *Fi, int *df1, int *df2)
{
	int i, nfac, factor_count[9];
	float MSB, MSW, SSB = 0, SSW = 0, SSI = 0, SSD, mY = 0;
	float mY_fac[9];
	float mean_row[3], mean_col[3];

	for(i = 0; i < 9; i++)
	{
		factor_count[i] = 0;
		mY_fac[i] = 0;
	}
	
	// calculate class means
	
	for(i = 0; i < n; i++)
	{
		mY += Y[i];
		factor_count[(int)X[i]]++;
		mY_fac[(int)X[i]] += Y[i];
	}

	mY = mY / n;


	for(i = 0; i < 3; i++)
	{
		mean_col[i] = (mY_fac[i]+mY_fac[i+3]+mY_fac[i+6]) / (factor_count[i]+factor_count[i+3]+factor_count[i+6]);
		mean_row[i] = (mY_fac[3*i]+mY_fac[3*i+1]+mY_fac[3*i+2]) / (factor_count[3*i]+factor_count[3*i+1]+factor_count[3*i+2]);
	}


	// calculate SSB
	nfac = 0;
	for(i = 0; i < 9; i++)
	{
		if(factor_count[i] > 0)
		{
			nfac++;
			mY_fac[i] /= factor_count[i];
			SSB += factor_count[i] * pow(mY_fac[i] - mY, 2);
			SSI += factor_count[i] * pow(mY_fac[i] - mean_row[(int)(i/3)] - mean_col[i%3] + mY, 2);
		}
	}

	// calculate SSW
	for (i = 0; i < n; i++)
	{
		SSW += pow(Y[i] - mY, 2);
	}

	*df1 = nfac - 1;
	*df2 = n - nfac;
	SSD = SSW - SSB;

	MSB = SSB / *df1;
	MSW = SSD / *df2;
	*F = MSB / MSW;
	*Fi = (SSI/4) / MSW;
}


void permute(int nid, ped *dat)
{
	int randnum, i;
	ped temp;
	char filename[50];
	FILE *out;
	for(i = nid; i > 1; i--)
	{
		randnum = rand() % i;
		temp = dat[randnum];
		dat[randnum] = dat[i - 1];
		dat[i - 1] = temp;
	}
	sprintf(filename,"perm%d.txt",UPERM);
	out = fopen(filename,"w");
	fprintf(out,"Family ID Phenotype\n");
	for(i = 0; i < nid; i++)
	{
		fprintf(out,"%s %s %f\n",dat[i].family, dat[i].id, dat[i].phen);
	}
	fclose(out);
}

// As in epiGPU, the scan is breaking down into chromosome x chromosome scans

int standard_scan(char *geno, int nid, int nsnp, int nchr, map *genmap, chromosome *chrstat, ped *dat, char *filename)
{
	char fileomp[25], line[256];
	FILE *dump, *out;
	int i, j, nthreads, tothits = 0;
	
	#pragma omp parallel private(dump,fileomp)
	{
    nthreads = omp_get_num_threads();
		sprintf(fileomp,"%s%d",filename,omp_get_thread_num());
		dump = fopen(fileomp,"w");
		fclose(dump);
	}

	for(i = 0; i < nchr; i++)
	{
		for(j = 0; j <= i; j++)
		{
			if(UTEST == 'f')
				tothits += squareomp(geno, dat, nid, i, j, chrstat, genmap, filename);
			else
				tothits += squareompi(geno, dat, nid, i, j, chrstat, genmap, filename);
		}
	}
	
	out = fopen(filename,"w");
	for(i = 0; i < nthreads; i++)
	{
	    sprintf(fileomp,"%s%d",filename,i);
	    dump = fopen(fileomp,"r");
	    while(fgets(line, sizeof(line), dump) != NULL)
    	{
			fputs(line, out);
		}
		fclose(dump);
		remove(fileomp);
	}
	
	return tothits;
}

int squareomp(char *geno, ped *dat, int nid, int chr1, int chr2, chromosome *chrstat, map *genmap, char *filename)
{
	int nthreads, tid;
	table *result;
	int i,ii,j,jj,k,df1,df2,hitcount,burncount,nclean;
	float F,*Y;
	char fileomp[25],*X;
	char *timetaken;
	FILE *dump;

	time_t oval,nval;
	double tres;
	oval = time(NULL);
	timetaken = (char *)malloc(sizeof(char) * 50);
	

	hitcount = 0;
	ii = chrstat[chr1].chrsize + chrstat[chr1].chrstart;

	#pragma omp parallel \
	default(shared) \
	private(i,j,k,tid,fileomp,X,Y,burncount,nclean,F,df1,df2,result,dump) \
	reduction(+: hitcount)
	{
		nthreads = omp_get_num_threads();
		tid = omp_get_thread_num();	
		burncount = 0; hitcount = 0;
	
		if(tid == 0)
		{
			if(chr1 == chr2)
			{
				printf("Scanning chromosome %s vs chromosome %s\n\tDimension %d x %d | triangle | %d threads\n\t", chrstat[chr1].chrname, chrstat[chr2].chrname, chrstat[chr1].chrsize, chrstat[chr2].chrsize,nthreads);
			} else {
				printf("Scanning chromosome %s vs chromosome %s\n\tDimension %d x %d | square | %d threads\n\t", chrstat[chr1].chrname, chrstat[chr2].chrname, chrstat[chr1].chrsize, chrstat[chr2].chrsize,nthreads);
			}
		}
    
		// each thread generates a different filename and opens it
		sprintf(fileomp,"%s%d",filename,tid);
		dump = fopen(fileomp,"a");
	
		// each thread allocates its own private variables
		Y = malloc(nid*sizeof(float));
		X = malloc(nid*sizeof(char));
		result = malloc(MAX*sizeof(table));
	
	
		// enter scan loop
	
		for(i = chrstat[chr1].chrstart; i < ii; i++)
		{			
			if(chr1 == chr2)
			{
				jj = i;
			} else {
				// specify chromosome2 scan range
				jj = chrstat[chr2].chrsize + chrstat[chr2].chrstart;
			}
	
			#pragma omp for	schedule(guided, 10)
			for(j = chrstat[chr2].chrstart; j < jj; j++)
			{
				nclean = 0;
				for(k = 0; k < nid; k++)
				{
					if((geno[(k+nid*i)] != 3) && (geno[(k+nid*j)] != 3))
					{
						X[nclean] = 3*geno[(i*nid + k)] + geno[(j*nid + k)];
						Y[nclean++] = dat[k].phen;
					}
				}
				ftest8df(X,Y,nclean,&F,&df1,&df2);
				if(F > UTHRESHOLD1)
				{
					result[burncount].snp1 = i;
					result[burncount].snp2 = j;
					result[burncount].df1 = df1;
					result[burncount].df2 = df2;
					result[burncount++].F = F;
					if(burncount == MAX)
					{
						for(k = 0; k < burncount; k++)
						{	
							fprintf(dump,"%d\t%d\t%d\t%d\t%f\n", result[k].snp1,result[k].snp2,result[k].df1,result[k].df2,result[k].F);
						}
						hitcount += burncount;
						burncount = 0;
					}
				}
			}
		}

		if(burncount > 0)
		{
			for(k = 0; k < burncount; k++)
			{
				fprintf(dump, "%d\t%d\t%d\t%d\t%f\n", result[k].snp1, result[k].snp2, result[k].df1, result[k].df2, result[k].F);
			}
			hitcount += burncount;
		}
		fclose(dump);
		free(X);
		free(Y);
		free(result);
	}
	
	nval = time(NULL);
	timediff(&tres, nval, oval);
	readtime(timetaken,tres);
	printf("\n\t%d interactions discovered in %s\n\n",hitcount,timetaken);
	free(timetaken);
	return hitcount;
}

int squareompi(char *geno, ped *dat, int nid, int chr1, int chr2, chromosome *chrstat, map *genmap, char *filename)
{
	int nthreads, tid;
	table *result;
	int i,ii,j,jj,k,df1,df2,hitcount,burncount,nclean;
	float F,Fi,*Y;
	char fileomp[25],*X;
	char *timetaken;
	FILE *dump;

	time_t oval,nval;
	double tres;
	oval = time(NULL);
	timetaken = (char *)malloc(sizeof(char) * 50);
	

	hitcount = 0;
	ii = chrstat[chr1].chrsize + chrstat[chr1].chrstart;

	#pragma omp parallel \
	default(shared) \
	private(i,j,k,tid,fileomp,X,Y,burncount,nclean,F,Fi,df1,df2,result,dump) \
	reduction(+: hitcount)
	{
		nthreads = omp_get_num_threads();
		tid = omp_get_thread_num();	
		burncount = 0; hitcount = 0;
	
		if(tid == 0)
		{
			if(chr1 == chr2)
			{
				printf("Scanning chromosome %s vs chromosome %s\n\tDimension %d x %d | triangle | %d threads\n\t", chrstat[chr1].chrname, chrstat[chr2].chrname, chrstat[chr1].chrsize, chrstat[chr2].chrsize,nthreads);
			} else {
				printf("Scanning chromosome %s vs chromosome %s\n\tDimension %d x %d | square | %d threads\n\t", chrstat[chr1].chrname, chrstat[chr2].chrname, chrstat[chr1].chrsize, chrstat[chr2].chrsize,nthreads);
			}
		}
    
		// each thread generates a different filename and opens it
		sprintf(fileomp,"%s%d",filename,tid);
		dump = fopen(fileomp,"a");
	
		// each thread allocates its own private variables
		Y = malloc(nid*sizeof(float));
		X = malloc(nid*sizeof(char));
		result = malloc(MAX*sizeof(table));
	
	
		// enter scan loop
	
		for(i = chrstat[chr1].chrstart; i < ii; i++)
		{			
			if(chr1 == chr2)
			{
				jj = i;
			} else {
				// specify chromosome2 scan range
				jj = chrstat[chr2].chrsize + chrstat[chr2].chrstart;
			}
	
			#pragma omp for	schedule(guided, 10)
			for(j = chrstat[chr2].chrstart; j < jj; j++)
			{
				nclean = 0;
				for(k = 0; k < nid; k++)
				{
					if((geno[(k+nid*i)] != 3) && (geno[(k+nid*j)] != 3))
					{
						X[nclean] = 3*geno[(i*nid + k)] + geno[(j*nid + k)];
						Y[nclean++] = dat[k].phen;
					}
				}
				ftest4df(X,Y,nclean,&F,&Fi,&df1,&df2);
				if(F > UTHRESHOLD1 || Fi > UTHRESHOLD2)
				{
					result[burncount].snp1 = i;
					result[burncount].snp2 = j;
					result[burncount].df1 = df1;
					result[burncount].df2 = df2;
					result[burncount].F = F;
					result[burncount++].F2 = Fi;
					if(burncount == MAX)
					{
						for(k = 0; k < burncount; k++)
						{	
							fprintf(dump,"%d\t%d\t%d\t%d\t%f\t%f\n", result[k].snp1,result[k].snp2,result[k].df1,result[k].df2,result[k].F, result[k].F2);
						}
						hitcount += burncount;
						burncount = 0;
					}
				}
			}
		}

		if(burncount > 0)
		{
			for(k = 0; k < burncount; k++)
			{
				fprintf(dump, "%d\t%d\t%d\t%d\t%f\t%f\n", result[k].snp1, result[k].snp2, result[k].df1, result[k].df2, result[k].F, result[k].F2);
			}
			hitcount += burncount;
		}
		fclose(dump);
		free(X);
		free(Y);
		free(result);
	}
	
	nval = time(NULL);
	timediff(&tres, nval, oval);
	readtime(timetaken,tres);
	printf("\n\t%d interactions discovered in %s\n\n",hitcount,timetaken);
	free(timetaken);
	return hitcount;
}
