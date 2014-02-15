/*****************************************************
* This code is distributed under a Non-commercial use 
* license.  For details see LICENSE.  Use of this
* code must be properly attributed to its author
* Gregory E. Sims provided that its use or derivative 
* use is non-commercial in nature.  Proper attribution        
* can be made by citing:
*
* Sims GE, et al (2009) Alignment-free genome 
* comparison with feature frequency profiles (FFP) and 
* optimal resolutions. Proc. Natl. Acad. Sci. USA.
* 106, 2677-82.
*
* Gregory E. Sims (C) 2010-2012
*
*****************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <math.h>
#include <errno.h>
#include <unistd.h>
#include <string.h>
#include "utils.h"
#include "vstring.h"
#include "sighandle.h"
#include "../config.h"


char PROG_NAME[FILENAME_MAX];/**< Program name */
#define TAXANAMELEN 50 /**< Maximum length string that can be used in a phylip format taxa name*/
#define DEFAULT_PRECISION 2 /**< Default precision for floating point output */
#define DEFAULT_NORM 2 /**< Default Norm for the Euclidean Distance Function */
#define STR_BUFF 255

int jsd(FILE * fp, double **D);
int jsdr(FILE * fp, double **D);
int euclidean_dist(FILE * fp, double **D);
int cosine_dist(FILE * fp, double **D);
int manhattan_dist(FILE * fp, double **D);
int chebyshev_dist(FILE * fp, double **D);
int jaccard_dist(FILE * fp, double **D);
int tanimoto_dist(FILE * fp, double **D);
int dice_dist(FILE * fp, double **D);
int antidice_dist(FILE * fp, double **D);
int hamming_dist(FILE * fp, double **D);
int evolution_dist(FILE * fp, double **D);
int yule_dist(FILE * fp, double **D);
int hamann_dist(FILE * fp, double **D);
int russel_dist(FILE * fp, double **D);
int sneath_dist(FILE * fp, double **D);
int matching_dist(FILE * fp, double **D);
int euclidean2_dist(FILE * fp, double **D);
int ochiai_dist(FILE * fp, double **D);
int canberra_dist(FILE * fp, double **D);
int anderberg_dist(FILE * fp, double **D);
int phi_dist(FILE * fp, double **D);
int gower_dist(FILE * fp, double **D);
int kulczynski_dist(FILE * fp, double **D);
float pearsons(double *x, double *y, int length);
int pearson_matrix(FILE * fp, double **D);
void printMatrix(double *D, int n);
void printInfile(double *D, int n);
void printLine(double *D, int n);

char usage_str[] = "Usage: %s [OPTION] vector ... \n\
Calculates a distance/divergence matrix from a columnar FFP.\n\n\
Given no options, the defeault behavior of the program is to\n\
generate a symmetric distance matrix.  The default distance used\n\
is the Jensen Shannon Divergence (JSD)\n\
\t-p FILE, --phylip=FILE\tPrint phylip output\n\
\t-d INT, --precision=INT\tSpecify decimal precision of matrix\n\
\t-r INT, --row=INT\tCalculate the INTth row of a JSD matrix\n\
\t-e, --euclid\t\tEuclidean Distance\n\
\t-E, --euclid2\t\tSquared Euclidean distance\n\
\t-n, --normval\t\tNorm val for -e, Default is 2\n\
\t-c, --cosine\t\tCosine distance [-s]\n\
\t-m, --manhattan\t\tManhattan distance\n\
\t-b, --canberra\t\tCanberra distance\n\
\t-M, --matching\t\tMatching distance [-s]\n\
\t-R, --pearson\t\tPearson Correlation distance [-s]\n\
\t-C, --chebyshev\t\tChebyshev distnace\n\
\t-j, --jaccard\t\tJaccard distance [-s]\n\
\t-t, --tanimoto\t\tRogers-Tanimoto distance [-s]\n\
\t-D, --dice\t\tDice distance [-s]\n\
\t-N, --antidice\t\tAnti-dice distance [-s]\n\
\t-S, --sneath\t\tSneath-Sokal distance [-s]\n\
\t-H, --hamming\t\tHamming distance\n\
\t-L, --evol\t\tEvolutionary distance\n\
\t-a, --hamman\t\tHamman distance [-s]\n\
\t-P, --phi\t\tPearson Phi distance [-s]\n\
\t-B, --anderberg\t\tAnderberg distance [-s]\n\
\t-g, --gower\t\tGower distance[-s]\n\
\t-u, --russel\t\tRussel-rao distance [-s]\n\
\t-y, --yule\t\tYule distance [-s]\n\
\t-o, --ochiai\t\tOchiai distance [-s]\n\
\t-k, --kulczynski\tKulczynski distance [-s]\n\
\t-s, --similarity\tForce similarity matrix (-koyrgapaSNDTjM)\n\
\t-q, --quiet\tSuppress warning messages\n\
\t-h, --help\t\tThis message\n\
\t-v, --version\t\tPrint version\n\n\
Copyright (c) %s\n\
%s\n\
Contact %s\n";


enum dist_modes { jensen_shannon, euclidean, cosine, manhattan, pearson_r,
    chebyshev, jaccard, tanimoto, dice, hamming, yule,
    russel, matching, hamann, antidice, sneath, ochiai,
    euclidean2, canberra, anderberg, phi, gower, kulczynski,evolution
};
enum matrix_modes { similarity, distance };

char matrix_mode = distance;
int precision = DEFAULT_PRECISION;
float euclidean_norm = DEFAULT_NORM;
				 /**< Precision for floating point output */
char **taxaNames;  /**< Array to store Taxa names for phylip infile output */


char rflag = 0;
	      /**< -l Option for calculating specific line in JSD */
int rflagN = 0;
char qFlag = 0;

int main(int argc, char **argv)
{
    FILE *fp;
    FILE *pp;
    int opt;
    unsigned rows = 0;
    char dist_mode = jensen_shannon;
    char pflag = 0;
    char dflag = 0;
    int dvalue = 0;
    char *pvalue = NULL;
    int i;
    double *D = NULL;
    int option_index = 0;
    char buffer[STR_BUFF];

    static struct option long_options[] = {
	{"phylip", required_argument, 0, 'p'},
	{"precision", required_argument, 0, 'd'},
	{"euclid", no_argument, 0, 'e'},
	{"euclid2", no_argument, 0, 'E'},
	{"canberra", no_argument, 0, 'b'},
	{"normval", required_argument, 0, 'n'},
	{"cosine", no_argument, 0, 'c'},
	{"manhattan", no_argument, 0, 'm'},
	{"pearson", no_argument, 0, 'R'},
	{"chebyshev", no_argument, 0, 'C'},
	{"jaccard", no_argument, 0, 'j'},
	{"tanimoto", no_argument, 0, 't'},
	{"dice", no_argument, 0, 'D'},
	{"hamming", no_argument, 0, 'H'},
	{"evol", no_argument, 0, 'L'},
	{"yule", no_argument, 0, 'y'},
	{"russel", no_argument, 0, 'u'},
	{"matching", no_argument, 0, 'M'},
	{"hamman", no_argument, 0, 'a'},
	{"antidice", no_argument, 0, 'N'},
	{"sneath", no_argument, 0, 'S'},
	{"ochiai", no_argument, 0, 'o'},
	{"version", no_argument, 0, 'v'},
	{"row", required_argument, 0, 'r'},
	{"anderberg", no_argument, 0, 'B'},
	{"phi", no_argument, 0, 'P'},
	{"gower", no_argument, 0, 'g'},
	{"kulczynski", no_argument, 0, 'k'},
	{"similarity", no_argument, 0, 's'},
	{"help", no_argument, 0, 'h'},
	{"quiet", no_argument, 0, 'q'},
	{0, 0, 0, 0}
    };

    initSignalHandlers();

    strcpy(PROG_NAME,basename( argv[0] ));

    while ((opt = getopt_long(argc, argv, "abp:d:ghkevr:cmBERCDHMNSPsn:ojtyuqL",
			      long_options, &option_index)) != -1)
	switch (opt) {
	case 'p':
	    pflag = 1;
	    pvalue = optarg;
	    break;
	case 'd':
	    dflag = 1;
	    dvalue = atoi(optarg);
	    break;
	case 'q':
	    qFlag = 1;
	    break;
	case 'c':
	    dist_mode = cosine;
	    break;
	case 'e':
	    dist_mode = euclidean;
	    break;
	case 'E':
	    dist_mode = euclidean2;
	    break;
	case 'b':
	    dist_mode = canberra;
	    break;
	case 'n':
	    euclidean_norm = atof(optarg);
	    if (euclidean_norm <= 1)
		fatal_msg("Norm value must be greater than 1\n");
	    break;
	case 'm':
	    dist_mode = manhattan;
	    break;
	case 'R':
	    dist_mode = pearson_r;
	    break;
	case 'C':
	    dist_mode = chebyshev;
	    break;
	case 'j':
	    dist_mode = jaccard;
	    break;
	case 't':
	    dist_mode = tanimoto;
	    break;
	case 'D':
	    dist_mode = dice;
	    break;
	case 'a':
	    dist_mode = hamann;
	    break;
	case 'H':
	    dist_mode = hamming;
	    break;
	case 'L':
	    dist_mode = evolution;
	    break;
	case 'S':
	    dist_mode = sneath;
	    break;
	case 'N':
	    dist_mode = antidice;
	    break;
	case 'y':
	    dist_mode = yule;
	    break;
	case 'u':
	    dist_mode = russel;
	    break;
	case 'M':
	    dist_mode = matching;
	    break;
	case 'o':
	    dist_mode = ochiai;
	    break;
	case 'B':
	    dist_mode = anderberg;
	    break;
	case 'P':
	    dist_mode = phi;
	    break;
	case 'g':
	    dist_mode = gower;
	    break;
	case 'k':
	    dist_mode = kulczynski;
	    break;
	case 'h':
	    printUsageStr();
	    exit(EXIT_SUCCESS);
	    break;
	case 'v':
	    printVersion();
	    exit(EXIT_SUCCESS);
	    break;
	case 'r':
	    rflag = 1;
	    rflagN = atoi(optarg) - 1;
	    break;
	case 's':
	    matrix_mode = similarity;
	    break;
	default:
	    printErrorUsageStr();
	    break;
	}


    if (dflag)
	precision = dvalue;

// process file arguments 
    argv += optind;

    do {
	fp = stdin;
	if (*argv) {
	    if (!strcmp(*argv, "-"))
		fp = stdin;
	    else if ((fp = fopen(*argv, "r")) == NULL)
		fatal_msg("%s: %s.\n", *argv,strerror(errno));

	    if ( isDirectory(*argv) ) 
	    	fatal_msg("%s: %s\n",*argv, strerror(EISDIR));

	    argv++;
	} else if (isatty(STDIN_FILENO))
	    printErrorUsageStr();

	// must confirm seekabilitiy

	if (!isRegularFile(fp)) {
	    fp = convertPipeToFile(fp);
	}


	switch (dist_mode) {
	case euclidean:
	    rows = euclidean_dist(fp, &D);
	    break;
	case euclidean2:
	    rows = euclidean2_dist(fp, &D);
	    break;
	case cosine:
	    rows = cosine_dist(fp, &D);
	    break;
	case manhattan:
	    rows = manhattan_dist(fp, &D);
	    break;
	case pearson_r:
	    rows = pearson_matrix(fp, &D);
	    break;
	case chebyshev:
	    rows = chebyshev_dist(fp, &D);
	    break;
	case jaccard:
	    rows = jaccard_dist(fp, &D);
	    break;
	case tanimoto:
	    rows = tanimoto_dist(fp, &D);
	    break;
	case dice:
	    rows = dice_dist(fp, &D);
	    break;
	case hamming:
	    rows = hamming_dist(fp, &D);
	    break;
	case evolution:
	    rows = evolution_dist(fp, &D);
	    break;
	case yule:
	    rows = yule_dist(fp, &D);
	    break;
	case russel:
	    rows = russel_dist(fp, &D);
	    break;
	case hamann:
	    rows = hamann_dist(fp, &D);
	    break;
	case antidice:
	    rows = antidice_dist(fp, &D);
	    break;
	case sneath:
	    rows = sneath_dist(fp, &D);
	    break;
	case ochiai:
	    rows = ochiai_dist(fp, &D);
	    break;
	case canberra:
	    rows = canberra_dist(fp, &D);
	    break;
	case anderberg:
	    rows = anderberg_dist(fp, &D);
	    break;
	case phi:
	    rows = phi_dist(fp, &D);
	    break;
	case gower:
	    rows = gower_dist(fp, &D);
	    break;
	case kulczynski:
	    rows = kulczynski_dist(fp, &D);
	    break;
	case matching:
	    rows = matching_dist(fp, &D);
	    break;
	case jensen_shannon:
	    if (rflag)
		rows = jsdr(fp, &D);
	    else
		rows = jsd(fp, &D);
	    break;
	}


	if (pflag)		// if phylip format requested
	{
	    if ((pp = fopen(pvalue, "r")) == NULL) {
		fprintf(stderr, "Error opening file %s", pvalue);
		exit(1);
	    }

	    taxaNames = (char **) malloc(sizeof(char *) * rows);
	    for (i = 0; i < rows; i++) {
		taxaNames[i] =
		    (char *) malloc(sizeof(char) * (TAXANAMELEN + 1));
		if (!fscanf(pp, "%s", buffer)) 
			fatal_msg("%: Read zero items.",pvalue);
		if (strlen(buffer) > TAXANAMELEN) {
		    if (!qFlag) 	
		    warn_msg("Taxaname: %s greater than %d. Truncating.\n",
			    buffer, TAXANAMELEN);
		    buffer[10] = '\0';
		}
		strcpy(taxaNames[i], buffer);
	    }
	}

	if (rflag)
	    printLine(D, rows);
	else if (pflag)
	    printInfile(D, rows);
	else
	    printMatrix(D, rows);

	free(D);


	if (fp != stdin)
	    fclose(fp);

    } while (*argv);

    return EXIT_SUCCESS;
}


/**
 * Calculates a Jensen Shannon Divergence of an FFP matrix
 * using a single line.
 *
 * The FFP is read from fp, and the JSD calculation is
 * perfomed for all pairs of FFP.  The form of the FFP 
 * must be columnar row normalized data.
 *
 * @param fp A file pointer to a row normalized columnar FFP
 * @param D Points to a dynamically allocated array of double precision floats.
 * @return Returns the number of rows read.
 */




int jsdr(FILE * fp, double **D)
{
    double *a;
    double val;
    double m;
    char ch[2];
    int tst = 0;
    long pos;
    int i;
    double ha;
    double hb;
    unsigned cols = 0;
    int row;
    unsigned colsize = 10000;
    int rowsize = 10;

    a = (double *) malloc(sizeof(double) * colsize);
    *D = (double *) malloc(sizeof(double) * rowsize);



    // Read length of first line

    while (tst != 2) {
	if (cols >= colsize)	// if vector is too large reallocate
	{
	    colsize += 10000;
	    if ((a = (double *) realloc(a, sizeof(double) * colsize)) == NULL) {
		fprintf(stderr, "Out of memory!\n");
		exit(ENOMEM);
	    }
	}
	tst = fscanf(fp, "%lf%[\n\r]", &a[cols++], ch);
    }

    pos = ftell(fp);

    // If line = 0 then good to go.  else jump forward to the right line and re-read
    if (rflagN > 0) {
	fseek(fp, pos * rflagN, SEEK_SET);
	//Now read that line;
	//check for EOF in case an invalid line was specified   
	for (i = 0; i < cols; i++)
	    tst = fscanf(fp, "%lf\n", &a[i]);
	// Possibly check for fscanf errors.    
    }

    clearerr(fp);
    fseek(fp, 0, SEEK_SET);
    row = 0;

    while (!feof(fp)) {

/*
	if (row == rflagN) {
		// last line
		(*D)[row++] = 0;
		fseek(fp,pos,SEEK_CUR);
		continue;
	}*/

	while (row + 1 > rowsize) {
	    rowsize *= 2;
	    if ((*D = (double *) realloc(*D, sizeof(double) * rowsize)) == NULL) {
		fprintf(stderr, "Out of memory!\n");
		exit(ENOMEM);
	    }
	}


	ha = hb = 0;
	for (i = 0; i < cols; i++) {
	    if (!fscanf(fp, "%lf\n", &val))
		fatal_msg("Read zero items.");
	    m = (val + a[i]) / 2.0;
	    if (!a[i] && !val) {
		continue;
	    } else if (!a[i]) {
		hb -= val;
		continue;
	    } else if (!val) {
		hb -= a[i];
		continue;
	    } else {
		ha += -a[i] * log2(m / a[i]);
		hb += -val * log2(m / val);
	    }
	}
	(*D)[row++] = fabs(0.5 * ha + 0.5 * hb);
    }
    free(a);
    return (row);
}





/**
 * Print out a distance matrix
 *
 * Prints out the distance matrix
 * stored in D with the precision specified
 *
 * @param n Dimensions of the distance matrix
 * @param D Points to a dynamically allocated array of double precision floats.
 * @return None
 * @see precision
 */


void printMatrix(double *D, int n)
{
    int i, j, k = 0, l = 0;

    for (i = 0; i < n; i++) {
	k += i;
	for (j = 0, l = 0; j < n; j++, l += j)
	    printf("%.*e ", precision,
		   (i > j ? D[j * n + i - l] : D[i * n + j - k]));
	printf("\n");

    }
}



/**
 * Print out a distance matrix
 *
 * Prints out the distance matrix
 * stored in D with the precision specified
 *
 * @param n Dimensions of the distance matrix
 * @param D Points to a dynamically allocated array of double precision floats.
 * @return None
 * @see precision
 */


void printLine(double *D, int n)
{
    int i;
    for (i = 0; i < n; i++)
	printf("%.*e ", precision, D[i]);
    printf("\n");

}



/**
 * Print a phylip format infile
 *
 * Prints out the distance matrix
 * stored in D, along with the taxa 
 * names stored in the global variable
 * taxaNames.  The precision used will
 * be as specified in the precision variable.
 *
 * @param n Dimensions of the distance matrix
 * @param D Points to a dynamically allocated array of double precision floats.
 * @return None
 * @see taxaNames
 * @see precision
 */




void printInfile(double *D, int n)
{
    int i, j, k, l;
    printf("%d\n", n);

    for (i = 0, k = 0; i < n; i++, k += i) {
	printf("%-*s", TAXANAMELEN, taxaNames[i]);
	for (j = 0, l = 0; j < n; j++, l += j)
	    printf("%.*e ", precision,
		   (i > j ? D[j * n + i - l] : D[i * n + j - k]));
	printf("\n");

    }



}



/**
 * Calculates pearson correlation coeficient
 * 
 * This function calculates a pairwise correlation
 * This is a generic function, and hasn't been optimized
 * to calculate correlations with FFPs quickly. 
 *
 * @param x A double array
 * @param y A double array
 * @param length The length of the array
 * @retrun Returns the correlation coefficient
 * @todo Incorporate code into pearson_dist so that the stdev
 *       of x isn't continually recalculated
 */

float pearsons(double *x, double *y, int length)
{
    int i;
    double N = (double) length;
    double sum_sq_x = 0.0;
    double sum_sq_y = 0.0;
    double sum_coproduct = 0.0;
    double mean_x = x[0];
    double mean_y = y[0];
    double sweep, delta_x, delta_y, index;
    for (i = 1; i < N; i++) {
	index = (double) i;
	sweep = (index - 1.0) / index;
	delta_x = x[i] - mean_x;
	delta_y = y[i] - mean_y;
	sum_sq_x += delta_x * delta_x * sweep;
	sum_sq_y += delta_y * delta_y * sweep;
	sum_coproduct += delta_x * delta_y * sweep;
	mean_x += delta_x / index;
	mean_y += delta_y / index;
    }

    if (sum_sq_x <= 0 || isnan(sum_sq_x) || isinf(sum_sq_x))
	return 0.0;

    double pop_sd_x = sqrt(sum_sq_x / N);
    double pop_sd_y = sqrt(sum_sq_y / N);
    double cov_x_y = sum_coproduct / N;

    return (float) cov_x_y / (pop_sd_x * pop_sd_y);
}

/**
 * This definition below is a macro that eliminates some repetitive code
 * below.  I chose the option of lazy coding, rather than find the best
 * way to implement the distance metrics with a smaller binary file.
 * The template is used to generate actual functions which follow the
 * template format
 */



#define DISTANCE_TEMPLATE(FUNC_NAME,DEFS,ALLOC, DIAGONAL,CALC_1,CALC_2,REALLOC) \
					\
int FUNC_NAME (FILE * fp, double **D)	\
{					\
    double *a;				\
    DEFS;				\
    double dist;			\
    char ch[2];				\
    int tst;				\
    long pos;				\
    int i, k;				\
    int kMax;				\
    int rMax;				\
    unsigned cols;			\
    int row;				\
    unsigned colsize = 10000;		\
    int rowsize = 10;			\
		                        \
    ALLOC;				\
    a = (double *) malloc(sizeof(double) * colsize);	\
    *D = (double *) malloc(sizeof(double) * rowsize);	\
							\
    k = 0;						\
    kMax = 1;						\
    rMax = 0;						\
    while (k < kMax) {					\
	tst = 0;					\
	row = 0;					\
	cols = 0;					\
	while (tst != 2) {				\
	    if (cols >= colsize)			\
	    {						\
		colsize += 10000;			\
		if ((a =				\
		     (double *) realloc(a, sizeof(double) * colsize)) == NULL) { \
		    fprintf(stderr, "Out of memory!\n");			 \
		    exit(ENOMEM);						 \
		}								 \
		REALLOC \
			\
	    }									 \
										 \
	    tst = fscanf(fp, "%lf%[\n\r]", &a[cols++], ch);		 \
	}									 \
	row++;									 \
	pos = ftell(fp);							 \
	(*D)[k++] = DIAGONAL;							\
	while (!feof(fp)) {							\
	    while (k + cols > rowsize ) 					\
	    {									\
		rowsize *= 2;							\
		if ((*D = 							\
		   (double *) realloc(*D, sizeof(double) * rowsize)) ==NULL) {	\
			fprintf(stderr,"Out of memory!\n");			\
			exit(ENOMEM);						\
		}								\
	    }									\
	    dist=0;								\
	    for (i = 0; i < cols; i++) {					\
		CALC_1;								\
	    }									\
	    row++;								\
	    CALC_2;								\
	}									\
										\
	if (row > rMax) {							\
	    rMax = row;								\
	    kMax = (row * row - row + 2 * row) / 2;				\
	}									\
										\
	clearerr(fp);								\
	fseek(fp, pos, SEEK_SET);						\
    }										\
										\
    free(a);									\
    return (rMax);								\
}





/**
 * Calculates a Jensen Shannon Divergence of an FFP matrix
 *
 * The FFP is read from fp, and the JSD calculation is
 * perfomed for all pairs of FFP.  The form of the FFP 
 * must be columnar row normalized data.
 *
 * @param fp A file pointer to a row normalized columnar FFP
 * @param D Points to a dynamically allocated array of double precision floats.
 * @return Returns the number of rows read.
 */

DISTANCE_TEMPLATE(jsd, double val;
		  double m;
		  double ha = 0;
		  double hb = 0,,
		  0,
		  if (!fscanf(fp, "%lf\n", &val)) { fatal_msg("Read zero items"); }
		  m = (val + a[i]) / 2.0; if (!a[i] && !val) {
		  continue;}

		  else
		  if (!a[i]) {
		  hb -= val; continue;}

		  else
		  if (!val) {
		  hb -= a[i]; continue;}

		  else {
		  ha += -a[i] * log2(m / a[i]); hb += -val * log2(m / val);}

		  , (*D)[k++] = fabs(0.5 * ha + 0.5 * hb); ha = 0; hb = 0;,)


/**
 * Calculates a Jaccard Distance matrix of an FFP matrix
 *
 * The FFP is read from fp, and the Jaccard distance metric is
 * calculated for all pairs of FFP.  The form of the FFP 
 * must be a columnar matrix, but not necessarily row normalized.
 * Each row is considered as a bit string.
 * This distance represents the number of features shared in 
 * common divided by the number of columns.
 * The -s option affects the output, producing a Jaccard similarity
 * marix rather than a distance matrix.
 * @todo If both measures are all zeros, this is undefined.  All zero vectors
 *       will have a similarity of 1, and zero to non-zero a similarity of 0.
 *
 * @param fp A file pointer to a columnar FFP
 * @param D Points to a dynamically allocated array of double precision floats.
 * @return Returns the number of rows read.
 */


DISTANCE_TEMPLATE(jaccard_dist,
	      double val,,
	      (matrix_mode == similarity),
	      if (!fscanf(fp, "%lf\n", &val)) { fatal_msg("Read zero items"); }
	      dist += (val && a[i]) ? 1 : 0, if (matrix_mode == similarity)
	      (*D)[k++] = dist / cols;
	      else
	      (*D)[k++] = 1 - dist / cols;,)



/**
 * Calculates a Tanimoto Distance matrix of an FFP matrix
 *
 * The FFP is read from fp, and the Jaccard distance metric is
 * calculated for all pairs of FFP.  The form of the FFP 
 * must be a columnar matrix, but not necessarily row normalized.
 * Each row is considered as a bit string.
 * This distance represents the number of features shared in 
 * common divided by the number of columns.
 * The -s option affects the output, producing a Jaccard similarity
 * marix rather than a distance matrix.
 *
 * @param fp A file pointer to a columnar FFP
 * @param D Points to a dynamically allocated array of double precision floats.
 * @return Returns the number of rows read.
 */

DISTANCE_TEMPLATE(tanimoto_dist, double val;
	      double amag = 0;
	      double bmag = 0,,
	      (matrix_mode == similarity),
	      if (!fscanf(fp, "%lf\n", &val)) { fatal_msg("Read zero items"); }
	      bmag += (val != 0);
	      amag += (a[i] != 0);
	      dist += (val && a[i]), if (matrix_mode == similarity)
	      (*D)[k++] = dist / (bmag + amag - dist);
	      else
	      (*D)[k++] = 1 - dist / (bmag + amag - dist); amag = bmag = 0;,)


/**
 * Calculates a Chebyshev distance matrix of an FFP matrix
 *
 * The FFP is read from fp, and the Chebyshev distance metric is
 * calculated for all pairs of FFP.  The form of the FFP 
 * must be a columnar matrix, but not necessarily row normalized.
 * This distance represented the the greatest difference among
 * vectors along any coordinate dimension.
 *
 * @param fp A file pointer to a columnar FFP
 * @param D Points to a dynamically allocated array of double precision floats.
 * @return Returns the number of rows read.
 */


DISTANCE_TEMPLATE(chebyshev_dist,
	      double val,,
	      0,
	      if (!fscanf(fp, "%lf\n", &val)) { fatal_msg("Read zero items"); }
	      if (fabs(val - a[i]) > dist)
	      dist = fabs(val - a[i]), (*D)[k++] = dist;,)



/**
 * Calculates a pearson correlation matrix of an FFP matrix
 *
 * The FFP is read from fp, and the correlation calculation is
 * perfomed for all pairs of FFP.  The form of the FFP 
 * must be a columnar matrix, but not necessarily row normalized.
 * In other words, raw frequencies can be used, for example
 * output directly from ffpcol.  If the -s option is specified
 * Matrix D will a similarity matrix otherwise it will be calculated
 * as 1-R_squared.
 *
 * @param fp A file pointer to a columnar FFP
 * @param D Points to a dynamically allocated array of double precision floats.
 * @return Returns the number of rows read.
 */

    DISTANCE_TEMPLATE(pearson_matrix,
		      double *b,
		      b = (double *) malloc(sizeof(double) * colsize),
		      (matrix_mode == similarity),
	      	      if (!fscanf(fp, "%lf\n", &b[i])) { fatal_msg("Read zero items"); },
		      if (matrix_mode == similarity)
		      (*D)[k++] = pearsons(a, b, cols);
		      else
		      (*D)[k++] = 1 - pow(pearsons(a, b, cols), 2);,
		      if ((b =  (double *) realloc(b, sizeof(double) * colsize)) == NULL) { 
		    		fprintf(stderr, "Out of memory!\n");			 
		    		exit(ENOMEM);						 
		      })								 


/**
 * Calculates a Manhattan Distance of an FFP matrix
 *
 * The FFP is read from fp, and the distance calculation is
 * perfomed for all pairs of FFP.  The form of the FFP 
 * must be a columnar matrix, but not necessarily row normalized.
 * In other words, raw frequencies can be used, for example
 * output directly from ffpcol.
 *
 * Distance Calculation is: sum(abs(ai - bi))
 *
 * @param fp A file pointer to a columnar FFP
 * @param D Points to a dynamically allocated array of double precision floats.
 * @return Returns the number of rows read.
 */

	DISTANCE_TEMPLATE(manhattan_dist,
			  double val,,
			  0,
	      	      	  if (!fscanf(fp, "%lf\n", &val)) { fatal_msg("Read zero items"); }
			  dist += fabs(val - a[i]), (*D)[k++] = dist;,)



/**
 * Calculates a Cosine Distance of an FFP matrix
 *
 * The FFP is read from fp, and the distance calculation is
 * perfomed for all pairs of FFP.  The form of the FFP 
 * must be columnar data. Row normalization is not required.
 *
 * D(A,B)=1-A.B/|A|/|B|
 *
 * @param fp A file pointer to a columnar FFP
 * @param D Points to a dynamically allocated array of double precision floats.
 * @return Returns the number of rows read.
 */
    DISTANCE_TEMPLATE(cosine_dist, double val;
		      double amag = 0;
		      double bmag = 0,,
		      (matrix_mode == similarity),
	      	      if (!fscanf(fp, "%lf\n", &val)) { fatal_msg("Read zero items"); }
		      bmag += val * val;
		      amag += a[i] * a[i];
		      dist += val * a[i],
		      bmag = sqrt(bmag);
		      amag = sqrt(amag); if (matrix_mode == similarity)
		      (*D)[k++] = dist / bmag / amag;
		      else
		      (*D)[k++] = 1 - dist / bmag / amag; amag = bmag = 0;,)


/**
 * Calculates a Euclidean Distance of an FFP matrix
 *
 * The FFP is read from fp, and the distance calculation is
 * perfomed for all pairs of FFP.  The form of the FFP 
 * must be a columnar matrix, but not necessarily row normalized.
 * In other words, raw frequencies can be used, for example
 * output directly from ffpcol.
 *
 * @param fp A file pointer to a columnar FFP
 * @param D Points to a dynamically allocated array of double precision floats.
 * @return Returns the number of rows read.
 */


    DISTANCE_TEMPLATE(euclidean_dist,
		      double val,,
		      0,
	      	      if (!fscanf(fp, "%lf\n", &val)) { fatal_msg("Read zero items"); }
		      dist += pow((val - a[i]), euclidean_norm),
		      (*D)[k++] = pow(dist, 1 / euclidean_norm);,)


/**
 * Calculates a Euclidean-squared Distance of an FFP matrix
 *
 * The FFP is read from fp, and the distance calculation is
 * perfomed for all pairs of FFP.  The form of the FFP 
 * must be a columnar matrix, but not necessarily row normalized.
 * In other words, raw frequencies can be used, for example
 * output directly from ffpcol.
 *
 * @param fp A file pointer to a columnar FFP
 * @param D Points to a dynamically allocated array of double precision floats.
 * @return Returns the number of rows read.
 */
DISTANCE_TEMPLATE(euclidean2_dist,
		      double val,,
		      0,
	      	      if (!fscanf(fp, "%lf\n", &val)) { fatal_msg("Read zero items"); }
		  dist += (val - a[i]) * (val - a[i]);
		  , (*D)[k++] = dist;,)

/**
 * Calculates a canberra Distance of an FFP matrix
 *
 * The FFP is read from fp, and the distance calculation is
 * perfomed for all pairs of FFP.  The form of the FFP 
 * must be a columnar matrix, but not necessarily row normalized.
 * In other words, raw frequencies can be used, for example
 * output directly from ffpcol.
 *
 * @todo if either amag or bmag is zero then undefined, 0 to 0 is 0;  0 to 1 is  1
 * @param fp A file pointer to a columnar FFP
 * @param D Points to a dynamically allocated array of double precision floats.
 * @return Returns the number of rows read.
 */
    DISTANCE_TEMPLATE(canberra_dist, double val;
		      double amag = 0;
		      double bmag = 0;
		      ,,
		      0,
	      	      if (!fscanf(fp, "%lf\n", &val)) { fatal_msg("Read zero items"); }
		      amag += val * val;
		      bmag += a[i] * a[i]; dist += fabs(val - a[i]);
		      ,
		      amag = sqrt(amag);
		      bmag = sqrt(bmag);
		      (*D)[k++] = dist / amag / bmag; amag = 0;
		      bmag = 0;,)


/**
 * Calculates a Dice Distance of an FFP matrix
 *
 * The FFP is read from fp, and the distance calculation is
 * perfomed for all pairs of FFP.  The form of the FFP 
 * must be a columnar matrix, but not necessarily row normalized.
 * In other words, raw frequencies can be used, for example
 * output directly from ffpcol.
 *
 * @param fp A file pointer to a columnar FFP
 * @param D Points to a dynamically allocated array of double precision floats.
 * @return Returns the number of rows read.
 */
    DISTANCE_TEMPLATE(dice_dist, double val;
		      double amag = 0;
		      double bmag = 0,,
		      (matrix_mode == similarity),
	      	      if (!fscanf(fp, "%lf\n", &val)) { fatal_msg("Read zero items"); }
		      bmag += (val != 0);
		      amag += (a[i] != 0);
		      dist += (val && a[i]), if (matrix_mode == similarity)
		      (*D)[k++] = 2 * dist / (bmag + amag);
		      else
		      (*D)[k++] = 1 - 2 * dist / (bmag + amag);
		      amag = bmag = 0;,)

/**
 * Calculates a bitwise hamming distance of an FFP matrix
 *
 *
 * The FFP is read from fp, and the distance calculation is
 * perfomed for all pairs of FFP.  The form of the FFP 
 * must be a columnar matrix, but not necessarily row normalized.
 * In other words, raw frequencies can be used, for example
 * output directly from ffpcol.
 *
 * @param fp A file pointer to a columnar FFP
 * @param D Points to a dynamically allocated array of double precision floats.
 * @return Returns the number of rows read.
 */


    DISTANCE_TEMPLATE(hamming_dist,
		      double val,,
		      0,
	      	      if (!fscanf(fp, "%lf\n", &val)) { fatal_msg("Read zero items"); }
		      dist += !(val && a[i]), (*D)[k++] = dist;,)


/**
 * Calculates the Evolutionary distance used in Ecoli paper.
 *
 *
 * The FFP is read from fp, and the distance calculation is
 * perfomed for all pairs of FFP.  The form of the FFP 
 * must be a columnar matrix, but not necessarily row normalized.
 * In other words, raw frequencies can be used, for example
 * output directly from ffpcol.
 *
 * @param fp A file pointer to a columnar FFP
 * @param D Points to a dynamically allocated array of double precision floats.
 * @return Returns the number of rows read.
 */


    DISTANCE_TEMPLATE(evolution_dist,
		      double val,,
		      0,
	      	      if (!fscanf(fp, "%lf\n", &val)) { fatal_msg("Read zero items"); }
		      dist += !(val == a[i]), (*D)[k++] = dist;,)



/**
 * Calculates a bitwise yule distance of an FFP matrix
 *
 *
 * The FFP is read from fp, and the distance calculation is
 * perfomed for all pairs of FFP.  The form of the FFP 
 * must be a columnar matrix, but not necessarily row normalized.
 * In other words, raw frequencies can be used, for example
 * output directly from ffpcol.
 *
 * @param fp A file pointer to a columnar FFP
 * @param D Points to a dynamically allocated array of double precision floats.
 * @return Returns the number of rows read.
 */
    DISTANCE_TEMPLATE(yule_dist, double val;
		      double b = 0;
		      double c = 0;
		      double d = 0,,
		      (matrix_mode == similarity),
	      	      if (!fscanf(fp, "%lf\n", &val)) { fatal_msg("Read zero items"); }
		      dist += (val && a[i]);
		      b += (!val && !a[i]);
		      c += (!val && a[i]); d += (val && !a[i]);
		      , if (matrix_mode == similarity)
		      (*D)[k++] = (dist * b - c * d) / (dist * b + c * d);
		      else
		      (*D)[k++] =
		      1 - pow((dist * b - c * d) / (dist * b + c * d), 2);
		      b = 0; c = 0; d = 0;,)



/**
 * Calculates a bitwise Russel/Rao distance of an FFP matrix
 *
 *
 * The FFP is read from fp, and the distance calculation is
 * perfomed for all pairs of FFP.  The form of the FFP 
 * must be a columnar matrix, but not necessarily row normalized.
 * In other words, raw frequencies can be used, for example
 * output directly from ffpcol.
 *
 * @param fp A file pointer to a columnar FFP
 * @param D Points to a dynamically allocated array of double precision floats.
 * @return Returns the number of rows read.
 */


    DISTANCE_TEMPLATE(russel_dist, double val;
		      double b = 0;
		      ,,
		      (matrix_mode == similarity),
	      	      if (!fscanf(fp, "%lf\n", &val)) { fatal_msg("Read zero items"); }
		      dist += (val && a[i]); b++, if (matrix_mode == similarity)
		      (*D)[k++] = dist / (b);
		      else
		      (*D)[k++] = 1 - dist / (b); b = 0;,)


/**
 * Calculates a bitwise matching distance of an FFP matrix
 *
 *
 * The FFP is read from fp, and the distance calculation is
 * perfomed for all pairs of FFP.  The form of the FFP 
 * must be a columnar matrix, but not necessarily row normalized.
 * In other words, raw frequencies can be used, for example
 * output directly from ffpcol.
 *
 * @param fp A file pointer to a columnar FFP
 * @param D Points to a dynamically allocated array of double precision floats.
 * @return Returns the number of rows read.
 */


    DISTANCE_TEMPLATE(matching_dist, double val;
		      double b = 0;
		      double c = 0;
		      ,,
		      (matrix_mode == similarity),
	      	      if (!fscanf(fp, "%lf\n", &val)) { fatal_msg("Read zero items"); }
		      dist += (val && a[i]);
		      b += (!val && !a[i]); c++, if (matrix_mode == similarity)
		      (*D)[k++] = (dist + b) / c;
		      else
		      (*D)[k++] = 1 - pow((dist + b) / c, 2); b = 0; c = 0;,)


/**
 * Calculates a bitwise Hamann distance of an FFP matrix
 *
 *
 * The FFP is read from fp, and the distance calculation is
 * perfomed for all pairs of FFP.  The form of the FFP 
 * must be a columnar matrix, but not necessarily row normalized.
 * In other words, raw frequencies can be used, for example
 * output directly from ffpcol.
 *
 * 
 * @param fp A file pointer to a columnar FFP
 * @param D Points to a dynamically allocated array of double precision floats.
 * @return Returns the number of rows read.
 */


    DISTANCE_TEMPLATE(hamann_dist, double val;
		      double b = 0;
		      double c = 0;
		      double d = 0,,
		      (matrix_mode == similarity),
	      	      if (!fscanf(fp, "%lf\n", &val)) { fatal_msg("Read zero items"); }
		      dist += (val && a[i]);
		      b += (!val && !a[i]);
		      c += (!val && a[i]);
		      d += (val && !a[i]), if (matrix_mode == similarity)
		      (*D)[k++] = ((dist + b) - (c + d)) / (dist + b + c + d);
		      else
		      (*D)[k++] =
		      1 - pow(((dist + b) - (c + d)) / (dist + b + c + d), 2);
		      b = 0; c = 0; d = 0;,)


/**
 * Calculates a bitwise antidice distance of an FFP matrix
 *
 *
 * The FFP is read from fp, and the distance calculation is
 * perfomed for all pairs of FFP.  The form of the FFP 
 * must be a columnar matrix, but not necessarily row normalized.
 * In other words, raw frequencies can be used, for example
 * output directly from ffpcol.
 *
 * @todo comparing two zero vectors is undefined and equals zero.
 * @param fp A file pointer to a columnar FFP
 * @param D Points to a dynamically allocated array of double precision floats.
 * @return Returns the number of rows read.
 */


    DISTANCE_TEMPLATE(antidice_dist, double val;
		      double b = 0;
		      double c = 0;
		      double d = 0,,
		      (matrix_mode == similarity),
	      	      if (!fscanf(fp, "%lf\n", &val)) { fatal_msg("Read zero items"); }
		      dist += (val && a[i]);
		      b += (!val && !a[i]);
		      c += (!val && a[i]);
		      d += (val && !a[i]), if (matrix_mode == similarity)
		      (*D)[k++] = dist / (dist + 2 * (c + d));
		      else
		      (*D)[k++] = 1 - dist / (dist + 2 * (c + d));
		      b = 0; c = 0; d = 0;,)


/**
 * Calculates a bitwise sneath distance of an FFP matrix
 *
 *
 * The FFP is read from fp, and the distance calculation is
 * perfomed for all pairs of FFP.  The form of the FFP 
 * must be a columnar matrix, but not necessarily row normalized.
 * In other words, raw frequencies can be used, for example
 * output directly from ffpcol.
 *
 * @param fp A file pointer to a columnar FFP
 * @param D Points to a dynamically allocated array of double precision floats.
 * @return Returns the number of rows read.
 */


    DISTANCE_TEMPLATE(sneath_dist, double val;
		      double b = 0;
		      double c = 0;
		      double d = 0,,
		      (matrix_mode == similarity),
	      	      if (!fscanf(fp, "%lf\n", &val)) { fatal_msg("Read zero items"); }
		      dist += (val && a[i]);
		      b += (!val && !a[i]);
		      c += (!val && a[i]);
		      d += (val && !a[i]), if (matrix_mode == similarity)
		      (*D)[k++] = 2 * (dist + b) / (2 * (dist + b) + (c + d));
		      else
		      (*D)[k++] =
		      1 - 2 * (dist + b) / (2 * (dist + b) + (c + d)); b = 0;
		      c = 0; d = 0;,)



/**
 * Calculates a bitwise ochiai distance of an FFP matrix
 *
 *
 * The FFP is read from fp, and the distance calculation is
 * perfomed for all pairs of FFP.  The form of the FFP 
 * must be a columnar matrix, but not necessarily row normalized.
 * In other words, raw frequencies can be used, for example
 * output directly from ffpcol.
 *
 * @todo similarity of zero to zero is one, zero to non-zero is zero
 * @param fp A file pointer to a columnar FFP
 * @param D Points to a dynamically allocated array of double precision floats.
 * @return Returns the number of rows read.
 */


    DISTANCE_TEMPLATE(ochiai_dist, double val;
		      double c = 0;
		      double d = 0,,
		      (matrix_mode == similarity),
	      	      if (!fscanf(fp, "%lf\n", &val)) { fatal_msg("Read zero items"); }
		      dist += (val && a[i]);
		      c += (!val && a[i]);
		      d += (val && !a[i]), if (matrix_mode == similarity)
		      (*D)[k++] = dist / sqrt((dist + c) * (dist + d));
		      else
		      (*D)[k++] = 1 - dist / sqrt((dist + c) * (dist + d));
		      c = 0; d = 0;,)




/**
 * Calculates a bitwise anderberg distance of an FFP matrix
 *
 *
 * The FFP is read from fp, and the distance calculation is
 * perfomed for all pairs of FFP.  The form of the FFP 
 * must be a columnar matrix, but not necessarily row normalized.
 * In other words, raw frequencies can be used, for example
 * output directly from ffpcol.
 *
 * @todo all zero or all 1 then 1; a+b, a+c,c+d, b+d is zero then 0.
 * @param fp A file pointer to a columnar FFP
 * @param D Points to a dynamically allocated array of double precision floats.
 * @return Returns the number of rows read.
 */


    DISTANCE_TEMPLATE(anderberg_dist, double val;
		      double b = 0;
		      double c = 0;
		      double d = 0,,
		      (matrix_mode == similarity),
	      	      if (!fscanf(fp, "%lf\n", &val)) { fatal_msg("Read zero items"); }
		      dist += (val && a[i]);
		      d += (!val && !a[i]);
		      b += (!val && a[i]);
		      c += (val && !a[i]), if (matrix_mode == similarity)
		      (*D)[k++] =
		      (dist / (dist + b) + dist / (dist + c) + d / (c + d) +
		       d / (b + d)) / 4;
		      else
		      (*D)[k++] =
		      1 -
		      ((dist / (dist + b) + dist / (dist + c) + d / (c + d) +
			d / (b + d)) / 4); b = 0; c = 0; d = 0;,)





/**
 * Calculates a bitwise pearson phi distance of an FFP matrix
 *
 *
 * The FFP is read from fp, and the distance calculation is
 * perfomed for all pairs of FFP.  The form of the FFP 
 * must be a columnar matrix, but not necessarily row normalized.
 * In other words, raw frequencies can be used, for example
 * output directly from ffpcol.
 *
 * @todo  b+c=0 then 1 a+d=0 then -1, ad-bc=0 then 0
 * @param fp A file pointer to a columnar FFP
 * @param D Points to a dynamically allocated array of double precision floats.
 * @return Returns the number of rows read.
 */


    DISTANCE_TEMPLATE(phi_dist, double val;
		      double b = 0;
		      double c = 0;
		      double d = 0,,
		      (matrix_mode == similarity),
	      	      if (!fscanf(fp, "%lf\n", &val)) { fatal_msg("Read zero items"); }
		      dist += (val && a[i]);
		      d += (!val && !a[i]);
		      b += (!val && a[i]);
		      c += (val && !a[i]), if (matrix_mode == similarity)
		      (*D)[k++] =
		      (dist * d -
		       b * c) / sqrt((dist + b) * (dist + c) * (d + b) * (d +
									  c));
		      else
		      (*D)[k++] =
		      1 -
		      pow((dist * d -
			   b * c) / sqrt((dist + b) * (dist + c) * (d +
								    b) * (d +
									  c)),
			  2); b = 0; c = 0; d = 0;,)


/**
 * Calculates a bitwise gower distance of an FFP matrix
 *
 *
 * The FFP is read from fp, and the distance calculation is
 * perfomed for all pairs of FFP.  The form of the FFP 
 * must be a columnar matrix, but not necessarily row normalized.
 * In other words, raw frequencies can be used, for example
 * output directly from ffpcol.
 *
 * @todo  ad=0 then 0, all zero or all 1 then 1
 * @param fp A file pointer to a columnar FFP
 * @param D Points to a dynamically allocated array of double precision floats.
 * @return Returns the number of rows read.
 */


    DISTANCE_TEMPLATE(gower_dist, double val;
		      double b = 0;
		      double c = 0;
		      double d = 0,,
		      (matrix_mode == similarity),
	      	      if (!fscanf(fp, "%lf\n", &val)) { fatal_msg("Read zero items"); }
		      dist += (val && a[i]);
		      d += (!val && !a[i]);
		      b += (!val && a[i]);
		      c += (val && !a[i]), if (matrix_mode == similarity)
		      (*D)[k++] =
		      dist * d / sqrt((dist + b) * (dist + c) * (d + b) *
				      (d + c));
		      else
		      (*D)[k++] =
		      1 -
		      dist * d / sqrt((dist + b) * (dist + c) * (d + b) *
				      (d + c)); b = 0; c = 0; d = 0;,)


/**
 * Calculates a bitwise gower distance of an FFP matrix
 *
 *
 * The FFP is read from fp, and the distance calculation is
 * perfomed for all pairs of FFP.  The form of the FFP 
 * must be a columnar matrix, but not necessarily row normalized.
 * In other words, raw frequencies can be used, for example
 * output directly from ffpcol.
 *
 * @todo  zero to zero 1 , zero to nonzero=0
 * @param fp A file pointer to a columnar FFP
 * @param D Points to a dynamically allocated array of double precision floats.
 * @return Returns the number of rows read.
 */


    DISTANCE_TEMPLATE(kulczynski_dist, double val;
		      double b = 0;
		      double c = 0;
		      double d = 0,,
		      (matrix_mode == similarity),
	      	      if (!fscanf(fp, "%lf\n", &val)) { fatal_msg("Read zero items"); }
		      dist += (val && a[i]);
		      d += (!val && !a[i]);
		      b += (!val && a[i]);
		      c += (val && !a[i]), if (matrix_mode == similarity)
		      (*D)[k++] = (dist / (dist + b) + dist / (dist + c)) / 2;
		      else
		      (*D)[k++] =
		      1 - (dist / (dist + b) + dist / (dist + c)) / 2; b = 0;
		      c = 0; d = 0;,)



// gower
// mixed continuous/ binary 
//
// Renkonen
// percentage similarity
// sum all i min(fi1,fi2)
// bray curtis dissimilarity
// S10+S11 + S01S11 - 2*C
// ~~~~~~~~~~~~~~~~~~~~~
//        S10+S11 + S01S11
//
// C=sum all i min(fi1,fi2);
