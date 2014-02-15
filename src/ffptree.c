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
#include <time.h>
#include <stdlib.h>
#include <float.h>
#include <stdbool.h>
#include <getopt.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <limits.h>
#include <regex.h>
#include <errno.h>
#include "utils.h"
#include "vstring.h"
#include "sighandle.h"
#include "../config.h"

#define FNMLNGTH        200  /* length of array to store a file name */
#define NMLNGTH         50   /* number of characters in species name    */
#define MAXNCH          60   /* must be greater than or equal to NMLNGTH */
#define MAX_COL_NEWICK 55    /* Maximum characters output per line in outtree */
                             /* The last line seems to ignore this limit */
#define SMALL_TREE_TH 10     /* Threshold for drawing tree with -- vs - or "  " vs " " */
#define MAX_COL_TREE 55      /* Maximum number of chracters per line for screen */
#define TREE_SCALE 0.43429448222   // Float width for Newick tree 
#define NON_ZERO 0.000000001 /* Maximum allowable difference between upper and 
                                lower triangle to allow for matrix symmetricity */
#define DOWN            2       /* For drawing trees: vertical distance between branches */
#define OVER          60      /* maximum width all branches of tree on screen */
#define LEFT_MARGIN 0.5       /* Margin for setting */
#define SHORT_DASH "-"
#define LONG_DASH "--"
#define SHORT_SPACE " "
#define LONG_SPACE "  "


typedef double * DBLVECTOR;    // Move to a .h file
typedef int * INTVECTOR;

/* NODE data structure */

typedef struct node {
  struct node *next, *back;
  int  index;
  double xcoord, ycoord;
  int ymin, ymax;             /* used by printree()  */
  double v;                             
  bool tip;
} NODE;

typedef NODE **pointarray;

/* TREE data structure */

typedef struct tree {

  /* An array of pointers to nodes. Each tip node and ring of nodes has a
   * unique index starting from one. The nodep array contains pointers to each
   * one, starting from 0. In the case of internal nodes, the entries in nodep
   * point to the rootward node in the group. Since the trees are otherwise
   * entirely symmetrical, except at the root, this is the only way to resolve
   * parent, child, and sibling relationships.
   *
   * Indices in range [0, txn) point to tips, while indices [txn, nonodes)
   * point to internal nodes
   */
  NODE ** nodep;
  NODE  * root;   // For rooted trees.  Points to internal node w/o back ptr.                 
  NODE  * start;  // For unrooted trees.  Points to outgroup node.                  

} TREE;




/* function prototypes */
void allocrest(void);
void doinit(void);
void getinput(void);
void describe(NODE *, double);
void summarize(void);
void jointree(void);
void maketree(void);
void freerest(void);
int readNumTaxa( void );
NODE ** allocTree(int nonodes);
void chkTxnNumEq(int ith);
void setupTree(TREE *a, int nonodes);
void freetree(NODE **treenode, int nonodes);
void connect(NODE *p, NODE *q);
void inputdata(bool, bool, bool,bool, DBLVECTOR *, INTVECTOR *);
void printree(NODE *);
void treeout(NODE *, int *,   NODE *);
void treeoutr(NODE *, int *, TREE *);
void readName(int );
double randum(INTVECTOR);
void coordinates(NODE *, double, int *, double *,NODE *);
void drawline(int i, double scale, NODE *start);
ssize_t getline(char **lineptr, size_t *n, FILE *stream);
void shuffle(int * a,int n);

char PROG_NAME[FILENAME_MAX];



char usage_str[] = "Usage: %s [OPTIONS]... [FILE]... \n\
This program generates trees in NEWICK format using either the\n\
neighbor joining or UPGMA methods. Given no options the default\n\
behavior is to generate a neighbor joining tree.\n\n\
\t-n, --upgma\t\tBuild a UPGMA tree.\n\
\t-j[N], --jumble=[N]\tJumble input order. Optional seed, N\n\
\t-m[N], --multiple=[N]\tUse up to N sets in input. Optional arg, N\n\
\t-l, --lower\t\tUse lower input matrix triangle.\n\
\t-l, --upper\t\tUse upper input matrix triangle. [Non functional]\n\
\t-t, --print-tree\tDisable human readable tree (stderr)\n\
\t-p, --progress\t\tDisable tree-build progress (stderr)\n\
\t-d, --print-data\tEnable input summary (stderr)\n\
\t-y, --symmetrize\tSymmetrize assymetric distance matrix\n\
\t-o I, --outgroup=I\tSet outgroup. Required species num, I\n\
\t-O FILE, --out=FILE\tWrite newick tree to file. (stdout).\n\
\t-P FILE, --out-prg=FILE\tRedirect stderr to file. (stderr).\n\
\t-w W, --precision=W\tSpecify float precision in tree [W=8].\n\
\t-q, --quiet\t\tSuppress output to stderr\n\
\t-h, --help\t\tThis message.\n\
\t-v, --version\t\tVersion Info.\n\n\
Copyright (c) %s\n\
%s\n\
Contact %s\n";

bool njoin  = true;       //Option -n
bool jumble = false;      //Option -j
unsigned seed=1;
int outgrno = 0;          //Option -o Arg
bool outgropt = false;    //Option -o
bool lower = false;       //Option -l
bool replicates = false;  //Option -s
bool symmetrize = false;  //option -y
bool trout = true;        //Option -w
bool upper = false;       //option -u
bool treeprint = true;    //option -t
bool progress = true;     //option -p
bool printdata = false;   //option -d
bool mulsets = false;     //Option -m
int  datasets = INT_MAX;        //Option -m Arg
bool quiet = false;
int precision = 8;  //Option -w arg


FILE * infile;
FILE * outfile;
FILE * outtree;
int txn; // Number of taxa
char * buffer=NULL;

/**@todo eliminate some of these globals or at least package them
 *  in a passable structure */

char infilename[FNMLNGTH], outfilename[FNMLNGTH], outtreename[FNMLNGTH];
int numnodes, col, datasets, ith;
DBLVECTOR *x;
INTVECTOR *reps;
TREE curtree;
int *taxaorder;
char ** name;                     /* taxa names */

NODE **cluster;  //used in maketree


int main(int argc, char *argv[])
{  /* main program */
    int opt;
    int option_index = 0;
   // FILE * infile;

    static struct option long_options[] = {
	{"multiple", optional_argument, 0, 'm'},
	{"upgma", no_argument, 0, 'n'},
	{"jumble", optional_argument, 0, 'j'},
	{"outgroup", required_argument, 0, 'o'},
	{"lower", no_argument, 0, 'l'},
	{"upper", no_argument, 0, 'u'},
	{"print-tree", no_argument, 0, 't'},
	{"progress", no_argument, 0, 'p'},
	{"print-data", no_argument, 0, 'd'},
	{"symmetrize", no_argument, 0, 'y'},
	{"out", required_argument, 0, 'O'},
	{"out-prg", required_argument, 0, 'P'},
	{"precision", required_argument, 0, 'w'},
	{"quiet", no_argument, 0, 'q'},
	{"help", no_argument, 0, 'h'},
	{"version", no_argument, 0, 'v'},
	{0, 0, 0, 0}
    };

  initSignalHandlers();

  strcpy(PROG_NAME,basename( argv[0] ));

  outfile=stderr;
  outtree=stdout;
  strcpy(outfilename,"stderr");
  strcpy(outtreename,"stdout");
  // add option for user readable tree.

  while ((opt = getopt_long(argc, argv, "m::nj::o:lutpdi:O:P:qhvw:y",
			      long_options, &option_index)) != -1)
	switch (opt) {
	case 'm':
	    if (optarg)  
	    	datasets=atoi(optarg);
    	    break;	    
	case 'n':
	    njoin=!njoin; // Toggle
    	    break;	    
	case 'o':
	    outgropt=!outgropt;
	    outgrno=atoi(optarg)-1;
    	    break;	    
        case 'j':
	    jumble=!jumble;
	    if (optarg)
		seed=(unsigned)atoi(optarg);
	    else
		seed=(unsigned)time(NULL) * getpid();
	break;
	case 'l':
	    lower=!lower;
    	    break;
	case 'u':
	    upper=!upper;
    	    break;
	case 't':  
	    treeprint=!treeprint;
    	    break;
	//case 'w':  // could change to O
	 //   trout=!trout;
    	 //   break;
	case 'd':
	    printdata=!printdata;
    	    break;
	case 'p':
	    progress=!progress;
    	    break;
	case 'q':
	    progress=!progress;
	    treeprint=!treeprint;
	break; 
        case 'O':
		strcpy(outtreename,optarg);
		if ((outtree=fopen(optarg,"w")) == NULL)
			fatal_msg("%s: Failed to open.",optarg);	
	break;	
        case 'P':
		strcpy(outfilename,optarg);
		if ((outfile=fopen(optarg,"w")) == NULL)
			fatal_msg("%s: Failed to open.",optarg);	
	break;
	case 'y':
		symmetrize=!symmetrize;
	break;
	case 'w':
		precision=atoi(optarg);
	break;
	case 'v':
	    printVersion();
	    exit(EXIT_SUCCESS);
	    break;
	case 'h':
	    printUsageStr();
	    exit(EXIT_SUCCESS);
	    break;
	default:
	    printErrorUsageStr();
	    exit(EXIT_FAILURE);
	    break;
	}

  if (outgrno < 0 )  //must check again to make sure less than txn
	fatal_msg("%ld: Outgroup number must be at least 1\n",outgrno);


  if (mulsets && datasets < 2) 
	fatal_msg("%ld: Number of sets must be greater than 1\n",datasets);
	  
  if (jumble)
   	srand(seed);


  argv+=optind;


  do {
	infile = stdin;
	if (*argv) {
	    if (!strcmp(*argv, "-"))
		infile = stdin;
	    else if ((infile = fopen(*argv, "r")) == NULL)
		fatal_msg("%s: %s.\n", *argv,strerror(errno));

	    if ( isDirectory(*argv) ) 
	    	fatal_msg("%s: %s\n",*argv, strerror(EISDIR));

	    argv++;
	} else if (isatty(STDIN_FILENO))
	    printErrorUsageStr();

  doinit(); // This allocates memory for the tree.

  ith = 1;
  do {
    if (progress && mulsets )
        fprintf(outfile,"Data set # %d:\n",ith);
    if (ith != 1)
        chkTxnNumEq(ith); // Confirm # of taxa is the same as txn
    maketree();
    ith++;
    ungetc(fgetc(infile),infile);
  } while ( !feof(infile) && ith <= datasets );
  freerest();
  freetree(curtree.nodep, numnodes+1); 



	if (infile != stdin)
	    fclose(infile);

    }
    while (*argv);


  free(buffer);

  fclose(outfile);
  fclose(outtree);
  return EXIT_SUCCESS;
}



void allocrest()
{
  int i;

  x = (DBLVECTOR *)chkcalloc(sizeof(DBLVECTOR),txn);
  reps = (INTVECTOR *)chkcalloc(sizeof(INTVECTOR),txn);
  name = (char **)chkcalloc(sizeof(char*),txn); 
  for (i = 0; i < txn; i++) {
    x[i] = (DBLVECTOR)chkcalloc(sizeof(double),txn);
    reps[i] = (INTVECTOR)chkcalloc(sizeof(int),txn);
    name[i]=(char *)chkcalloc(sizeof(char),MAXNCH);
  }
	  
  taxaorder = (int *)chkcalloc(sizeof(int),txn);
  cluster = (NODE **)chkcalloc(sizeof(NODE *),txn);
}  


void freerest()
{
  int i;

  for (i = 0; i < txn; i++) {
    free(x[i]);
    free(reps[i]);
  }

  free(x);
  free(reps);
  free(name);
  free(taxaorder);
  free(cluster);
}  


void doinit()
{
  /* initializes variables */
  //NODE *p;
  txn=readNumTaxa();
  if (progress)
  	fprintf(outfile, "%d Taxa\n", txn);
  numnodes = 2 * txn - 2;  // number of nodes in an unrooted bifurcating tree is always 2n-2
  numnodes += (njoin ? 0 : 1);  // Add 1 node (for root) if the upgma method
  curtree.nodep =  allocTree(numnodes+1); 
  //This code is unnecessary. It is deallocating memory which just got allocated.
  //Why not just not allocate in the first place.
  //This makes no sense....
  //create a circular connected list
  //p = curtree.nodep[numnodes]->next;   // Assign p the address pointed to by the Nth element next potr 
 // curtree.nodep[numnodes]->next = curtree.nodep[numnodes]; //Make the Nth element point to itself
 // free(p->next);                                           // Free curtree.nodep[numnodes]->next->next
 // free(p);                                                 // Free curtree.nodep[numnodes]->next
  allocrest();
  
} 


//nice use of recursion.

/* print out information for one branch */
void describe(NODE *p, double height)
{
  NODE *q;

  q = p->back;
  if (njoin)
    fprintf(outfile, "%d\t", q->index - txn);
  else
    fprintf(outfile, "%d\t", q->index - txn);
  if (p->tip) {
      fprintf(outfile,"%10s\t",name[p->index]);	  
  } else {
      fprintf(outfile, "%10d\t", p->index - txn);
  }
  if (njoin)
    fprintf(outfile, "%.4e\n", q->v);
  else
    fprintf(outfile, "%6.4e\t%6.4e\n", q->v, q->v+height);
  if (!p->tip) {
    describe(p->next->back, height+q->v);
    describe(p->next->next->back, height+q->v);
  }
}  


/* print out branch lengths etc. */
void summarize()
{
  putc('\n', outfile);
  if (njoin) {
    if (outgropt)
      fprintf(outfile, "Rearranged with outgroup at root.\n");
    fprintf(outfile, "Neighbor joining is an unrooted method.\n");
  }

  fprintf(outfile, "\n\tNodes\n");
  fprintf(outfile, "------------------\n");
  if (njoin) {
    fprintf(outfile, "i\t\t j\tLength\n");
    fprintf(outfile, "-----------------------------------\n");
  } else {
    fprintf(outfile, "i\t\tj\tLength\t\tRoot_Length\n");
    fprintf(outfile, "--------------------------------------------------\n");
  }
  describe(curtree.start->next->back, 0.0);
  describe(curtree.start->next->next->back, 0.0);
  if (njoin)
    describe(curtree.start->back, 0.0);
  fprintf(outfile, "\n\n");
}  /* summarize */



#define printNodeLabel(a) fputc( (a) ? 'N':'T',outfile)  

  /* calculate the tree */
void jointree()
{
  int c, nextnode, mini=0, minj=0, i, j, ia, ja, ii, jj, nude, cycles;
  double otu, q, qmin, dio, djo, bi, bj, bk, dmin=0, da;
  int el[3];
  DBLVECTOR av;
  INTVECTOR oc;

  double *R;   
  R = (double *)chkmalloc(sizeof(double),txn);

  if (symmetrize) 
   for (i = 0; i < txn - 1; i++) 
     for (j = i + 1; j < txn; j++) {
       da = (x[i][j] + x[j][i]) / 2.0;
       x[i][j] = da;
       x[j][i] = da;
     }

 if (progress) { 
 	fprintf(outfile,"Cycle\tType\ti\tLength\t        Type\tj\tLength\n");
 }	fprintf(outfile,"----------------------------------------------------------------\n");

  // First initialization 
  otu = txn - 2.0;
  nextnode = txn;
  av = (DBLVECTOR)chkmalloc(sizeof(double),txn);
  oc = (INTVECTOR)chkmalloc(sizeof(int),txn);

  for (i = 0; i < txn; i++) {
    av[i] = 0.0;
    oc[i] = 1;
  }

  // Enter the main cycle 
  if (njoin)
    cycles = txn - 3;
  else
    cycles = txn - 1;
  for (c = 0; c < cycles; c++) {
    for (j = 1; j < txn; j++) {
      for (i = 0; i < j; i++)
        x[j][i] = x[i][j];
    }
    qmin = DBL_MAX;

  // Compute Row sum of observable taxonomic units (otu)
  // If group has been joined use the aggregate cluster.

    if (njoin) {     
      for (i = 0; i < txn; i++)
        R[i] = 0.0;

      for (ja = 0; ja < txn; ja++) {
        jj = taxaorder[ja];
        if (cluster[jj] != NULL) 
          for (ia = 0; ia < ja; ia++) {
            ii = taxaorder[ia];
            if (cluster[ii] != NULL) {
              R[ii] += x[ii][jj];
              R[jj] += x[ii][jj];
            }
          }
      }
    }

    // Compute Q matrix
    for (ja = 0; ja < txn; ja++) {
      jj = taxaorder[ja];
      if (cluster[jj] != NULL) {
        for (ia = 0; ia < ja; ia++) {
          ii = taxaorder[ia];
          if (cluster[ii] != NULL) {
            if (njoin) {
              q = otu * x[ii][jj] - R[ii] - R[jj];
            } else
              q = x[ii][jj];
            if (q < qmin) {
              qmin = q;
              mini = ii;
              minj = jj;
            }
          }
        }
      }
    }
    
    // compute lengths and print 
    if (njoin) {
      dio = 0.0;
      djo = 0.0;
      for (i = 0; i < txn; i++) {
        dio += x[i][mini];
        djo += x[i][minj];
      }
      dmin = x[mini][minj];
      dio = (dio - dmin) / otu;
      djo = (djo - dmin) / otu;
      bi = (dmin + dio - djo) * 0.5;
      bj = dmin - bi;
      bi -= av[mini];
      bj -= av[minj];
    } else {
      bi = x[mini][minj] / 2.0 - av[mini];
      bj = x[mini][minj] / 2.0 - av[minj];
      av[mini] += bi;
    }
    if (progress) {
      fprintf(outfile,"%d\t", cycles - c );
      if (njoin)
        printNodeLabel(av[mini] > 0.0 );
      else
        printNodeLabel(oc[mini] > 1.0);
      fprintf(outfile,"\t%d\t%.2e\t", mini+1, bi);
      if (njoin)
        printNodeLabel(av[minj] > 0.0);
      else
        printNodeLabel(oc[minj] > 1.0);
      fprintf(outfile,"\t%d\t%.2e\n", minj+1, bj);
    }
    connect(curtree.nodep[nextnode]->next, cluster[mini]);
    connect(curtree.nodep[nextnode]->next->next, cluster[minj]);
    cluster[mini]->v = bi;
    cluster[minj]->v = bj;
    cluster[mini]->back->v = bi;
    cluster[minj]->back->v = bj;
    cluster[mini] = curtree.nodep[nextnode];
    cluster[minj] = NULL;
    nextnode++;
    if (njoin)
      av[mini] = dmin * 0.5;
    
    // re-initialization 
    otu -= 1.0;
    for (j = 0; j < txn; j++) {
      if (cluster[j] != NULL) {
        if (njoin) {
          da = (x[mini][j] + x[minj][j]) * 0.5;
          if (mini - j < 0)
            x[mini][j] = da;
          if (mini - j > 0)
            x[j][mini] = da;
        } else {
          da = x[mini][j] * oc[mini] + x[minj][j] * oc[minj];
          da /= oc[mini] + oc[minj];
          x[mini][j] = da;
          x[j][mini] = da;
        }
      }
    }
    for (j = 0; j < txn; j++) {
      x[minj][j] = 0.0;
      x[j][minj] = 0.0;
    }
    oc[mini] += oc[minj];
  }
  // Final cycle 
  nude = 0;
  for (i = 0; i < txn; i++) {
    if (cluster[i] != NULL) {
      el[nude] = i; 
      nude++;
    }
  }
  if (!njoin) {
    curtree.start = cluster[el[0]];
    curtree.start->back = NULL;
    free(av);
    free(oc);
    return;
  }
  bi = (x[el[0]][el[1]] + x[el[0]][el[2]] - x[el[1]]
        [el[2]]) * 0.5;
  bj = x[el[0]][el[1]] - bi;
  bk = x[el[0]][el[2]] - bi;
  bi -= av[el[0]];
  bj -= av[el[1]];
  bk -= av[el[2]];
  if (progress) {
    fprintf(outfile,"0\t");
    printNodeLabel(av[el[0]] > 0.0);
    fprintf(outfile,"\t%d\t%.2e\t", el[0]+1, bi);
    printNodeLabel(av[el[1]] > 0.0);
    fprintf(outfile,"\t%d\t%.2e\t", el[1]+1, bj);
    printNodeLabel(av[el[2]] > 0.0);
    fprintf(outfile,"\t%d\t%.2e\n", el[2]+1, bk);
  }
  connect(curtree.nodep[nextnode], cluster[el[0]]);
  connect(curtree.nodep[nextnode]->next, cluster[el[1]]);
  connect(curtree.nodep[nextnode]->next->next, cluster[el[2]]);
  cluster[el[0]]->v = bi;
  cluster[el[1]]->v = bj;
  cluster[el[2]]->v = bk;
  cluster[el[0]]->back->v = bi;
  cluster[el[1]]->back->v = bj;
  cluster[el[2]]->back->v = bk;
  curtree.start = cluster[el[0]]->back;
  free(av);
  free(oc);
  free(R);
}  


 

  /* Build the tree */
void maketree()
{
  int i;

  inputdata(replicates, printdata, lower, upper, x, reps);
  if (njoin && (txn < 3)) 
    fatal_msg("\nMust have at least 3 taxa.\n");
  
  if (progress)
    fprintf(outfile,"\n");

  if (ith == 1)
    setupTree(&curtree, numnodes + 1);

  for (i = 0; i < txn; i++)
    taxaorder[i] = i;

  if (jumble)
    shuffle(taxaorder,txn);

  for (i = 0; i < txn; i++)
    cluster[i] = curtree.nodep[i];
  jointree();
  if (njoin)
    curtree.start = curtree.nodep[outgrno]->back;
  if (treeprint)
  	printree(curtree.start);
  if (treeprint)
    summarize();
  if (trout) {
    col = 0;
    if (njoin)
      treeout(curtree.start, &col, curtree.start);
    else
      curtree.root = curtree.start,
      treeoutr(curtree.start,&col,&curtree);
  }
} 


int readNumTaxa( void )
{
  int numTaxa;
  char newline[2];
  /* read species number */
  if (fscanf(infile, "%d%[\n\r]", &numTaxa,newline) != 2 || numTaxa <= 0) 
	fatal_msg("Can not read the number of species in input\n");
  return numTaxa;
} 




/* Allocate taxa nodes and internal nodes 
* treenode is an array of pointers to nodes 
* Taxa nodes are stored from 0 to txn-1 and
* internal nodes stored from txn to nonodes-1 
*
* structure looks like this:
*
* [0] ->   Node  
* [1]
*   
*  Each pointer points to a connected list of three
*  nodes.
*/



NODE ** allocTree(int n)
{
  NODE ** treenode;
  int i, j;
  NODE *p, *q;

  treenode = (NODE **)chkmalloc(sizeof(NODE *),n); 
  for (i = 0; i < txn; i++)
   treenode[i] = (NODE *)chkmalloc(sizeof(NODE),1);

  // For each internal tree node create a circular 
  // connected list or 'ring' of three nodes.
  for (i = txn; i < n; i++) {
    q = NULL;
    for (j = 0; j < 3; j++) {
      p = (NODE *)chkmalloc(sizeof(NODE),1);
      p->next = q;
      q = p;
    }
    p->next->next->next = p;
    treenode[i] = p;
  }
return treenode;  
} 



void chkTxnNumEq(int ith)
{
  /* check if txn is same as the first set in other data sets */
  int curtxn;
  char c[3];

  if (fscanf(infile, "%d%[\n]", &curtxn,c) != 2) 
    fatal_msg("Set: %d: Unable to read taxa number.\n",ith);
  if (curtxn != txn) 
    fatal_msg("Set: %d: Inconsist taxa number.\n",ith);	  
} 

/* initialize a tree */
void setupTree(TREE *a, int n)
{
  int i=0;
  NODE *p;

  for (i = 0; i < n; i++) {
    a->nodep[i]->back = NULL;
    a->nodep[i]->tip = (i < txn);
    a->nodep[i]->index = i;
    a->nodep[i]->v = 0.0;
    if (i >= txn) {
      p = a->nodep[i]->next;
      while (p != a->nodep[i]) {
        p->back = NULL;
        p->tip = false;
        p->index = i;
        p = p->next;
      }
    }
  }
  a->start = a->nodep[0];
  a->root = NULL;
} 

void freetree(NODE **treenode, int n)
{
  int i;
  NODE *p, *q;

  for (i = 0; i < txn; i++)
    free(treenode[i]);
  for (i = txn; i < n; i++) {
    p = treenode[i];
    q = p->next;
    while(q != p) {
        NODE * r = q;
        q = q->next;
        free(r);
    }
    free(p);
  }
  free(treenode);
} 


 /* connect two nodes */
void connect(NODE *p, NODE *q)
{
  p->back = q;
  q->back = p;
}  



  /* read in distance matrix */

void inputdata(bool replicates, bool printdata, bool lower,
                        bool upper, DBLVECTOR *x, INTVECTOR *reps)
{
  int i=0, j=0, k=0, columns=0;
  bool skipit=false, skipother=false;
  char c[3];

  if (replicates)
    columns = 4;
  else
    columns = 6;
  if (printdata) {
    fprintf(outfile, "\nName                       Distances");
    if (replicates)
      fprintf(outfile, " (replicates)");
    fprintf(outfile, "\n----                       ---------");
    if (replicates)
      fprintf(outfile, "-------------");
    fprintf(outfile, "\n\n");
  }
  for (i = 0; i < txn; i++) {
    x[i][i] = 0.0;
    readName(i);
    for (j = 0; j < txn; j++) {
      skipit = ((lower && j + 1 >= i + 1) || (upper && j + 1 <= i + 1));
      skipother = ((lower && i + 1 >= j + 1) || (upper && i + 1 <= j + 1));
      if (!skipit) {
        if (fscanf(infile, "%lf%*[ ]%[\n]", &x[i][j],c) < 1)  
          fatal_msg("The infile is of the wrong type\n");
        if (replicates) { // decide how replicates are handled.
          if (fscanf(infile, "%d", &reps[i][j]) != 1) 
            fatal_msg("The infile is of the wrong type\n");
        } else
          reps[i][j] = 1;
      }
      if (!skipit && skipother) {
          x[j][i] = x[i][j];
          reps[j][i] = reps[i][j];
      }
      if ((i == j) && (fabs(x[i][j]) > NON_ZERO)) 
       fatal_msg("Diagonal of row %d from input matrix is not zero.", i+1);
     
      if (!symmetrize)  
	      if ((j < i) && (fabs(x[i][j]-x[j][i]) > NON_ZERO))
       	 	fatal_msg("Matrix is assymetric: (%d,%d) not equal to (%d,%d).\n",
               		i+1, j+1, j+1, i+1);
    }
  }
  if (!printdata)
    return;
  for (i = 0; i < txn; i++) {
    for (j = 0; j < NMLNGTH; j++)
      putc(name[i][j], outfile);
    putc(' ', outfile);
    for (j = 0; j < txn; j++) {
      fprintf(outfile, "%10.5f", x[i][j]);
      if (replicates)
        fprintf(outfile, " (%3d)", reps[i][j]);
      if (j % columns == 0 && j < txn) {
        putc('\n', outfile);
        for (k = 0; k < NMLNGTH + 1; k++)
          putc(' ', outfile);
      }
    }
    putc('\n', outfile);
  }
  putc('\n', outfile);
} 


/***********************************
*
* shuffle()
* Performs a Fischer-Yates shuffle
* Durstenfeld, Richard (July 1964). 
* Algorithm 235: Random permutation. 
* Communications of the ACM 7, 420.
*
*************************************/

void shuffle(int * a,int n) {
	int k,tmp;
	while (n > 1) {
		n--;
		k = rand() % n;
		tmp=a[k];	
		a[k]=a[n];
		a[n]=tmp;	
	}	
}




  /* prints out diagram of the tree */
void printree(NODE *start)
{
  int i,tipy;
  double scale,tipmax;

  putc('\n', outfile);
  tipy = 1;
  tipmax = 0.0;
  coordinates(start, 0.0, &tipy, &tipmax, start);
  scale = 1.0 / (int)(tipmax + 1.0);
  for (i = 0; i <= (tipy - DOWN); i++)
    drawline(i, scale, start );
  putc('\n', outfile);
}  



/*UPGMA write out file with representation of final tree. */
   
void treeoutr(NODE *p, int *col, TREE *curtree)
{
  char * cptr;


  if (p->tip) {
    //replace spaces in name with underscores.
    while ((cptr=strchr(name[p->index],' ')) != NULL ) 
	    *cptr='_';

    fprintf(outtree,"%s",name[p->index]); 

    (*col) += strlen(name[p->index]);
  } else {
    putc('(', outtree);
    (*col)++;
    treeoutr(p->next->back,col,curtree);
    putc(',', outtree);
    (*col)++;
    if ((*col) > MAX_COL_TREE ) {
      putc('\n', outtree);
      (*col) = 0;
    }
    treeoutr(p->next->next->back,col,curtree);
    putc(')', outtree);
    (*col)++;
  }
  if (p == curtree->root)
    fprintf(outtree, ";\n");
  else {
    fprintf(outtree, ":%.*e", precision,p->v);
    *col += precision+5;
  }
} 


/*Neighbor joining: write out rerpesentation of tree */
void treeout(NODE *p, int *col, NODE *start)
{
  char * cptr;

  if (p->tip) {
    //replace spaces in name with underscores.
    while ((cptr=strchr(name[p->index],' ')) != NULL ) 
  	    *cptr='_';
    fprintf(outtree,"%s",name[p->index]);

    *col += strlen(name[p->index]);
  } else {
    putc('(', outtree);
    (*col)++;
    treeout(p->next->back, col, start);
    putc(',', outtree);
    (*col)++;
    if (*col > MAX_COL_NEWICK) {
      putc('\n', outtree);
      *col = 0;
    }
    treeout(p->next->next->back, col, start);
    if (p == start && njoin) {
      putc(',', outtree);
      if (*col > MAX_COL_NEWICK) {
        putc('\n', outtree);
        *col = 0;
      }
      treeout(p->back, col, start);
    }
    putc(')', outtree);
    (*col)++;
  }
  if (p == start)
    fprintf(outtree, ";\n");
  else {
    fprintf(outtree, ":%.*e", precision,p->v);
    *col += precision+5;
  }
}  

/* read in taxa name */
void readName(int i)
{
	regex_t    re;
	regmatch_t pm;
	char * c;

	if (regcomp(&re, "[]:;(),\n[]", 0) != 0) 
		fatal_msg("Failed to compile regular expression");

	if (!fread(name[i],sizeof(char),NMLNGTH,infile)) 
		fatal_msg("Read zero items.");

	if (regexec(&re, name[i], (size_t) 1, &pm, 0) == 0 ) 
		fatal_msg("Unexpected character: %c in taxa %s\n",name[i][pm.rm_so],name[i]);

	//trim terminating whitespace
	while ((c=strrchr(name[i],' ')) != NULL) 
		*c='\0';

	regfree(&re);

} 


/* establishes coordinates of nodes */
//Uses recursion
void coordinates(NODE *p, double lengthsum, int *tipy, double *tipmax,
                        NODE *start )
{
  NODE *q, *first, *last;

  if (p->tip) {
    p->xcoord = (int)(OVER * lengthsum + LEFT_MARGIN);
    p->ycoord = *tipy;
    p->ymin = *tipy;
    p->ymax = *tipy;
    (*tipy) += DOWN;
    if (lengthsum > *tipmax)
      *tipmax = lengthsum;
    return;
  }
  q = p->next;
  do {
    if (q->back)
      coordinates(q->back, lengthsum + q->v, tipy,tipmax, start);
    q = q->next;
  } while ((p == start || p != q) && (p != start || p->next != q));
  first = p->next->back;
  q = p;
  while (q->next != p && q->next->back)  /* is this right ? */
    q = q->next;
  last = q->back;
  p->xcoord = (int)(OVER * lengthsum + LEFT_MARGIN);
  if (p == start && p->back) 
    p->ycoord = p->next->next->back->ycoord;
  else
    p->ycoord = (first->ycoord + last->ycoord) / 2;
  p->ymin = first->ymin;
  p->ymax = last->ymax;
}  




  /* draws one row of the tree diagram by moving up tree */
void drawline(int i, double scale, NODE *start)
{
  NODE *p, *q;
  int n=0, j=0;
  bool extra=false, trif=false;
  NODE *r, *first =NULL, *last =NULL;
  bool done=false;

  p = start;
  q = start;
  extra = false;
  trif = false;
  if (i == (int)p->ycoord && p == start) {  /* display the root */
    if (!njoin) {
      if (p->index - txn >= SMALL_TREE_TH)
        fprintf(outfile, SHORT_DASH);
      else
        fprintf(outfile, LONG_DASH);
    }
    else {
      if (p->index - txn >= SMALL_TREE_TH)
        fprintf(outfile, SHORT_SPACE);
      else
        fprintf(outfile, LONG_SPACE);
    }
    if (p->index - txn >= SMALL_TREE_TH)
      fprintf(outfile, "%2d", p->index - txn + 1);
    else
      fprintf(outfile, "%d", p->index - txn + 1);
    extra = true;
    trif = true;
  } else
    fprintf(outfile, LONG_SPACE);
  do {
    if (!p->tip) { /* internal nodes */
      r = p->next;
      /* r->back here is going to the same node. */
      do {
        if (!r->back) {
          r = r->next;
          continue;
        }
        if (i >= r->back->ymin && i <= r->back->ymax) {
          q = r->back;
          break;
        }
        r = r->next;
      } while (!((p != start && r == p) || (p == start && r == p->next)));
      first = p->next->back;
      r = p;
      while (r->next != p) 
        r = r->next;
      last = r->back;
      if (njoin && (p == start))
        last = p->back;
    } /* end internal node case... */
    /* draw the line: */
    done = (p->tip || p == q);
    n = (int)(scale * (q->xcoord - p->xcoord) + LEFT_MARGIN);
    if (!q->tip) {
      if ((n < 3) && (q->index - txn >= SMALL_TREE_TH))
        n = 3;
      if ((n < 2) && (q->index - txn < SMALL_TREE_TH))
        n = 2;
    }
    if (extra) {
      n--;
      extra = false;
    }
    if ((int)q->ycoord == i && !done) {
      if (p->ycoord != q->ycoord)
        putc('+', outfile);
      if (trif) {
        n++;
        trif = false;
      } 
      if (!q->tip) {
        for (j = 2; j < n; j++)
          putc('-', outfile);
        if (q->index - txn >= SMALL_TREE_TH)
          fprintf(outfile, "%2d", q->index - txn + 1);
        else
          fprintf(outfile, "-%d", q->index - txn + 1);
        extra = true;
      } else {
        for (j = 1; j < n; j++)
          putc('-', outfile);
      }
    } else if (!p->tip) {
      if ((int)last->ycoord > i && (int)first->ycoord < i
           && i != (int)p->ycoord) {
        putc('|', outfile);
        for (j = 1; j < n; j++)
          putc(' ', outfile);
      } else {
        for (j = 0; j < n; j++)
          putc(' ', outfile);
        trif = false;
      }
    }
    if (q != p)
      p = q;
  } while (!done);
  if ((int)p->ycoord == i && p->tip) {
    fprintf(outfile,"%s",name[p->index]);	  
  }
  putc('\n', outfile);
} 


