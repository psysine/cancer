//#include <iostream>
#include <string>
#include <vector>
//#include <stack>
#include <math.h>
#include <error.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <limits.h>
using namespace std;
#define rep(i,n) for(int i = 0; i < (n); i++)
#define rap(i,n,m) for(int i = (n); i < (m); i++)
#define iter(it, v) for(typeof((v).begin()) it = (v).begin(); it != (v).end(); ++it)
#define MP(x,y) make_pair(x,y)
#define F first
#define S second
#define all(v) (v).begin(), (v).end() //ex: sort(all(v));
#define mset(v, byte) memset(v, byte, sizeof(v))
typedef unsigned int uint;
typedef unsigned long long ull;
typedef long long ll;
//typedef pair<int,int> pii;

string fname, region;
const int maxreads = 100000;

double phred2prob[255];

void precompute() { //called from main
  rep(i,255) phred2prob[i] = pow(10,-i/10.);
  }

void init(string cmd, string fname, string region, int *nfiles, vector<FILE*> *files) {
  FILE *fi = fopen(fname.c_str(), "r");
  if(!fi) error(1, errno, "Can't open list of bam files");
  char *samtools;
  fscanf(fi, "%a[^\n]%*c", &samtools);
  fscanf(fi, "%d", nfiles);
  if(*nfiles < 1) error(1, 0, "Need at least one file");

  rep(f, *nfiles) {
    char *bamfile = 0;
    if(fscanf(fi, "%as", &bamfile) != 1)
      error(1, 0, "Not enough filenames");
    FILE *file = popen((string(samtools)+" "+cmd+" -r "+region+" \""+bamfile+'"').c_str(), "r");
    if(file <= 0)
      error(1, errno, "Can't run samtools");
    files->push_back(file);
    free(bamfile);
    }
  free(samtools), samtools = 0;
  fclose(fi);
  }

void findcover(int argc, char **argv) {
  if(argc != 2) printf("Usage: findcover <mindepth> <minlength>\n"
"The output will be a list of regions that are at least <minlength> long and in which each base pair is covered by at least <mindepth> reads.\n"
"The first column is the starting posision of each region, the second column is the ending position and the third column is the length of the region.\n"), exit(1);
  int mindepth = atoi(argv[0]), minlen = atoi(argv[1]);
  if(mindepth <= 0) error(1, 0, "mindepth must be greater than zero");
  int nfiles;
  vector <FILE*> files;
  init("depth", fname, region, &nfiles, &files);

  int poss[nfiles]; //last read position
  mset(poss, 0);
  
  bool in = 0, end = 0;
  int startpos = 0;
  for(int p = 1; !end; p++) {
    bool allcovered = 1;
    rep(f, nfiles) {
      int depth = 0;
      if(poss[f] < p && fscanf(files[f], "%*s %d %d", poss+f, &depth) != 2)
        poss[f] = INT_MAX, end = 1;
      assert(poss[f] >= p);
      bool covered = poss[f] == p && depth >= mindepth;
      if(!covered) allcovered = 0;
      }

    if(!in && allcovered) startpos = p, in = 1;
    if(in && !allcovered) {
      int endpos = p-1;
      if(endpos - startpos + 1 >= minlen) printf("%12d%12d%10d\n", startpos, endpos, endpos-startpos+1);
      in = 0;
      }
    }
  }


//for parsing mpilup output:
struct parsebase_state {
  int nfiles, p;
  vector<int> poss;
  char **reads1, **quals1;
  vector<FILE*> files;
  };
void parsebase_init(parsebase_state *s, int nfiles, vector<FILE*> files) {
  s->p = 0;
  s->poss.resize(nfiles, 0); //last read position
  s->nfiles = nfiles;
  s->reads1 = new char *[nfiles];
  rep(f, nfiles) s->reads1[f] = new char [maxreads];
  s->quals1 = new char *[nfiles];
  rep(f, nfiles) s->quals1[f] = new char [maxreads];
  s->files = files;
  }
void parsebase_clear(parsebase_state *s) {
  rep(f, s->nfiles) delete[] s->reads1[f];
  delete[] s->reads1;
  rep(f, s->nfiles) delete[] s->quals1[f];
  delete[] s->quals1;
  }
bool parsebase_get(parsebase_state *s,
                   int atgc[][4], //counts of A/T/G/C
                   int atgc_r[][4], //counts of A/T/G/C only in reverse strand
                   int reads[][maxreads],
                   int quals[][maxreads],
                   int nreads[],
                   int *pos //position read
                   ) { //returns 1 if successful
  s->p++;
  rep(f, s->nfiles) {
    if(s->poss[f] < s->p) {
      if(fscanf(s->files[f], "%*s %d %*s %d", &s->poss[f], nreads+f) != 2)
        s->poss[f] = INT_MAX;
      else {
	if(nreads[f] > 0) {
	  if(fscanf(s->files[f], "%99999s %99999s%*[^\n]%*c", s->reads1[f], s->quals1[f]) != 2) {
	    s->poss[f] = INT_MAX;
	    assert(0); //remove later (it shouldn't be possible to end up here)
	    }
	  }
	else s->reads1[f][0] = 0, s->quals1[f][0] = 0, fscanf(s->files[f], "%*[^\n]%*c");
	}
      }
    }
  s->p = INT_MAX;
  //check which position to focus on in this iteration of the outer loop
  rep(f, s->nfiles) if(s->poss[f] < s->p) s->p = s->poss[f];
  if(s->p == INT_MAX) return 0;

  *pos = s->p;

  rep(f, s->nfiles) {
    assert(s->p <= s->poss[f]); //remove later
    rep(j,4) atgc[f][j] = atgc_r[f][j] = 0;
    int i = 0;
    int slen = strlen(s->reads1[f]), num_asterisk = 0;
    rep(k, slen) { //convert mpileup output into numerical data
      int indel_len;
      switch(s->reads1[f][k]) {
	case 'a': atgc_r[f][0]++;
	case 'A': reads[f][i] = 0, atgc[f][0]++, quals[f][i] = s->quals1[f][i]-33, i++;
	break;
	case 't': atgc_r[f][1]++;
	case 'T': reads[f][i] = 1, atgc[f][1]++, quals[f][i] = s->quals1[f][i]-33, i++;
	break;
	case 'g': atgc_r[f][2]++;
	case 'G': reads[f][i] = 2, atgc[f][2]++, quals[f][i] = s->quals1[f][i]-33, i++;
	break;
	case 'c': atgc_r[f][3]++;
	case 'C': reads[f][i] = 3, atgc[f][3]++, quals[f][i] = s->quals1[f][i]-33, i++;
	break;
	case '^': k++; //ignore mapping quality
	break;
	case '+': case '-': indel_len = atoi(s->reads1[f]+k+1), k++; //we ignore indels
	  while(isdigit(s->reads1[f][k])) k++;
	  k += indel_len - 1; //the -1 is because the for loop will increase it by 1
	break; 
	case '$': break;
	case '*': nreads[f]--, num_asterisk++;
	break;
	default: printf("WARNING: strange character \'%c\': p: %d, f: %d\n", s->reads1[f][k], s->p, f), exit(1);
	}
      }

    if(nreads[f] != i) printf("WARNING: nreads[f] != i: p: %d, f: %d\n", s->p, f), exit(1);
    if(s->reads1[f][slen] != 0 || s->quals1[f][i+num_asterisk] != 0) //DELETE THIS LATER
      printf("reads1[f][nreads+1]!= 0||quals1[f][i]!=0, pos: %d, file: %d\n", s->p, f);
    }
  return 1;
  }

//a normal distribution is represented by a double[2], first mean, then variance
int ncmp(const void *a, const void *b) { //compare means of normals's (for sorting, in ascending order, that's why things are upside-down, see man qsort)
  const double *aa = (const double *) a, *bb = (const double *) b;
  if(aa[0] > bb[0]) return -1;
  return aa[0] < bb[0];
  }
double dist(double a[2], double b[2]) {return fabs(a[0]-b[0]) + 0/*fabs(a[1]-b[1])*/;} //CHANGE THIS
double dist3(double a[3][2], double b[3][2]) {double d = 0; rep(i,3) d += dist(a[i], b[i]); return d;}

void dostuff(int argc, char **argv) {
  if(argc != 1) printf("Usage: dostuff <args to mpileup>\n"), exit(1);
  int nfiles;
  vector <FILE*> files;
  init(string("mpileup ")+argv[0], fname, region, &nfiles, &files);

  //things related to clustering:
  //reference: http://www.cs.princeton.edu/courses/archive/fall08/cos436/Duda/C/sk_means.htm
  const int maxclusters = 300;
  const double thres = 0.1;
  double means[nfiles][maxclusters][3][2]; //store the three nucleotides that we have the most of, sorted by how much we have
  int counts[nfiles][maxclusters];
  int nclusters[nfiles];
  mset(nclusters, 0);

  //progress printing:
  int every_n_pos = 100000;
  int next_pos_to_print = 100000;

  parsebase_state s;
  parsebase_init(&s, nfiles, files);
  int atgc[nfiles][4];
  int atgc_r[nfiles][4];
  int reads[nfiles][maxreads];
  int quals[nfiles][maxreads];
  int nreads[nfiles];
  int pos;
  while(parsebase_get(&s, atgc, atgc_r, reads, quals, nreads, &pos)) {
    rep(f, nfiles) {
      //approximate variance for atgc
      double natgc[4][2]; //approximate prior normal distribution of amount of dna
      rep(i, 4) natgc[i][0] = natgc[i][1] = atgc[f][i];
      rep(i, nreads[f]) rep(j, 4) {
	double errprob = phred2prob[quals[f][i]];
	natgc[j][1] += reads[f][i] == j ? errprob : errprob/2; //heuristic
	}

      //sort normal dists by mean
      qsort(natgc, 4, 2*sizeof(double), ncmp);

      //might be better to cluster the logarighms of the parameters
      rep(i, 3) rep(j, 2) natgc[i][j] = log(natgc[i][j]+10);

      //do clustering
      if(!nclusters[f]) {
        memcpy(means[f][0], natgc, 3*2*sizeof(double));
        counts[f][0] = 1;
        nclusters[f]++;
        }
      else {
	int closest = 0;
	double mindist = 1e100;
	rep(i, nclusters[f]) {
          double d = dist3(means[f][i], natgc);
          if(d < mindist) closest = i, mindist = d;
          }
        if(mindist > thres) {
	  memcpy(means[f][nclusters[f]], natgc, 3*2*sizeof(double));
	  counts[f][nclusters[f]] = 1;
	  nclusters[f]++;
          }
        else {
          double c = ++counts[f][closest];
          rep(i,3) rep(j,2) means[f][closest][i][j] = (1-1/c)*means[f][closest][i][j] + 1/c*natgc[i][j];
          }
        }
      if(nclusters[f] == maxclusters) {
        rep(c, nclusters[f]) {
          printf("%03d  %6d", f, counts[f][c]);
          rep(i, 3) rep(j, 2) printf("  %9.3f", exp(means[f][c][i][j])-10);
          printf("\n");
          }
        nclusters[f] = 0;
        }
      }
    //print progress:
    if(pos > next_pos_to_print) {
      next_pos_to_print += every_n_pos;
      fprintf(stderr, "pos: %d\n", pos);
      rep(f, nfiles) fprintf(stderr, "%03d, %6d clusters\n", f, nclusters[f]);
      fprintf(stderr, "\n");
      }
    }
  rep(f, nfiles) fprintf(stderr, "%03d, %8d clusters\n", f, nclusters[f]);
  parsebase_clear(&s);
  }

void compare(int argc, char **argv) {
  if(argc != 1) printf("Usage: compare <args to mpileup>\n"
"This tool lets you compare the counts of the nucleotide types for ranges of positions, over several samples.\n"
"If you don't want to give any extra arguments to mpileup, write '\"\"'\n"), exit(1);
  int nfiles;
  vector <FILE*> files;
  init(string("mpileup ")+argv[0], fname, region, &nfiles, &files);
  parsebase_state s;
  parsebase_init(&s, nfiles, files);
  int atgc[nfiles][4], atgc_r[nfiles][4], reads[nfiles][maxreads], quals[nfiles][maxreads], nreads[nfiles], pos;
  while(parsebase_get(&s, atgc, atgc_r, reads, quals, nreads, &pos)) {
    printf("Position %d:\n", pos);
    rep(f, nfiles) {
      printf("%03d ", f);
      rep(i, 4)
        printf("\t%d", atgc[f][i]);
      printf("\n");
      }
    }
  parsebase_clear(&s);
  }

const double cutoff = 0.99, //only print if the probability of the most likely genotype is less than this or if genotype is heterozygous
             a = 0.01; //prior prob. of heterozygosity

void genotype(int, char**) {
  error(1, 0, "genotype is temporarilty unavailable, because of rewriting the code.");
  }
//some compilers complain if you store this big things on the stack
//char reads1[maxfiles][maxreads];
//char quals1[maxfiles][maxreads];
//int atgc[maxfiles][4]; //TODO: think about weighting these based on read quality
//int atgc_r[maxfiles][4]; //reverse read counts
/*
void genotype2(int argc, char **argv) {
  if(argc != 3) printf("Usage: genotype <cutoff> <mindepth> \"<arguments to samtools>\"\n"
"Positions where at least one sample has a pobability less than <cutoff> for all of AA, TT, GG or CC will be shown.\n"
"Only samples where a postion is covered by at least <mindepth> reads will be used to determine if that position is interestin.\n"
"If you don't want any to send any extra arguments to samtools, write \"\".\n"
"The output is under development, please see the source code.\n"), exit(1);
  init(string("mpileup ")+argv[2]);
  double P_AB_X[maxfiles][10]; //P(AB|X)
  precompute();
  double cutoff = atof(argv[0]);
  int mindepth = atoi(argv[1]);
  int poss[maxfiles]; //last read position
  mset(poss, 0);
  bool end = 0; //remove?
  int p = 0;
  int nreads[maxfiles];
  while(!end) {
    bool interesting = 0;
    p++;
    rep(f, nfile) {
      if(poss[f] < p) {
        if(fscanf(files[f], "%*s %d %*s %d", poss+f, nreads+f) != 2)
          poss[f] = INT_MAX;
	else {
          if(nreads > 0) {
	    if(fscanf(files[f], "%99999s %99999s%*[^\n]%*c", reads1[f], quals1[f]) != 2) {
	      poss[f] = INT_MAX;
	      assert(0); //remove later (it shouldn't be possible to end up here)
	      }
	    }
	  else reads1[f][0] = 0, quals1[f][0] = 0, fscanf(files[f], "%*[^\n]%*c");
	  }
        }
      }
    p = INT_MAX;
    //check which position to focus on in this iteration of the outer loop
    rep(f, nfile) if(poss[f] < p) p = poss[f];
    if(p == INT_MAX) break;

    rep(f, nfile) {
      assert(p <= poss[f]); //remove later
      rep(j,4) atgc[f][j] = atgc_r[f][j] = 0;
      if(p < poss[f]) nreads[f] = 0;
      int i = 0;
      int reads[100000];
      int quals[100000];
      int slen = strlen(reads1[f]), num_asterisk = 0;
      rep(k, slen) { //convert mpileup output into numerical data
        int indel_len;
	switch(reads1[f][k]) {
	  case 'a': atgc_r[f][0]++;
          case 'A': reads[i] = 0, atgc[f][0]++, quals[i] = quals1[f][i]-33, i++;
	  break;
	  case 't': atgc_r[f][1]++;
          case 'T': reads[i] = 1, atgc[f][1]++, quals[i] = quals1[f][i]-33, i++;
	  break;
	  case 'g': atgc_r[f][2]++;
          case 'G': reads[i] = 2, atgc[f][2]++, quals[i] = quals1[f][i]-33, i++;
	  break;
	  case 'c': atgc_r[f][3]++;
          case 'C': reads[i] = 3, atgc[f][3]++, quals[i] = quals1[f][i]-33, i++;
	  break;
	  case '^': k++; //ignore mapping quality
	  break;
	  case '+': case '-': indel_len = atoi(reads1[f]+k+1), k++; //we ignore indels
            while(isdigit(reads1[f][k])) k++;
            k += indel_len - 1; //the -1 is because the for loop will increase it by 1
	  break; 
          case '$': break;
          case '*': nreads[f]--, num_asterisk++;
          break;
          default: printf("WARNING: strange character \'%c\': p: %d, f: %d\n", reads1[f][k], p, f), exit(1);
	  }
	}

      if(nreads[f] != i) printf("WARNING: nreads[f] != i: p: %d, f: %d\n", p, f), exit(1);
      if(reads1[f][slen] != 0 || quals1[f][i+num_asterisk] != 0) //DELETE THIS LATER
        printf("reads1[f][nreads+1]!= 0||quals1[f][i]!=0, pos: %d, file: %d\n", p, f);

      //calculate genotype probabilities
      //order: AA TT GG CC AT AG AC TG TC GC,  ATGC = 0123
      int c[10][2] = {{0,0},{1,1},{2,2},{3,3},{0,1},{0,2},{0,3},{1,2},{1,3},{2,3}};
      double l_P_AB[10]; //log(P(AB)), i.e prior prob. of the genotypes
      double l_P_Xi_AB[3][255]; //log(P(Xi=xi|AB))
      if(s->p < s->poss[f] || !nreads[f]) {
	rep(i, 10) P_AB_X[f][i] = 0;
	continue;
	}
      static int inited = 0;
      if(!inited) {
        rep(i,10) l_P_AB[i] = c[i][0]==c[i][1] ? log((1-a)/4) : log(a/6); //"a" is prior prob of heterozygosity
        rep(j,255) l_P_Xi_AB[0][j] = log(phred2prob[j]/3);
        rep(j,255) l_P_Xi_AB[1][j] = log(.5*(1-phred2prob[j]) + .5*phred2prob[j]/3);
        rep(j,255) l_P_Xi_AB[2][j] = log(1-phred2prob[j]);
        inited = 1;
        }
      double l_P_AB_X[10];
      rep(i, 10) {
	double l_P_X_AB = 0; //P(X|AB)
	rep(k, nreads[f]) {
	//P(Xi=xi|Yi=A) = P(Yi=A|Xi=xi)*P(Xi=xi)/P(Yi=A) = P(Yi=A|Xi=xi)*.25/.25 = P(Yi=A|Xi=xi)
	//P(Xi=...) ........... (write down how to derive this stuff)
          //this is the old way of computing this:
	  //P_Xi_AB = .5*(reads[k]==c[i][0]?1-quals[k]:quals[k]/3) + .5*(reads[k]==c[i][1]?1-quals[k]:quals[k]/3);
	  //l_P_X_AB += log(P_Xi_AB);
          //this is the new way, using precomputed values:
          int nsame = (reads[k]==c[i][0]) + (reads[k]==c[i][1]);
          l_P_X_AB += l_P_Xi_AB[nsame][quals[k]];
	  }
	l_P_AB_X[i] = l_P_X_AB + l_P_AB[i];
	}
      double m = -1e100;
      rep(i,10) if(l_P_AB_X[i] > m) m = l_P_AB_X[i];
      rep(i,10) P_AB_X[f][i] = exp(l_P_AB_X[i] - m);
      double denom = 0;
      rep(i,10) denom += P_AB_X[f][i];
      rep(i,10) P_AB_X[f][i] /= denom;

      //check if interesting:
      m = 0;
      rep(i, 4) m = P_AB_X[f][i] > m ? P_AB_X[f][i] : m;
      if(m < cutoff && nreads[f] >= mindepth) interesting = 1;
      }

    if(interesting) {
      static bool first = 1; //(this way of doing it depends on this function only being called once)
      if(first) first = 0;
      else printf("\n");
      printf("Position %d:\n", p);
      rep(f, nfile) {
        printf("%03d", f);
        rep(i, 10) printf(" %7.5f", P_AB_X[f][i]);
        rep(i, 4) printf("\t%d", atgc[f][i]);
        printf("\t");
        rep(i, 4) printf("\t%d", atgc_r[f][i]);
        printf("\n");
        }
      }
    }
  }
*/

void usage(char **argv) {
  printf("Usage:\n%s <filename> <region> <command> [command options].\n"
	 "Available commands: compare, findcover, genotype, dostuff.\n"
         "<region> can be something like \"chr1\" to get one whole region, or something like \"chr1:123-234\" to get a certain range of positions within this region.\n"
	 "The first line of the file with name filename contains the path of the samtools command.\n"
	 "The second line contains an integer n specifying the number of files to be opened.\n"
	 "Each of the following n lines contains the path of a .bam-file.\n"
	 "It is ok if there are empty lines *between* the filenames, but there can't be spaces *inside* the filnames.\n", argv[0]);
  }

int main(int argc, char **argv) {
  if(argc < 4)
    usage(argv), exit(1);
  fname = argv[1];
  region = argv[2];
  const char *command = argv[3];
  argc -= 4, argv += 4;
  precompute();
  if(!strcmp(command, "compare"))
    compare(argc, argv);
  else if(!strcmp(command, "findcover"))
    findcover(argc, argv);
  else if(!strcmp(command, "genotype"))
    genotype(argc, argv);
  else if(!strcmp(command, "dostuff"))
    dostuff(argc, argv);
  else usage(argv); 
  }
