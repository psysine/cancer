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

int sumi(int *a, int n) {
  int sum = 0;
  rep(i,n) sum += a[i];
  return sum;
  }
double sumd(double *a, int n) {
  double sum = 0;
  rep(i,n) sum += a[i];
  return sum;
  }

string fname, region;
const int maxreads = 100000;

double phred2prob[255];

int prob2phred(double prob) {
  return lrint(-10*log10(prob));
  }

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
    if(region != "all")
      cmd += " -r " + region;
    FILE *file = popen((string(samtools)+" "+cmd+" \""+bamfile+'"').c_str(), "r");
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
char *chrs[22] = {"chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9",
                  "chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17",
                  "chr18","chr19","chr20","chr21","chr22"};
int lastchr = 0; //for performance, to avoid calling strcmp too many times
int chrtoi(const char *c) {
  if(!strcmp(c, chrs[lastchr])) return lastchr;
  rep(i, 22) if(!strcmp(c, chrs[i])) return lastchr = i;
  return -1;
  }

struct parsebase_state {
  int nfiles, c, p;
  vector<int> chrs, poss, nreads;
  char **reads1, **quals1;
  vector<FILE*> files;
  };
void parsebase_init(parsebase_state *s, int nfiles, vector<FILE*> files) {
  s->c = INT_MAX, s->p = INT_MAX;
  s->chrs.resize(nfiles, -1); //last read chromosome
  s->poss.resize(nfiles, 0); //last read position
  s->nreads.resize(nfiles); //alleged number of reads
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
void parsebase_parse(const char *reads1, const char *quals1, int nreads, int *atgc, int *atgc_r, int *reads, int *quals, int *nreads2) {
  rep(j,4) atgc[j] = atgc_r[j] = 0;
  int i = 0;
  int slen = strlen(reads1), num_asterisk = 0;
  rep(k, slen) { //convert mpileup output into numerical data
    int indel_len;
    switch(reads1[k]) {
      case 'a': atgc_r[0]++;
      case 'A': reads[i] = 0, atgc[0]++, quals[i] = quals1[i]-33, i++;
      break;
      case 't': atgc_r[1]++;
      case 'T': reads[i] = 1, atgc[1]++, quals[i] = quals1[i]-33, i++;
      break;
      case 'g': atgc_r[2]++;
      case 'G': reads[i] = 2, atgc[2]++, quals[i] = quals1[i]-33, i++;
      break;
      case 'c': atgc_r[3]++;
      case 'C': reads[i] = 3, atgc[3]++, quals[i] = quals1[i]-33, i++;
      break;
      case '^': k++; //ignore mapping quality
      break;
      case '+': case '-': indel_len = atoi(reads1+k+1), k++; //we ignore indels
	while(isdigit(reads1[k])) k++;
	k += indel_len - 1; //the -1 is because the for loop will increase it by 1
      break; 
      case '$': break;
      case '*': nreads--, num_asterisk++;
      break;
      default: printf("WARNING: strange character \'%c\'\n", reads1[k]), exit(1);
      }
    }
  *nreads2 = nreads;
  if(nreads != i)
    printf("WARNING: nreads != i"), *(int*)0=0, exit(1);
  if(reads1[slen] != 0 || quals1[i+num_asterisk] != 0) //DELETE THIS LATER
    printf("reads1[nreads+1]!= 0||quals1[i]!=0"), exit(1);
  }
void parsebase_read(FILE *fi, int *pos, int *chr, int *nreads, char *reads, char *quals) 
  { // read lines from mplieup input until we get one from an interesting chromosome or until there is no more input
    // interesting chromosomes are specified in char *chrs[] above
  assert(maxreads == 100000);
  char c[1000];
  *chr = -1;
  while(*chr == -1) {
    if(fscanf(fi, "%999s %d %*s %d", c, pos, nreads) != 3)
      *chr = *pos = INT_MAX, *nreads = 0;
    else {
      *chr = chrtoi(c); //somewhat inefficient
      if(*nreads > 0) {
	if(fscanf(fi, "%99999s %99999s%*[^\n]%*c", reads, quals) != 2) {
	  assert(0); // (it shouldn't be possible to end up here)
	  *chr = *pos = INT_MAX, *nreads = 0; //but this is what we'd to if it were
	  }
	}
      else reads[0] = 0, quals[0] = 0, fscanf(fi, "%*[^\n]%*c");
      }
    }
  }
bool parsebase_get(parsebase_state *s,
                   int atgc[][4], //counts of A/T/G/C
                   int atgc_r[][4], //counts of A/T/G/C only in reverse strand
                   int reads[][maxreads],
                   int quals[][maxreads],
                   int nreads[],
                   int *chr, //chromosome read
                   int *pos //position read
                   ) { //returns 1 if successful
  rep(f, s->nfiles)
    if(s->chrs[f] < s->c || s->chrs[f] == s->c && s->poss[f] < s->p)
      parsebase_read(s->files[f], &s->poss[f], &s->chrs[f], &s->nreads[f], s->reads1[f], s->quals1[f]);

//if the minimum of the chrs is c and (the minimun of the poss for which the chr is c) is p
//  do parsing for c,p, then increment p
//if the minimum of the chrs is INT_MAX
//  return 0

  //check which chromosome to focus on
  s->c = INT_MAX;
  rep(f, s->nfiles) if(s->chrs[f] < s->c) s->c = s->chrs[f];
  if(s->c == INT_MAX) return 0;

  //check which position to focus on
  s->p = INT_MAX;
  rep(f, s->nfiles) if(s->chrs[f] == s->c && s->poss[f] < s->p) s->p = s->poss[f];
  assert(s->p != INT_MAX);

  *chr = s->c, *pos = s->p;

  rep(f, s->nfiles) {
    assert(s->p <= s->poss[f] && s->c <= s->chrs[f]);
    if(s->chrs[f] == *chr && s->poss[f] == *pos)
      parsebase_parse(s->reads1[f], s->quals1[f], s->nreads[f], atgc[f], atgc_r[f], reads[f], quals[f], nreads+f);
    else {
      rep(i,4) atgc[f][i] = atgc_r[f][i] = 0;
      nreads[f] = 0;
      }
    assert(sumi(atgc[f], 4) == nreads[f]);
    }
  s->p++;
  return 1;
  }

const double max_ratio = 1.1; //maximum allowed (biggest belonging to bin)/(smallest belonging to bin)
int maxcount[4] = {15000, 200, 50, 10};
int dims[5];
const int maxc = 20000;
int count2bin[maxc];
int nquals = 60;
struct postype {
  double sums[5];
  int count;
  };
postype ***counts;
int totsize;

int dimcoefs[4];
void initdimcoefs() { //called from abg_init
  dimcoefs[3] = dims[4];
  dimcoefs[2] = dimcoefs[3]*dims[3];
  dimcoefs[1] = dimcoefs[2]*dims[2];
  dimcoefs[0] = dimcoefs[1]*dims[1];
  }
int ind(int a[5]) {
  return a[0]*dimcoefs[0] + a[1]*dimcoefs[1] + a[2]*dimcoefs[2] + a[3]*dimcoefs[3] + a[4];
  }


void abg_init(int nfiles) {
  count2bin[0] = 0;
  count2bin[1] = 1;
  int current_bin = 1;
  int lastbinmin = 1;
  for(int i = 2; i < maxc; i++)  {
    if(double(i)/lastbinmin > max_ratio)
      current_bin++, lastbinmin = i;
    count2bin[i] = current_bin;
    }
  int maxbin[4];
  for(int d = 0; d < 4; d++) {
   maxbin[d] = count2bin[maxcount[d]];
   rep(i, maxc) if(count2bin[i] == maxbin[d]) maxcount[d] = i; //adjust maxcount
   dims[d] = maxbin[d] + 1;
   }
  dims[4] = nquals;
  totsize = 1;
  rep(d, 5) totsize *= dims[d];
  rep(d, 5) fprintf(stderr, "dims[%d] = %d\n", d, dims[d]);
  const int ptrsize = sizeof(int32_t*);
  initdimcoefs();
  counts = new postype**[nfiles];
  rep(f, nfiles) {
    counts[f] = new postype*[totsize];
    memset(counts[f], 0, ptrsize*totsize);
    }
  }

void abg(int argc, char **argv) {
  if(argc != 1) fprintf(stderr, "Usage: abg <args to mpileup>\n"
"This function does some yet undefined stuff.\n"
"If you don't want to give any extra arguments to mpileup, write '\"\"'\n"), exit(1);
  int nfiles;
  vector <FILE*> files;
  init(string("mpileup ")+argv[0], fname, region, &nfiles, &files);
  abg_init(nfiles);
  parsebase_state s;
  parsebase_init(&s, nfiles, files);
  int atgc[nfiles][4], atgc_r[nfiles][4], reads[nfiles][maxreads], quals[nfiles][maxreads], nreads[nfiles], chr, pos;
  while(parsebase_get(&s, atgc, atgc_r, reads, quals, nreads, &chr, &pos)) {
    if(! (pos % 100000)) fprintf(stderr, "chr: %d, pos: %d\n", chr, pos);
    rep(f, nfiles) {
      double avgqual = 0;
      rep(r, nreads[f]) avgqual += phred2prob[quals[f][r]];
      avgqual /= nreads[f];
      int qual = prob2phred(avgqual);
      rep(i, 4) for(int j = i+1; j < 4; j++) if(atgc[f][j] > atgc[f][i]) swap(atgc[f][j], atgc[f][i]);
      assert(atgc[f][0] >= atgc[f][1] && atgc[f][1] >= atgc[f][2] && atgc[f][2] >= atgc[f][3]);
      assert(atgc[f][0] < maxc);
      int bins[5];
      rep(i, 4) bins[i] = count2bin[atgc[f][i]];
      bins[4] = qual;
      int depth = sumi(atgc[f], 4);
      if(depth >= 3 && qual >= 5) {
        if(counts[f][ind(bins)]) {
          postype *ppt = counts[f][ind(bins)];
          rep(i, 4) ppt->sums[i] += atgc[f][i];
          ppt->sums[4] += avgqual;
          ppt->count++;
          }
        else {
          postype *ppt = new postype;
          rep(i, 4) ppt->sums[i] = atgc[f][i];
          ppt->sums[4] = avgqual;
          ppt->count = 1;
          counts[f][ind(bins)] = ppt;
          }
        }
      }
    }
  rep(f, nfiles) {
    int i[5] = {0}, mins[5], maxs[5], nnonzero = 0;
    rep(d, 5) mins[d] = INT_MAX, maxs[d] = INT_MIN;
    for(i[0] = 0; i[0] < dims[0]; i[0]++)
      for(i[1] = 0; i[1] < dims[1]; i[1]++)
	for(i[2] = 0; i[2] < dims[2]; i[2]++)
	  for(i[3] = 0; i[3] < dims[3]; i[3]++)
	    for(i[4] = 0; i[4] < dims[4]; i[4]++) {
	      if(counts[f][ind(i)]) {
                nnonzero++;
                rep(d, 5) {
                  if(i[d] < mins[d]) mins[d] = i[d];
                  if(i[d] > maxs[d]) maxs[d] = i[d];
                  }
                postype *ppt = counts[f][ind(i)];
                if(nfiles > 1) printf("%5d ", f);
                printf("%10d", ppt->count);
                rep(j, 5) printf(" %10f", ppt->sums[j]/ppt->count);
                printf("\n");
                }
	      }
    fprintf(stderr, "file: %d\n", f);
    rep(d, 5)
      fprintf(stderr, "dim %d: min=%d, max=%d\n", d, mins[d], maxs[d]);
    fprintf(stderr, "total nonzero: %d\n", nnonzero);
    } 
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
  int atgc[nfiles][4], atgc_r[nfiles][4], reads[nfiles][maxreads], quals[nfiles][maxreads], nreads[nfiles], chr, pos;
  while(parsebase_get(&s, atgc, atgc_r, reads, quals, nreads, &chr, &pos)) {
    printf("chr %d, pos %d:\n", chr, pos);
    rep(f, nfiles) {
      printf("%03d ", f);
      rep(i, 4)
        printf("\t%d", atgc[f][i]);
      printf("\n");
      }
    }
  parsebase_clear(&s);
  }

void usage(char **argv) {
  printf("Usage:\n%s <filename> <region> <command> [command options].\n"
	 "Available commands: compare, findcover, abg.\n"
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
  else if(!strcmp(command, "abg"))
    abg(argc, argv);
  else usage(argv); 
  }
