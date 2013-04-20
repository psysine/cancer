#include <algorithm>
#include <string>
#include <vector>
#include <map>
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
#define repu(i,n) for(unsigned i = 0; i < (n); i++)
#define rap(i,n,m) for(int i = (n); i < (m); i++)
#define rapu(i,n,m) for(unsigned i = (n); i < (m); i++)
#define iter(it, v) for(typeof((v).begin()) it = (v).begin(); it != (v).end(); ++it)
#define MP(x,y) make_pair(x,y)
#define F first
#define S second
#define all(v) (v).begin(), (v).end() //ex: sort(all(v));
#define mset(v, byte) memset(v, byte, sizeof(v))
typedef unsigned int uint;
typedef unsigned long long ull;
typedef long long ll;
typedef pair<int,int> pii;
typedef pair<double,double> pdd;

template<typename T> T sum(T *a, int n) { T sum = 0; rep(i,n) sum += a[i]; return sum; }
template<typename T> T max(T *a, int n) { T m = a[0]; rap(i,1,n) if(a[i] > m) m = a[i]; return m; }
template<typename T> T min(T *a, int n) { T m = a[0]; rap(i,1,n) if(a[i] < m) m = a[i]; return m; }
template<typename T> T sq(T a) { return a*a; }

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
  if(argc != 3) printf("Usage: findcover <mindepth> <minlength> <arguments to \'samtools depth\'>\n"
"The output will be a list of regions that are at least <minlength> long and in which each base pair is covered by at least <mindepth> reads.\n"
"The first column is the starting posision of each region, the second column is the ending position and the third column is the length of the region.\n"
"If you don't want to give any extra arguments to mpileup, write '\"\"'\n"), exit(1);
  int mindepth = atoi(argv[0]), minlen = atoi(argv[1]);
  if(mindepth <= 0) error(1, 0, "mindepth must be greater than zero");
  int nfiles;
  vector <FILE*> files;
  init(string("depth ")+argv[2], fname, region, &nfiles, &files);

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


const int nchrs = 23;
const char *chrs[nchrs] = {"chrM", "chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9",
                           "chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17",
                           "chr18","chr19","chr20","chr21","chr22"};
int lastchr = 0; //for performance, to avoid calling strcmp too many times
int chrtoi(const char *c) {
  if(!strcmp(c, chrs[lastchr])) return lastchr;
  rep(i, nchrs) if(!strcmp(c, chrs[i])) return lastchr = i;
  return -1;
  }

void readranges(FILE *fi, vector<pair<int, int> > vs[nchrs]) {
  char c[10000];
  int s, e;
  while(fscanf(fi, "%9999s %d %d", c, &s, &e) == 3) {
    int chrn = chrtoi(c);
    if(chrn < 0) continue;
    vs[chrn].push_back(make_pair(s, e));
    }
  if(fclose(fi)) error(1, errno, "couldn't close protein regions file");
  rep(i, nchrs) {
    sort(vs[i].begin(), vs[i].end());
    }
  }

//isinprotein_* are used by parsebase_*
bool use_isinprotein;
vector<pair<int,int> > pchrs[nchrs];
void merge_overlaps(vector<pair<int,int> > v) {
  unsigned i = 0;
  while(i < v.size()-1)
    if(v[i].S >= v[i+1].F-1) {
      v[i].S = max(v[i].S, v[i+1].S);
      v.erase(v.begin()+i+1);
      }
    else i++;
  }
void isinprotein_init(string fname) {
  FILE *fi = fopen(fname.c_str(), "r");
  if(fi) {
    fprintf(stderr, "Found list of protein coding regions: %s, will filter based on this.\n", fname.c_str());
    use_isinprotein = 1;
    }
  else {
    fprintf(stderr, "Did NOT find list of protein coding regions: %s, will NOT filter based on this.\n", fname.c_str());
    use_isinprotein = 0;
    return;
    }
  readranges(fi, pchrs);
  rep(i, nchrs) merge_overlaps(pchrs[i]);
  }
bool isinprotein(int chr, int pos) { //assumes pos to be increasing from call to call, unless chr changes. may only be called for chr \in [0,nchrs)
  static int c = -1, i = -1;
  if(c != chr) c = chr, i = 0;
  while(i < int(pchrs[chr].size()) && pos > pchrs[chr][i].S) i++;
  return i < int(pchrs[chr].size()) && pos >= pchrs[chr][i].F;
  }


//###########################//
//##### PARSEBASE STUFF #####//
//###########################//

//the parsebase_*-functions are for parsing mpilup output
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
  isinprotein_init("proteinranges");
  }
void parsebase_clear(parsebase_state *s) {
  rep(f, s->nfiles) delete[] s->reads1[f];
  delete[] s->reads1;
  rep(f, s->nfiles) delete[] s->quals1[f];
  delete[] s->quals1;
  }
void parsebase_parse(const char *reads1, const char *quals1, int nreads, int *atgc, int *atgc_r, int *reads, int *quals, int *nstartend, int *nreads2) {
  rep(j,4) atgc[j] = atgc_r[j] = 0;
  rep(j,2) nstartend[j] = 0;
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
      case '^': nstartend[0]++, k++; //ignore mapping quality
      break;
      case '+': case '-': indel_len = atoi(reads1+k+1), k++; //we ignore indels
	while(isdigit(reads1[k])) k++;
	k += indel_len - 1; //the -1 is because the for loop will increase it by 1
      break; 
      case '$': nstartend[1]++;
      break;
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
                   int nstartend[][2], //number of reads starting/ending at postition (i.e., number of ^'s and $'s in mpileup output)
                   int *chr, //chromosome read
                   int *pos //position read
                   ) { //returns 1 if successful
  rep(f, s->nfiles)
    if(s->chrs[f] < s->c || (s->chrs[f] == s->c && s->poss[f] < s->p))
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
  s->p++; //for next time

  if(use_isinprotein && !isinprotein(*chr, *pos))
    return parsebase_get(s, atgc, atgc_r, reads, quals, nreads, nstartend, chr, pos);

  rep(f, s->nfiles) {
    assert(*pos <= s->poss[f] && *chr <= s->chrs[f]);
    if(s->chrs[f] == *chr && s->poss[f] == *pos)
      parsebase_parse(s->reads1[f], s->quals1[f], s->nreads[f], atgc[f], atgc_r[f], reads[f], quals[f], nstartend[f], nreads+f);
    else {
      rep(i,4) atgc[f][i] = atgc_r[f][i] = 0;
      nreads[f] = 0;
      }
    assert(sum(atgc[f], 4) == nreads[f]);
    }
  return 1;
  }

//###########################//
//####### ABG STUFF #########//
//###########################//
const double max_ratio = 1.1; //maximum allowed (biggest belonging to bin)/(smallest belonging to bin)
int dims[5];
const int maxc = 1000000;
int count2bin[maxc];
void count2bin_init() {
  count2bin[0] = 0;
  count2bin[1] = 1;
  int current_bin = 1;
  int lastbinmin = 1;
  for(int i = 2; i < maxc; i++)  {
    if(double(i)/lastbinmin > max_ratio)
      current_bin++, lastbinmin = i;
    count2bin[i] = current_bin;
    }
  }
struct posbin {
  int b[5];
  };
bool posbincmp(const posbin &l, const posbin &r) {
  rep(i, 5) 
    if(l.b[i] < r.b[i]) return 1;
    else if(l.b[i] > r.b[i]) return 0;
  return 0;
  }
struct poscluster {
  poscluster() {rep(i,5) sums[i] = 0; count = 0;}
  double sums[5];
  int count;
  };

typedef map<posbin, poscluster, bool(*)(const posbin &,const posbin &)> posmap_t;
posmap_t **clusterings;

void abg(int argc, char **argv) {
  if(argc != 1) fprintf(stderr, "Usage: abg <args to mpileup>\n"
"This function does some yet undefined stuff.\n"
"If you don't want to give any extra arguments to mpileup, write '\"\"'\n"), exit(1);
  int nfiles;
  vector <FILE*> files;
  init(string("mpileup ")+argv[0], fname, region, &nfiles, &files);
  count2bin_init();
  clusterings = new posmap_t*[nfiles];
  rep(f, nfiles) {
    clusterings[f] = new posmap_t(posbincmp);
    }
  parsebase_state s;
  parsebase_init(&s, nfiles, files);
  int atgc[nfiles][4], atgc_r[nfiles][4], reads[nfiles][maxreads], quals[nfiles][maxreads], nreads[nfiles], nstartend[nfiles][2], chr, pos;
  while(parsebase_get(&s, atgc, atgc_r, reads, quals, nreads, nstartend, &chr, &pos)) {
    if(! (pos % 1000000)) fprintf(stderr, "chr: %d, pos: %d\n", chr, pos);
    rep(f, nfiles) {
      double avgqual = 0;
      rep(r, nreads[f]) avgqual += phred2prob[quals[f][r]];
      avgqual /= nreads[f];
      int qual = prob2phred(avgqual);
      rep(i, 4) for(int j = i+1; j < 4; j++) if(atgc[f][j] > atgc[f][i]) swap(atgc[f][j], atgc[f][i]);
//      assert(atgc[f][0] >= atgc[f][1] && atgc[f][1] >= atgc[f][2] && atgc[f][2] >= atgc[f][3]);
//      assert(atgc[f][0] < maxc);
      if(!atgc[f][0]) continue; //we're not interested in positions with no reads
      posbin pb;
      rep(i, 4) pb.b[i] = count2bin[atgc[f][i]];
      pb.b[4] = qual;
      poscluster &pc = (*clusterings[f])[pb];
      rep(i, 4) pc.sums[i] += atgc[f][i];
      pc.sums[4] += avgqual;
      pc.count++;
      }
    }
  rep(f, nfiles) {
    int mins[5], maxs[5];
    rep(d, 5) mins[d] = INT_MAX, maxs[d] = INT_MIN;
    iter(it, *clusterings[f]) {
      const posbin &pb = it->first;
      const poscluster &pc = it->second;
      rep(d, 5) {
	if(pb.b[d] < mins[d]) mins[d] = pb.b[d];
	if(pb.b[d] > maxs[d]) maxs[d] = pb.b[d];
	}
      if(nfiles > 1) printf("%5d ", f);
      printf("%12d", pc.count);
      rep(j, 4) printf(" %8.2f", pc.sums[j]/pc.count);
      printf(" %13.8f", pc.sums[4]/pc.count);
      printf("\n");
      }
    fprintf(stderr, "file: %d\n", f);
    rep(d, 5)
      fprintf(stderr, "dim %d: min=%d, max=%d\n", d, mins[d], maxs[d]);
    fprintf(stderr, "total nonzero: %d\n", int(clusterings[f]->size()));
    } 
  parsebase_clear(&s);
  }



//###########################//
//###### PEARSON STUFF ######//
//###########################//

double pearsonstat_mx2(int m, int N, int a[][2], int *rowsums, int *colsums) {
  double s = 0;
  rep(i, m) rep(j, 2) s += double(sq(a[i][j]))/rowsums[i]/colsums[j];
  return N*(s-1);
  }
double pearsonmean_mx2(int m, int N) {
  return (m-1.)*N/(N-1);
  }
double pearsonvar_mx2(int m, int N, int *rowsums, int *colsums) {
  int Nsq = sq(N), msq = sq(m);
  double T = double(N)/colsums[0]/colsums[1];
  double S = 0;
  rep(i, m) S += 1./rowsums[i];
  double NS = N*S;
  double var = double(Nsq)/sq(N-1)/(N-2)/(N-3) * (2*(m-1)*Nsq + 2*(2*msq+m)*N - 6*msq - 6*(N-1)*NS + (N-1)*T*((N+1)*NS - (msq + 2*(m-1))*N + msq - 2*m));
  if(var < -1e-12) fprintf(stderr, "m %d, N %d, Nsq %d, T %g, S %g\n", m, N, Nsq, T, S);
  return var;
  }
double effectsize_est(int N) {
  return N;
  }
pdd test_atgc(int *atgc_no, int *atgc_ca) {
  bool posind[4];
  rep(i, 4) posind[i] = atgc_no[i] + atgc_ca[i] > 0;
  int atgc[4][2];
  int nrows = 0;
  const int ncols = 2;
  rep(i, 4)
    if(posind[i])
      atgc[nrows][0] = atgc_no[i], atgc[nrows][1] = atgc_ca[i], nrows++;
  if(nrows < 2) return pdd(0, 0);
  int rowsums[nrows], colsums[ncols] = {};
  rep(i, nrows) rowsums[i] = sum(atgc[i], ncols);
  rep(i, ncols) rep(j, nrows) colsums[i] += atgc[j][i];
  int N = sum(colsums, ncols);
  if(N < 4 || min(colsums, 2) == 0) return pdd(0, 0);
  double chisq = pearsonstat_mx2(nrows, N, atgc, rowsums, colsums);
  double mean = pearsonmean_mx2(nrows, N);
  double var = pearsonvar_mx2(nrows, N, rowsums, colsums);
  if(var <= -1e-12) {
    fprintf(stderr, "\nvar: %g\n", var);
    rep(i, nrows) fprintf(stderr, "%3d %3d | %3d\n", atgc[i][0], atgc[i][1], rowsums[i]);
    fprintf(stderr, "%3d %3d\n", colsums[0], colsums[1]);
    }
  double centralized = chisq-mean;
  double effectsize = effectsize_est(N);
  return pdd(centralized*effectsize, var*sq(effectsize));
  }


//###########################//
//###### PERGENE STUFF ######//
//###########################//

struct pergenedata { //one of these for each gene and each normal/cancer pair
  pergenedata() : sum(), var(), nstartend_no(), nstartend_ca() {}
  double sum, var;
  long long nstartend_no, nstartend_ca;
  };
struct pergenestuff {
  pergenestuff(pii ra, int npairs) : range(ra), touched(), pgd(npairs) {}
  pii range;
  bool touched;
  vector<pergenedata> pgd;
  };
struct genestore_state { int lastchr, lastpos, gb, npairs; } gss;
void printheader() {
  printf("pair_number\tchromosome\tblock\tgene_start\t"
         "gene_end\tchisq_sum\tchisq_stddev\tchisq_z\t"
         "count_no\tcount_ca\tcount_ratio\n");
  }
void printgene(const pergenestuff &pgs, int chr, int gb) {
  rep(p, gss.npairs) {
    const pergenedata &gd = pgs.pgd[p];
    printf("%02d %5s %8d %9d %9d %15g %15g %15g %15lld %15lld %15g\n",
	   p, chrs[chr], gb, pgs.range.F,
           pgs.range.S, gd.sum, sqrt(gd.var), gd.sum/sqrt(gd.var),
           gd.nstartend_no, gd.nstartend_ca, double(gd.nstartend_no)/gd.nstartend_ca);
    }
  }
void printgeneblock(const vector<pergenestuff> &pgs, int chr, int gb) {
  repu(g, pgs.size()) {
    if(pgs[g].touched)
      printgene(pgs[g], chr, gb);
    }
  }
vector<pair<pii, vector<pergenestuff> > > genedata[nchrs];
void genestore_init(int npairs, const char *fname) {
  FILE *fi = fopen(fname, "r");
  if(!fi) error(1, errno, "Couldn't open list of genes");
  vector<pii> generanges[nchrs];
  readranges(fi, generanges);
  rep(chr, nchrs) {
    repu(i, generanges[chr].size()) {
      if(genedata[chr].empty() || generanges[chr][i].F > genedata[chr].back().F.S) {
        genedata[chr].push_back(make_pair(generanges[chr][i], vector<pergenestuff>()));
        genedata[chr].back().S.push_back(pergenestuff(generanges[chr][i], npairs));
        }
      else { //it overlaps with previous block
        genedata[chr].back().F.S = max(genedata[chr].back().F.S, generanges[chr][i].S);
        genedata[chr].back().S.push_back(pergenestuff(generanges[chr][i], npairs));
        }
      }
    }
  gss.lastchr = -1, gss.lastpos = -1, gss.gb = 0, gss.npairs = npairs;
  printheader();
  }
int genestore_advance(int chr, int pos) {
  int &lastchr = gss.lastchr, &lastpos = gss.lastpos, &gb = gss.gb;
  if(chr < lastchr || (chr == lastchr && pos < lastpos)) error(1, 0, "genestore called backwards: chr < lastchr || (chr == lastchr && pos < lastpos)");
  if(chr > lastchr) {
    if(lastchr != -1)
      rapu(gbb, gb, genedata[lastchr].size())
        printgeneblock(genedata[lastchr][gbb].S, lastchr, gbb);
    lastchr = chr, lastpos = -1, gb = 0;
    }
  if(pos > lastpos) {
    while(uint(gb) < genedata[chr].size() && pos > genedata[chr][gb].F.S)
      printgeneblock(genedata[chr][gb].S, chr, gb), gb++;
    lastpos = pos;
    }
  return gb;
  }
void genestore(int pa, int npairs, int chr, int pos, int *atgc_no, int *atgc_ca, int nstartend_no, int nstartend_ca, int nreads_no, int nreads_ca) {
  int gb = genestore_advance(chr, pos);
  if(uint(gb) == genedata[chr].size() || pos < genedata[chr][gb].F.F) return;
assert(pos >= genedata[chr][gb].F.F && pos <= genedata[chr][gb].F.S);
  vector<pergenestuff> &vpgs = genedata[chr][gb].S;
  pdd st = pdd(0, 0);
  if(nreads_no >= 20 && nreads_ca >= 20)
    st = test_atgc(atgc_no, atgc_ca);
  repu(i, vpgs.size())
    if(pos >= vpgs[i].range.F && pos <= vpgs[i].range.S) {
      vpgs[i].touched = 1;
      pergenedata &p = vpgs[i].pgd[pa];
      p.sum += st.F, p.var += st.S;
      p.nstartend_no += nstartend_no, p.nstartend_ca += nstartend_ca;
      }
  }
void genestore_finish() {
  printgeneblock(genedata[lastchr][gss.gb].S, gss.lastchr, gss.gb);
  }

void pergene(int argc, char **argv) {
  if(argc != 1) printf("Usage: pergene <args to mpileup>\n"
"This tool compares normal and cancer samples, gene per gene.\n"
"The bam-files listed in <filename> should be even in number,\n"
"alternating between normal and cancer samples.\n"
"If you don't want to give any extra arguments to mpileup, write '\"\"'\n"), exit(1);
  int nfiles;
  vector <FILE*> files;
  init(string("mpileup ")+argv[0], fname, region, &nfiles, &files);
  if(nfiles % 2) error(1, 0, "must specify even number of bam-files.");
  int npairs = nfiles / 2;
  parsebase_state s;
  parsebase_init(&s, nfiles, files);
  genestore_init(npairs, "proteinranges");
  int atgc[nfiles][4], atgc_r[nfiles][4], reads[nfiles][maxreads], quals[nfiles][maxreads], nreads[nfiles], nstartend[nfiles][2], chr, pos;
  while(parsebase_get(&s, atgc, atgc_r, reads, quals, nreads, nstartend, &chr, &pos))
    rep(pa, npairs) {
      int no = 2*pa, ca = 2*pa+1;
      genestore(pa, npairs, chr, pos, atgc[no], atgc[ca], sum(nstartend[no], 2), sum(nstartend[ca], 2), nreads[no], nreads[ca]);
      }
//  rep(c, nchrs)
//    repu(gb, genedata[c].size())
//      printgeneblock(genedata[c][gb].S, gb);
  genestore_finish();
  parsebase_clear(&s);
  }

//###########################//
//### ISINTERESTING STUFF ###//
//###########################//

bool isinteresting(int pa, int npairs, int chr, int pos, int *atgc_no, int *atgc_ca, int nstartend_no, int nstartend_ca, int nreads_no, int nreads_ca) {
  pdd st = test_atgc(atgc_no, atgc_ca);
  return st.F < -1e-12 || st.F > 1e-12;
  }

void interesting(int argc, char **argv) {
  if(argc != 1) printf("Usage: interesting <args to mpileup>\n"
"This tool lets you show positions that are pairwisely interesting.\n"
"If you don't want to give any extra arguments to mpileup, write '\"\"'\n"), exit(1);
  int nfiles;
  vector <FILE*> files;
  init(string("mpileup ")+argv[0], fname, region, &nfiles, &files);
  if(nfiles % 2) error(1, 0, "must specify even number of bam-files.");
  int npairs = nfiles / 2;
  parsebase_state s;
  parsebase_init(&s, nfiles, files);
  int atgc[nfiles][4], atgc_r[nfiles][4], reads[nfiles][maxreads], quals[nfiles][maxreads], nreads[nfiles], nstartend[nfiles][2], chr, pos;
  while(parsebase_get(&s, atgc, atgc_r, reads, quals, nreads, nstartend, &chr, &pos)) {
    int ninteresting = 0;
    rep(pa, npairs) {
      int no = 2*pa, ca = 2*pa+1;
      ninteresting += isinteresting(pa, npairs, chr, pos, atgc[no], atgc[ca], sum(nstartend[no], 2), sum(nstartend[ca], 2), nreads[no], nreads[ca]);
      }
    if(ninteresting > 0) {
      printf("chr %d, pos %d:\n", chr, pos);
      rep(f, nfiles) {
	printf("%03d ", f);
	rep(i, 4)
	  printf("\t%d", atgc[f][i]);
        if(! f % 2) {
          pdd st = test_atgc(atgc[f], atgc[f+1]);
          printf(" %15g %15g", st.F, sqrt(st.S));
          }
	printf("\n");
	}
      }
    }
  parsebase_clear(&s);
  }

//###########################//
//###### COMPARE STUFF ######//
//###########################//

void compare(int argc, char **argv) {
  if(argc != 1) printf("Usage: compare <args to mpileup>\n"
"This tool lets you compare the counts of the nucleotide types for ranges of positions, over several samples.\n"
"If you don't want to give any extra arguments to mpileup, write '\"\"'\n"), exit(1);
  int nfiles;
  vector <FILE*> files;
  init(string("mpileup ")+argv[0], fname, region, &nfiles, &files);
  parsebase_state s;
  parsebase_init(&s, nfiles, files);
  int atgc[nfiles][4], atgc_r[nfiles][4], reads[nfiles][maxreads], quals[nfiles][maxreads], nreads[nfiles], nstartend[nfiles][2], chr, pos;
  while(parsebase_get(&s, atgc, atgc_r, reads, quals, nreads, nstartend, &chr, &pos)) {
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


//###########################//
//####### MAIN STUFF ########//
//###########################//

void usage(char **argv) {
  printf("Usage:\n%s <filename> <region> <command> [command options].\n"
	 "Available commands: compare, findcover, abg, pergene.\n"
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
  if(!strcmp(command, "interesting"))
    interesting(argc, argv);
  else if(!strcmp(command, "findcover"))
    findcover(argc, argv);
  else if(!strcmp(command, "abg"))
    abg(argc, argv);
  else if(!strcmp(command, "pergene"))
    pergene(argc, argv);
  else usage(argv); 
  }
