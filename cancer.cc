//#include <iostream>
#include <string>
//#include <vector>
//#include <stack>
//#include <queue>
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

const int maxfiles = 99;

const char *fname, *region;
char reads[maxfiles][100000];
FILE *files[maxfiles];
int nfile;
bool readfail[maxfiles];

void init(const char *cmd) {
  FILE *fi = fopen(fname, "r");
  if(!fi) error(1, errno, "Can't open list of bam files");
  char *samtools;
  fscanf(fi, "%a[^\n]%*c", &samtools);
  fscanf(fi, "%d", &nfile);
  if(nfile < 1) error(1, 0, "Need at least one file");
  if(nfile > maxfiles) error(1, 0, "Can't open more than %d bam files", maxfiles);

  rep(f, nfile) {
    char *bamfile = 0;
    if(fscanf(fi, "%as", &bamfile) != 1)
      error(1, 0, "Not enough filenames");
    files[f] = popen((string(samtools)+" "+cmd+" -r "+region+" \""+bamfile+'"').c_str(), "r");
    if(files[f] <= 0)
      error(1, errno, "Can't run samtools");
    free(bamfile);
    }
  free(samtools), samtools = 0;
  fclose(fi);
  }

void compare(int argc, char **argv) {
  if(argc != 2) error(1, 0, "Usage: compare <startpos> <endpos>");
  init("mpileup");
  int startpos = atoi(argv[0]), endpos = atoi(argv[1]);
  if(startpos <= 0) error(1, 0, "Starting position must be greater than zero");
  int poss[maxfiles]; //last read position
  mset(poss, 0);
  for(int p = startpos; p <= endpos; p++) {
    printf("Position %d\n", p);
    rep(f, nfile) {
      while(poss[f] < p)
	if(fscanf(files[f], "%*s %d %*s %*d %99999s%*[^\n]%*c", poss+f, reads[f]) != 2)
          poss[f] = INT_MAX;

      if(p < poss[f])
        printf("%02d\t0\t0\t0\t0\n", f);
      else {
	int slen = strlen(reads[f]);
	int A = 0, T = 0, G = 0, C = 0;
	rep(k, slen) {
	  switch(reads[f][k]) {
	    case 'a': case 'A': A++;
	    break;
	    case 't': case 'T': T++;
	    break;
	    case 'g': case 'G': G++;
	    break;
	    case 'c': case 'C': C++;
	    break;
	    }
	  }
	printf("%02d\t%d\t%d\t%d\t%d\n", f, A, T, G, C);
        }
      }
    if(p != endpos) printf("\n");
    }
  }

void findcover(int argc, char **argv) {
  if(argc != 2) error(1, 0, "Usage: findcover <mindepth> <minlength>");
  int mindepth = atoi(argv[0]), minlen = atoi(argv[1]);
  if(mindepth <= 0) error(1, 0, "mindepth must be greater than zero");
  init("depth");

  int poss[maxfiles]; //last read position
  mset(poss, 0);
  
  bool in = 0, end = 0;
  int startpos = 0;
  for(int p = 1; !end; p++) {
    bool allcovered = 1;
    rep(f, nfile) {
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

void usage(char **argv) {
  printf("Usage:\n%s <filename> <region> <command> [command options].\n"
	 "Available commands: compare, findcovers.\n"
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
  if(!strcmp(argv[3], "compare"))
    compare(argc-4, argv+4);
  else if(!strcmp(argv[3], "findcover"))
    findcover(argc-4, argv+4);
  else usage(argv); 
  }
