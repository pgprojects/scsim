
/*************************** scSim.cpp *********************** ply 2016-7-11 *
*                                                                            *
*  Simulation of single cell evolution: diploid ACGT                         *
*                                                                            *
*****************************************************************************/

#include <mpi.h>
#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <math.h>
#include <algorithm>
#include <time.h>
#include "nr3.h"
#include "ran.h"
#include "gamma.h"
#include "deviates.h"

using namespace std;

void PopGrowth(int&,int,int&,int,double,double,int*,int**,int);

void TypGrowth(int,int,int&,int&,double,double,int*,int**,int);

void SamplingWOreplacement(int,int,int*,int);

double CalculateMean(double *, int);

double CalcVar(double, double *, int);

int main (int argc, char *argv[]){

    /* MPI related variables */
    int nprocs,myrank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD,&myrank);

    Ran myran(time(0)+myrank); 

    string s1;
    stringstream out1;
    out1<<argv[1];
    s1=out1.str();
    int plateau = atoi(s1.c_str());

    string s2;
    stringstream out2;
    out2<<argv[2];
    s2=out2.str();
    double br = atof(s2.c_str());

    double dr = br;


    int i, j, k, ngen, tgen=48+plateau;
    int ss=100;

    int siteNo;
    double mr=0.0000001; 
    int **treesampled = new int* [tgen];
    for(k=0; k<tgen; k++) treesampled[k]= new int[ss];
    int *nopopg = new int [tgen];

    int myquote = (int)myrank/10;
    int myremainder = myrank % 10;
    if(myremainder==0){

      int ngen=0;
      int nopop=0;

      int **tree = new int* [tgen];
      for(k=0; k<tgen; k++) tree[k]= new int[5000000];
      int *pop0 = new int[tgen];
      int *pop1 = new int[tgen];
      int *pop2 = new int[tgen];


      TypGrowth(5,plateau,nopop,ngen,br,dr,pop0,tree,myran.int64());


    for(k=0; k<tgen; k++) for(i=0; i<ss; i++) treesampled[k][i]=-1;
    
    int *samplingO = new int[ss];
    int *sampling = new int[ss];
    SamplingWOreplacement(nopop,ss,sampling,myran.int64());
    int *prevsampling = new int[ss];
    for(i=0; i<ss; i++) prevsampling[i]=-1;

    int totalcount =ss;
    nopopg[tgen-1]=ss;
    for(k=tgen; k>0; k--){
      for(i=0; i<nopopg[k-1]; i++) samplingO[i]=sampling[i];
      int coalescence = 0;
      for(i=0; i<nopopg[k-1]; i++){
	prevsampling[i]=tree[k-1][samplingO[i]];
	if(prevsampling[i]==0){
	  treesampled[k-1][i]=0;
	  sampling[i-coalescence]=samplingO[i];
	} else {
          if(prevsampling[i]<0){
            treesampled[k-1][i]=0;
            sampling[i-coalescence]=samplingO[i]-prevsampling[i];
          } else {
	    int checkprev=0;
	    for(j=0; j<i; j++) if(prevsampling[j]>0 && prevsampling[j]<100000 && prevsampling[j]==prevsampling[i]) checkprev=j+1;
	    if(checkprev==0){
	      treesampled[k-1][i]=i+1-coalescence;
	      sampling[i-coalescence]=tree[k-1][samplingO[i]]-1;
	    } else {
	      treesampled[k-1][i]=treesampled[k-1][checkprev-1];
	      totalcount--;
	      coalescence++;
	    }
	  }
	}
	
      }

      if(k>1){
	nopopg[k-2]=totalcount;
      } 

    }


    for(i=0; i<tgen; i++) delete [] tree[i];
    delete[] tree;
    delete[] prevsampling;
    delete[] samplingO;
    delete[] sampling;

    delete[] pop0;
    delete[] pop1;
    delete[] pop2;

    } // myrank==0

    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Comm simGroup;
    int newprocs, newrank;
    MPI_Comm_split(MPI_COMM_WORLD,myquote,myremainder, &simGroup);
    MPI_Comm_size(simGroup, &newprocs);
    MPI_Comm_rank(simGroup, &newrank);

    for(k=0; k<tgen; k++){
      MPI_Bcast(treesampled[k], ss, MPI_UNSIGNED,0,simGroup);
    }
    MPI_Bcast(nopopg, tgen, MPI_UNSIGNED,0,simGroup);
    


    string filename;
    stringstream infile;
    infile<<argv[newrank+3];
    filename=infile.str();
    ifstream fin;
    fin.open(filename.c_str());

    string line, sequence1, sequence2;
    getline(fin, line);
    sequence1 = line;

    while(getline(fin,line)){
      sequence1.append(line);
    }

    siteNo = sequence1.length();
    sequence2 = sequence1; 
    
    char acgt[4] = {'A','C','G','T'};

    string osequence1 = sequence1;
    string osequence2 = sequence2;
    


    int *gen = new int[siteNo];
    int *nmut = new int[siteNo];
    for(j=0; j<siteNo; j++){
      gen[j]=0;
      nmut[j]=0;
    }

    ngen = 1;

    string *samples = new string[2*ss];
    string *samplesO = new string[2*ss];
    samples[0] = osequence1;
    samples[1] = osequence2;
    samples[2] = osequence1;
    samples[3] = osequence2;
    int temp, pick1;;
    for(i=0; i<4; i++){
      Binomialdev binom(siteNo,mr,myran.int64());
      for (j=0; j<binom.dev(); j++){
        temp = myran.int64() % siteNo;
        pick1 = myran.int64() % 4;
        samples[i][temp]=acgt[pick1];
        nmut[temp]++;
        if(gen[temp]==0) gen[temp]=ngen;
      }
    }

    int *prev = new int[ss];
    for(k=1; k<tgen; k++){
      ngen++;
      for (i=0; i<nopopg[k-1]; i++){
        samplesO[2*i]=samples[2*i];
        samplesO[2*i+1]=samples[2*i+1];
      }
      int nocoal =0;
      for (i=0; i<nopopg[k]; i++){
        prev[i]=treesampled[k][i];
        for(j=0; j<i; j++) if (prev[j]>0 && prev[j]==prev[i]) nocoal++;
        if(treesampled[k][i]>0){
          samples[2*i] = samplesO[2*(treesampled[k][i]-1)];
          samples[2*i+1] = samplesO[2*(treesampled[k][i]-1)+1];
          Binomialdev binom1(siteNo,mr,myran.int64());
          for (j=0; j<binom1.dev(); j++){
            temp = myran.int64() % siteNo;
            pick1 = myran.int64() % 4;
            samples[2*i][temp]=acgt[pick1];
            nmut[temp]++;
            if(gen[temp]==0) gen[temp]=ngen;
          }
          Binomialdev binom2(siteNo,mr,myran.int64());
          for (j=0; j<binom2.dev(); j++){
            temp = myran.int64() % siteNo;
            pick1 = myran.int64() % 4;
            samples[2*i+1][temp]=acgt[pick1];
            nmut[temp]++;
            if(gen[temp]==0) gen[temp]=ngen;
          }
        } else {
          samples[2*i] = samplesO[2*(i-nocoal)];
          samples[2*i+1] = samplesO[2*(i-nocoal)+1];
        }

      }
    }



    string genofilename ("genotypes..txt");
    string sg,sgr;
    stringstream rgout,rgrank;
    rgout<<myquote;
    rgrank<<newrank;
    sg=rgout.str();
    sgr=rgrank.str();
    genofilename.insert(10,sgr);
    genofilename.insert(9,sg);
    ofstream fgeno;
    fgeno.open(genofilename.c_str());


    if(newrank==0){
      fgeno << "#chr\tpos\tref\talt";
      for(i=0; i<ss; i++) fgeno << "\tind" << i+1;
      fgeno << "\n";
    }

    int *genotypes = new int [ss];

    for(j=0; j<siteNo; j++){
      int geno = 0;
      for(i=0; i<ss; i++){
        if(samples[2*i][j]!=osequence1[j]) geno++;
        if(samples[2*i+1][j]!=osequence2[j]) geno++;
      }


      int test = 0;
      if(geno>0){
        fgeno << myrank+1 << "\t" << j+1 << "\t" << osequence1[j];
        for(i=0; i<ss; i++){
          int a1 = 0;
          int a2 = 0;
          if(samples[2*i][j]!=osequence1[j]){
            if(test==0){
              fgeno << "\t" << samples[2*i][j];
              test=1;
            }
            a1 = 1;
          }
          if(samples[2*i+1][j]!=osequence2[j]){
            if(test==0){
              fgeno << "\t" << samples[2*i][j];
              test=1;
            }
            a2 = 1;
          }

          genotypes[i] = a1 + a2;

        }

        for(i=0; i<ss; i++) fgeno << "\t" << genotypes[i];
        fgeno << "\n";

      }
    }

    fgeno.close();


    delete [] genotypes;
    delete [] prev;
    delete [] nopopg;
    for(i=0; i<tgen; i++) delete [] treesampled[i];
    delete[] treesampled;

    delete [] gen;
    delete [] nmut;
    delete [] samples;
    delete [] samplesO;

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    return 0;

} // main

	  
  double CalculateMean(double *array, int max)

    { 

        double sum = 0; 

        for(int i = 0; i < max; i++) 

            sum += array[i]; 

        return (sum/max); 

    } 

  double CalcVar(double m, double *array, int max) 

    { 

        double temp = 0; 

        for(int i = 0; i < max; i++) 

        { 

             temp += (array[i] - m) * (array[i] - m) ; 

        } 

        return temp / (max-1); 

    } 

 
void SamplingWOreplacement
(
 int populationSize,    // size of set sampling from
 int sampleSize,        // size of each sample
 int * samples,  // output
 int rseed
 )
{

  Ran myrans(rseed);
  // Use Knuth's 
  int& n = sampleSize;
  int& N = populationSize;

  int t = 0; 
  int m = 0; 
  double u;

  while (m < n)
    {
      u = myrans.doub(); 

      if ( (N - t)*u >= n - m )
        {
          t++;
        }
      else
        {
          samples[m] = t;
          t++; m++;
        }
    }
}



void PopGrowth 
(
 int &popSize, 
 int ps,   
 int &gen, 
 int gens, 
 double br,
 double dr,
 int *popgn,
 int ** tree,
 int prseed
 )
{
  int i,j,k, nopop;
  Ran myranp(prseed);

    int finalPop = popSize;
    int psmax = finalPop;
    for(k=0; k<gens; k++){
      double temp=(double)finalPop*(1+br-dr);
      finalPop = round(temp);
      if(finalPop>psmax) psmax=finalPop;
    }
    int *birth = new int[psmax];
    int *death = new int[psmax];

    nopop = popSize;


    for(k=0; k<gens; k++){
      int nob = (int)(nopop*br);
      int nod = (int)(nopop*dr);
      SamplingWOreplacement(nopop,nob,birth,myranp.int64());
      SamplingWOreplacement(nopop,nod,death,myranp.int64());
      if(nob==nod){
        for(i=ps; i<ps+nopop; i++) tree[gen][i]=0;
        for (i=0; i<nod; i++){
          tree[gen][ps+birth[i]]=ps+birth[i]+1;
          tree[gen][ps+death[i]]=ps+birth[i]+1;
        }
      } else if(nob>nod){
        for(i=ps; i<ps+nopop; i++) tree[gen][i]=0;
        for (i=0; i<nod; i++){
          tree[gen][ps+birth[i]]=ps+birth[i]+1;
          tree[gen][ps+death[i]]=ps+birth[i]+1;
        }
        for (i=0; i<nob-nod; i++){
          tree[gen][ps+birth[i]]=ps+birth[i]+1;
          tree[gen][ps+nopop+i]=ps+birth[i]+1;
        }
        nopop+=nob-nod;
      } else {
        for(i=ps; i<ps+nopop-nod+nob; i++) tree[gen][i]=0;
        for (i=0; i<nob; i++){
          tree[gen][ps+birth[i]]=ps+birth[i]+1;
          tree[gen][ps+death[i]]=ps+birth[i]+1;
        }
        int curd = ps+death[nob];
        for (i=1; i<=nod-nob; i++){
          int next = ps+death[nob+i];
          if(i==nod-nob) next = nopop-1;
          for(j=curd; j<next-i; j++){
            if(tree[gen][j]==0) tree[gen][j]=-i;
          }
          curd = next-i;
        }
        for(j=curd; j<nopop; j++){
          if(tree[gen][j]>0){
            tree[gen][curd]=tree[gen][j];
            tree[gen][j]=0;
            curd++;
          }
        }
        nopop-=nod-nob;
      } 
      popgn[gen]=nopop;
      gen++;
    }
    delete [] birth;
    delete [] death;



  popSize = nopop;

}


void TypGrowth 
(
 int kgen,    
 int plateau, 
 int &tps,  
 int &tgen, 
 double br1,
 double br2,
 int *popn,
 int ** tree,
 int prseed
 )
{
  
  Ran myrang(prseed);
  int i,j,k,npop=1;


  PopGrowth(npop,tps,tgen,2*kgen,1,0,popn,tree,myrang.int64());
  

  PopGrowth(npop,tps,tgen,kgen,0.5,0,popn,tree,myrang.int64());


  int adjgen2=abs((int)(10*br1-5));
  double br1t = 0.5;
  double br1s = 0.1;
  double dr1=0;

  if(br1<=0.5){
    br1s=-0.1;
    for(k=0; k<adjgen2; k++){
      br1t += br1s;
      PopGrowth(npop,tps,tgen,1,br1t,0,popn,tree,myrang.int64());
    }
  } else {
    for(k=0; k<adjgen2; k++){
      br1t += br1s;
      PopGrowth(npop,tps,tgen,1,br1t,dr1,popn,tree,myrang.int64());
    }
  }


  int adjgen1=(int)(10*br1-dr1*10);
  for(k=0; k<adjgen1-1; k++){
    dr1+=0.1;
    if(dr1<br1t){
      PopGrowth(npop,tps,tgen,1,br1t,dr1,popn,tree,myrang.int64());
    }
  }


  double adjarr[5] = {0.06,0.04,0.03,0.02,0.01};
  if(br1>=0.1){
    for(k=0; k<5; k++){
      PopGrowth(npop,tps,tgen,1,br1t,br1-adjarr[k],popn,tree,myrang.int64());
    }
  }


  PopGrowth(npop,tps,tgen,plateau,br1,br1,popn,tree,myrang.int64());


  PopGrowth(npop,tps,tgen,2,br2,br2-0.05,popn,tree,myrang.int64());
  PopGrowth(npop,tps,tgen,22,br2,br2-0.1,popn,tree,myrang.int64());
  tps+=npop;

}

