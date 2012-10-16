#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <sstream>
#include <ctime>
#include <sys/timeb.h>  
using namespace std;

#include "qpalma_dp.h"

void initialize_quality_plifs(penalty_struct* qplifs, double* params) {

   printf("begin of initialize_quality_plifs \n");
   double limits[10] = {-5.0, 0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0};

   int idx, pos;
   int ctr = 36;

   for(idx=0;idx<30;idx++) {
      penalty_struct new_plif;
      init_penalty_struct(new_plif);
   
      new_plif.min_len = -5;
      new_plif.max_len = 40;
      new_plif.len = 10;

      new_plif.limits    = (double*) malloc(sizeof(double)*10);
      if (new_plif.limits == NULL)
         printf("malloc (limits)\n");

      new_plif.penalties = (double*) malloc(sizeof(double)*10);
      if (new_plif.penalties == NULL)
         printf("malloc (penalties)\n");

      for(pos=0;pos<10;pos++) {
         new_plif.limits[pos]     = limits[pos];
         new_plif.penalties[pos]  = params[ctr];
         ctr++;
      }

      qplifs[idx] = new_plif;
   }
   
   printf("end of initialize_quality_plifs \n");
   return;
}


int main(int argc, char** argv){

   printf("Starting test_fm\n");
   int idx;
   char* alphabet = (char*) malloc(sizeof(char)*4);

   alphabet[0] = 'a';
   alphabet[1] = 'c';
   alphabet[2] = 'g';
   alphabet[3] = 't';

   int est_len =  36;
   int dna_len = 10000;
   int nr_paths = 2;

   int nr_donor_sites = 0 ;
   int d_len = dna_len;

   printf("begin of donor/acceptor stuff\n");
   double* donor     = (double*) malloc(sizeof(double)*dna_len);
   double* acceptor  = (double*) malloc(sizeof(double)*dna_len);

   for(idx=0;idx<dna_len;idx++) {
      donor[idx]     = drand48();
      acceptor[idx]  = drand48();
   }

   for (int ii=0; ii<d_len; ii++) {
    if(isfinite(donor[ii])) {
      nr_donor_sites++ ;
    }
   }

   int* donor_sites = new int[nr_donor_sites];
   int donor_idx = 0 ;
   for (int ii=0; ii<d_len; ii++) {
    if(isfinite(donor[ii])) {
      donor_sites[donor_idx] = ii+1 ;
      donor_idx++ ;
    }  
   }

   printf("end of donor/acceptor stuff\n");

   Pre_score** matrices = new Pre_score*[nr_paths];
   for (int z=0; z<nr_paths; z++)
      matrices[z] = new Pre_score[dna_len * est_len];

   char* est = (char*) malloc(sizeof(char)*est_len);
   char* dna = (char*) malloc(sizeof(char)*dna_len);

   for(idx=0;idx<est_len;idx++)
      est[idx] = alphabet[lrand48()%4];

   for(idx=0;idx<dna_len;idx++)
      dna[idx] = alphabet[lrand48()%4];
      
   double* prb = (double*) malloc(sizeof(double)*est_len);

   for(idx=0;idx<est_len;idx++)
      prb[idx] = lrand48() % 45;

   double* matchmatrix = (double*) malloc(sizeof(double)*6);
   matchmatrix[0] = 0;
   matchmatrix[1] = -1;
   matchmatrix[2] = -1;
   matchmatrix[3] = -1;
   matchmatrix[4] = -1;
   matchmatrix[5] = 0;

   REAL limits_array[10] = {19.999999999999993, 33.362010744001154, 55.651188044142465, 92.831776672255558, 154.85273653622525, 258.30993300297655, 430.88693800637651, 718.76273276092525, 1198.968500637882, 1999.9999999999982};
   REAL pen_array[10] = {0.024414829630412388, 0.025836097651978609, 0.17473437525993854, 0.46902269254060858, -0.14383509580738071, -0.20766155250201096, -0.32734767088988909, -0.18507220858585502, -0.097862557903448902, -0.034886887426478781};

   penalty_struct functions;

   init_penalty_struct(functions);

   functions.min_len = 20;
   functions.max_len = 2000;
   functions.len = 10;
   functions.limits = &limits_array[0];
   functions.penalties = &pen_array[0];

   penalty_struct* qualityScores = (penalty_struct*) malloc(sizeof(penalty_struct)*30);

   int num_params = 336;
   double* params = (double*) malloc(sizeof(double)*num_params);

   double parameter;
   int ctr = 0;

   string s;
   ifstream infile("param.txt");

   if (infile.fail()) 
     printf("Could not open file param.txt for reading.\n");

   while (getline(infile, s)) {
      if (s.length() == 0 || s[0] == '#') continue;

      istringstream iss(s);
      iss >> parameter;
      params[ctr] = parameter;
      ctr++;
   }

   infile.close();
   
   initialize_quality_plifs(qualityScores,params);

   bool remove_duplicate_scores = true;
   int* max_score_positions = new int[nr_paths*2];
   mode currentMode = USE_QUALITY_SCORES;

   //time_t startTime;
   //time_t endTime;

   printf("calling fill_matrix\n");

   //time(&startTime);

   for(idx=0;idx<1;idx++)
      fill_matrix(nr_paths, matrices, est_len, dna_len, est, dna, prb, &functions, matchmatrix, qualityScores, donor, acceptor, remove_duplicate_scores, nr_donor_sites, donor_sites, max_score_positions);

   //time(&endTime);

   //printf ("Scan time: %6.3f\n", difftime(endTime,startTime));

   printf("End of test_fm \n");
}
