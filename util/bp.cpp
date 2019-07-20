// BatchProcess.cpp
// by Scott Campbell campbel7@hughes.net

/******************************************************************************
*                                                                             *
*                                                                             *
*                                                                             *
*******************************************************************************/
#include <ctype.h>     // for isdigit in read_tle
#include "date.h"      // for Date
using namespace std;

// subroutines
int read_tle(char *file);
void get_obs(char *file, char iod_line[][81], double odata[]);
void sort(char iod_line[][81], double odata[]);
inline char *s_in(char *prompt, char *buffer);

// global variables
int nobs;
char ssn[5], name[81], tle1[81], tle2[81], buf[81];
double tle, nn, odata[36];
Date t1;  // now

/////////////////// MAIN //////////////////////////////////////////////////////
main(int argc, char *argv[])
{
  int i, j;
  double jd;
  char iod_line[200][81], file[81], file_in[81];
  FILE *fp;

  for(i = 1; i < argc; i++)
  {
     sprintf(file, argv[i]);

     // read epoch and mean motion
     read_tle(file);

     // load obs less than 15 days old, /* except for low drag sats? */
     get_obs(file, iod_line, odata);

     // sort and delete duplicates
     sort(iod_line, odata);


     // delete all but most recent 36 obs
     if(nobs > 36)
     {
       for(j = 0; j < 36; j++) iod_line[j] = iod_line[nobs - 36 + j];
       nobs = 36;
     }


     // write to file
     fp = fopen(file, "w");
     fprintf(fp, "%s", name);
     fprintf(fp, "%s", tle1);
     fprintf(fp, "%s\n", tle2);
     for(j = 0; j < nobs; j++)
     {
       fprintf(fp, "%s", iod_line[j]);
     }
     fclose(fp);

     // create and run time batch file
     fp = fopen("time.inp", "w");
     fprintf(fp, "%s\n\n\n\n\n\n\n", "t");
     fprintf(fp, "%s\n", "w");
     fprintf(fp, "%s\n", "u");
     fprintf(fp, "%s\n", "q");
     fprintf(fp, "%s\n", "q");
     fclose(fp);

     // process sat file
     if(nobs == 0) continue;
     if(nobs > 0)                // do this for all cases with obs
     {
       // create and run time batch file
       fp = fopen("time.inp", "w");
       fprintf(fp, "%s\n\n\n\n\n\n\n", "t");
       fprintf(fp, "%s\n", "w");
       fprintf(fp, "%s\n", "u");
       fprintf(fp, "%s\n", "q");
       fprintf(fp, "%s\n", "q");
       fclose(fp);
       sprintf(file_in, "satfit %s ba<time.inp", file);
       system(file_in);
       sprintf(file_in, "SATFIT %s.txt", ssn);   // check for outliers
       system(file_in);

       // create and run elcor batch file for anomaly and node
       fp = fopen("batch.inp", "w");
       fprintf(fp, "%s\n", "s");
       fprintf(fp, "%d %d\n", 0, 4);
       fprintf(fp, "%s\n", "d");
       fprintf(fp, "%s\n", "w");
       fprintf(fp, "%s\n", "u");
       fprintf(fp, "%s\n", "q");
       fprintf(fp, "%s\n", "q");
       fclose(fp);
       sprintf(file_in, "elcord %s ba<batch.inp", file);
       system(file_in);
       sprintf(file_in, "DEL batch.inp");
       system(file_in);
     }
     if(nobs > 3)
     {
       // create and run bstar batch file
       fp = fopen("bstar.inp", "w");
       fprintf(fp, "%s\n\n\n\n\n\n\n", "s");
       fprintf(fp, "%s\n", "q");
       if(nobs > 7)
       {
          fprintf(fp, "%s\n", "b");
          fprintf(fp, "%s\n\n\n", "a");
          fprintf(fp, "%s\n", "q");
       }
       fprintf(fp, "%s\n\n\n\n", "s");
       fprintf(fp, "%s\n", "q");
       if(nobs > 7)
       {
          fprintf(fp, "%s\n", "b");
          fprintf(fp, "%s\n\n\n", "a");
          fprintf(fp, "%s\n", "q");
       }
       fprintf(fp, "%s\n\n\n\n", "s");
       fprintf(fp, "%s\n", "q");
       if(nobs > 7)
       {
          fprintf(fp, "%s\n", "b");
          fprintf(fp, "%s\n\n\n", "a");
          fprintf(fp, "%s\n", "q");
       }
       fprintf(fp, "%s\n\n\n\n", "s");
       fprintf(fp, "%s\n", "q");
       if(nobs > 7)
       {
          fprintf(fp, "%s\n", "b");
          fprintf(fp, "%s\n\n\n", "a");
          fprintf(fp, "%s\n", "q");
       }
       fprintf(fp, "%s\n\n\n\n", "z");
       fprintf(fp, "%s\n", "q");
       if(nobs > 7)
       {
          fprintf(fp, "%s\n", "b");
          fprintf(fp, "%s\n\n\n", "a");
          fprintf(fp, "%s\n", "q");
       }
       fprintf(fp, "%s\n\n\n\n", "z");
       fprintf(fp, "%s\n", "q");
       if(nobs > 7)
       {
          fprintf(fp, "%s\n", "b");
          fprintf(fp, "%s\n\n\n", "a");
          fprintf(fp, "%s\n", "q");
       }
       fprintf(fp, "%s\n", "w");
       fprintf(fp, "%s\n", "u");
       fprintf(fp, "%s\n", "q");
       fprintf(fp, "%s\n", "q");
       fclose(fp);
       sprintf(file_in, "satfit %s ba<bstar.inp", file);
       system(file_in);
       sprintf(file_in, "DEL bstar.inp");
       system(file_in);
     }
     if(nobs > 9 && nn > 3)
     {
       // create and run elcor batch file
       fp = fopen("elcor.inp", "w");
       fprintf(fp, "%s\n", "s");
       fprintf(fp, "%d %d\n", 0, 4);
       fprintf(fp, "%s\n", "d");
       fprintf(fp, "%s\n", "s");
       fprintf(fp, "%d %d\n", 3, 4);
       fprintf(fp, "%s\n", "d");
       fprintf(fp, "%s\n", "s");
       fprintf(fp, "%d %d %d\n", 0, 3, 4);
       fprintf(fp, "%s\n", "d");
       fprintf(fp, "%s\n", "s");
       fprintf(fp, "%d %d\n", 5, 6);
       fprintf(fp, "%s\n", "d");
       fprintf(fp, "%s\n", "s");
       fprintf(fp, "%d %d %d %d\n", 0, 3, 4, 5);
       fprintf(fp, "%s\n", "d");
       fprintf(fp, "%s\n", "d");
       fprintf(fp, "%s\n", "s");
       fprintf(fp, "%d %d %d %d\n", 3, 4, 5, 6);
       fprintf(fp, "%s\n", "d");
       fprintf(fp, "%s\n", "d");
       fprintf(fp, "%s\n", "s");
       fprintf(fp, "%d %d %d %d %d\n", 0, 3, 4, 5, 6);
       fprintf(fp, "%s\n", "d");
       fprintf(fp, "%s\n", "d");
       fprintf(fp, "%s\n", "s");
       fprintf(fp, "%d %d\n", 0, 1);
       fprintf(fp, "%s\n", "d");
       fprintf(fp, "%s\n", "s");
       fprintf(fp, "%d\n", 2);
       fprintf(fp, "%s\n", "d");
       fprintf(fp, "%s\n", "w");
       fprintf(fp, "%s\n", "u");
       fprintf(fp, "%s\n", "q");
       fprintf(fp, "%s\n", "q");
       fclose(fp);
       sprintf(file_in, "elcord %s ba<elcor.inp", file);
       system(file_in);
       sprintf(file_in, "DEL elcor.inp");
       system(file_in);
     }
     // do this for all cases with obs > 0
     // create and run last batch file, re-adjust time
     fp = fopen("last.inp", "w");
     fprintf(fp, "%s\n", "a");
     fprintf(fp, "%s\n", "l");
     fprintf(fp, "%s\n", "q");
     fprintf(fp, "%s\n", "w");
     fprintf(fp, "%s\n", "u");
     fprintf(fp, "%s\n", "q");
     fprintf(fp, "%s\n", "q");
     fclose(fp);
     sprintf(file_in, "satfit %s ba<last.inp", file);
     system(file_in);
     sprintf(file_in, "DEL last.inp");
     system(file_in);
     sprintf(file_in, "satfit %s ba<time.inp", file);
     system(file_in);
     sprintf(file_in, "DEL time.inp");
     system(file_in);
     sprintf(file_in, "SATFIT %s.txt", ssn);  // check results
     system(file_in);
  }
}    // end main

////////////////// FUNCTIONS //////////////////////////////////////////////////

inline char *s_in(char *prompt, char *buffer)
{
  printf("%s", prompt);
  if (feof(stdin)) exit(0);
  gets(buffer);
  return buffer;
}

// smart tle reader, reads first tle in file, no conversion of elements
int read_tle(char *file)
{
  FILE *fpo;
  if ((fpo = fopen(file, "r")) == NULL) return(1);
  tle1[0] = '0';
  do
  {
    if(*tle1 == '1' && strlen(tle1) == 70)
    {
      fgets(tle2, 80, fpo);
      if(*tle2 == '2' && strlen(tle2) == 70)
      {
        fclose(fpo);
        sscanf(tle1 + 2, "%5s", &ssn);
        tle = atof(tle1 + 18);      // epoch in tle format
        buf[0] = tle2[63];
        tle2[63] = '\0';
        nn = atof(tle2 + 52);       // mean motion, revolutions/day
        tle2[63] = buf[0];          // restore line2[63]
        return(0);
      }
    }
    else    //  if(*tle1 != '1')
    {
      sprintf(name, tle1);
    }
  }while(fgets(tle1, 80, fpo));
  return(1);   //  return here if no tle in file
}

// load IOD lines from input file_in
void get_obs(char *file, char iod_line[][81], double odata[])
{
  int year, month, day, hour, min;
  double sec;
  char sbuf[5];
  FILE *fpo;
  if ((fpo = fopen(file, "r")) == NULL)
  {
     printf("Can't open input file\n");
     s_in("[exit]", buf);
     return;
  }
  nobs = 0;

  // find observation entries and write to iod_line matrix
  while(fgets(iod_line[nobs], 80, fpo))
  {
    if(strlen(iod_line[nobs]) > 58 && strlen(iod_line[nobs]) != 70 &&
       isdigit(iod_line[nobs][0]))
    {
      sscanf(iod_line[nobs] + 23, "%4d %2d %2d %2d %2d %5s",
                              &year, &month, &day, &hour, &min, &sbuf);
      sec = atof(sbuf);
      sec /= pow(10, strlen(sbuf) - 2);

      Date t2(year, month, day, hour, min, sec);
      odata[nobs] = t2.jd;
      // leo sats +15  deep sats +30
      if((nn <= 14 && t1.jd < t2.jd + 15) ||
         (nn <   3 && t1.jd < t2.jd + 22) ||
         (nn >  14 && t1.jd < t2.jd + 10)) nobs++;
    }
  }
  fclose(fpo);
  return;
}

void sort(char iod_line[][81], double odata[])
{
  int i, j;
  double k;
  sort:
  j = 0;
  for(i = 0; i < nobs - 1; i++)
  {
    if(odata[i] > odata[i + 1])    // switch low/high
    {
      k = odata[i];
      odata[i] = odata[i + 1];
      odata[i + 1] = k;
      buf = iod_line[i];
      iod_line[i] = iod_line[i + 1];
      iod_line[i + 1] = buf;
      j++;
    }
    if(odata[i] == odata[i + 1])   // remove duplicates
    {
      for(int s = i; s < nobs - 1; s++)
      {
         odata[s] = odata[s + 1];
         iod_line[s] = iod_line[s + 1];
      }
      nobs--;
    }
  }
  if(j > 0)  goto sort;
}

