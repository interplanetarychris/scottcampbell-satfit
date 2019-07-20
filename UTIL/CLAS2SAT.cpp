// CLAS2SAT.cpp
// by Scott Campbell campbel7@the-i.net

/******************************************************************************
*                                                                             *
*                                                                             *
*                                                                             *
*******************************************************************************/
#include <cstdlib>
#include <cstdio>
using namespace std;

char name[81], iod_line[300][81];
char tle1[81], tle2[81];
int nobs;

inline char *s_in(char *prompt, char *buffer);
void print_file(char *file, char iod_line[][81]);
void get_obs(char *file, char iod_line[][81]);

/////////////////// MAIN //////////////////////////////////////////////////////
main(int argc, char *argv[])
{
  char file_in[20], file_out[20];
  char buf[81], buf1[81];
  char ssn[5];
  double date_old, date_new;

  if(argc == 2)
    sprintf(file_in, argv[1]);
  else          // default input files
    sprintf(file_in, "classfd.tle");

  // load tle from input file_in

  FILE *fp, *fpi;

  if ((fp = fopen(file_in, "r")) == NULL)
  {
    printf("can't open %s\n", file_in);
    s_in("[Exit]", buf);
    exit(0);
  }

  // classified file tle reader
  tle1[0] = '0';
  do
  {
    if(*tle1 == '1')
    {
      fgets(tle2, 80, fp);
      if(*tle2 == '2')         // found a tle in the classified file
      {
        // get date
        date_new = atof(tle1 + 18);      // classified epoch

        // scan for NORAD Catalog Number
        sscanf(tle2 + 2, "%5s", &ssn);

        // file_out from file_in ssn
        sprintf(file_out, "%5s.txt", ssn);

        // open file_out or go to next classified tle if file doesn't exist
        if((fpi = fopen(file_out, "r")) == NULL) continue;

        // read tle line 1 and line 2
        date_old = 0.;
        buf1[0] = '0';
        do
        {
          if(*buf1 == '1' && strlen(buf1) == 70)
          {
            date_old = atof(buf1 + 18);
            break;
          }
        }while(fgets(buf1, 80, fpi));
        fclose(fpi);

        // printf("\ndate_new = %f date_old = %f\n", date_new, date_old);
        // compare dates
        if(date_new > date_old + 1)
        {
           get_obs(file_out, iod_line);
           print_file(file_out, iod_line);
        }
      }  // if line2
    }    // if line1
    else  sprintf(name, "%s", tle1);
  }while(fgets(tle1, 80, fp));
  fclose(fp);
  system("sats2history");
  // s_in("[Done]", buf);
}    // end main

////////////////// FUNCTIONS //////////////////////////////////////////////////

inline char *s_in(char *prompt, char *buffer)
{
  printf("%s", prompt);
  if (feof(stdin)) exit(0);
  gets(buffer);
  return buffer;
}

void print_file(char *file, char iod_line[][81])
{
   FILE *fpo;
   fpo = fopen(file, "w");
   fprintf(fpo, "%s", name);
   fprintf(fpo, "%s", tle1);
   fprintf(fpo, "%s\n", tle2);
   for(int i = 0; i < nobs; i++)
   {
      fprintf(fpo, "%s", iod_line[i]);
   }
   fclose(fpo);
}

// load IOD lines from input file_in
void get_obs(char *file, char iod_line[][81])
{
  FILE *fpo;
  fpo = fopen(file, "r");

  nobs = 0;

  // find observation entries and write to iod_line matrix
  while(fgets(iod_line[nobs], 80, fpo))
  if(strlen(iod_line[nobs]) > 58 && strlen(iod_line[nobs]) != 70) nobs++;

  fclose(fpo);
  return;
}

