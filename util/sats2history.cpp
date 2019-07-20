// SATS2HISTORY.cpp
// by Scott Campbell campbel7@the-i.net

/******************************************************************************
*                                                                             *
*                                                                             *
*                                                                             *
*******************************************************************************/
#include <cstdlib>
#include <cstdio>
using namespace std;

int read_tle(char *file);
inline char *s_in(char *prompt, char *buffer);
char buf1[81], buf2[81];

/////////////////// MAIN //////////////////////////////////////////////////////
main(int argc, char *argv[])
{
  char file_in[20], file_out[20];
  char tle[81];
  char ssn[5], buf[81];
  double date_old, date_new;

  FILE *fp, *fpi, *fpo;

  // open ssr.id
  fp = fopen("ssr.id", "r");

  while(fgets(buf, 80, fp))   // get a line from ssr.id
  {
     // read the ssn from ssr.id
     sscanf(buf, "%5s", &ssn);
     // identify the file by that name (ssn.txt)
     sprintf(file_in, "%5s.txt", ssn);

     // read tle line 1 and line 2
     // if read sat file fails then create files and go to next ssr.id line
     if ((fpi = fopen(file_in, "r")) == NULL)
     {
       fclose(fpi);
       // create empty sat file
       sprintf(file_out, "%05s.txt", ssn);
       fpo = fopen(file_out, "w");
       fprintf(fpo, "\0");
       fclose(fpo);
       // create empty history file
       sprintf(file_out, "history/%05s.txt", ssn);
       fpo = fopen(file_out, "w");
       fprintf(fpo, "\0");
       fclose(fpo);
       // go to next ssr.id line
       continue;
     }

     // read tle line 1 and line 2
     // if read fails then go to next ssr.id line
     if(read_tle(file_in)) continue;

     // get date_new (date in the ssn.txt file)
     date_new = atof(buf1 + 18);      // epoch in tle format at this point

     // open history file
     sprintf(file_out, "history/%05s.txt", ssn);
     if ((fpi = fopen(file_out, "r")) == NULL)
     {
       date_old = 0.;
       goto write;
     }

     // read second to last line in history file
     while(fgets(tle, 80, fpi))
     // get date_old
     if(*tle == '1') date_old = atof(tle + 18);
     fclose(fpi);

     write:
     // printf("\ndate_new = %f date_old = %f\n", date_new, date_old);
     // compare dates
     if(date_new > date_old)
     {
       //  write to history if new date is later than old date
       fpo = fopen(file_out, "a");
       fprintf(fpo, "\n");
       fprintf(fpo, "%s", buf1);
       fprintf(fpo, "%s", buf2);
       fclose(fpo);
     }
  }
  fclose(fp);
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
  FILE *fp;
  if ((fp = fopen(file, "r")) == NULL)
  {
    // printf("can't open %s\n", file);
    // s_in("", buf1);
    return(1);
  }

  buf1[0] = '0';
  do
  {
    if(*buf1 == '1' && strlen(buf1) == 70)
    {
      fgets(buf2, 80, fp);
      if(*buf2 == '2' && strlen(buf2) == 70)
      {
        fclose(fp);
        return(0);
      }
    }
  }while(fgets(buf1, 80, fp));
  // if no tle then return
  return(1);
}

