// SATS2FILE.cpp
// by Scott Campbell campbel7@the-i.net

/******************************************************************************
*                                                                             *
*                                                                             *
*                                                                             *
*******************************************************************************/
#include <cstdlib>
#include <cstdio>
#include <time.h>
using namespace std;

int read_tle(char *file);
inline char *s_in(char *prompt, char *buffer);
char  name[81], buf1[81], buf2[81];

/////////////////// MAIN //////////////////////////////////////////////////////
main(int argc, char *argv[])
{
  char file_in[81];
  char ssn[5], buf[81];
  FILE *fp, *fpo;

  // open ssr.id
  fp = fopen("ssr.id", "r");

  // open classified file
  fpo = fopen("classfd.sat", "w");

  while(fgets(buf, 80, fp))   // get a line from ssr.id
  {
     // read the ssn from ssr.id
     sscanf(buf, "%5s", &ssn);
     // identify the file by that name (ssn.txt)
     sprintf(file_in, "%5s.txt", ssn);

     // read tle name, line 1 and line 2
     // if read fails then go to next ssr.id line
     if(read_tle(file_in)) continue;

     //  write to classified file
     fprintf(fpo, "%s", name);
     fprintf(fpo, "%s", buf1);
     fprintf(fpo, "%s", buf2);
  }
  fclose(fp);
  fclose(fpo);
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
  if ((fp = fopen(file, "r")) == NULL) return(1);
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
    else    //  if(*buf1 != '1')
    {
      sprintf(name, "%s", buf1);
    }
  }while(fgets(buf1, 80, fp));
  return(1);
}

