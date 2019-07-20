// MERGE.CPP
#include <cstdlib>
#include <cstdio>
using namespace std;

/*
Merges two sorted files where the newer tles are found in two_file.

merge one_file two_file out_file
*/

//////////// DECLARE FUNCTIONS ////////////////////////////////////////////////

int n;
inline char *s_in(char *prompt, char *buffer);
void sort(char big_file[][81]);

/////////////////// MAIN //////////////////////////////////////////////////////
int main(int argc, char *argv[])
{
  int i;
  char buf11[81], buf12[81], buf13[81], buf21[81], buf22[81], buf23[81];
  char one_file[81], two_file[81], out_file[81];
  char big_file[50000][81];

  if(argc == 4)
  {
    sprintf(one_file, argv[1]);
    sprintf(two_file, argv[2]);
    sprintf(out_file, argv[3]);
  }
  else
  {
    printf("file name arguments error");
    s_in("[exit]", buf11);
    exit(0);
  }

  // read all TLE's into big_file, three line format, no echo
  FILE *fp1, *fp2, *fpo;
  // load one_file
  if ((fp1 = fopen(one_file, "r")) == NULL)
  {
    printf("Can't open one_file\n");
    s_in("[exit]", buf11);
    exit(0);
  }
  // load two_file
  if ((fp2 = fopen(two_file, "r")) == NULL)
  {
     printf("Can't open two_file\n");
     s_in("[exit]", buf11);
     exit(0);
  }

  // start by reading one_file and two_file into buffers
  fgets(buf11, 80, fp1);        // name
  fgets(buf12, 80, fp1);        // line1
  fgets(buf13, 80, fp1);        // line2
  fgets(buf21, 80, fp2);        // name
  fgets(buf22, 80, fp2);        // line1
  fgets(buf23, 80, fp2);        // line2

  // which ssn is smaller?
  n = 0;
  do
  {
    // if buf1 = buf2, load buf1
    // and replace both buf1 and buf2
    if(atoi(buf13 + 1) == atoi(buf23 + 1))
    {
      big_file[n]     = buf11;
      big_file[n + 1] = buf12;
      big_file[n + 2] = buf13;
      n += 3;
      // replace buf1 with one_file unless it's empty, in which case
      if(fgets(buf11, 80, fp1))
      {
        fgets(buf12, 80, fp1);        // line1
        fgets(buf13, 80, fp1);        // line2
      }
      // replace buf1 with two_file unless it's empty, in which case
      else if(fgets(buf11, 80, fp2))
      {
        fgets(buf12, 80, fp2);        // line1
        fgets(buf13, 80, fp2);        // line2
      }
      // we are done
      else break;
      // replace buf2 with two_file unless it's empty, in which case
      if(fgets(buf21, 80, fp2))
      {
        fgets(buf22, 80, fp2);        // line1
        fgets(buf23, 80, fp2);        // line2
        continue;
      }
      // replace buf2 with one_file unless it's empty, in which case
      else if(fgets(buf21, 80, fp1))
      {
        fgets(buf22, 80, fp1);        // line1
        fgets(buf23, 80, fp1);        // line2
        continue;
      }
      // we are done
      else break;
    }
    // if tle1 < tle2 then load tle1 into big_file
    // and replace buf1 from one_file unless one_file is empty
    if(atoi(buf13 + 1) < atoi(buf23 + 1))
    {
      big_file[n]     = buf11;
      big_file[n + 1] = buf12;
      big_file[n + 2] = buf13;
      n += 3;
      // replace from one_file unless it's empty, in which case
      if(fgets(buf11, 80, fp1))
      {
        fgets(buf12, 80, fp1);        // line1
        fgets(buf13, 80, fp1);        // line2
        continue;
      }
      // replace from two_file unless it's empty, in which case
      if(fgets(buf11, 80, fp2))
      {
        fgets(buf12, 80, fp2);        // line1
        fgets(buf13, 80, fp2);        // line2
        continue;
      }
      // load the last file
      else
      {
        big_file[n]     = buf21;
        big_file[n + 1] = buf22;
        big_file[n + 2] = buf23;
        n += 3;
        break;  // and we are done
      }
    }
    // if tle2 < tle1 then load tle1 into big_file
    // and replace buf2 from two_file unless two_file is empty
    if(atoi(buf23 + 1) < atoi(buf13 + 1))
    {
      big_file[n]     = buf21;
      big_file[n + 1] = buf22;
      big_file[n + 2] = buf23;
      n += 3;
      // replace from two_file unless it's empty, in which case
      if(fgets(buf21, 80, fp2))
      {
        fgets(buf22, 80, fp2);        // line1
        fgets(buf23, 80, fp2);        // line2
        continue;
      }
      // replace from one_file unless it's empty, in which case
      if(fgets(buf21, 80, fp1))
      {
        fgets(buf22, 80, fp1);        // line1
        fgets(buf23, 80, fp1);        // line2
        continue;
      }
      // load the last file
      else
      {
        big_file[n]     = buf11;
        big_file[n + 1] = buf12;
        big_file[n + 2] = buf13;
        n += 3;
        break;  // and we are done
      }
    }
  }while(n < 50000);

  fclose(fp1);
  fclose(fp2);

  // sort and remove duplicates
  sort(big_file);

  // output to out_file
  fpo = fopen(out_file, "w");
  for(i = 0; i < n+6; i++) fprintf(fpo, "%s", big_file[i]);
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

void sort(char big_file[][81])
{
  int i, j;
  char buf1[81], buf2[81], buf3[81];
  double ssn1, ssn2;

  sort:
  j = 0;
  for(i = 0; i <= n; i = i + 3)
  {
    ssn1 = atoi(big_file[i + 1] + 2);
    ssn2 = atoi(big_file[i + 4] + 2);

    if(ssn1 > ssn2)    // switch low/high
    {
      buf1 = big_file[i];
      buf2 = big_file[i + 1];
      buf3 = big_file[i + 2];
      big_file[i]     = big_file[i + 3];
      big_file[i + 1] = big_file[i + 4];
      big_file[i + 2] = big_file[i + 5];
      big_file[i + 3] = buf1;
      big_file[i + 4] = buf2;
      big_file[i + 5] = buf3;
      j++;
    }
    if(ssn1 == ssn2)   // remove duplicates
    {
      for(int s = i; s <= n ; s = s + 3)
      {
        big_file[s]     = big_file[s + 3];
        big_file[s + 1] = big_file[s + 4];
        big_file[s + 2] = big_file[s + 5];
      }
      n = n - 3;
    }
  }
  if(j > 0) goto sort;
}

