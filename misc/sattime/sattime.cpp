// SatTime.cpp
// by Scott Campbell campbel7@the-i.net


#include <cstdlib>  // for atof
#include "date.h"   // for Date, Locate
using namespace std;

//////////// declare functions ////////////////////////////////////////////////

inline char *s_in(char *prompt, char *buffer)    // input utility
{
  printf("%s", prompt);
  if (feof(stdin)) exit(0);
  gets(buffer);
  return buffer;
}

/////////////////// MAIN //////////////////////////////////////////////////////
int main()
{
  char buf[25];
  s_in("Enter Calendar Date [Y/N] : ", buf);
  printf("\n");
  if(buf[0] == 'y')
  {
     Date t1;
     printf("\n");
     t1.input();
     printf("\n\nJulian Date   : %f", t1.jd);
     printf("\nTLE Date      : %014.8f", t1.tle);
     printf("\nSidereal Date : %f", t1.thetag);
     printf("\nDay of Year   : %d", t1.doy);
  }
  else
  {
     s_in("Julian or TLE Date : ", buf);
     double time = atof(buf);
     Date t1(time);
     printf("\n");
     t1.print();
     printf("\n\nJulian Date   : %f", t1.jd);
     printf("\nTLE Date      : %014.8f", t1.tle);
     printf("\nSidereal Date : %f", t1.thetag);
     printf("\nDay of Year   : %d", t1.doy);
  }
  s_in("\n\n[Exit]", buf);
  exit(0);
}

////////////////// end main ///////////////////////////////////////////////////

