
#include "elcor.h"   /* global header */

/*
file stations.in:

0001 MM   30.334   -97.761    160.    5105 Crestway Dr.
0002 MM  30.3138  -97.8661    280.    Bee Caves Rsrch Ctr
2420 RE  55.9486   -3.1386     40.    Russell Eberst
2675 DB  52.1358   -2.3264     70.    David Brierley (2675)
2018 PW  51.0945   -1.1188    150.    Peter Wakelin (2018)
*/

/* Load file of sites */
void getsites(void)
{
   FILE *fp;
   char inp_str[81];

   if((fp = fopen("stations.in", "r")) == NULL)
   {
      printf("can't open file stations.in\n");
      s_in("", buf);
      exit(2);
   }

   num_sites = 0;

   while(fgets(inp_str, 80, fp))
   {
      /* printf("%s", inp_str); */

      if (num_sites >= MAXSITES) {
         printf("sites file too many sites error (limit 50)\n");
         s_in("", buf);
         exit(2);
      }

      if (sscanf(inp_str, "%d %2s %lf %lf %lf",
              &sitenum[num_sites], &siteabbr[num_sites],
              &xlat[num_sites], &xlong[num_sites],
              &xhgt[num_sites]) < 5) {
         printf("sites file error in site number lat long height\n");
         s_in("", buf);
         exit(2);
      }
      /* printf("%d ->%s<-\n",
              sitenum[num_sites], siteabbr[num_sites]); */
      /* printf("%d <- %.3f\n",
              sitenum[num_sites], xlat[num_sites]); */

      /* Increment */
      num_sites++;
   }
   fclose(fp);
}
