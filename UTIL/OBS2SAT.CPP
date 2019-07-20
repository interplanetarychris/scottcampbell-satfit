// OBS2SAT.cpp
// by Scott Campbell campbel7@hughes.net

/******************************************************************************
*                                                                             *
*                                                                             *
*                                                                             *
*******************************************************************************/
#include <cstdlib>
#include <cstdio>
#include <ctype.h>
using namespace std;


/////////////////// MAIN //////////////////////////////////////////////////////
main(int argc, char *argv[])
{
  char file_in[20], file_out[20], buf[81];
  char ssn[6], obs[3000][160], iod_line[3000][81];
  int sy, sn, nobs = 0;

  // UK variables
  char line[81];
  char ra[7], dc[7];
  int i, site = 0;
  int type, sd, dy, dn, ep;
  char date[15], dx[1], dl[1];
  char sufx[] = "ABCDEFGHJKLMNPQ";

  // RDE variables
  int yy, mm, dd, ymd, ss, hms;
  char xa[3], xb[5], ya[3], yb[5];


  if(argc == 2)
    sprintf(file_in, argv[1]);
  else          // default input files
    sprintf(file_in, "obs.in");

  // load obs from input file_in

  FILE *fp, *fpw, *fpx;

  if ((fp = fopen(file_in, "r")) == NULL)
  {
    printf("can't open %s\n", file_in);
    system("PAUSE");
    exit(0);
  }

  // obs reader

  while(fgets(obs[nobs], 160, fp))
  {
    sscanf(obs[0], "%4d", &sy);
    if(sy == 2420) goto rde;

    if(isdigit(obs[nobs][0]) &&
       isdigit(obs[nobs][1]) &&
       isdigit(obs[nobs][2]) &&
       isdigit(obs[nobs][3]) &&
       isdigit(obs[nobs][4])    ) nobs++;
  }
  fclose(fp);

  if(isdigit(obs[0][5])) goto uk;

  // IOD OBS
  // scan ssr.id for NORAD Catalog Number

  fpx = fopen("obs.out", "a");
  for(int i = 0; i < nobs; i++)
  {
     fp = fopen("ssr.id", "r");
     sscanf(obs[i], "%5s", &ssn);
     sprintf(file_out, "%5s.txt", ssn);
     sn = atoi(ssn);
     while(fgets(buf, 20, fp))
     {
        sscanf(buf, "%5d", &sy);
        if(sn == sy)
        {
           fclose(fp);
           sprintf(file_out, "%5s.txt", ssn);
           fp = fopen(file_out, "a");
           printf("%s", obs[i]);
           fprintf(fp, "%s", obs[i]);
           fprintf(fpx, "%s", obs[i]);
           fclose(fp);
           break;
        }
     }
  }
  fclose(fpx);
  printf("\nDone\n");
  system("PAUSE");
  exit(0);
  // END IOD OBS

  // BEGIN UK OBS
  uk:
  for(i = 0; i < nobs; i++)
  {
    sscanf(obs[i], "%2d", &sy);
    sscanf(obs[i] + 2, "%3d", &sn);
    sscanf(obs[i] + 5, "%2d", &sd);
    sd -= 1;
    sscanf(sufx + sd, "%1s", &dl);
    sscanf(obs[i] + 7, "%4d", &site);
    sscanf(obs[i] + 11, "%s", &date);
    sscanf(obs[i] + 33, "%1d", &type);
    sscanf(obs[i] + 34, "%s", &ra);
    sscanf(obs[i] + 42, "%s", &dc);
    sscanf(obs[i] + 54, "%1d", &ep);
    sprintf(iod_line[i], "%02d %.3d%1s   %4d G ",
                           sy,   sn,dl,  site);
    if(!isdigit(date[13]))
    {
      sprintf(iod_line[i], "%s20%-13s00 17 %1d%1d %-6s0",
              iod_line[i],    date,     type,ep,  ra);
    }
    else if(isdigit(date[14]))
    {
      sprintf(iod_line[i], "%s20%-14s 17 %1d%1d %-6s0",
              iod_line[i],    date,     type,ep,  ra);
    }
    else
    {
      sprintf(iod_line[i], "%s20%-14s0 17 %1d%1d %-6s0",
              iod_line[i],    date,     type,ep,  ra);
    }
    if(!isdigit(dc[5]))
    {
      sprintf(dc, "%-5s00", dc);
    }
    if(!isdigit(dc[6]))
    {
      sprintf(dc, "%-6s0", dc);
    }
    sprintf(iod_line[i], "%s%-6s 57", iod_line[i], dc);
  }


  // check for existence of ssr.id file
  if ((fp = fopen("ssr.id", "r")) == NULL)
  {
    printf("can't open %s\n", "ssr.id");
    system("PAUSE");
    exit(0);
  }
  fclose(fp);


  // scan ssr.id for NORAD Catalog Number

  fpx = fopen("obs.out", "a");
  for(int i = 0; i < nobs; i++)
  {
     sscanf(iod_line[i], "%2d", &dy);
     sscanf(iod_line[i] + 3, "%3d", &dn);
     sscanf(iod_line[i] + 6, "%1s", &dl);
     fp = fopen("ssr.id", "r");
     while(fgets(buf, 20, fp))
     {
        sscanf(buf, "%5s", &ssn);
        sscanf(buf + 6, "%2d", &sy);
        sscanf(buf + 9, "%3d", &sn);
        sscanf(buf + 12, "%1s", &dx);

        if(sy == dy && sn == dn && dx[0] == dl[0])
        {
           printf("%s", ssn);
           printf(" %s\n", iod_line[i]);
           sprintf(file_out, "%5s.txt", ssn);
           fpw = fopen(file_out, "a");
           fprintf(fpw, "%s", ssn);
           fprintf(fpw, " %s\n", iod_line[i]);
           fprintf(fpx, "%s", ssn);
           fprintf(fpx, " %s\n", iod_line[i]);
           fclose(fpw);
           break;
        }
     }
     fclose(fp);
  }
  fclose(fpx);
  printf("\nDone\n");
  system("PAUSE");
  exit(0);

  // END UK OBS


  // BEGIN RDE OBS
  rde:
  sscanf(obs[0] + 5, "%2d", &yy);
  sscanf(obs[0] + 7, "%2d", &mm);
  yy += 2000;                      // change 01/01/2100


  while(fgets(line, 80, fp))       // find observation entries
  {
    if(strlen(line) < 4 && atoi(line) > 0 && atoi(line) < 32)
    {
       sscanf(line, "%2d", &dd);
       ymd = yy * 10000 + mm * 100 + dd;
       fgets(line, 80, fp);
    }
    if(strlen(line) >= 44)
    {
      sscanf(line, "%2d", &sy);
      sscanf(line + 2, "%3d", &sn);
      sscanf(line + 5, "%2d", &sd);
      sd -= 1;
      sscanf(sufx + sd, "%1s", &dl);
      sscanf(line + 8, "%6d", &hms);
      sscanf(line + 15, "%2d", &ss);
      sscanf(line + 18, "%6s", &ra);
      sscanf(line + 24, "%7s", &dc);
      hms = hms * 100 + ss;
      sprintf(iod_line[nobs], "%02d %.3d%1s   2420 G %8d%.8d0 17 14 %6s0%7s 57",
                                 sy, sn, dl,         ymd, hms,       ra, dc);
      nobs++;
    }
    if(atoi(line) == 999) break;
  }
  fclose(fp);

  // scan ssr.id for NORAD Catalog Number

  fpx = fopen("obs.out", "a");
  for(int i = 0; i < nobs; i++)
  {
     sscanf(iod_line[i], "%2s", &xa);
     sd = atoi(xa);
     sscanf(iod_line[i] + 3, "%4s", &xb);
     dn = atoi(xb);
     fp = fopen("ssr.id", "r");
     while(fgets(buf, 80, fp))
     {
        sscanf(buf, "%5s", &ssn);
        sscanf(buf + 6, "%2s", &ya);
        sy = atoi(ya);
        sscanf(buf + 9, "%4s", &yb);
        sn = atoi(yb);
        if(sd == sy && dn == sn && xb[3] == yb[3])
        {
           printf("%s", ssn);
           printf(" %s\n", iod_line[i]);
           sprintf(file_out, "%5s.txt", ssn);
           fpw = fopen(file_out, "a");
           fprintf(fpw, "%s", ssn);
           fprintf(fpw, " %s\n", iod_line[i]);
           fprintf(fpx, "%s", ssn);
           fprintf(fpx, " %s\n", iod_line[i]);
           fclose(fpw);
           break;
        }
     }
     fclose(fp);
  }
  fclose(fpx);
  printf("\nDone\n");
  system("PAUSE");
  exit(0);

  // END RDE OBS


}

////////////////// end main ///////////////////////////////////////////////////

