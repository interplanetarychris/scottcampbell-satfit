// dpfit.cpp
#include <ctype.h>      // for isdigit in read_tle
#include <string.h>     // for command line input
#include "orbit.h"      // for Date, Locate, Satellite
using namespace std;

///////////// DECLARE GLOBAL VARIABLES /////////////////////////////////////////

double
   tle,      // epoch of elements in tle format
   ii,       // inclination, degrees
   om,       // right ascension of ascending node, degrees
   ec,       // eccentricity
   ww,       // argument of the perigee, degrees
   uu,       // longitude
   ma,       // mean anomaly, degrees
   nn,       // mean motion, revolutions/day
   c2,       // internal drag term
   bstar;    // BSTAR drag term

// please note that these variables are in TLE format, i.e. degrees, whereas
// their counterparts inside SGP4 are in radians.

int dopp = 0;               // doppler min/mx flag
int pass = 0;               // sim pass flag
int Z_flag = 0;             // AOS/LOS included
int off = 0, ofs;           // frequency offset
int batch = 0, out = 0;
int nobs, aobs, site;       // number of observations, site code
double la, lo, hh;          // observer parameters
double s_data[50][3];       // observer site data
double tol[50];
char file[81];              // input file string declared
char buf[81];               // input buffer
char name[81], tle1[81], tle2[81];    // original tle data lines
char line1[81], line2[81];            // read and write tle buffers
char iod_line[1000][81];
char aos_los[100][81];
double sum, nsum, xsum;     // sum of squares of residuals

// declare variables for search option
double e, ik, ok, ek, wk, mk, nk, bk, theta, ra, dc;

double
  astep,   amin,    amax,
  istep,   imin,    imax,
  ostep,   omin,    omax,
  estep,   emin,    emax,    esize,
  wstep,   wmin,    wmax,    wsize,
  nstep,   nmin,    nmax,    nsize,
  bstep,   bmin,    bmax,
           mmax,    mmin;

//////////// DECLARE FUNCTIONS ////////////////////////////////////////////////

void get_obs(char *file, char iod_line[][81]);
void read_obs(char iod_line[][81],
            double odata[][4],
            double rd[][3]);
void read_dopp_obs(char iod_line[][81],
            double odata[][4],
            double rd[][3]);
void sort(char iod_line[][81],
            double odata[][4],
            double rd[][3]);
void read_tle(char *file);
void read_site(char line[]);
void print_fit(Satellite sat,  double rd[][3], double odata[][4]);
double find_sum(Satellite sat,  double rd[][3], double odata[][4]);
void anomaly(Satellite sat,  double rd[][3], double odata[][4]);
void motion(Satellite sat,  double rd[][3], double odata[][4]);
void perigee(Satellite sat,  double rd[][3], double odata[][4]);
void node(Satellite sat,  double rd[][3], double odata[][4]);
void align(Satellite sat,  double rd[][3], double odata[][4]);
double find_dopp(Satellite sat,  double rd[][3], double odata[][4]);
void do_sim(Satellite sat);
void list(Satellite sat);
void longitude(void);
double ww2ma(double wx);
void write_tle(char *file);
inline double rtw(double ao, double ac);
inline char *s_in(char *prompt, char *buffer);
inline double mod(double var, double div);

/////////////////// MAIN //////////////////////////////////////////////////////

int main(int argc, char *argv[])           // command line args for read data
{
  char file_in[81];
  int i, j, ct;                            // counter and check variables
  double max, min, minsum;
  char srch[] = {'W'};                     // initialize wide search
  FILE *fp;

  if(argc == 1)
    sprintf(file, "dp.obs");               // default input/output file
  if(argc == 2)                            // command line file or default?
    sprintf(file, argv[1]);                // command line input/output file
  if (argc == 3 && strcmp(argv[2],"ba") == 0)
  {
    batch = 1;
    sprintf(file, argv[1]);
  }
  if(strchr(file, '\\'))
    sprintf(file, "%s", strrchr(file, '\\')+1);

  // read first tle from input file
  // orbit parameters are passed from read_tle as global variables
  // echo tle to screen
  read_tle(file);
  sprintf(tle1, "%s", line1);  // copy original line1 to tle1
  sprintf(tle2, "%s", line2);  // copy original line2 to tle2

  // Create Satellite variables, sat and save_sat, from TLE
  Satellite sat(tle, ii, om, ec, ww, ma, nn, bstar);
  Satellite save_sat(tle, ii, om, ec, ww, ma, nn, bstar);

  // calculate uu, degrees, for search
  longitude();

  // find observation lines in input file
  // number of observations = nobs (global)
  // store observations to iod_line matrix (up to 200 observations)
  get_obs(file, iod_line);
  printf("%d Observtions Found\n",nobs);

  double rd[aobs+nobs][3];     // topocentric vectors to observer positions
  double odata[aobs+nobs][4];  // observational data: jd, fo, df, obscode

  // read all observations from input
  // Observed topocentric vectors for each observer position.
  // Equatorial coordinates.
  if(nobs > 0)
  {
    if(dopp) read_dopp_obs(iod_line, odata, rd);
    else read_obs(iod_line, odata, rd);
    sort(iod_line, odata, rd);
  }

  // print dates
  Date t2;
  printf("\nTODAY   : %d", t2.doy);
  if(nobs > 0)
  {
    Date t3(odata[nobs - 1][0]);
    printf("\nLAST OB : %d\n", t3.doy);
    if(!dopp) goto offset;
  }

accept:
   // Accept a new command

   printf("\nEnter command\n");
   printf(" (I)ncl  (N)ode  (X)ntrcty  (P)erigee  (A)nomaly   (M)otion  (B)star\n");
   printf(" (S)tep  (O)fft  (T)ime  (F)it (V)iew  (U)Simulate (W)rite   (Q)uit\n");

   s_in(": ", buf);

   if ((*buf & 0x5f) == 'G') goto fit_out;           // hidden
   if ((*buf & 0x5f) == 'E') goto edit_file;         // hidden
   if ((*buf & 0x5f) == 'H') goto history;           // hidden
   if ((*buf & 0x5f) == 'Y') goto edit_history;      // hidden
   if ((*buf & 0x5f) == 'R') goto remove;            // hidden
   if ((*buf & 0x5f) == 'D') goto discover;          // hidden
   if ((*buf & 0x5f) == 'S') goto step;
   if ((*buf & 0x5f) == 'O')
   {
      if(nobs > 0) goto offset;
      else printf("\nNO OBS TO OFFSET\n");
   }
   if ((*buf & 0x5f) == 'I')
   {
      sat.print_el();
      goto incl;
   }
   if ((*buf & 0x5f) == 'N')
   {
      sat.print_el();
      goto node;
   }
   if ((*buf & 0x5f) == 'X')
   {
      sat.print_el();
      goto xntrcty;
   }
   if ((*buf & 0x5f) == 'P')
   {
      sat.print_el();
      goto perigee;
   }
   if ((*buf & 0x5f) == 'A')
   {
      sat.print_el();
      goto anomaly;
   }
   if ((*buf & 0x5f) == 'M')
   {
      sat.print_el();
      goto motion;
   }
   if ((*buf & 0x5f) == 'B')
   {
      sat.print_el();
      goto bstar;
   }
   if ((*buf & 0x5f) == 'T') goto time;
   if ((*buf & 0x5f) == 'F') goto fit;
   if ((*buf & 0x5f) == 'V') goto viewobs;
   if ((*buf & 0x5f) == 'W') goto write_el;
   if ((*buf & 0x5f) == 'U') goto simulate;
   if ((*buf & 0x5f) == 'Q') exit(0);

   // try again
   goto accept;

fit_out:
   out = 1;
   print_fit(sat, rd, odata);
   out = 0;
   goto accept;

edit_file:
   sprintf(file, "%s", file);
   if(ssn < 10000) file[0] = '0';
   system(file);
   goto accept;

edit_history:
   sprintf(file_in, "history%c%5d.txt", '\\', ssn);
   if(ssn < 10000) file_in[8] = '0';
   system(file_in);
   goto accept;

history:
   sprintf(file_in, "history/%5d.txt", ssn);
   if(ssn < 10000) file_in[8] = '0';
   if ((fp = fopen(file_in, "r")) == NULL)
   {
     printf("NO HISTORY FOR THIS SAT");
     goto write_el;
   }
   while(fgets(buf, 80, fp))   // get a line from ssn.txt
     printf("%s", buf);
   fclose(fp);
   printf("\nCURRENT TLE");
   write_tle("");
   goto accept;

viewobs:
   printf("\n");
   for(i = 0; i < nobs; i++)
   {
      printf("(%2d) %s",i + 1, iod_line[i]);
      if (mod(i+1, 24) == 0 && i > 0 && i+1 < nobs) s_in("[more]", buf);
   }
   if(mod(nobs, 24) == 21) printf("\n");
   if (mod(nobs, 24) > 19 || mod(nobs, 24) == 0) s_in("[continue]", buf);
   printf("\n");
   goto accept;

remove:
   s_in("\nRemove Line Number : ", buf);
   j = atoi(buf);
   for(i = j; i < aobs+nobs; i++)
   {
      rd[i - 1] = rd[i];
      odata[i - 1] = odata[i];
      iod_line[i - 1] = iod_line[i];
      s_data[i - 1] = s_data[i];
      tol[i - 1] = tol[i];
   }
   if(j <= nobs) nobs--;
   else aobs--;

fit:
   out = 0;
   print_fit(sat, rd, odata);
   sat.print_el();
   goto accept;

discover:       // partition search

   printf("\n(W)ide  (N)arrow  [%c", srch[0]);
   if (strlen(s_in("]: ", buf))) srch[0] = buf[0];
   if ((*srch & 0x5f) == 'W')
   {
      emin = 0;             emax = .2;            esize = 80;
      wmin = 0;             wmax = 360;           wsize = 144;
      nmin = 11;            nmax = 15;            nsize = 80;
   }
   if ((*srch & 0x5f) == 'N')
   {
      emin = ec * .9;       emax = ec * 1.1;      esize = 20;
      wmin = ww - 2;        wmax = ww + 2;        wsize = 20;
      nmin = nn - .1;       nmax = nn + .1;       nsize = 20;
   }
   if ((*srch & 0x5f) == 'Z')
   {
      emin = ec - estep;    emax = ec + estep;    esize = 20;
      wmin = ww - wstep;    wmax = ww + wstep;    wsize = 20;
      nmin = nn - nstep;    nmax = nn + nstep;    nsize = 20;
   }

   minsum = 1.e6;
   printf("\nemax [%.7f", emax);
   if (strlen(s_in("]: ", buf)))
     emax = atof(buf);
   printf("emin [%.7f", emin);
   if (strlen(s_in("]: ", buf)))
     emin = atof(buf);
   estep = (emax - emin) / esize;

   printf("\nwmax [%.4f", wmax);
   if (strlen(s_in("]: ", buf)))
     wmax = atof(buf);
   printf("wmin [%.4f", wmin);
   if (strlen(s_in("]: ", buf)))
     wmin = atof(buf);
   wstep = (wmax - wmin) / wsize;

   printf("\nnmax [%.8f", nmax);
   if (strlen(s_in("]: ", buf)))
     nmax = atof(buf);
   printf("nmin [%.8f", nmin);
   if (strlen(s_in("]: ", buf)))
     nmin = atof(buf);
   nstep = (nmax - nmin) / nsize;

   for(wk = wmin; wk < wmax; wk += wstep)
   {
      theta = (uu - wk) * de2ra;
      for(ek = emin; ek < emax; ek += estep)
      {
         e = acos((ek + cos(theta)) / (1 + ek * cos(theta)));
         if(theta > pi) e = 2 * pi - e;
         mk = e - ek * sin(e);
         mk = mk / de2ra;
         for(nk = nmin; nk < nmax; nk += nstep)
         {
            Satellite satx(sat.jd, ii, om, ek, wk, mk, nk, bstar);
            // establish the computed doppler at jdo with no perturbations
            sum = find_sum(satx, rd, odata);
            if(sum < minsum)
            {
               minsum = sum;
               ec = ek;
               ww = wk;
               ma = mk;
               nn = nk;
            }
         }   // for nk
      }   // for ek
   }   // for wk
   ww = mod(ww, 360);
   ma = mod(ma, 360);
   sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar);

   // update mean_anomaly
   anomaly(sat, rd, odata);
   sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar);

   // print new elements
   print_fit(sat, rd, odata);
   sat.print_el();

   longitude();

   srch[0] = 'Z';
   goto accept;

step:       // partition search within limits set below

   // update mean_anomaly
   anomaly(sat, rd, odata);
   motion(sat, rd, odata);
   sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar);

   emax = ec * 1.1;
   emin = ec * .9;
   wmax = ww + 2;
   wmin = ww - 2;
   imax = ii;            // no change on first pass
   imin = ii;            //
   omax = om;            //
   omin = om;            //

   printf("\nPress Q to Quit    :\n\n");
   while(1)
   {
      perigee(sat, rd, odata);
      sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar);
      node(sat, rd, odata);
      sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar);

      // set new limits
      emax = 1.01 * ec;
      emin = .99 * ec;
      wmax = ww + .5;
      wmin = ww - .5;
      imax = ii + .5;
      imin = ii - .5;
      omax = om + .5;
      omin = om - .5;

      printf("sum%12.5f", sum);
      s_in("    : ", buf);
      if ((*buf & 0x5f) == 'Q') break;
   }


   print_fit(sat, rd, odata);
   sat.print_el();       // print new elements

   srch[0] = 'N';
   goto accept;

incl:

   printf("\n(A)uto  (ii)  (Q)uit  ");
   s_in(": ", buf);

   if (isdigit(buf[0]))
   {
     ii = atof(buf);
     sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar);
     print_fit(sat, rd, odata);
     sat.print_el();       // print new elements
   }
   if ((*buf & 0x5f) == 'A')
   {
     // partition search

     imax = ii + 2;
     imin = ii - 2;
     omax = om + 2;
     omin = om - 2;

     printf("\nimax [%.4f", imax);
     if (strlen(s_in("]: ", buf)))
       imax = atof(buf);
     printf("imin [%.4f", imin);
     if (strlen(s_in("]: ", buf)))
       imin = atof(buf);

     printf("\nomax [%.4f", omax);
     if (strlen(s_in("]: ", buf)))
       omax = atof(buf);
     printf("omin [%.4f", omin);
     if (strlen(s_in("]: ", buf)))
       omin = atof(buf);

     node(sat, rd, odata);

     sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar);
     print_fit(sat, rd, odata);
     sat.print_el();       // print new elements

     srch[0] = 'N';
   }
   if ((*buf & 0x5f) == 'Q') goto accept;
   goto incl;

node:

   printf("\n(A)uto  (om)  (Q)uit  ");
   s_in(": ", buf);

   if (isdigit(buf[0]))
   {
     om = atof(buf);
     sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar);
     print_fit(sat, rd, odata);
     sat.print_el();       // print new elements
   }
   if ((*buf & 0x5f) == 'A')
   {
     Satellite satx = sat;

     // partition search
     omax = om + 2;
     omin = om - 2;

     printf("\nomax [%.4f", omax);
     if (strlen(s_in("]: ", buf)))
       omax = atof(buf);
     printf("omin [%.4f", omin);
     if (strlen(s_in("]: ", buf)))
       omin = atof(buf);

     minsum = 1.e6;
     while((omax - omin) > 1e-5)
     {
       ostep = fabs(rtw(omax, omin) / 20);
       for(ok = omin; ok < omax; ok += ostep)
       {
          satx.delta_el(sat.jd, ii, ok, ec, ww, ma, nn, bstar);
          // establish the computed ra, dc, at jdo with no perturbations
          sum = find_sum(satx, rd, odata);
          if(sum < minsum)
          {
             minsum = sum;
             om = ok;
          }
       }   // for ok
       satx.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar);

       omin = om - ostep;
       omax = om + ostep;
     }
     om = mod(om, 360);

     sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar);
     print_fit(sat, rd, odata);
     sat.print_el();       // print new elements

     srch[0] = 'N';
   }
   if ((*buf & 0x5f) == 'Q') goto accept;
   goto node;

xntrcty:

   printf("\n(S)earch  (A)uto  (ec)  (Q)uit  ");
   s_in(": ", buf);

   emax = ec * 1.1;
   emin = ec * .9;
   minsum = 1.E6;;

   if (isdigit(buf[0]))
   {
     ec = atof(buf);
     sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar);
     print_fit(sat, rd, odata);
     sat.print_el();       // print new elements
   }
   if ((*buf & 0x5f) == 'A')
   {
     printf("\nemax [%.7f", emax);
     if (strlen(s_in("]: ", buf)))
       emax = atof(buf);
     printf("emin [%.7f", emin);
     if (strlen(s_in("]: ", buf)))
       emin = atof(buf);

     while((emax - emin) > 1.e-8)
     {
        estep = (emax - emin) / 20;
        for(ek = emin; ek < emax; ek += estep)
        {
           sat.delta_el(sat.jd, ii, om, ek, ww, ma, nn, bstar);
           // establish the computed ra, dc, at jdo with no perturbations
           sum = find_sum(sat, rd, odata);
           if(sum < minsum)
           {
              minsum = sum;
              ec = ek;
           }
        }   // for ek
        sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar);
        emin = ec - estep;
        emax = ec + estep;
     }

     sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar);
     print_fit(sat, rd, odata);
     sat.print_el();       // print new elements
     srch[0] = 'N';
   }
   if ((*buf & 0x5f) == 'S')
   {
     printf("\nemax [%.7f", emax);
     if (strlen(s_in("]: ", buf)))
       emax = atof(buf);
     printf("emin [%.7f", emin);
     if (strlen(s_in("]: ", buf)))
       emin = atof(buf);
     printf("estep [%.0f", estep);
     if (strlen(s_in("]: ", buf)))
       estep = atof(buf);
     estep = (emax - emin) / estep;

     ek = ec;
     printf("\neccentricity    sum");
     for(ec = emin; ec < emax + estep; ec += estep)
     {
       sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar);
       sum = find_sum(sat, rd, odata);
       printf("\n%.7f     %7.4f", ec, sum);
     }
     printf("\n");
     ec = ek;              // restore
     sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar);
   }
   if ((*buf & 0x5f) == 'Q') goto accept;
   goto xntrcty;

perigee:

   printf("\n(S)earch  (A)uto  (ww)  (Q)uit  ");
   s_in(": ", buf);

   if (isdigit(buf[0]))
   {
     longitude();
     ww = atof(buf);
     ma = ww2ma(ww);
     sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar);
     print_fit(sat, rd, odata);
     sat.print_el();       // print new elements
   }
   if ((*buf & 0x5f) == 'A')
   {
     double wx = ww + 1;
     while(fabs(wx - ww) > 1e-4)
     {
       wx = ww;
       printf("\n%8.4f  %7.4f", ww, ec);
       wmax = ww + 2;
       wmin = ww - 2;
       emax = ec * 1.01;
       emin = ec * .99;
       perigee(sat, rd, odata);
       sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar);
     }
     print_fit(sat, rd, odata);
     sat.print_el();       // print new elements
     srch[0] = 'N';
   }
   if ((*buf & 0x5f) == 'S')
   {
     longitude();
     wmax = ww + 2;
     wmin = ww - 2;
     wstep = 20;

     printf("\nwmax [%.4f", wmax);
     if (strlen(s_in("]: ", buf)))
       wmax = atof(buf);
     printf("wmin [%.4f", wmin);
     if (strlen(s_in("]: ", buf)))
       wmin = atof(buf);
     printf("wstep [%.0f", wstep);
     if (strlen(s_in("]: ", buf)))
       wstep = atof(buf);
     wstep = (wmax - wmin) / wstep;

     wk = ww;
     mk = ma;
     printf("\nperigee        sum");
     for(ww = wmin; ww < wmax + wstep; ww += wstep)
     {
       ma = ww2ma(ww);
       sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar);
       sum = find_sum(sat, rd, odata);
       printf("\n%8.4f     %7.4f", mod(ww, 360), sum);
     }
     printf("\n");
     ww = wk;
     ma = mk;
     sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar);
   }
   if ((*buf & 0x5f) == 'Q') goto accept;
   goto perigee;

offset:

   printf("\n(offset)  (A)uto  (Q)uit  ");
   s_in(": ", buf);

   if (isdigit(buf[1]))
   {
     ofs = atoi(buf);
     for(int j = 0; j < nobs; j++)
     {
       odata[j][2] += ofs;
     }
     print_fit(sat, rd, odata);
     off += ofs;
     printf("\ntotal offset :  %5d\n", off);
     sat.print_el();       // print new elements
   }
   if ((*buf & 0x5f) == 'A')
   {
     double zum;
     int step;

     // determine direction
     zum = find_sum(sat, rd, odata);   // sum at zero offset
     for(int j = 0; j < nobs; j++)
       odata[j][2] += 100.;
     ofs = 100;
     sum = find_sum(sat, rd, odata);   // +100
     if(sum < zum)
     {
       step = 100;
       goto loop;
     }
     for(int j = 0; j < nobs; j++)
       odata[j][2] -= 200.;
     ofs = -100;
     sum = find_sum(sat, rd, odata);   // -100
     if(sum < zum)
     {
       step = -100;
       goto loop;
     }
     step = 10;
     loop:
     do
     {
       zum = sum;
       for(int j = 0; j < nobs; j++)
         odata[j][2] += step;
       ofs += step;
       sum = find_sum(sat, rd, odata);
     }while(sum < zum);
     if(abs(step) != 1)
     {
       step = (int)(step / -10);
       goto loop;
     }
     // back up one
     for(int j = 0; j < nobs; j++)
       odata[j][2] -= step;
     ofs -= step;
     off += ofs;
     printf("\ntotal offset :  %5d\n", off);
     s_in("", buf);
     print_fit(sat, rd, odata);
     sat.print_el();       // print new elements
   }
   if ((*buf & 0x5f) == 'Q') goto accept;
   goto offset;

e_longitude:

   // amax and amin are used as temporary variables
   longitude();
   amax = uu * de2ra;
   amin = sin(sat.xincl/de2ra) * sin(amax);
   amin *= amin;
   amin = sqrt(1 - amin);
   amin = (acos(cos(amax) / amin)) / de2ra;
   if(fmod(uu, 360) > 180) amin = 360 - amin;
   if(sat.xincl/de2ra < 90)
        amax = fmod(sat.xnodeo/de2ra + amin - sat.thetag + 360, 360.);
   else amax = fmod(sat.xnodeo/de2ra - amin - sat.thetag + 720, 360.);
   printf("\nE Longitude = %8.4f\n", amax);

anomaly:

   printf("\n(S)earch  (A)uto  (ma)  (L)ong  (Q)uit  ");
   s_in(": ", buf);

   amax = 360.;
   amin = 0.;
   astep = 20;

   if (isdigit(buf[0]))
   {
     ma = atof(buf);
     longitude();
     sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar);
     print_fit(sat, rd, odata);
     sat.print_el();       // print new elements
     goto e_longitude;
   }
   if ((*buf & 0x5f) == 'A')
   {
     anomaly(sat, rd, odata);
     sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar);
     print_fit(sat, rd, odata);
     sat.print_el();       // print new elements
     goto e_longitude;

     srch[0] = 'N';
   }
   if ((*buf & 0x5f) == 'S')
   {
     printf("\namax [%.7f", amax);
     if (strlen(s_in("]: ", buf)))
       amax = atof(buf);
     printf("amin [%.7f", amin);
     if (strlen(s_in("]: ", buf)))
       amin = atof(buf);
     printf("astep [%.0f", astep);
     if (strlen(s_in("]: ", buf)))
       astep = atof(buf);
     astep = (amax - amin) / astep;

     mk = ma;
     printf("\nanomaly        sum");
     for(ma = amin; ma < amax + astep; ma += astep)
     {
       sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar);
       sum = find_sum(sat, rd, odata);
       printf("\n%8.4f     %7.4f", ma, sum);
     }
     printf("\n");
     ma = mk;              // restore
     sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar);
   }
   if ((*buf & 0x5f) == 'L') goto e_longitude;
   if ((*buf & 0x5f) == 'Q') goto accept;
   goto anomaly;

motion:

   printf("\n(A)uto  (nn)  (Q)uit  ");
   s_in(": ", buf);

   if (isdigit(buf[0]))
   {
     nn = atof(buf);
     sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar);
     print_fit(sat, rd, odata);
     sat.print_el();       // print new elements
   }
   if ((*buf & 0x5f) == 'A')
   {
     // update mean motion, no limits
     motion(sat, rd, odata);
     sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar);
     print_fit(sat, rd, odata);
     sat.print_el();       // print new elements
   }
   if ((*buf & 0x5f) == 'Q') goto accept;
   goto motion;

bstar:

   printf("\n(A)uto  (b*)  (Q)uit  ");
   s_in(": ", buf);

   if (isdigit(buf[0]))
   {
     bstar = atof(buf);
     sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar);
     print_fit(sat, rd, odata);
     sat.print_el();       // print new elements
   }
   if ((*buf & 0x5f) == 'A')
   {
     minsum = 1.E6;
     // update Bstar within limits
     bmax = bstar * 1.1;
     bmin = bstar * .9;
     if(bstar < 0)
     {
       bmax = bmin;
       bmin = bstar * 1.1;
     }
     printf("\nbmax [%.8f", bmax);
     if (strlen(s_in("]: ", buf)))
       bmax = atof(buf);
     printf("bmin [%.8f", bmin);
     if (strlen(s_in("]: ", buf)))
       bmin = atof(buf);

     while((bmax - bmin) > 1.e-9)
     {
        bstep = (bmax - bmin) / 20;
        for(bk = bmin; bk < bmax; bk += bstep)
        {
           sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bk);
           // establish the computed ra, dc, at jdo with no perturbations
           sum = find_sum(sat, rd, odata);
           if(sum < minsum)
           {
              minsum = sum;
              bstar = bk;
           }
        }   // for bk
        sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar);
        bmin = bstar - bstep;
        bmax = bstar + bstep;
     }
     sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar);
     print_fit(sat, rd, odata);
     sat.print_el();       // print new elements
   }
   if ((*buf & 0x5f) == 'Q') goto accept;
   goto bstar;

simulate:

   printf("\n(S)im  (P)ass  (F)ile  (Q)uit  ");
   s_in(": ", buf);

   // toggle
   if ((*buf & 0x5f) == 'F')
   {
     if(out == 0)
     {
       out = 1;
       printf("PRINT TO FILE AND SCREEN");
     }
     else
     {
       out = 0;
       printf("PRINT TO SCREEN ONLY");
     }
   }
   if ((*buf & 0x5f) == 'S')
   {
     pass = 0;
     do_sim(sat);
     sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar);
   }
   if ((*buf & 0x5f) == 'P')
   {
     pass = 1;
     do_sim(sat);
     sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar);
   }
   if ((*buf & 0x5f) == 'Q') goto accept;
   goto simulate;

write_el:
   tle = sat.tle;
   ii = sat.xincl / de2ra;
   om = sat.xnodeo / de2ra;
   ec = sat.eo;
   ww = sat.omegao / de2ra;
   ma = sat.xmo / de2ra;
   nn = sat.xno / nocon;
   c2 = sat.c2;
   bstar = sat.bstar;

   printf("\n(U)pdate (V)iew (A)ppend (O)riginal (R)estore (Q)uit");
   s_in(": ", buf);

   if ((*buf & 0x5f) == 'A') write_tle(file);
   if ((*buf & 0x5f) == 'V') sat.print_el();
   if ((*buf & 0x5f) == 'O')
   {
      // view original
      save_sat.print_el();       // print original elements
   }
   if ((*buf & 0x5f) == 'R')
   {
      // replace TLE in file with original
      FILE *fp;
      fp = fopen(file, "w");
      fprintf(fp, "%s", name);
      fprintf(fp, "%s", tle1);
      fprintf(fp, "%s\n", tle2);
      for(int i = 0; i < nobs+aobs; i++)
      {
         fprintf(fp, "%s", iod_line[i]);
      }
      fclose(fp);
      // replace working elements with original
      sat = save_sat;
      sat.print_el();            // print original elements

   }
   if ((*buf & 0x5f) == 'U')
   {
      write_tle("");             // udates lines and prints to screen
      //  print_file
      FILE *fp;
      fp = fopen(file, "w");
      fprintf(fp, "%s", name);
      fprintf(fp, "%s\n", line1);
      fprintf(fp, "%s\n\n", line2);
      for(int i = 0; i < nobs+aobs; i++)
      {
         fprintf(fp, "%s", iod_line[i]);
      }
      fclose(fp);
   }
   if ((*buf & 0x5f) == 'Q') goto accept;
   goto write_el;

time:
   double z1, z2, z3, xns, sinio, beta, time1, time2, time3;
   ec = sat.eo;
   ww = sat.omegao / de2ra;
   nn = sat.xno / nocon;
   xns = 2160 * sat.bstar * sat.c2 * sat.xno / nocon;
   if(nobs > 0) time2 = odata[nobs - 1][0];
   else
   {
      t2.now();
      time2 = t2.jd;
   }

   sat.delta_t(time2);
   z2 = sat.rr[2];
   do
   {
      time1 = time2;
      z1 = z2;
      time2 -= .01;
      sat.delta_t(time2);
      z2 = sat.rr[2];
   }while(z1 < 0 || z2 > 0);

   while(time1 - time2 > 1e-9)
   {
      time3 = (time1 + time2) / 2;
      sat.delta_t(time3);
      z3 = sat.rr[2];
      if(z3 < 0) time2 = time3;
      else time1 = time3;
   }
   Date t1(time2);
   t1.input();
   sat.delta_t(t1.jd);
   sat.rv2el(sat.rr, sat.vv);

   if(sat.xno / nocon < 6.4)
   {
      double rr2[3], vv2[3];
      rr2[0] = sat.rr[0];          // save state vectors
      rr2[1] = sat.rr[1];
      rr2[2] = sat.rr[2];
      vv2[0] = sat.vv[0];
      vv2[1] = sat.vv[1];
      vv2[2] = sat.vv[2];

      double
      ik = 2*sat.xincl,
      ok = 2*sat.xnodeo,
      ek = 2*sat.eo,
      wk = 2*sat.omegao,
      mk = 2*sat.xmo,
      nk = 2*sat.xno;

      sat.sxp4(0.0);              // sdp4 propagation
      sat.rv2el(sat.rr, sat.vv);  // sdp4 elements

      sat.xincl  = ik - sat.xincl;
      sat.xnodeo = ok - sat.xnodeo;
      sat.eo     = ek - sat.eo;
      sat.omegao = wk - sat.omegao;
      sat.xmo    = mk - sat.xmo;
      sat.xno    = nk - sat.xno;

      sat.rr[0] = rr2[0];     // reset state vectors from sdp4 propagation
      sat.rr[1] = rr2[1];
      sat.rr[2] = rr2[2];
      sat.vv[0] = vv2[0];
      sat.vv[1] = vv2[1];
      sat.vv[2] = vv2[2];
   }

   tle = t1.tle;
   sat.jd = t1.jd;

   ii = sat.xincl / de2ra;
   om = mod(sat.xnodeo / de2ra, 360);
   ec = sat.eo;
   ww = mod(sat.omegao / de2ra, 360);
   ma = mod(sat.xmo / de2ra, 360);
   nn = sat.xno / nocon;
   bstar = xns / (2160 * sat.c2 * nn);
   sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar);
   sat.print_el();       // print new elements
   sum = find_sum(sat, rd, odata);
   printf("sum%12.5f\n", sum);

   longitude();

   goto accept;

}    // end main

////////////////// FUNCTIONS //////////////////////////////////////////////////

/* Converts quasi scientific notation to double.  Format:  SnSn where S is
either '-' or '\0' and n represents a sequence of 1 or more digits.  An
implied decimal point exists to the left of the left n.  Total length 27
chars max. Subroutine of read_tle*/
double sci(char *string)
{
   char buf[30], *bufptr;

   bufptr = buf;

   if (string[1] == '\0')
      return 0.;

   /* get significand */
   if (*string == '-')
      *bufptr++ = '-';
   *bufptr = '.';
   while (isdigit(*++bufptr = *++string))
      ;   /* copy significand */
   *bufptr++ = 'E';

   /* get exponent */
   if (*string == '\0')   /* no sign on exponent */
      ++string;

   strcpy(bufptr, string);      /* copy exponent */
   return atof(buf);
}

// smart tle reader, reads first tle in file, no conversion of elements
void read_tle(char *file)
{
  FILE *fp;
  if ((fp = fopen(file, "r")) == NULL)
  {
    printf("\ncan't open %s\n", file);
    s_in("[exit]", buf);
    exit(0);
  }
  line1[0] = '0';
  int i = 0;            // counter to determine if name string found
  do
  {
    if(*line1 == '1' && strlen(line1) == 70)
    {
      fgets(line2, 80, fp);
      if(*line2 == '2' && strlen(line2) == 70)
      {
        if(i)
        {
           printf("\n%s", name);    // print name, string above line1
        }
        fclose(fp);

        // first data line
        printf("%s", line1);         // print first line on screen
        ssn = atoi(line1 + 2);        // satellite number
        sscanf(line1 + 9, "%6s", &desig);      // international designation
        tle = atof(line1 + 18);      // epoch in tle format at this point
        bstar = sci(line1 + 53);

        // second data line
        printf("%s\n", line2);        // print second line on screen
        ii = atof(line2 + 8);         // inclination, degrees
        om = atof(line2 + 17);        // ascending node, degrees
          line2[25] = '.';            // add decimal point to eccentricity
        ec = atof(line2 + 25);        // eccentricity
        ww = atof(line2 + 34);        // perigee, degrees
        ww = mod(ww, 360);
        ma = atof(line2 + 43);        // mean anomaly, degrees
        // Make sure mean motion is null-terminated, since rev. no.
        // may immediately follow.
        buf[0] = line2[63];
        line2[63] = '\0';
        nn = atof(line2 + 52);        // mean motion, revolutions/day
        line2[63] = buf[0];           // restore line2[63]
        return;
      }
    }
    else    //  if(*line1 != '1')
    {
      i++;
      sprintf(name, line1);
    }
  }while(fgets(line1, 80, fp));

  printf("\nNo TLE found in file %s\n", file);
  s_in("[exit]", buf);
  exit(0);
}

// load IOD lines from input file_in
void get_obs(char *file, char iod_line[][81])
{
  FILE *fp;
  int i;
  // raw doppler data
  if ((fp = fopen(file, "r")) == NULL)
  {
     printf("Can't open input file\n");
     s_in("[exit]", buf);
     exit(0);
  }
  nobs = 0;

  // find observation entries and write to iod_line matrix
  while(fgets(iod_line[nobs], 80, fp))
  {
    if(strlen(iod_line[nobs]) > 53 && strlen(iod_line[nobs]) != 70 &&
       isdigit(iod_line[nobs][0]))
    {
      nobs++;
      if(nobs > 1000)
      {
         printf("Exceeding 1000 Observation Maximum \n");
         s_in("[continue]", buf);
         return;
      }
    }
  }
  fclose(fp);

  aobs = 0;
  // doppler min max type data
  if ((fp = fopen(file, "r")) == NULL)
  {
     printf("Can't open input file\n");
     s_in("[exit]", buf);
     exit(0);
  }

  // find observation entries and write to iod_line matrix
  while(fgets(iod_line[nobs], 80, fp))
  {
    if(strlen(iod_line[nobs]) > 31 && (iod_line[nobs][15] == ':'))
    {
      dopp = 1;

      // radio aos/los type data
      // find observation entries and write to aos_los matrix
      if(iod_line[nobs][32] == 'O')
      {
         sprintf(aos_los[aobs], "%34s", iod_line[nobs]);
         aobs++;
      }
      else nobs++;
    }
  }
  // write aos/los on the bottom of the iod_line matrix
  for(i = 0; i < aobs; i++) sprintf(iod_line[nobs+i], "%33s", aos_los[i]);

  fclose(fp);
  return;
}

// read site data, subroutine of read_obs
void read_site(char line[])
{
  int num;
  char* fl;
  char inp_str[81];

  if(dopp) sscanf(line+6, "%4d", &site);
  else sscanf(line, "%4d", &site);    // read site number from iod_line
  // open stations.in file to read site data for each observation
  FILE *fp;
  if((fp = fopen("stations.in", "r")) != NULL)
  {
    while(fgets(inp_str, 80, fp))
    {
      // global la, lo, hh
      sscanf(inp_str, "%d %s %lf %lf %lf",&num, &fl, &la, &lo, &hh);
      if(num == site)
      {
        fclose(fp);
        return;
      }
    } // end while
    printf("\nno site data for %d\n", site);
    s_in("[exit]", buf);
    fclose(fp);
    exit(0);
  } // end if
  else
  {
    printf("\nno stations.in file found\n");
    s_in("[exit]", buf);
    exit(0);;
  }
}

void sort(char iod_line[][81],
            double odata[][4],
            double rd[][3])
{
  int i, j;
  double x, k[3], m[4];
  sort:
  j = 0;
  for(i = 0; i < nobs - 1; i++)
  {
    if(odata[i][0] > odata[i + 1][0])    // switch low/high
    {
      k = rd[i];
      rd[i] = rd[i + 1];
      rd[i + 1] = k;
      m = odata[i];
      odata[i] = odata[i + 1];
      odata[i + 1] = m;
      buf = iod_line[i];
      iod_line[i] = iod_line[i + 1];
      iod_line[i + 1] = buf;
      k = s_data[i];
      s_data[i] = s_data[i + 1];
      s_data[i + 1] = k;
      x = tol[i];
      tol[i] = tol[i + 1];
      tol[i + 1] = x;
      j++;
    }
    if(odata[i][0] == odata[i + 1][0])   // remove duplicates
    {
      for(int s = i; s < aobs+nobs - 1; s++)
      {
         rd[s] = rd[s + 1];
         odata[s] = odata[s + 1];
         iod_line[s] = iod_line[s + 1];
         s_data[s] = s_data[s + 1];
         tol[s] = tol[s + 1];
      }
      nobs--;
    }
  }
  for(i = nobs; i < aobs+nobs - 1; i++)
  {
    if(odata[i][0] > odata[i + 1][0])    // switch low/high
    {
      k = rd[i];
      rd[i] = rd[i + 1];
      rd[i + 1] = k;
      m = odata[i];
      odata[i] = odata[i + 1];
      odata[i + 1] = m;
      buf = iod_line[i];
      iod_line[i] = iod_line[i + 1];
      iod_line[i + 1] = buf;
      k = s_data[i];
      s_data[i] = s_data[i + 1];
      s_data[i + 1] = k;
      x = tol[i];
      tol[i] = tol[i + 1];
      tol[i + 1] = x;
      j++;
    }
  }
  if(j > 0)  goto sort;
}

// decodes the iod_line data
void read_obs(char iod_line[][81],
            double odata[][4],
            double rd[][3])
{

  int year, month, day, hour, min;
  double sec, fo, df, sf = 0.;

  for(int i = 0; i < nobs; i++)
  {
    // 0         1         2         3         4         5
    // 0123456789012345678901234567890123456789012345678901234
    // 0078,90027,20060704120431,01,FM,2232500000,2232487503,5

    sscanf(iod_line[i], "%4d", &site);
    sscanf(iod_line[i] + 11, "%4d %2d %2d %2d %2d %2lf",
          &year, &month, &day, &hour, &min, &sec);
    sscanf(iod_line[i] + 32, "%10lf", &fo);
    sscanf(iod_line[i] + 43, "%10lf", &df);

    Date t1(year, month, day, hour, min, sec);

    odata[i][0] = t1.jd;                   // julian date
    odata[i][1] = df;                      // received frequency
    odata[i][2] = df - fo;                 // doppler shift
    sf = (i*sf + df)/(i+1);                // running average
    odata[i][3] = site;                    // observer code

    read_site(iod_line[i]);
    Locate rre(t1.jd, la, lo, hh);
    rd[i] =  rre.rre;
  } // end for

  fo = sf;
  sf = fo / 1.E6;                          // for MHz display

  // manual input override data
  printf("(A)ccept average or ENTER tx freq (MHz) [%8.3f", sf);
  s_in("]: ", buf);
  if((*buf & 0x5f) == 'A')
  {
    for(int i = 0; i < nobs; i++)
    {
      odata[i][2] = odata[i][1] - fo;      // use average nominal tx freq
    }
  }
  else if(isdigit(buf[1]))
  {
    fo = atof(buf);
    fo *= 1.E6;
    for(int i = 0; i < nobs; i++)
    {
      odata[i][2] = odata[i][1] - fo;      // manually input nominal tx freq
    }
  }
  else
  {
    return;                                // use observer nominal tx freq
  }
}

// decodes the iod_line data
void read_dopp_obs(char iod_line[][81],
                   double odata[][4],
                   double rd[][3])
{

  int year, month, day, hour, min, i;
  double sec, pm;

  for(i = 0; i < nobs+aobs; i++)
  {
    // 0         1         2         3
    // 012345678901234567890123456789012
    // 90025 0070 2007:10:01:22:23:17 05

    sscanf(iod_line[i] + 6, "%4d", &site);
    sscanf(iod_line[i] + 11, "%4d:%2d:%2d:%2d:%2d:%2lf",
          &year, &month, &day, &hour, &min, &sec);

    Date t1(year, month, day, hour, min, sec);

    odata[i][0] = t1.jd;                   // julian date
    if(iod_line[i][31] == 'A') odata[i][1] = 0;
    else if(iod_line[i][31] == 'L') odata[i][1] = 1;
    else sscanf(iod_line[i] + 31, "%2lf", &pm);
    odata[i][3] = site;                    // observer code
    tol[i] = pm;                           // observation tolerance

    read_site(iod_line[i]);
    s_data[i][0] = la;
    s_data[i][1] = lo;
    s_data[i][2] = hh;
    Locate rre(t1.jd, la, lo, hh);
    rd[i] =  rre.rre;
  } // end for

  return;

}

// write TLE to output file and to screen
void write_tle(char *file_out)
{
  char ec_string[9];
  char xns_string[12];
  char bstar_string[13];
  char bstar_fract[7];
  char bstar_exp[3];

  sprintf(ec_string, "%.7f", ec);
  ec_string[0] = ec_string[2];
  ec_string[1] = ec_string[3];
  ec_string[2] = ec_string[4];
  ec_string[3] = ec_string[5];
  ec_string[4] = ec_string[6];
  ec_string[5] = ec_string[7];
  ec_string[6] = ec_string[8];

  sprintf(bstar_string, "%12.4e", bstar*10);
  bstar_fract[0] = bstar_string[0]; // sign
  bstar_fract[1] = bstar_string[1];
  bstar_fract[2] = bstar_string[3];
  bstar_fract[3] = bstar_string[4];
  bstar_fract[4] = bstar_string[5];
  bstar_fract[5] = bstar_string[6];
  bstar_fract[6] = '\0';
  bstar_exp[0] = bstar_string[8];
  bstar_exp[1] = bstar_string[11];
  bstar_exp[2] = '\0';

  double xns = 2160 * bstar * nn * c2;

  sprintf(xns_string, "%.8f", xns);
  if(xns_string[0] == '-')
  {
    xns_string[1] = xns_string[2];
    xns_string[2] = xns_string[3];
    xns_string[3] = xns_string[4];
    xns_string[4] = xns_string[5];
    xns_string[5] = xns_string[6];
    xns_string[6] = xns_string[7];
    xns_string[7] = xns_string[8];
    xns_string[8] = xns_string[9];
    xns_string[9] = xns_string[10];
    xns_string[10] = xns_string[11];
  }

  sprintf(line1, "1 %05dU %-8s %014.8f %10s  00000-0 %6s%2s 0    00"
           ,ssn, desig, tle, xns_string, bstar_fract, bstar_exp);
  ccksum(line1);
  sprintf(line2, "2 %05d %8.4lf %8.4lf %.7s %8.4lf %8.4lf %11.8lf    00",
           ssn, ii, om, ec_string, ww, ma, nn);
  ccksum(line2);

  if(strlen(file_out) > 0)
  {
    FILE *fp;
    fp = fopen(file_out, "a");
    fprintf(fp, "\n%s\n", line1);
    fprintf(fp, "%s\n", line2);
    fclose(fp);
  }

  printf("\n%s\n", line1);
  printf("%s\n", line2);
}

inline char *s_in(char *prompt, char *buffer)
{
  printf("%s", prompt);
  if (feof(stdin)) exit(0);
  gets(buffer);
  if (batch) printf("%s\n", buffer);
  return buffer;
}

// round the world
inline double rtw(double ao, double ac)
{
  if(fabs(ao - ac) > 180)
  {
    if(ao < 180) ao += 360;
    else ac += 360;
  }
  return ao - ac;
}

// set true longitude
void longitude(void)
{
  ma *= de2ra;
  e = ma;
  ek = 0.;
  while(fabs(e - ek) > 1.e-6)
  {
     ek = e;
     e = ma + ec * sin(e);
  }
  ma /= de2ra;
  theta = (ec - cos(e)) / (ec * cos(e) - 1);
  theta = acos(theta);
  if(e > pi) theta = 2 * pi - theta;
  uu = ww + theta / de2ra;
}
// find mean anomaly from true longitude and perigee
double ww2ma(double wx)
{
   theta = mod(uu - wx, 360) * de2ra;
   e = acos((ec + cos(theta)) / (1 + ec * cos(theta)));
   if(theta > pi) e = 2 * pi - e;
   ma = e - ec * sin(e);
   ma /= de2ra;
   return ma;
}

// box search, no limits
void anomaly(Satellite sat,  double rd[][3], double odata[][4])
{
   double min, max, step;
   step = .1;
   mk = ma;
   sum = find_sum(sat, rd, odata);
   do
   {
      min = mk;
      max = mk;
      nsum:
      min = mk - step;
      sat.delta_el(sat.jd, ii, om, ec, ww, min, nn, bstar);
      nsum = find_sum(sat, rd, odata);
      if(nsum < sum)
      {
         mk = min;
         sum = nsum;
         goto nsum;
      }
      xsum:
      max = mk + step;
      sat.delta_el(sat.jd, ii, om, ec, ww, max, nn, bstar);
      xsum = find_sum(sat, rd, odata);
      if(xsum < sum)
      {
         mk = max;
         sum = xsum;
         goto xsum;
      }
      step /= 2;
   }while(fabs(max - min) > 1e-5);
   ma = mod(mk, 360);
}

// box search, no limits
void motion(Satellite sat,  double rd[][3], double odata[][4])
{
   double min, max, step;
   step = .1;
   nk = nn;
   sum = find_sum(sat, rd, odata);
   do
   {
      min = nk;
      max = nk;
      nsum:
      min = nk - step;
      sat.delta_el(sat.jd, ii, om, ec, ww, ma, min, bstar);
      nsum = find_sum(sat, rd, odata);
      if(nsum < sum)
      {
         nk = min;
         sum = nsum;
         goto nsum;
      }
      xsum:
      max = nk + step;
      sat.delta_el(sat.jd, ii, om, ec, ww, ma, max, bstar);
      xsum = find_sum(sat, rd, odata);
      if(xsum < sum)
      {
         nk = max;
         sum = xsum;
         goto xsum;
      }
      step /= 2;
   }while(fabs(max - min) > 1e-10);
   nn = nk;
}

// partition search on node and inclination within global limits
void node(Satellite sat,  double rd[][3], double odata[][4])
{
   double minsum = 1.e6;
   while((imax - imin) > 1.e-5)
   {
      istep = (imax - imin) / 20;
      ostep = fabs(rtw(omax, omin) / 20);
      for(ik = imin; ik < imax; ik += istep)
      {
         for(ok = omin; ok < omax; ok += ostep)
         {
            Satellite satx(tle, ik, ok, ec, ww, ma, nn, bstar);
            // establish the computed ra, dc, at jdo with no perturbations
            sum = find_sum(satx, rd, odata);
            if(sum < minsum)
            {
               minsum = sum;
               ii = ik;
               om = ok;
            }
         }   // for ok
      }   // for ik
      sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar);
      imin = ii - istep;
      imax = ii + istep;
      omin = om - ostep;
      omax = om + ostep;
   }
   om = mod(om, 360);
}

// partition search on perigee and eccentricity within global limits
void perigee(Satellite sat,  double rd[][3], double odata[][4])
{
   Satellite satx = sat;
   double minsum = 1e6;
   while((wmax - wmin) > 1e-5)
   {
      estep = (emax - emin) / 20;
      wstep = (wmax - wmin) / 20;
      for(wk = wmin; wk < wmax; wk += wstep)
      {
         theta = (uu - wk) * de2ra;
         for(ek = emin; ek < emax; ek += estep)
         {
            e = acos((ek + cos(theta)) / (1 + ek * cos(theta)));
            if(theta > pi) e = 2 * pi - e;
            mk = e - ek * sin(e);
            mk = mk / de2ra;

            satx.delta_el(sat.jd, ii, om, ek, wk, mk, nn, bstar);
            // establish the computed ra, dc, at jdo with no perturbations
            sum = find_sum(satx, rd, odata);
            if(sum < minsum)
            {
               minsum = sum;
               ec = ek;
               ww = wk;
               ma = mk;
            }
         }   // for ek
      }   // for wk
      sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar);

      emin = ec - estep;
      emax = ec + estep;
      wmin = ww - wstep;
      wmax = ww + wstep;
   }
                
   // update mean_anomaly
   anomaly(sat, rd, odata);
   sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar);

   // update mean_motion
   motion(sat, rd, odata);
   sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar);

   // calculate uu, degrees
   longitude();

   ww = mod(ww, 360);
   ma = mod(ma, 360);
}

void do_sim(Satellite sat)
{
   int n_obs, year, month, day, hour, min, sec;
   double dy, time1, time2, topt, tmin, tmax, step, freq;

   double U, Uopt, Umin, Umax, cmp, sv, ev, nrr, nvv, az, el, asp, range,
          ll[3], tempv[3], temp[3], nv[3], zz[3] = {0, 0, 1};

   // copy sat
   Satellite satx = sat;

   s_in("\n               STA : ", buf);
   read_site(buf);

   s_in("          FREQ MHz : ", buf);
   freq = atof(buf);

   reset:

   s_in("START YYYYMMDDHHMM : ", buf);
   sscanf(buf, "%4d %2d %2d %2d %2d", &year, &month, &day, &hour, &min);
   Date t1(year, month, day, hour, min, 0.01);
   time1 = t1.jd;
   s_in("STOP  YYYYMMDDHHMM : ", buf);
   sscanf(buf, "%4d %2d %2d %2d %2d", &year, &month, &day, &hour, &min);
   Date t2(year, month, day, hour, min, 0.01);
   time2 = t2.jd;
   s_in("         STEP MMSS : ", buf);
   sscanf(buf, "%2d %2d", &min, &sec);
   step = (min + sec / 60.) / 1440.;

   zoom:
   n_obs = ((int)((time2 - time1) / step)) + 2;
   double srd[n_obs][3], sdata[n_obs];

   n_obs = 0;
   for(dy = time1; dy <= time2; dy += step)
   {
      // load arrays
      sdata[n_obs] = dy;                  // julian date
      Locate rre(dy, la, lo, hh);
      srd[n_obs] = rre.rre;               // observer location
      n_obs++;
   }

   while(1)
   {
     if(out)
     {
       FILE *fp;
       fp = fopen(file, "a");
       fprintf(fp, "\n\n  TIME       AZ    EL     RNG          Freq\n");
       fclose(fp);
     }
     printf("\n  TIME       AZ     EL     RNG          Freq\n");
     for(int j = 0; j < n_obs; j++)
     {

       // advance satellite position
       satx.delta_t(sdata[j]);
       nrr = norm(satx.rr);
       nvv = norm(satx.vv);
       vmadd(satx.rr, srd[j], temp, -1);

       // convert ll to unit vector
       el = norm(temp);
       range = el * 6378.135;
       smult(1/el, temp, ll);

       //  computing elevation
       el = acos(dot(srd[j], ll) / norm(srd[j]));
       el /= de2ra;
       el = 90 - el;

       // computing aspect
       asp = acos(dot(ll, satx.vv) / nvv);
       asp /= de2ra;
       asp = 180 - asp;

       // computing azimuth
       cross(srd[j], zz, tempv);
       cross(tempv, srd[j], nv);
       cross(srd[j], ll, tempv);
       cross(tempv, srd[j], temp);
       az = acos(dot(nv, temp) / (norm(nv)*norm(temp)));
       cross(temp, nv, tempv);
       if(dot(tempv, srd[j]) < 0) az = 2*pi - az;
       az /= de2ra;

       // computed doppler and velocity
       cross(zz, srd[j], tempv);
       smult(.004351409367, tempv, temp);
       sv = dot(satx.vv, ll);
       ev = dot(temp, ll);
       U  = ev - sv;
       cmp = sqrt((2820.188631 + U) / (2820.188631 - U));
       cmp = (cmp - 1.) * freq;
       U  = U * 6378135./60.;
       if(j == 0)
       {
         Umin = U;
         Umax = U;
         tmin = sdata[j];
         tmax = sdata[j];
       }
       if(j > 0 && U < Umin)
       {
         Umin = U;
         tmin = sdata[j];
       }
       if(j > 0 && U > Umax)
       {
         Umax = U;
         tmax = sdata[j];
       }
       topt = sdata[j];

       // time
       Date t3(sdata[j]);
       if(!pass || el > 0)
       {
         printf("%02d:%02d:%02d  %5.1f  %5.1f  %6.0f  %12.5f\n",
               (int)t3.hr, (int)t3.mn, (int)t3.ss,
               az, el, range, freq+cmp);
         // print fit to file
         if(out)
         {
           FILE *fp;
           fp = fopen(file, "a");
           fprintf(fp, "%02d:%02d:%02d  %5.1f  %5.1f  %6.0f  %12.5f\n",
                      (int)t3.hr, (int)t3.mn, (int)t3.ss,
                      az, el, range, freq+cmp);
           fclose(fp);
         }
       }
     }
     out = 0;

     printf("\n(A)nomaly  (M)otion  (Z)oom  (R)eset  (Q)uit");
     s_in(": ", buf);
     if ((*buf & 0x5f) == 'A')
     {
       s_in(": ", buf);
       if (isdigit(buf[0]))
       {
         ma = atof(buf);
         satx.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar);
         satx.delta_t(sat.jd);
         satx.print_el();       // print new elements
       }
     }
     if ((*buf & 0x5f) == 'M')
     {
       s_in(": ", buf);
       if (isdigit(buf[0]))
       {
         nn = atof(buf);
         satx.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar);
         satx.delta_t(sat.jd);
         satx.print_el();       // print new elements
       }
     }
     if ((*buf & 0x5f) == 'Z')
     {
       if(tmin > time1 && tmin < topt) topt = tmin;   // not an endpoint
       else topt = tmax;
       time1 = topt - step;
       time2 = topt + step;
       step /= 5.;
       goto zoom;
     }
     if ((*buf & 0x5f) == 'R')
     {
       printf("\n");
       goto reset;
     }
     if ((*buf & 0x5f) == 'Q') break;
   } // end while(1)
}

void print_fit(Satellite sat,  double rd[][3], double odata[][4])
{
  double U, sv, ev, nrr, nvv, az, el, asp, cmp, Diff, zum = 0,
         ll[3], tempv[3], temp[3], nv[3], zz[3] = {0, 0, 1};

  char type[3];
  // copy sat
  Satellite satx = sat;

  if(dopp)
  {
     sum = find_dopp(satx, rd, odata);
     if(out)
     {
         FILE *fp;
         fp = fopen(file, "a");
         fprintf(fp, "\n\n       sta  observed              predicted             Diff");
         fclose(fp);
     }
     printf("\n       sta  observed              predicted             Diff");
     for(int j = 0; j < nobs; j++)
     {
        Date t1(odata[j][0]);
        Date t2(odata[j][1]);
        Diff = (odata[j][0] - odata[j][1]) * 86400.;
        // remove the "non observation"
        if(fabs(Diff) > 1943)
           printf("\n(%2d)  %04.0f  %4.0f:%02.0f:%02.0f:%02.0f:%02.0f:%02.0f   none",
             j + 1, odata[j][3], t1.yy, t1.mm, t1.dd, t1.hr, t1.mn, t1.ss);
        else
           printf("\n(%2d)  %04.0f  %4.0f:%02.0f:%02.0f:%02.0f:%02.0f:%02.0f   %4.0f:%02.0f:%02.0f:%02.0f:%02.0f:%02.0f  %+5.0f",
             j + 1, odata[j][3], t1.yy, t1.mm, t1.dd, t1.hr, t1.mn, t1.ss,
                                 t2.yy, t2.mm, t2.dd, t2.hr, t2.mn, t2.ss, Diff);

        // print fit to file
        if(out)
        {
           FILE *fp;
           fp = fopen(file, "a");
           if(fabs(Diff) > 1943)
              fprintf(fp, "\n(%2d)  %04.0f  %4.0f:%02.0f:%02.0f:%02.0f:%02.0f:%02.0f   none",
                 j + 1, odata[j][3], t1.yy, t1.mm, t1.dd, t1.hr, t1.mn, t1.ss);
           else
              fprintf(fp, "\n(%2d)  %04.0f  %4.0f:%02.0f:%02.0f:%02.0f:%02.0f:%02.0f   %4.0f:%02.0f:%02.0f:%02.0f:%02.0f:%02.0f  %+5.0f",
                 j + 1, odata[j][3], t1.yy, t1.mm, t1.dd, t1.hr, t1.mn, t1.ss,
                                     t2.yy, t2.mm, t2.dd, t2.hr, t2.mn, t2.ss, Diff);
           fclose(fp);
        }
     }

     // aos/los data
     if(aobs && out)
     {
         FILE *fp;
         fp = fopen(file, "a");
         fprintf(fp, "\n\n       sta  observed               az   el  type");
         fclose(fp);
     }
     if(aobs) printf("\n       sta  observed               az   el  type");
     for(int j = nobs; j < nobs+aobs; j++)
     {
        Date t1(odata[j][0]);

        // advance satellite position
        satx.delta_t(odata[j][0]);
        nrr = norm(satx.rr);
        nvv = norm(satx.vv);
        vmadd(satx.rr, rd[j], temp, -1);

        // convert ll to unit vector
        el = norm(temp);
        smult(1/el, temp, ll);

        //  computing elevation
        el = acos(dot(rd[j], ll) / norm(rd[j]));
        el /= de2ra;
        el = 90 - el;

        // computing azimuth
        cross(rd[j], zz, tempv);
        cross(tempv, rd[j], nv);
        cross(rd[j], ll, tempv);
        cross(tempv, rd[j], temp);
        az = acos(dot(nv, temp) / (norm(nv)*norm(temp)));
        cross(temp, nv, tempv);
        if(dot(tempv, rd[j]) < 0) az = 2*pi - az;
        az /= de2ra;

        // type of obs
        if(odata[j][1] == 0) sprintf(type, "AOS");
        else sprintf(type, "LOS");

        printf("\n(%2d)  %04.0f  %4.0f:%02.0f:%02.0f:%02.0f:%02.0f:%02.0f   %3.0f  %+3.0f   %3s",
          j + 1, odata[j][3], t1.yy, t1.mm, t1.dd, t1.hr, t1.mn, t1.ss, az, el, type);
        // print fit to file
        if(out)
        {
           FILE *fp;
           fp = fopen(file, "a");
           fprintf(fp, "\n(%2d)  %04.0f  %4.0f:%02.0f:%02.0f:%02.0f:%02.0f:%02.0f   %3.0f  %+3.0f   %3s",
              j + 1, odata[j][3], t1.yy, t1.mm, t1.dd, t1.hr, t1.mn, t1.ss, az, el, type);
           fclose(fp);
        }
     }

     printf("\n");
     printf("\nsum%12.5f\n", sum);
     return;
  }

  if(out)
  {
      FILE *fp;
      fp = fopen(file, "a");
      fprintf(fp, "\n\n          TIME     STA   AZ      EL      ASP        OBS      CMP     Diff\n");
      fclose(fp);
  }

  printf("\n          TIME     STA   AZ      EL      ASP        OBS      CMP     Diff\n");
  for(int j = 0; j < nobs; j++)
  {

     // advance satellite position
     satx.delta_t(odata[j][0]);
     nrr = norm(satx.rr);
     nvv = norm(satx.vv);
     vmadd(satx.rr, rd[j], temp, -1);

     // convert ll to unit vector
     el = norm(temp);
     smult(1/el, temp, ll);

     //  computing elevation
     el = acos(dot(rd[j], ll) / norm(rd[j]));
     el /= de2ra;
     el = 90 - el;

     // computing aspect
     asp = acos(dot(ll, satx.vv) / nvv);
     asp /= de2ra;
     asp = 180 - asp;

     // computing azimuth
     cross(rd[j], zz, tempv);
     cross(tempv, rd[j], nv);
     cross(rd[j], ll, tempv);
     cross(tempv, rd[j], temp);
     az = acos(dot(nv, temp) / (norm(nv)*norm(temp)));
     cross(temp, nv, tempv);
     if(dot(tempv, rd[j]) < 0) az = 2*pi - az;
     az /= de2ra;

     // computed doppler and %diff
     cross(zz, rd[j], tempv);
     smult(.004351409367, tempv, temp);
     sv = dot(satx.vv, ll);
     ev = dot(temp, ll);
     U  = ev - sv;
     cmp = sqrt((2820.188631 + U) / (2820.188631 - U));
     cmp = (cmp - 1.) * (odata[j][1] - odata[j][2]);
     Diff = odata[j][2] - cmp;

     // sum percentage difference error squared
     zum += Diff*Diff;

     // time
     Date t1(odata[j][0]);

     printf("(%3d)  %02d%02d:%02d  %04d  %6.2f  %5.2f  %6.2f   %6.0f   %6.0f   %6.0f\n",
             j + 1, (int)t1.hr, (int)t1.mn, (int)t1.ss,
             (int)odata[j][3], az, el, asp, odata[j][2], cmp, Diff);

     // print fit to file
     if(out)
     {
        FILE *fp;
        fp = fopen(file, "a");
        fprintf(fp, "(%3d)  %02d%02d:%02d  %04d  %6.2f  %5.2f  %6.2f   %6.0f   %6.0f   %6.0f\n",
                     j + 1, (int)t1.hr, (int)t1.mn, (int)t1.ss,
                     (int)odata[j][3], az, el, asp, odata[j][2], cmp, Diff);
        fclose(fp);
     }
  }
  printf("\nsum%12.5f\n", sqrt(zum / nobs));
}

double find_dopp(Satellite sat,  double rd[][3], double odata[][4])
{
   int k, n_obs;
   double dy, time1, time2, topt, tmin, tmax, step;

   double nrr, U, Uopt, Umin, Umax, sv, ev, el, Diff, zum = 0.,
          ll[3], tempv[3], temp[3], hiz[3], zz[3] = {0, 0, 1};

   // copy sat
   Satellite satx = sat;

   for(int i = 0; i < nobs; i++)
   {
      // initial partition
      time1 = odata[i][0] - .02;
      time2 = odata[i][0] + .02;
      step  = 6./1440.;             // 6 min steps
      k = 0;

      ///////////////////////
      zoom:
      k++;
      n_obs = ((int)((time2 - time1) / step)) + 2;
      double srd[n_obs][3], sdata[n_obs];

      n_obs = 0;
      for(dy = time1; dy <= time2; dy += step)
      {
         // load arrays
         sdata[n_obs] = dy;                  // julian date
         Locate rre(dy, s_data[i][0], s_data[i][1], s_data[i][2]);
         srd[n_obs] = rre.rre;               // observer location
         n_obs++;
      }

      for(int j = 0; j < n_obs; j++)
      {

         // advance satellite position
         satx.delta_t(sdata[j]);
         vmadd(satx.rr, srd[j], temp, -1);

         // convert ll to unit vector
         el = norm(temp);
         smult(1/el, temp, ll);

         // computed doppler and velocity
         cross(zz, srd[j], tempv);
         smult(.004351409367, tempv, temp);
         sv = dot(satx.vv, ll);
         ev = dot(temp, ll);
         U  = ev - sv;
         if(j == 0)
         {
            Umin = U;
            Umax = U;
            tmin = sdata[j];
            tmax = sdata[j];
         }
         if(j > 0 && U < Umin)
         {
            Umin = U;
            tmin = sdata[j];
         }
         if(j > 0 && U > Umax)
         {
            Umax = U;
            tmax = sdata[j];
         }
         topt = sdata[j];

      }  // end for j

      if (k < 5)
      {
         if(tmin > time1 && tmin < topt) topt = tmin;   // not an endpoint
         else topt = tmax;
         time1 = topt - step;
         time2 = topt + step;
         step /= 6.;
         goto zoom;
      }
      odata[i][1] = topt;
      Diff = (odata[i][0] - odata[i][1]) * 86400.;
      if(fabs(Diff) < tol[i]) Diff = 0.;
      else Diff = Diff - tol[i];
      zum += Diff*Diff;

   } // end for i

   for(int j = nobs; j < nobs+aobs; j++)
   {
      Date t1(odata[j][0]);

      // advance satellite position
      satx.delta_t(odata[j][0]);
      // line of sight vector, temp
      vmadd(satx.rr, rd[j], temp, -1);

      // convert ll to unit vector
      el = norm(temp);              // magnitude of ll
      smult(1/el, temp, ll);

      // observer angular velocity vector, temp
      cross(zz, rd[j], tempv);
      ev = norm(tempv);
      smult(.004351409367/ev, tempv, temp);

      // sat angular velocity vector
      nrr = norm(satx.rr);
      smult(1/nrr, satx.vv, zz);

      // direction toward horizon from satellite
      cross(rd[j], ll, tempv);
      ev = norm(rd[j]);
      smult(1/ev, tempv, tempv);
      cross(tempv, ll, hiz);

      // sat angular velocity toward horizon - observer
      Diff = dot(zz, hiz) - dot(temp, ll);

      //  computed elevation, radians
      ev = acos(dot(rd[j], ll) / norm(rd[j]));
      el = pi/2 - ev;

      // minutes from motionless observer horizon
      Diff = el/Diff;
      // seconds
      Diff *= 60;

      // convert to degrees for display
      el /= de2ra;

      zum += Diff*Diff;

    }

   return sqrt(zum / (aobs+nobs));
}

double find_sum(Satellite sat,  double rd[][3], double odata[][4])
{
  double U, sv, ev, el, cmp, Diff, zum = 0,
         ll[3], tempv[3], temp[3], zz[3] = {0, 0, 1};

  // copy sat
  Satellite satx = sat;

  if(dopp) return find_dopp(satx, rd, odata);

  for(int j = 0; j < nobs; j++)
  {

     // advance satellite position
     satx.delta_t(odata[j][0]);
     vmadd(satx.rr, rd[j], temp, -1);

     // convert ll to unit vector
     el = norm(temp);
     smult(1/el, temp, ll);

     // computed doppler
     cross(zz, rd[j], tempv);
     smult(.004351409367, tempv, temp);
     sv = dot(satx.vv, ll);
     ev = dot(temp, ll);
     U  = ev - sv;
     cmp = sqrt((2820.188631 + U) / (2820.188631 - U));
     cmp = (cmp - 1.) * (odata[j][1] - odata[j][2]);
     Diff = fabs(odata[j][2] - cmp);

     // sum position error squared
     zum += Diff*Diff;
  }
  return sqrt(zum / nobs);
}

