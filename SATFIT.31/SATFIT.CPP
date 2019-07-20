// satfit.cpp               // retain perturbation values at epoch
#include <ctype.h>          // for isdigit in read_tle
#include <string.h>         // for command line input
#include "orbit_gray.h"     // for Date, Locate, Satellite
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
   pgee,     // perigee
   agee,     // apogee
   bstar;    // BSTAR drag term

int
   xi = 0, xe = 0, xw = 0, xn = 0;   // flags

// please note that these variables are in TLE format, i.e. degrees, whereas
// their counterparts inside SGP4 are in radians.

double la, lo, hh;          // observer parameters
int first = 1, batch = 0, out = 0;
int nobs, len;              // number of observations
char file[81];              // input file string declared
char buf[81];               // input buffer
char name[81], tle1[81], tle2[81];    // original tle data lines
char line1[81], line2[81];  // read and write tle buffers
double sum, nsum, xsum;     // sum of squares of residuals
double xrms;
double zero, tleh;          // history parameters
FILE *fp, *fpo;

// declare variables for search and history options
double e, ik, ok, ek, wk, mk, nk, bk, theta, minsum, ra, dc;

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
            double ll[][3],
            double rd[][3]);
void sort(char iod_line[][81],
            double odata[][4],
            double ll[][3],
            double rd[][3]);
void read_tle(char *file);
void print_fit(Satellite sat,  double rd[][3],
              double ll[][3],  double odata[][4]);
double find_rms(Satellite sat,  double rd[][3],
               double ll[][3],  double odata[][4]);
void anomaly(Satellite sat,  double rd[][3],
            double ll[][3],  double odata[][4]);
void motion(Satellite sat,  double rd[][3],
           double ll[][3],  double odata[][4]);
void perigee(Satellite sat,  double rd[][3],
            double ll[][3],  double odata[][4]);
void node(Satellite sat,  double rd[][3],
         double ll[][3],  double odata[][4]);
void align(Satellite sat,  double rd[][3],
         double ll[][3],   double odata[][4]);
void rref(double m[][7], double b[]);
void zrll(Satellite sat, double* rd, double& ra, double& dc);
void diff_el(Satellite sat,  double rd[][3],
         double ll[][3],   double odata[][4]);
void longitude(void);
double ww2ma(double wx);
void write_tle(char *file);
void so2r(double r, double rd[], double ll[], double rr[]);
inline double rtw(double ao, double ac);
inline char *s_in(char *prompt, char *buffer);

/////////////////// MAIN //////////////////////////////////////////////////////

int main(int argc, char *argv[])           // command line args for read data
{
  char iod_line[200][81], file_in[81], file_out[81];
  int i, j;                                // counters
  double max, min, rms;
  char srch[] = {'W'};                     // initialize wide search

  if(argc == 1)
    sprintf(file, "sat.txt");              // default input/output file
  if(argc == 2)                            // command line file or default?
    sprintf(file, argv[1]);                // command line input/output file
  if (argc == 3 && strcmp(argv[2],"ba") == 0)
  {
    batch = 1;
    sprintf(file, argv[1]);
  }
  if(strchr(file, '\\'))
    sprintf(file, "%s", strrchr(file, '\\')+1);

  restart:

  // find observation lines in input file
  // number of observations = nobs (global)
  // store observations to iod_line matrix (up to 200 observations)
  get_obs(file, iod_line);

  double ll[nobs+1][3];       // line of sight vectors; observer -> satellite
  double rd[nobs+1][3];       // topocentric vectors to observer positions
  double odata[nobs+1][4];    // observational data: jd, ra, dc, obscode

  // storage for maneuver function
  int nobsx;
  char iod_linex[200][81];
  double llx[nobs][3];
  double rdx[nobs][3];
  double odatax[nobs][4];


  // read first tle from input file
  // orbit parameters are passed from read_tle as global variables
  // echo tle to screen
  read_tle(file);
  sprintf(tle1, "%s", line1);  // copy original line1 to tle1
  sprintf(tle2, "%s", line2);  // copy original line2 to tle2

  // Create Satellite variable, sat, from TLE
  Satellite sat(tle, ii, om, ec, ww, ma, nn, bstar);
  Satellite save_sat(tle, ii, om, ec, ww, ma, nn, bstar);
  // maneuvering sat
  Satellite satm(tle, ii, om, ec, ww, ma, nn, bstar);

  // calculate uu, degrees, for search
  longitude();

  // read all observations from input
  // Observed topocentric vectors for each observer position.
  // Equatorial coordinates.
  read_obs(iod_line, odata, ll, rd);
  sort(iod_line, odata, ll, rd);
  sum = find_rms(sat, rd, ll, odata);
  printf("%d Observtions Found\n",nobs);

  // print dates
  Date t2;
  printf("\nTODAY   : %d", t2.doy);
  if(nobs > 0)
  {
    Date t3(odata[nobs - 1][0]);
    printf("\nLAST OB : %d\n", t3.doy);
  }

accept:
   // Accept a new command

   printf("\nEnter command\n");
   printf(" (I)ncl  (N)ode  (X)ntrcty  (P)erigee   (A)nomaly  (M)otion  (B)star\n");
   printf(" (S)tep   I(D)   (T)ime  (F)it  (V)iew  (R)emove   (W)rite   (Q)uit\n");

   s_in(": ", buf);

   if ((*buf & 0x5f) == 'G') goto fit_out;           // hidden
   if ((*buf & 0x5f) == 'E') goto edit_file;         // hidden
   if ((*buf & 0x5f) == 'H') goto history;           // hidden
   if ((*buf & 0x5f) == 'Y') goto edit_history;      // hidden
   if ((*buf & 0x5f) == 'C') goto elcor;             // hidden
   if ((*buf & 0x5f) == 'O') goto discover;          // hidden
   if ((*buf & 0x5f) == 'U') goto maneuver;          // hidden
   if ((*buf & 0x5f) == 'S') goto step;
   if ((*buf & 0x5f) == 'Z') goto step;
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
   if ((*buf & 0x5f) == 'D') goto id;
   if ((*buf & 0x5f) == 'T') goto time;
   if ((*buf & 0x5f) == 'F') goto fit;
   if ((*buf & 0x5f) == 'V') goto viewobs;
   if ((*buf & 0x5f) == 'R') goto remove;
   if ((*buf & 0x5f) == 'W') goto write_el;
   if ((*buf & 0x5f) == 'Q') exit(0);

   // try again
   goto accept;

fit_out:
   out = 1;
   print_fit(sat, rd, ll, odata);
   out = 0;
   goto accept;

elcor:
   sprintf(file_in, "elcord %s", file);
   system(file_in);
   goto restart;

edit_file:
   sprintf(file_in, "%s", file);
   system(file_in);
   goto restart;

edit_history:
   sprintf(file_in, "history%c%05d.txt", '\\', ssn);
   system(file_in);
   goto accept;

viewobs:
   printf("\n");
   for(i = 0; i < nobs; i++)
      printf("(%2d) %s",i + 1, iod_line[i]);
   printf("\n");
   goto accept;

remove:
   s_in("\nRemove Line Number : ", buf);
   j = atoi(buf);
   for(i = j; i < nobs; i++)
   {
      rd[i - 1] = rd[i];
      ll[i - 1] = ll[i];
      odata[i - 1] = odata[i];
      iod_line[i - 1] = iod_line[i];
   }
   nobs--;

fit:
   out = 0;
   print_fit(sat, rd, ll, odata);
   sat.print_el();
   goto accept;

id:
   // requires SATID in same folder
   sprintf(buf, "satid %s", file);
   system(buf);
   goto restart;

history:

   printf("\n(A)dd  (G)raphs  (E)dit  (Q)uit  ");
   s_in(": ", buf);

   if ((*buf & 0x5f) == 'Q')
   {
     write_tle("");
     goto accept;
   }

   sprintf(file_in, "history/%05d.txt", ssn);
   if ((fp = fopen(file_in, "r")) == NULL)
   {
     fclose(fp);
     // create history file
     fpo = fopen(file_in, "w");
     //  write current TLE to history
     fprintf(fpo, "\n");
     fprintf(fpo, "%s", line1);
     fprintf(fpo, "%s", line2);
     fclose(fpo);
     printf("\nHISTORY created with current TLE\n");
   }
   fclose(fp);

   if ((*buf & 0x5f) == 'G')
   {
     i = 1;
     fp = fopen(file_in, "r");

     sprintf(file_out, "trendA.csv");
     fpo = fopen(file_out, "w");
     fprintf(fpo, "time,inc,omega,ecc,perigee,motion,YYddd,,,%05d\n", ssn);

     while(fgets(line1, 80, fp))
     {
       if(*line1 == '1' && strlen(line1) == 70)
       {
         fgets(line2, 80, fp);
         if(*line2 == '2' && strlen(line2) == 70)
         {
           // first data line
           tleh = atof(line1 + 18);      // epoch in tle format at this point
           Date t1(tleh);
           tleh = t1.jd;
           if(i)
           {
             zero = tleh;
             i = 0;
           }
           tleh -= zero;

           // second data line
           ik = atof(line2 + 8);         // inclination, degrees
           ok = atof(line2 + 17);        // ascending node, degrees
             line2[25] = '.';            // add decimal point to eccentricity
           ek = atof(line2 + 25);        // eccentricity
           wk = atof(line2 + 34);        // perigee, degrees
           mk = atof(line2 + 43);        // mean anomaly, degrees
           // Make sure mean motion is null-terminated, since rev. no.
           // may immediately follow.
           line2[63] = '\0';
           nk = atof(line2 + 52);        // mean motion, revolutions/day
         }
         fprintf(fpo, "%f,%f,%f,%f,%f,%f,%.5s\n", tleh, ik, ok, ek, wk, nk, line1 + 18);
         theta = tleh;
       }
     }
     Date t1(tle);
     tleh = t1.jd;
     tleh -= zero;
     if(tleh > theta)
        fprintf(fpo, "%f,%f,%f,%f,%f,%f,%5.0f\n", tleh, ii, om, ec, ww, nn, tle);
     fclose(fpo);
     sprintf(file_out, "trendA.xls");
     system(file_out);

     fclose(fp);
     // sprintf(file_out, "del trendA.csv");
     // system(file_out);
   }
   if ((*buf & 0x5f) == 'E')
   {
     sprintf(file_in, "history%c%05d.txt", '\\', ssn);
     system(file_in);
   }
   if((*buf & 0x5f) == 'A')
   {
     fclose(fp);
     //  write current TLE to history
     write_tle(file_in);
     printf("\nCurrent TLE added to HISTORY\n");
   }

   goto history;

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
         e = acose((ek + cos(theta)) / (1 + ek * cos(theta)));
         if(theta > pi) e = 2 * pi - e;
         mk = e - ek * sin(e);
         mk = mk / de2ra;
         for(nk = nmin; nk < nmax; nk += nstep)
         {
            Satellite satx(sat.jd, ii, om, ek, wk, mk, nk, bstar);
            // establish the computed ra, dc, at jdo with no perturbations
            rms = find_rms(satx, rd, ll, odata);
            if(rms < sum)
            {
               sum = rms;       // global
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
   anomaly(sat, rd, ll, odata);
   sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar);

   // print new elements
   print_fit(sat, rd, ll, odata);
   sat.print_el();

   longitude();

   srch[0] = 'Z';
   goto accept;

step:       // partition search within limits set below

   // update mean_anomaly
   anomaly(sat, rd, ll, odata);
   sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar);

   emax = ec * 1.1;
   emin = ec * .9;
   wmax = ww + 2;
   wmin = ww - 2;
   imax = ii + 2;
   imin = ii - 2;
   omax = om + 2;
   omin = om - 2;

   printf("\nPress Q to Quit    :\n\n");
   while(1)
   {
      perigee(sat, rd, ll, odata);
      sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar);
      node(sat, rd, ll, odata);
      sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar);
      if(nobs > 3 && (*buf & 0x5f) == 'Z')
      {
        diff_el(sat, rd, ll, odata);
        sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar);
      }

      // set new limits
      emax = 1.01 * ec;
      emin = .99 * ec;
      wmax = ww + .5;
      wmin = ww - .5;
      imax = ii + .5;
      imin = ii - .5;
      omax = om + .5;
      omin = om - .5;

      printf("rms%12.5f", sum);
      s_in("    : ", buf);
      if ((*buf & 0x5f) == 'Q') break;
   }

   print_fit(sat, rd, ll, odata);
   sat.print_el();       // print new elements

   srch[0] = 'N';
   goto accept;

incl:

   printf("\n(A)uto  (ii)  (Q)uit  ");
   s_in(": ", buf);

   if (isdigit(buf[0]))
   {
     ii = atof(buf);
     xi = 1;
     sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar);
     print_fit(sat, rd, ll, odata);
     sat.print_el();       // print new elements
   }
   if ((*buf & 0x5f) == 'A')
   {
     // partition search

     imax = ii + 2;
     imin = ii - 2;
     omax = om + 2;
     omin = om - 2;
     xi = 0;

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

     node(sat, rd, ll, odata);

     sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar);
     print_fit(sat, rd, ll, odata);
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
     print_fit(sat, rd, ll, odata);
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

     while((omax - omin) > 1e-5)
     {
       ostep = fabs(rtw(omax, omin) / 20);
       for(ok = omin; ok < omax; ok += ostep)
       {
          satx.delta_el(sat.jd, ii, ok, ec, ww, ma, nn, bstar);
          // establish the computed ra, dc, at jdo with no perturbations
          rms = find_rms(satx, rd, ll, odata);
          if(rms < sum)
          {
             sum = rms;      // global
             om = ok;
          }
       }   // for ok
       satx.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar);

       omin = om - ostep;
       omax = om + ostep;
     }
     om = mod(om, 360);

     sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar);
     print_fit(sat, rd, ll, odata);
     sat.print_el();       // print new elements

     srch[0] = 'N';
   }
   if ((*buf & 0x5f) == 'Q') goto accept;
   goto node;

xntrcty:

   printf("\n(S)earch  (A)uto  (ec)  (Q)uit  ");
   s_in(": ", buf);

   emax = ec * 1.05;
   emin = ec * .95;

   if (isdigit(buf[0]))
   {
     ec = atof(buf);
     xe = 1;
     sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar);
     print_fit(sat, rd, ll, odata);
     sat.print_el();       // print new elements
   }
   if ((*buf & 0x5f) == 'A')
   {
     xe = 0;
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
           rms = find_rms(sat, rd, ll, odata);
           if(rms < sum)
           {
              sum = rms;
              ec = ek;
           }
        }   // for ek
        sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar);
        emin = ec - estep;
        emax = ec + estep;
     }
 
     sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar);
     print_fit(sat, rd, ll, odata);
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
     printf("\neccentricity  sum");
     for(ec = emin; ec < emax + estep; ec += estep)
     {
       sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar);
       sum = find_rms(sat, rd, ll, odata);
       printf("\n%.7f     %7.4f", ec, sum);
     }
     printf("\n");
     ec = ek;
     sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar);
   }
   if ((*buf & 0x5f) == 'Q') goto accept;
   goto xntrcty;

perigee:

   pgee = (1-ec)*sat.aodp;
   agee = (1+ec)*sat.aodp;

   printf("\nPerigee = %f er", pgee);
   printf("     %d X %d km\n", (int)((pgee-1)*6378.135),
                               (int)((agee-1)*6378.135));

   printf("\n(S)earch  (A)uto  (ww)  (Q)uit  ");
   s_in(": ", buf);

   if (isdigit(buf[0]))
   {
     longitude();
     ww = atof(buf);
     ma = ww2ma(ww);
     xw = 1;
     sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar);
     print_fit(sat, rd, ll, odata);
     sat.print_el();       // print new elements
   }
   if ((*buf & 0x5f) == 'A')
   {
     xw = 0;
     double wx = ww + 1;
     while(fabs(wx - ww) > 1e-4)
     {
       wx = ww;
       printf("\n%8.4f  %.7f", ww, ec);
       wmax = ww + 1;
       wmin = ww - 1;
       emax = ec * 1.01;
       emin = ec * .99;
       perigee(sat, rd, ll, odata);
       sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar);
     }
     print_fit(sat, rd, ll, odata);
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
     printf("\nperigee      sum");
     for(ww = wmin; ww < wmax + wstep; ww += wstep)
     {
       ma = ww2ma(ww);
       sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar);
       sum = find_rms(sat, rd, ll, odata);
       printf("\n%8.4f  %7.4f", mod(ww, 360), sum);
     }
     printf("\n");
     ww = wk;
     ma = mk;
     sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar);
   }
   if ((*buf & 0x5f) == 'Q') goto accept;
   goto perigee;

   // amax and amin are used as temporary variables
   longitude();
   amax = uu * de2ra;
   amin = sin(sat.xincl/de2ra) * sin(amax);
   amin *= amin;
   amin = sqrt(1 - amin);
   amin = (acose(cos(amax) / amin)) / de2ra;
   if(fmod(uu, 360) > 180) amin = 360 - amin;
   if(sat.xincl/de2ra < 90)
        amax = fmod(sat.xnodeo/de2ra + amin - sat.thetag + 360, 360.);
   else amax = fmod(sat.xnodeo/de2ra - amin - sat.thetag + 720, 360.);
   printf("\nE Longitude = %8.4f\n", amax);

anomaly:

   if(nn < 1.5)
   {
     // amax and amin are used as temporary variables
     longitude();
     amax = uu * de2ra;
     amin = sin(sat.xincl/de2ra) * sin(amax);
     amin *= amin;
     amin = sqrt(1 - amin);
     amin = (acose(cos(amax) / amin)) / de2ra;
     if(fmod(uu, 360) > 180) amin = 360 - amin;
     if(sat.xincl/de2ra < 90)
          amax = fmod(sat.xnodeo/de2ra + amin - sat.thetag + 360, 360.);
     else amax = fmod(sat.xnodeo/de2ra - amin - sat.thetag + 720, 360.);
     printf("\nE Longitude = %8.4f\n", amax);
   }

   printf("\n(S)earch  (A)uto  (ma)  (L)ast  (Q)uit  ");
   s_in(": ", buf);

   amax = 360.;
   amin = 0.;
   astep = 20;

   if (isdigit(buf[0]))
   {
     ma = atof(buf);
     longitude();
     sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar);
     print_fit(sat, rd, ll, odata);
     sat.print_el();       // print new elements
   }
   if ((*buf & 0x5f) == 'A')
   {
     anomaly(sat, rd, ll, odata);
     sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar);
     print_fit(sat, rd, ll, odata);
     sat.print_el();       // print new elements
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
       sum = find_rms(sat, rd, ll, odata);
       printf("\n%8.4f     %7.4f", ma, sum);
     }
     printf("\n");
     ma = mk;              // restore
     sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar);
   }
   if ((*buf & 0x5f) == 'L')
   {
     align(sat, rd, ll, odata);
     sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar);
     longitude();
     print_fit(sat, rd, ll, odata);
     sat.print_el();       // print new elements
   }
   if ((*buf & 0x5f) == 'Q') goto accept;
   goto anomaly;

motion:

   printf("\n(A)uto  (nn)  (Q)uit  ");
   s_in(": ", buf);

   if (isdigit(buf[0]))
   {
     nn = atof(buf);
     xn = 1;
     sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar);
     print_fit(sat, rd, ll, odata);
     sat.print_el();       // print new elements
   }
   if ((*buf & 0x5f) == 'A')
   {
     xn = 0;
     // update mean motion, no limits
     motion(sat, rd, ll, odata);
     sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar);
     print_fit(sat, rd, ll, odata);
     sat.print_el();       // print new elements
   }
   if ((*buf & 0x5f) == 'Q') goto accept;
   goto motion;

bstar:

   printf("\n(A)uto  (b*)  (B)atch  (Q)uit  ");
   s_in(": ", buf);

   if (isdigit(buf[0]))
   {
     bstar = atof(buf);
     sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar);
     print_fit(sat, rd, ll, odata);
     sat.print_el();       // print new elements
   }
   if ((*buf & 0x5f) == 'A')
   {
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
           rms = find_rms(sat, rd, ll, odata);
           if(rms < sum)
           {
              sum = rms;
              bstar = bk;
           }
        }   // for bk
        sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar);
        bmin = bstar - bstep;
        bmax = bstar + bstep;
     }
 
     sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar);
     print_fit(sat, rd, ll, odata);
     sat.print_el();       // print new elements
   }
   if ((*buf & 0x5f) == 'B')
   {
     // create batch file
     fp = fopen("batch.inp", "w");
     fprintf(fp, "%s\n\n\n\n\n\n\n", "s");
     fprintf(fp, "%s\n", "q");
     fprintf(fp, "%s\n", "b");
     fprintf(fp, "%s\n\n\n", "a");
     fprintf(fp, "%s\n", "q");
     fprintf(fp, "%s\n\n\n\n", "s");
     fprintf(fp, "%s\n", "q");
     fprintf(fp, "%s\n", "b");
     fprintf(fp, "%s\n\n\n", "a");
     fprintf(fp, "%s\n", "q");
     fprintf(fp, "%s\n\n\n\n", "s");
     fprintf(fp, "%s\n", "q");
     fprintf(fp, "%s\n", "b");
     fprintf(fp, "%s\n\n\n", "a");
     fprintf(fp, "%s\n", "q");
     fprintf(fp, "%s\n\n\n\n", "s");
     fprintf(fp, "%s\n", "q");
     fprintf(fp, "%s\n", "b");
     fprintf(fp, "%s\n\n\n", "a");
     fprintf(fp, "%s\n", "q");
     fprintf(fp, "%s\n\n\n\n", "z");
     fprintf(fp, "%s\n", "q");
     fprintf(fp, "%s\n", "b");
     fprintf(fp, "%s\n\n\n", "a");
     fprintf(fp, "%s\n", "q");
     fprintf(fp, "%s\n\n\n\n", "z");
     fprintf(fp, "%s\n", "q");
     fprintf(fp, "%s\n", "b");
     fprintf(fp, "%s\n\n\n", "a");
     fprintf(fp, "%s\n", "q");
     fprintf(fp, "%s\n", "w");
     fprintf(fp, "%s\n", "u");
     fprintf(fp, "%s\n", "q");
     fprintf(fp, "%s\n", "q");
     fclose(fp);
     sprintf(file_in, "satfit %s ba<batch.inp", file);
     system(file_in);
     sprintf(file_in, "DEL batch.inp");
     system(file_in);
     sprintf(file_in, "satfit %s", file);
     system(file_in);
     exit(0);
   }
   if ((*buf & 0x5f) == 'Q') goto accept;

   goto bstar;

maneuver:

   printf("\nBoundary#  (P)erigee  (A)pogee  (O)b  (E)nergy  (R)estore  (Q)uit  ");
   s_in(": ", buf);
   if ((*buf & 0x5f) == 'Q')
   {
     tle = sat.tle;
     ii = sat.xincl / de2ra;
     om = mod(sat.xnodeo / de2ra, 360);
     ec = sat.eo;
     ww = mod(sat.omegao / de2ra, 360);
     ma = mod(sat.xmo / de2ra, 360);
     nn = sat.xno / nocon;
     bstar = sat.bstar;
     goto accept;
   }
   if (isdigit(buf[0]))
   {
     int p = atoi(buf);

     // store old obs
     for(i = 0; i < nobs; i++)
     {
       sprintf(iod_linex[i], "%s", iod_line[i]);
       llx[i] = ll[i];
       rdx[i] = rd[i];
       odatax[i] = odata[i];
     }
     for(i = p; i <= nobs; i++)
     {
       sprintf(iod_line[i-p], "%s", iod_line[i-1]);
       ll[i-p] = ll[i-1];
       rd[i-p] = rd[i-1];
       odata[i-p] = odata[i-1];
     }
     nobsx = nobs;
     nobs  = nobs - p + 1;

     out = 0;
     print_fit(sat, rd, ll, odata);
     sat.print_el();
     printf("\nperiod = %f days\n", nocon/sat.xno);
   }
   if ((*buf & 0x5f) == 'P')
   {
     double time;      // days

     printf("\n(P)revious  (N)ext  (Q)uit  ");
     s_in(": ", buf);
     if ((*buf & 0x5f) == 'Q') goto maneuver;
     if ((*buf & 0x5f) == 'P')
     {
       // if previously at perigee, back up one revolution
       if(ma < 0.1) time = satm.jd - nocon/satm.xno*(1 + satm.xmo/twopi);
       // if not previously at perigee, back up to perigee
       else         time = satm.jd - satm.xmo*nocon/(satm.xno*twopi);
     }
     // advance one revolution
     if ((*buf & 0x5f) == 'N')
     {
       // if previously at perigee, go forward one revolution
       if(ma > 359.9) time = satm.jd + nocon/satm.xno*(2 - satm.xmo/twopi);
       // if not previously at perigee, go forward to perigee
       else           time = satm.jd + nocon/satm.xno*(1 - satm.xmo/twopi);
     }
     // move to time and ma at perigee
     Date t2(time);
     satm.delta_t(t2.jd);
     satm.rv2el(satm.rr, satm.vv);
     satm.jd = t2.jd;
     ma = mod(satm.xmo / de2ra, 360);
     // refine perigee
     for(i = 0; i < 3; i++)
     {
       // go forward
       if(ma > 359.9) time = nocon/satm.xno*(satm.xmo/twopi - 1);
       // back up
       else time = satm.xmo*nocon/(satm.xno*twopi);
       Date t1(satm.jd - time);
       satm.delta_t(t1.jd);
       satm.jd = t1.jd;
       satm.rv2el(satm.rr, satm.vv);
       ma = mod(satm.xmo / de2ra, 360);
     }
     printf("\nPERIGEE");
     satm.print_el();       // print new elements

     // perigee residual
     double delr[3];
     time = sat.jd;                     // save sat epoch
     sat.delta_t(satm.jd);              // move sat to perigee
     vmadd(sat.rr, satm.rr, delr, -1);  // compare sat and satm perigees
     sat.delta_t(time);                 // restore sat epoch
     printf("\nperigee delta %5.0f\n", norm(delr)*6378.135);
   }
   if ((*buf & 0x5f) == 'A')
   {
     double time;      // days

     // time to travel from perigee to apogee, days
     printf("\n(P)revious  (N)ext  (Q)uit  ");
     s_in(": ", buf);
     if ((*buf & 0x5f) == 'Q') goto maneuver;
     if ((*buf & 0x5f) == 'P')
     {
       // if previously at or past apogee and past next perigee, back up to apogee
       if(ma < 180.1) time = satm.jd - 0.5*nocon/satm.xno*(1 + satm.xmo/pi);
       // if previously past apogee and before next perigee, back up to apogee
       else           time = satm.jd + 0.5*nocon/satm.xno*(1 - satm.xmo/pi);
     }
     // advance to apogee
     if ((*buf & 0x5f) == 'N')
     {
       // if previously at or past apogee and before perigee, go forward to apogee
       if(ma > 179.9) time = satm.jd + 0.5*nocon/satm.xno*(3 - satm.xmo/pi);
       // if previously past apogee and past next perigee, go forward to apogee
       else           time = satm.jd + 0.5*nocon/satm.xno*(1 - satm.xmo/pi);
     }
     // move time and satm.xmo to apogee
     Date t2(time);
     satm.delta_t(t2.jd);
     satm.jd = t2.jd;
     satm.rv2el(satm.rr, satm.vv);
     for(i = 0; i < 3; i++)
     {
       // loop to refine apogee, find when satm.xmo = pi
       Date t1(satm.jd + 0.5*nocon/satm.xno*(1 - satm.xmo/pi));
       satm.delta_t(t1.jd);
       satm.jd = t1.jd;
       satm.rv2el(satm.rr, satm.vv);
     }
     ma = mod(satm.xmo / de2ra, 360);
     printf("\nAPOGEE");
     satm.print_el();       // print new elements

     // apogee residual
     double delr[3];
     time = sat.jd;                     // save sat epoch
     sat.delta_t(satm.jd);              // move sat to apogee
     vmadd(sat.rr, satm.rr, delr, -1);  // compare sat and satm apogees
     sat.delta_t(time);                 // restore sat epoch
     printf("\napogee delta %5.0f\n", norm(delr)*6378.135);  // kilometers
   }
   if ((*buf & 0x5f) == 'O')
   {
     // move up one to make room for pseudo ob
     for(i = nobs-1; i >= 0; i--)
     {
       sprintf(iod_line[i+1], "%s", iod_line[i]);
       ll[i+1] = ll[i];
       rd[i+1] = rd[i];
       odata[i+1] = odata[i];
     }
     nobs++;
     // ll is unit vector in same direction as satm.rr
     smult(1./norm(satm.rr), satm.rr, ll[0]);
     // rd is unit vector (1 er) in same direction as satm.rr
     smult(1./norm(satm.rr), satm.rr, rd[0]);
     // odata holds epoch, ra, dc, obscode data
     double ra, dc;
     odata[0][0] = satm.jd;
     zrll(satm, rd[0], ra, dc);
     odata[0][1] = ra*de2ra;
     odata[0][2] = dc*de2ra;
     odata[0][3] = 0000;

     // print pseuso ob
     char desig[5];
     int norad, yy;
     double rm, dm;

/*
          1         2         3         4         5         6
0123456789012345678901234567890123456789012345678901234567890
36868 10 039A   2018 G 20100909203108300 17 25 2216290+034350 57
*/

    norad = atoi(iod_line[1]);
    yy    = atoi(iod_line[1] + 6);
    sscanf(iod_line[1] + 9, "%4s", &desig);
    desig[4] = '\0';
    Date t1(satm.jd);
    ra /= 15;
    rm = modf(ra, &rm)*60000;
    dm = fabs(modf(dc, &dm)*6000);
    printf("\nPSEUDO OB:");
    printf("\n%05d %02d %4s   0000 P %4.0f%02.0f%02.0f%02.0f%02.0f%05.0f 16 25 %02.0f%05.0f%+03.0f%04.0f\n",
            norad, yy, desig, t1.yy, t1.mm, t1.dd, t1.hr,
            t1.mn, t1.ss*1000, ra, rm, dc, dm);
    sprintf(iod_line[0], "%05d %02d %4s   0000 P %4.0f%02.0f%02.0f%02.0f%02.0f%05.0f 16 25 %02.0f%05.0f%+03.0f%04.0f\n",
            norad, yy, desig, t1.yy, t1.mm, t1.dd, t1.hr,
            t1.mn, t1.ss*1000, ra, rm, dc, dm);


   }
   if ((*buf & 0x5f) == 'E')
   {
     s_in("\ndelta specific energy(m^2/s^2): ", buf);
     if (isdigit(buf[0]))
     {
       double dE = atof(buf);
       dE /= 11300.168;   // er^2 / min^2
       double E2, VV, dV[3], vec[3];
       cross(satm.rr, satm.vv, vec);
       double mu = norm(vec);
       mu = mu*mu;
       mu = mu / satm.aodp;
       mu = mu / (1 - satm.eo*satm.eo);
       E2 = -0.5*mu / satm.aodp;
       VV = sqrt(2*(E2 + dE + mu/norm(satm.rr)));  // new velocity magnitude
       unitv(satm.vv, dV);         // velocity direction
       smult(VV, dV, vec);         // new velocity vector
       sat = satm;
       sat.rv2el(sat.rr, vec);     // State vectors to mean elements
       sat.print_el();             // print new elements
     }
   }
   if ((*buf & 0x5f) == 'R')
   {
     printf("\nRestore (O)bs  (E)lements  (Q)uit  ");
     s_in(": ", buf);
     // re-store old obs
     if ((*buf & 0x5f) == 'O')
     {
       nobs = nobsx;
       for(i = 0; i < nobs; i++)
       {
         sprintf(iod_line[i], "%s", iod_linex[i]);
         ll[i] = llx[i];
         rd[i] = rdx[i];
         odata[i] = odatax[i];
       }
     }
     // re-store old TLE
     if ((*buf & 0x5f) == 'E')
     {
       // replace working elements with original
       satm = save_sat;        // original elements maneuverd to a node
       sat  = save_sat;        // new fit elements
       sat.print_el();         // print original elements
     }
   }

   goto maneuver;

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
      save_sat.print_el();       // print new elements
   }
   if ((*buf & 0x5f) == 'R')
   {
      // replace TLE in file with original
      FILE *fp;
      fp = fopen(file, "w");
      fprintf(fp, "%s", name);
      fprintf(fp, "%s", tle1);
      fprintf(fp, "%s\n", tle2);
      for(int i = 0; i < nobs; i++)
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

      len = strlen(name);
      if (name[len-1] == '\n')
          name[len-1] = '\0';
      pgee = (1-ec)*sat.aodp;
      agee = (1+ec)*sat.aodp;
      sprintf(buf, "%d X %d km", (int)((pgee-1)*6378.135),
                                 (int)((agee-1)*6378.135));
      fp = fopen(file, "w");
      fprintf(fp, "%-50s%19s\n", name, buf);
      fprintf(fp, "%s\n", line1);
      fprintf(fp, "%s\n\n", line2);
      for(int i = 0; i < nobs; i++)
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
   sat.delta_t(t1.jd);             // advance to node
   sat.rv2el(sat.rr, sat.vv);      // sgp4 elements
   tle = t1.tle;
   sat.jd = t1.jd;

   ii = sat.xincl / de2ra;
   om = mod(sat.xnodeo / de2ra, 360);
   ec = sat.eo;
   ww = mod(sat.omegao / de2ra, 360);
   ma = mod(sat.xmo / de2ra, 360);
   nn = sat.xno / nocon;
   bstar = sat.bstar;
   sat.print_el();       // print new elements
   sum = find_rms(sat, rd, ll, odata);
   printf("\nrms%12.5f\n", sum);

   longitude();

   goto accept;

}    // end main

////////////////// FUNCTIONS //////////////////////////////////////////////////


//  asymptotic extrapolation
inline double asym(double a1, double a2, double a3)
{
  double b;
  if(fabs(a1 - a2) < 1.e-5) b = 0.;
  else b = (a2 - a3) / (a1 - a2);
  return (b * a2 - a3) / (b - 1.);
}

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
// if no tle found, calls elfind to create one
void read_tle(char *file)
{
  FILE *fp;
  reread:
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
           printf("\n%s", name);    // print name string above line1
        }
        fclose(fp);

        // first data line
        printf("%s", line1);         // print first line on screen
        ssn = atoi(line1 + 2);       // satellite number
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
      strncpy(name, line1, 50);
      len = strlen(name);
      name[len-1] = '\n';
    }
  }while(fgets(line1, 80, fp));

  if(nobs > 1)
  {
    char file_in[81];
    sprintf(file_in, "elfind %s", file);
    system(file_in);
    goto reread;
  }
  else exit(0);
}

// load IOD lines from input file_in
void get_obs(char *file, char iod_line[][81])
{
  FILE *fp;
  if ((fp = fopen(file, "r")) == NULL)
  {
     printf("Can't open %s\n", file);
     s_in("[exit]", buf);
     exit(0);
  }
  nobs = 0;

  // find observation entries and write to iod_line matrix
  while(fgets(iod_line[nobs], 80, fp))
  {
    if(strlen(iod_line[nobs]) > 58 && strlen(iod_line[nobs]) != 70 &&
       isdigit(iod_line[nobs][0]))
    {
      nobs++;
      if(nobs > 200)
      {
         printf("Exceeding 200 Observation Maximum \n");
         s_in("[continue]", buf);
         return;
      }
    }
  }
  fclose(fp);
  return;
}

// read site data, subroutine of read_obs
void read_site(char line[])
{
  int site, num;
  char* fl;
  char inp_str[81];
  sscanf(line + 16, "%4d", &site);    // read site number from iod_line
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
            double ll[][3],
            double rd[][3])
{
  int i, j;
  double k[3], m[4];
  sort:
  j = 0;
  for(i = 0; i < nobs - 1; i++)
  {
    if(odata[i][0] > odata[i + 1][0])    // switch low/high
    {
      k = rd[i];
      rd[i] = rd[i + 1];
      rd[i + 1] = k;
      k = ll[i];
      ll[i] = ll[i + 1];
      ll[i + 1] = k;
      m = odata[i];
      odata[i] = odata[i + 1];
      odata[i + 1] = m;
      buf = iod_line[i];
      iod_line[i] = iod_line[i + 1];
      iod_line[i + 1] = buf;
      j++;
    }
    // remove duplicates
    if(odata[i][0] == odata[i + 1][0] && odata[i][1] == odata[i + 1][1])
    {
      for(int s = i; s < nobs - 1; s++)
      {
         rd[s] = rd[s + 1];
         ll[s] = ll[s + 1];
         odata[s] = odata[s + 1];
         iod_line[s] = iod_line[s + 1];
      }
      nobs--;
    }
  }
  if(j > 0)  goto sort;
}

// decodes the iod_line data
void read_obs(char iod_line[][81],
            double odata[][4],
            double ll[][3],
            double rd[][3])
{
  int year, month, day, hour, min, sign, format, age, epoch, obscode;
  double sec, ra, mm, ss, dc, dd, ds;
  double t, a, b, c, dcx;
  char sn[1], sbuf[5], dbuf[2];

  for(int i = 0; i < nobs; i++)
  {
    sscanf(iod_line[i] + 16, "%4d", &obscode);
    sscanf(iod_line[i] + 23, "%4d %2d %2d %2d %2d %5s",
          &year, &month, &day, &hour, &min, &sbuf);
    sec = atof(sbuf);
    sec /= pow(10, strlen(sbuf) - 2);
    sscanf(iod_line[i] + 44, "%1d", &format);
    sscanf(iod_line[i] + 45, "%1d", &age);
    epoch = age;
    sscanf(iod_line[i] + 47, "%2lf %2lf %3lf %1s", &ra, &mm, &ss, &sn);
    sscanf(iod_line[i] + 55, "%2lf %2lf %2s", &dc, &dd, &dbuf);
    ds = atof(dbuf);
    if(strlen(dbuf) == 1) ds *= 10;
    sign = (sn[0] == '-') ? -1 : 1 ;

    switch(format)      // change to HHhhhh, DDdddd
    {
      /* Format 1: RA/DEC = HHMMSSs+DDMMSS */
      case 1 : ra += mm / 60 + ss / 36000;
               dc  = sign * (dc + dd / 60 + ds / 3600); break;
      /* Format 2: RA/DEC = HHMMmmm+DDMMmm */
      case 2 : ra += mm / 60 + ss / 60000;
               dc  = sign * (dc + dd / 60 + ds / 6000); break;
      /* Format 3: RA/DEC = HHMMmmm+DDdddd */
      case 3 : ra += mm / 60 + ss / 60000;
               dc  = sign * (dc + dd / 100 + ds / 10000); break;
      /* Format 7: RA/DEC = HHMMSSs+DDdddd */
      case 7 : ra += mm / 60 + ss / 36000;
               dc  = sign * (dc + dd / 100 + ds / 10000); break;
      default : printf("\nFormat not found\n");
                s_in("[exit]", buf);
                exit(0);
    }

    Date t1(year, month, day, hour, min, sec);
    // change from hours and degrees to radians
    ra *= 15;
    ra *= de2ra;
    dc *= de2ra;

    double csi = .0055878713278878,
           zet = .0055888307019922,
           the = .0048580335354883;

    if(epoch == 4)   // precess from B1950 to J2000
    {
      a = cos(dc) * sin(ra + csi);
      b = cos(the) * cos(dc) * cos(ra + csi)
          - sin(the) * sin(dc);
      c = sin(the) * cos(dc) * cos(ra + csi)
          + cos(the) * sin(dc);
      ra = atan(a / b);        // ra - zet
      if(b < 0) ra += pi;
      ra += zet;               // right ascension, radians
      ra += 1 / 30000;
      dc = asin(c);
      if(fabs(dc) > 1.4) dc = c / fabs(c) * acose(sqrt(a*a + b*b));
    }

    if(epoch != 0)
    {
      // precession from J2000
      t = (t1.jd - 2451545) / 36525;
      csi = (2306.2181 +  .30188 * t + .017998 *t*t) * t * de2ra / 3600;
      zet = (2306.2181 + 1.09468 * t + .018203 *t*t) * t * de2ra / 3600;
      the = (2004.3109 -  .42665 * t - .041833 *t*t) * t * de2ra / 3600;
      a = cos(dc) * sin(ra + csi);
      b = cos(the) * cos(dc) * cos(ra + csi)
          - sin(the) * sin(dc);
      c = sin(the) * cos(dc) * cos(ra + csi)
          + cos(the) * sin(dc);
      ra = atan(a / b);          // ra - zet
      if(b < 0) ra += pi;
      ra += zet;                 // right ascension, radians
      dc = asin(c);
      if(fabs(dc) > 1.4) dc = c / fabs(c) * acose(sqrt(a*a + b*b));
    }

    // line-of-sight vectors
    ll[i][0] = cos(dc) * cos(ra);
    ll[i][1] = cos(dc) * sin(ra);
    ll[i][2] = sin(dc);

    odata[i][0] = t1.jd;                   // julian date
    odata[i][1] = ra;                      // ra radians (observed)
    odata[i][2] = dc;                      // dc radians (observed)
    odata[i][3] = obscode;                 // observer code

    read_site(iod_line[i]);
    Locate rre(t1.jd, la, lo, hh);
    rd[i] =  rre.rre;
  } // end for
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
  while(fabs(e - ek) > 1e-6)
  {
     ek = e;
     e = ma + ec * sin(e);
  }
  ma /= de2ra;
  theta = (ec - cos(e)) / (ec * cos(e) - 1);
  theta = acose(theta);
  if(e > pi) theta = 2 * pi - theta;
  uu = ww + theta / de2ra;
}
// find mean anomaly from true longitude and perigee
double ww2ma(double wx)
{
   theta = mod(uu - wx, 360) * de2ra;
   e = acose((ec + cos(theta)) / (1 + ec * cos(theta)));
   if(theta > pi) e = 2 * pi - e;
   ma = e - ec * sin(e);
   ma /= de2ra;
   return ma;
}

//   gaussian elimination
void rref(double m[][7], double b[])    // scaled partial pivoting
{
  int i, j, k, ROW;
  double s[6], bin, mult;

////////////////// calculate permanent scale factors //////////////////////////

  for(i = 0; i < 6; i++)
  {
    s[i] = abs(m[i][0]);
    for(j = 1; j < 6; j++)
      if(s[i] < abs(m[i][j]))
      {
        s[i] = abs(m[i][j]);
      }
  }  // end for i

///////////////////// swap rows according to scale ////////////////////////////

  for(j = 0; j < 5; j++)
  {
    ROW = j;
    for(i = j + 1; i < 6; i++)
    {
      if(abs(m[ROW][j] / s[ROW]) < abs(m[i][j] / s[i]))
        ROW = i;
    }
    if(ROW != j)
    {
      for(k = j; k < 7; k++)     // swap rows
      {
        bin = m[j][k];
        m[j][k] = m[ROW][k];
        m[ROW][k] = bin;
      }
      bin = s[j];                 // swap scales
      s[j] = s[ROW];
      s[ROW] = bin;
    }  // end if
  }

///////////////////// forward elimination /////////////////////////////////////

  for(j = 0; j < 6; j++)
  {
    mult = m[j][j];
    if(mult != 0.)
    {
      for(k = j; k < 7; k++)
      {
        m[j][k] = m[j][k] / mult;
      }
    }
    for(i = j + 1; i < 6; i++)
    {
      mult = m[i][j];
      for(k = j + 1; k < 7; k++)
      {
        m[i][k] = m[i][k] - mult * m[j][k];
      }
      m[i][j] = 0;
    }  // end for i
  }  // end for j

///////////// test for singular matrix ///////////////////////////////////////

  bin = 1;
  for(i = 0; i < 6; i++)
  {
    bin *= m[i][i];
  }
  if(bin == 0)
  {
    printf("Singular matrix");
    exit(0);
  }

////////////////// back sustitution //////////////////////////////////////////

  b[5] = m[5][6] / m[5][5];
  for(i = 5; i >= 0; i--)
  {
    bin = 0;
    for(k = i + 1; k < 6; k++)
      bin = bin + m[i][k] * b[k];
    b[i] = (m[i][6] - bin) / m[i][i];
  }
}  // end rref

/* Calculates predicted direction angles, right ascension(ra) and
declination(dc), of the line of sight vector, rll, in equatorial coordinates
from an observer position, rd, to the satellite position, sat.rr.
A subroutine of diff_el.  Satellite object, sat, simplifies prediction */
void zrll(Satellite satx,  // perturbed Satellite object at TLE
          double* rd,      // input, topocentric vector to observer
          double& ra,      // calculated right ascension, degrees
          double& dc)      // calculated declination, degrees
{
  double rll[3];                 // declare line of sight vector
  vmadd(satx.rr, rd, rll, -1);   // rll = satx.rr - rd
  // return ra and dc, degrees
  dc = asin(rll[2] / norm(rll)) / de2ra;
  ra = atan(rll[1] / rll[0]) / de2ra;
  if(rll[0] < 0) ra += 180.;
  ra = mod(ra, 360);
}

//////////////////////////////////////////////////////

// differential correction of elements
void diff_el(Satellite sat,  double rd[][3],
            double ll[][3],  double odata[][4])
{

   int i, j, k, c = 0;
   double ra, dc;                 // right ascension, declination variables
   double delta, el;              // small change amount
   double mdata[150][7];          // define matrix for multiple regression
   double dat[2][7];              // define output matrix of zpde utility
   double b[6];                   // output deltas, sat position, velocity, b*
   double rv[6][7];               // normal matrix
   double rac, dcc, rms;

   loop:

   // begin the main loop of the program
   for(i = 0; i < nobs; i++)
   {
     Satellite satx = sat;
     Satellite satz = sat;             // differential sat
     satx.delta_t(odata[i][0]);        // relocate satellite at new time

     // first establish the computed ra, dc, at jdo with no perturbations
     zrll(satx, rd[i], ra, dc);        // output ra, dc, degrees
     rac = ra;                         // store computed ra and dc
     dcc = dc;

     // find the deltas and load into output matrix, dat
     dat[0][6] = rtw(odata[i][1]/de2ra, rac);      // store delta_ra
     dat[1][6] = odata[i][2]/de2ra - dcc;          // store delta_dc

     // 6 steps to build differential correction matrix
     j = 0;
     if(xi)
     {
        dat[0][j] = .001;
        dat[1][j] = .001;
     }
     else
     {
        delta = .001;                         // change
        el = satz.xincl;                      // store reference
        satz.xincl += delta;                  // delta element
        satz.delta_t(odata[i][0]);            // recalculate with perturbed element
        zrll(satz, rd[i], ra, dc);            // perturbed ra, dc
        satz.xincl = el;                      // restore reference
        dat[0][j] = rtw(ra, rac) / delta;     // perturbed - computed
        dat[1][j] = (dc - dcc) / delta;
     }

     j = 1;
     delta = .001;                         // change
     el = satz.xnodeo;                     // store reference
     satz.xnodeo += delta;                 // delta element
     satz.delta_t(odata[i][0]);            // recalculate with perturbed element
     zrll(satz, rd[i], ra, dc);            // perturbed ra, dc
     satz.xnodeo = el;                     // restore reference
     dat[0][j] = rtw(ra, rac) / delta;     // perturbed - computed
     dat[1][j] = (dc - dcc) / delta;

     j = 2;
     if(xe)
     {
        dat[0][j] = .00001;
        dat[1][j] = .00001;
     }
     else
     {
        delta = .0001;                        // change
        el = satz.eo;                         // store reference
        satz.eo += delta;                     // delta element
        satz.delta_t(odata[i][0]);            // recalculate with perturbed element
        zrll(satz, rd[i], ra, dc);            // perturbed ra, dc
        satz.eo = el;                         // restore reference
        dat[0][j] = rtw(ra, rac) / delta;     // perturbed - computed
        dat[1][j] = (dc - dcc) / delta;
     }

     j = 3;
     if(xw)
     {
        dat[0][j] = .001;
        dat[1][j] = .001;
     }
     else
     {
        delta = .001;                         // change
        el = satz.omegao;                     // store reference
        satz.omegao += delta;                 // delta element
        satz.delta_t(odata[i][0]);            // recalculate with perturbed element
        zrll(satz, rd[i], ra, dc);            // perturbed ra, dc
        satz.omegao = el;                     // restore reference
        dat[0][j] = rtw(ra, rac) / delta;     // perturbed - computed
        dat[1][j] = (dc - dcc) / delta;
     }

     j = 4;
     delta = .001;                         // change
     el = satz.xmo;                        // store reference
     satz.xmo += delta;                    // delta element
     satz.delta_t(odata[i][0]);            // recalculate with perturbed element
     zrll(satz, rd[i], ra, dc);            // perturbed ra, dc
     satz.xmo = el;                        // restore reference
     dat[0][j] = rtw(ra, rac) / delta;     // perturbed - computed
     dat[1][j] = (dc - dcc) / delta;

     j = 5;
     if(xn)
     {
        dat[0][j] = .000001;
        dat[1][j] = .000001;
     }
     else
     {
        delta = .00001;                       // change
        el = satz.xno;                        // store reference
        satz.xno += delta;                    // delta element
        satz.delta_t(odata[i][0]);            // recalculate with perturbed element
        zrll(satz, rd[i], ra, dc);            // perturbed ra, dc
        satz.xno = el;                        // restore reference
        dat[0][j] = rtw(ra, rac) / delta;     // perturbed - computed
        dat[1][j] = (dc - dcc) / delta;
     }

     mdata[2 * i]     = dat[0];        // numerical deltas transferred to
     mdata[2 * i + 1] = dat[1];        // multiple regresssion matrix
   }   // end for i

   // multiple regression
   for(j = 0; j < 6; j++)
   {
     for(k = 0; k < 7; k++)
     {
       rv[j][k] = 0;
       for(i = 0; i < nobs*2; i++)
       {
         rv[j][k] = rv[j][k] + mdata[i][j] * mdata[i][k];
       }
     }
   }

   rref(rv, b);

   Satellite saty = sat;
   // test update components with deltas
   saty.xincl  += b[0]*.1;
   saty.xnodeo += b[1]*.1;
   saty.eo     += b[2]*.1;
   saty.omegao += b[3]*.1;
   saty.xmo    += b[4]*.1;
   saty.xno    += b[5]*.1;
   rms = find_rms(saty, rd, ll, odata);
   if(rms < sum)
   {
     sum = rms;
     // update components with deltas
     sat.xincl  += b[0]*.1;
     sat.xnodeo += b[1]*.1;
     sat.eo     += b[2]*.1;
     sat.omegao += b[3]*.1;
     sat.xmo    += b[4]*.1;
     sat.xno    += b[5]*.1;
     c++;
     if(c < 20) goto loop;
   }
/*
   // display computed deltas for all 6 components and an rms
   for(i = 0; i < 6; i++)
   {
     printf("\n");
     printf("   %9.6f", b[i]);
   }
   printf("\n\nrms%9.6f\n", rms);
   s_in(": ", buf);
*/
   // global elements updated
   ii = sat.xincl/de2ra;
   om = sat.xnodeo/de2ra;
   ec = sat.eo;
   ww = sat.omegao/de2ra;
   ma = sat.xmo/de2ra;
   nn = sat.xno/nocon;
}

//////////////////////////////////////////////////////

// box search, no limits
void anomaly(Satellite sat,  double rd[][3],
            double ll[][3],  double odata[][4])
{
   double min, max, step;
   step = .1;
   mk = ma;
   sum = find_rms(sat, rd, ll, odata);
   do
   {
      min = mk;
      max = mk;
      nsum:
      min = mk - step;
      sat.delta_el(sat.jd, ii, om, ec, ww, min, nn, bstar);
      nsum = find_rms(sat, rd, ll, odata);
      if(nsum < sum)
      {
         mk = min;
         sum = nsum;
         goto nsum;
      }
      xsum:
      max = mk + step;
      sat.delta_el(sat.jd, ii, om, ec, ww, max, nn, bstar);
      xsum = find_rms(sat, rd, ll, odata);
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
void motion(Satellite sat,  double rd[][3],
           double ll[][3],  double odata[][4])
{
   double min, max, step;
   step = .1;
   nk = nn;
   sum = find_rms(sat, rd, ll, odata);
   do
   {
      min = nk;
      max = nk;
      nsum:
      min = nk - step;
      sat.delta_el(sat.jd, ii, om, ec, ww, ma, min, bstar);
      nsum = find_rms(sat, rd, ll, odata);
      if(nsum < sum)
      {
         nk = min;
         sum = nsum;
         goto nsum;
      }
      xsum:
      max = nk + step;
      sat.delta_el(sat.jd, ii, om, ec, ww, ma, max, bstar);
      xsum = find_rms(sat, rd, ll, odata);
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

// partition search on node and inclination within set limits
void node(Satellite sat,  double rd[][3],
         double ll[][3],  double odata[][4])
{
   double rms;
   while((imax - imin) > 1e-5)
   {
      istep = (imax - imin) / 20;
      ostep = fabs(rtw(omax, omin) / 20);
      for(ik = imin; ik < imax; ik += istep)
      {
         if(xi)
         {
            imin = ii;
            imax = ii;
            ik = ii;
            istep = 0;
         }
         for(ok = omin; ok < omax; ok += ostep)
         {
            Satellite satx(tle, ik, ok, ec, ww, ma, nn, bstar);
            // establish the computed ra, dc, at jdo with no perturbations
            rms = find_rms(satx, rd, ll, odata);
            if(rms < sum)
            {
               sum = rms;
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

// partition search on perigee and eccentricity
void perigee(Satellite sat,  double rd[][3],
            double ll[][3],  double odata[][4])
{
   double rms;
   if(ec > .1)
   {
      wmax = ww + .1;
      wmin = ww - .1;
      emax = ec * 1.01;
      emin = ec * .99;
   }
   while((wmax - wmin) > 1e-5)
   {
      estep = (emax - emin) / 20;
      wstep = (wmax - wmin) / 20;
      for(wk = wmin; wk < wmax; wk += wstep)
      {
         if(xw)
         {
            wmin = ww;
            wmax = ww;
            wk = ww;
            wstep = 0;
         }
         theta = (uu - wk) * de2ra;
         for(ek = emin; ek < emax; ek += estep)
         {
            if(xe)
            {
               emin = ec;
               emax = ec;
               ek = ec;
               estep = 0;
            }
            e = acose((ek + cos(theta)) / (1 + ek * cos(theta)));
            if(theta > pi) e = 2 * pi - e;
            mk = e - ek * sin(e);
            mk = mk / de2ra;

            Satellite satx(sat.jd, ii, om, ek, wk, mk, nn, bstar);
            // establish the computed ra, dc, at jdo with no perturbations
            rms = find_rms(satx, rd, ll, odata);
            if(rms < sum)
            {
               sum = rms;
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
   anomaly(sat, rd, ll, odata);
   sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar);

   // update mean_motion
   if(!xn)
   {
      motion(sat, rd, ll, odata);
      sat.delta_el(sat.jd, ii, om, ec, ww, ma, nn, bstar);
   }

   // calculate uu, degrees
   longitude();

   ww = mod(ww, 360);
   ma = mod(ma, 360);
}

/* find geocentric vector, rr, to satellite given the topocentric unit
vector, ll, from the observer to the satellite, the geocentric vector, rd,
to the observer, and the length, r, of the geocentric vector sat.rr. */
void so2r(double r, double rd[], double ll[], double rr[])
{
  double ang1, ang2, nrd, rho;

  nrd = norm(rd);
  ang1 = acose(dot(rd, ll) / nrd);
  if(ang1 < .001) rho = r - nrd;
  else
  {
    ang2 = asin(nrd * sin(ang1) / r);
    rho = r * sin(ang1 - ang2) / sin(ang1);
  }
  vmadd(rd, ll, rr, rho);
}

// sets deltaT error at last obs epoch to zero
void align(Satellite sat,  double rd[][3],
         double ll[][3],  double odata[][4])
{
   double mk, nrr, nvv, Perr, delt, rr[3], delr[3];
   int last = nobs - 1;

   Satellite satx = sat;        // copy sat

   do
   {
     // advance satellite position
     satx.delta_t(odata[last][0]);
     nrr = norm(satx.rr);
     nvv = norm(satx.vv);

     // observed geocentric position vector, rr
     so2r(nrr, rd[last], ll[last], rr);

     // position error in radians
     Perr = acose(dot(satx.rr, rr) / (nrr*nrr));

     // difference between computed and observed position vectors, er
     vmadd(satx.rr, rr, delr, -1);

     // magnitude of delta r in direction of v, radians
     delt = Perr * dot(satx.vv, delr) / (nvv * norm(delr));

     // degrees
     delt /= de2ra;

     if(first)
     {
       first = 0;
       ma = ma - 0.75*delt;
     }
     else ma = ma - delt/5.;
     satx.delta_el(satx.jd, ii, om, ec, ww, ma, nn, bstar);
   }while(fabs(delt) > 1.E-5);
}

void print_fit(Satellite sat,  double rd[][3],
              double ll[][3],  double odata[][4])
{
  double nrr, nvv, Perr, delt, xtrk, az, el, asp, alpha, sum = 0,
         tempv[3], temp[3], rr[3], nv[3], delr[3], zz[3] = {0, 0, 1};
  int yy, day, hh, mm, ss, sign;
  if(nobs == 0)
  {
     printf("\nno obs");
     return;
  }

  // copy sat
  Satellite satx = sat;

  if(out)
  {
      FILE *fp;
      fp = fopen(file, "a");
      fprintf(fp, "\n\n      STA  YYday HHMM:SSsss   AZ     EL     ASP     XTRK    deltaT   Perr\n");
      fclose(fp);
  }

  printf("\n      STA  YYday HHMM:SSsss   AZ     EL     ASP     XTRK    deltaT   Perr\n");
  for(int j = 0; j < nobs; j++)
  {
     // advance satellite position
     satx.delta_t(odata[j][0]);
     nrr = norm(satx.rr);
     nvv = norm(satx.vv);

     //  computing elevation
     el = acose(dot(rd[j], ll[j]) / norm(rd[j]));
     el /= de2ra;
     el = 90 - el;
     el = round(el, 1);

     // computing aspect
     asp = acose(dot(ll[j], satx.vv) / nvv);
     asp /= de2ra;
     asp = 180 - asp;
     asp = round(asp, 1);

     // computing azimuth
     cross(rd[j], zz, tempv);
     cross(tempv, rd[j], nv);
     cross(rd[j], ll[j], tempv);
     cross(tempv, rd[j], temp);
     az = acose(dot(nv, temp) / (norm(nv)*norm(temp)));
     cross(temp, nv, tempv);
     if(dot(tempv, rd[j]) < 0) az = 2*pi - az;
     if(norm(temp) == 0) az = -0.0;
     az /= de2ra;
     az  = round(az, 1);

     // observed satellite geocentric position vector, rr
     so2r(nrr, rd[j], ll[j], rr);

     // geocentric position error angle in radians
     Perr = acose(dot(satx.rr, rr) / (nrr*nrr));

     // difference between computed and observed position vectors, delr, in e.r.
     vmadd(satx.rr, rr, delr, -1);
     cross(satx.rr, satx.vv, temp);  // xtrk reference vector points left of track
     sign = dot(delr, temp) < 0 ? -1 : 1;

     // observer velocity vector
     cross(zz, rd[j], tempv);
     smult(.004351409367, tempv, temp);
     // observed satellite velocity
     vmadd(satx.vv, temp, tempv, -1);
     nvv = norm(tempv);

     // angle between delr vector and tempv vector, radians
     alpha = acose(dot(tempv, delr) / (nvv * norm(delr)));

     // magnitude of delr in direction of tempv, radians
     delt = atan(cos(alpha) * tan(Perr));   // geocentric range error

     // time error
     delt *= nrr / nvv;                     // delta r in min
     delt *= 60;                            // seconds
     delt  = round(delt, 2);

     // predicted topocentric coordinates (sat xyz - observer xyz)
     // new use of delr variable, predicted line of sight vector
     vmadd(satx.rr, rd[j], delr, -1);
     nrr = norm(delr);

     // convert to unit vector
     smult(1/nrr, delr, delr);

     // topocentric position error angle in radians
     Perr = acose(dot(delr, ll[j]));

     // cross track error, as component of topocentric Perr
     xtrk  = asin(sin(alpha) * sin(Perr));  // cross track magnitude, radians
     xtrk /= de2ra;                         // degrees
     xtrk *= sign;                          // left of track is positive
     xtrk  = round(xtrk, 2);

     // sum position error in squared degrees
     Perr /= de2ra;
     sum += Perr*Perr;

     yy = (int)satx.yy - 2000;
     day = (int)satx.doy;
     hh = (int)satx.hr;
     mm = (int)satx.mn;
     ss = (int)((satx.ss + .0001) * 1000);
     if(ss >= 60000)
     {
        ss = (int)mod(ss, 60000);
        mm += 1;
     }
     if(mm >= 60)
     {
        mm = (int)mod(mm, 60);
        hh += 1;
     }
     if(hh >= 24)
     {
        hh = (int)mod(hh, 24);
        day += 1;
     }


     printf("(%2d) %04d  %02d%03d %02d%02d:%05d  %5.1f  %5.1f  %5.1f  %6.2f   %6.2f  %7.3f\n",
     j + 1, (int)odata[j][3], yy, day, hh, mm, ss, az, el, asp, xtrk, delt, Perr);

     // print fit to file
     if(out)
     {
        FILE *fp;
        fp = fopen(file, "a");
        fprintf(fp, "(%2d) %04d  %02d%03d %02d%02d:%05d  %5.1f  %5.1f  %5.1f  %6.2f   %6.2f  %7.3f\n",
          j + 1, (int)odata[j][3], yy, day, hh, mm, ss, az, el, asp, xtrk, delt, Perr);
        fclose(fp);
     }
  }
  printf("\nrms%12.5f\n", sqrt(sum / nobs));
}

double find_rms(Satellite sat,  double rd[][3],
               double ll[][3],  double odata[][4])
{
   double nrr, Perr, delr[3], zum = 0;

   // copy sat
   Satellite satx = sat;

   for(int j = 0; j < nobs; j++)
   {

      // advance satellite position
      satx.delta_t(odata[j][0]);

      // predicted topocentric coordinates (sat xyz - observer xyz)
      vmadd(satx.rr, rd[j], delr, -1);
      nrr = norm(delr);       // range

      // convert to unit vector
      smult(1/nrr, delr, delr);

      // topocentric position error in degrees
      Perr = acose(dot(delr, ll[j])) / de2ra;

      // sum position error in squared degrees
      zum += Perr*Perr;
   }
   return sqrt(zum / nobs);
}

