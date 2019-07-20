// sst.h
// by Scott Campbell campbel7@hughes.net


/*************************** Class Satellite **********************************
*                                                                             *
*                            Two Constructors                                 *
*                                                                             *
*                                                                             *
*  Satellite sat(tle_date,                                                    *
*                inclination_degrees,                                         *
*                ascending_node_degrees,                                      *
*                eccentricity,                                                *
*                argument_of_perigee_degrees,                                 *
*                mean_anomaly_degrees,                                        *
*                mean_motion_revolutions_per_day,                             *
*                bstar)                                                       *
*     Eight argument constructor initializes satellite at epoch tle_date.     *
*                                                                             *
*                                                                             *
*  Satellite sat(tle_date, position_pointer, velocity_pointer, bstar)         *
*     Four argument constructor produces mean(SGP4)orbital elements           *
*                                                                             *
*                                                                             *
*                      Seven Member Functions                                 *
*                                                                             *
*  sat.delta_t(time);                                                         *
*  sat.delta_sst(time);                                                       *
*                Updates satellite elements to new position at time.  Time    *
*                is either Julian day or TLE date. Program decides which      *
*                format based on relative magnitude.                          *
*                                                                             *
*  sat.delta_el(epoch_Julian date,                                            *
*               inclination_degrees,                                          *
*               ascending_node_degrees,                                       *
*               eccentricity,                                                 *
*               argument_of_perigee_degrees,                                  *
*               mean_anomaly_degrees,                                         *
*               mean_motion_degrees_per_day,                                  *
*               bstar)                                                        *
*         Change elements for an existing satellite object.                   *
*                                                                             *
*  sat.rv2el(rr, vv);  State vectors to mean elements                         *
*                                                                             *
*  sat.print_el();     Output mean elements to screen                         *
*                                                                             *
*  sat.print_rv();     Output state vectors to screen                         *
*                                                                             *
*  sat.convert();      Mean elements to Equinocital elements                  *
*                                                                             *
*                                                                             *
*                                                                             *
*                        Twelve Output Values                                 *
*                                                                             *
*  sat.rr;  Topocentric position vector to satellite at epoch, earth radii.   *
*                                                                             *
*  sat.vv;  Satellite velocity vector at epoch, earth radii / minute.         *
*                                                                             *
*  sat.jd;  Julian date epoch of satellite elements.                          *
*                                                                             *
*  sat.thetag;  Greenwich Hour Angle, degrees.                                *
*                                                                             *
*  sat.xincl;   Mean inclination, radians.                                    *
*                                                                             *
*  sat.xnodeo;  Mean longitude of the ascending node, radians.                *
*                                                                             *
*  sat.eo;      Mean eccentricity.                                            *
*                                                                             *
*  sat.omegao;  Mean argument of perigee, radians.                            *
*                                                                             *
*  sat.xmo;     Mean anomaly, radians.                                        *
*                                                                             *
*  sat.xno;     Mean motion, radians/min.                                     *
*                                                                             *
*  sat.c2;      internal drag term                                            *
*                                                                             *
*  sat.bstar;   Drag term.                                                    *
*                                                                             *                                                         *                                                                            **                                                                             **                                                                             *                                                                             *
*******************************************************************************/

#include <cstdlib>
#include <complex>
#include "date.h"

using namespace std;

/* Modified Harris-Preister density tables */
/* densities in kg/km^3 */
static const double dens[60][3] = {
   100.,   4.974e+02,  4.974e+02,
   110.,   7.800e+01,  7.800e+01,
   120.,   2.490e+01,  2.400e+01,
   130.,   8.377e+00,  8.710e+00,
   140.,   3.899e+00,  4.059e+00,
   150.,   2.122e+00,  2.215e+00,
   160.,   1.263e+00,  1.344e+00,
   170.,   8.008e-01,  8.758e-01,
   180.,   5.283e-01,  6.010e-01,
   190.,   3.618e-01,  4.297e-01,
   200.,   2.557e-01,  3.162e-01,
   210.,   1.839e-01,  2.396e-01,
   220.,   1.341e-01,  1.853e-01,
   230.,   9.949e-02,  1.455e-01,
   240.,   7.488e-02,  1.157e-01,
   250.,   5.709e-02,  9.308e-02,
   260.,   4.403e-02,  7.555e-02,
   270.,   3.430e-02,  6.182e-02,
   280.,   2.697e-02,  5.095e-02,
   290.,   2.139e-02,  4.226e-02,
   300.,   1.708e-02,  3.526e-02,
   320.,   1.099e-02,  2.511e-02,
   340.,   7.214e-03,  1.819e-02,
   360.,   4.824e-03,  1.337e-02,
   380.,   3.274e-03,  9.955e-03,
   400.,   2.249e-03,  7.492e-03,
   420.,   1.558e-03,  5.684e-03,
   440.,   1.091e-03,  4.355e-03,
   460.,   7.701e-04,  3.362e-03,
   480.,   5.474e-04,  2.612e-03,
   500.,   3.916e-04,  2.042e-03,
   520.,   2.819e-04,  1.605e-03,
   540.,   2.042e-04,  1.267e-03,
   560.,   1.488e-04,  1.005e-03,
   580.,   1.092e-04,  7.997e-04,
   600.,   8.070e-05,  6.390e-04,
   620.,   6.012e-05,  5.123e-04,
   640.,   4.519e-05,  4.121e-04,
   660.,   3.430e-05,  3.325e-04,
   680.,   2.632e-05,  2.691e-04,
   700.,   2.043e-05,  2.185e-04,
   720.,   1.607e-05,  1.779e-04,
   740.,   1.281e-05,  1.452e-04,
   760.,   1.036e-05,  1.190e-04,
   780.,   8.496e-06,  9.776e-05,
   800.,   7.069e-06,  8.059e-05,
   850.,   4.800e-06,  5.500e-05,
   900.,   3.300e-06,  3.700e-05,
   950.,   2.450e-06,  2.400e-05,
   1000.,  1.900e-06,  1.700e-05,
   1100.,  1.180e-06,  8.700e-06,
   1200.,  7.500e-07,  4.800e-06,
   1300.,  5.300e-07,  3.200e-06,
   1400.,  4.100e-07,  2.000e-06,
   1500.,  2.900e-07,  1.350e-06,
   1600.,  2.000e-07,  9.500e-07,
   1700.,  1.600e-07,  7.700e-07,
   1800.,  1.200e-07,  6.300e-07,
   1900.,  9.600e-08,  5.200e-07,
   2000.,  7.300e-08,  4.400e-07};

/*##################### MATHEMATICAL CONSTANTS #####################*/

static const double
   nocon = 4.36332312998582e-3,   // 2pi / 1440
   de2ra = .0174532925199433,
   pi = 3.141592653589793238462643383279502884197,
   xkmper = 6378.1363,   // equatorial earth radius, km
   twopi = 2.*pi,
   two_thirds = 2./3.;

////////////// define Satellite class /////////////////////////////////////////

   char desig[9];              // international designation number
   int ssn;                    // spacetrack number

class Satellite
{
private:

/*##################### Degree of Model #####################*/

static const int
N;

/*####################### PHYSICAL CONSTANTS #######################*/

/* dimensions & gravity of earth */

static const double
xj3,
ck2,
ck4,
xke,    /* = (G*M)^(1/2)*(er/min)^(3/2) where G =
                               Newton's grav const, M = earth mass */
xkmper,   // equatorial earth radius, km

/* SGP4 density constants.  qoms2t = ((qo - so) / xkmper) ** 4,
s = 1 + so / xkmper, where qo = 120 and so = 78 */

qoms2t,
a3ovk2,
s;

/*####################### UNITS & CONVENTIONS #########################

Unless otherwise indicated, throughout this program
quantities are measured in the following units:

angles        radians
length        equatorial earth radii (1 unit = xkmper km)
velocity      equatorial earth radii / minute
South latitudes are negative.
East longitude is positive, west negative.


/*################### SATELLITE ORBITAL ELEMENTS ####################*/

  double theta,
         tsince;    // time since epoch in minutes

  double *rr2,   // auxillary position, velocity pointers
         *vv2;

  inline double delos(double s)
  {
    return s == 0. ? 1. : 0.0;
  }

  // analytic orbit propogators
  void sgp4(double tsince);
  void sdp4(double tsince);
  // vectors to mean SXP4 elements subroutine
  void rvel(double* rr1, double* vv1);

  // semi-analytic propagator
  double sst(double tsince, int NODE);

  void Vxx(void);
  void GdG(double eq[6]);
  void QdQ(void);
  void nKdK(void);
  void KdK(void);
  void zonal(double eq[6], double rate[6]);
  void srp(double ep, double eq[6], double rate[6]);
  void rates(double ep, double eq[6], double rate[6]);
  // Runge Kutta
  void rk4(double ep, double *eq, double dt);

  // equinocital orbit elements
  double
   I,   // retrograde factor
   F,   // eccentric longitude
   L,   // true longitude
   nn;  // mean motion

   // equinocital coefficients
   double A, B, CC, alf, bet, gam;
   // equinocital derivatives
   double  dUa, dUh, dUk, dUalf, dUbet, dUgam, dUlam;

   // initialization
   void init(void);

   // Sun Moon
   double Zsun[3], Zmoon[3];
   void sun(double ep);
   void moon(double ep);
   void sun_moon(double ep, double eq[6], double rate[6]);

   // drag
   double DG[6];
   void drag(double ep, double eq[6]);

public:

  // Satellite geocentric position and velocity from SXP4.
  double rr[3], vv[3];
  // SXP4 mean orbit elements and Julian date
  double
   jd,        // Julian date
   tle,       // epoch of elements in tle format
   thetag,    // Greenwich hour angle in radians
   xincl,     // inclination
   xnodeo,    // right ascension of ascending node
   eo,        // eccentricity
   omegao,    // argument of the perigee
   xmo,       // mean anomaly
   xno,       // mean motion, radians/min
   c2,        // internal drag term for bstar conversion and ndot
   Ksrp,      // solar radiation coefficient
   REF,
   Kdr,       // drag coefficient
   bstar;     // BSTAR drag term

  // equinocital orbit elements
  double eq[6];
  // equinocital deltas
  double DEQ[6];
  // Mean element deltas
  double DEL[6];

  double
   yy,       // year
   doy,      // day of year
   hr,       // hour
   mn,       // minute
   ss;       // second

  // drag
  // always need global Zsun vector
  double rho(double hgt);

  Satellite(double tm, double* rr1, double* vv1, double bsr);
  Satellite(double tm,  double ii, double om,  double ec,
            double ww,  double ma, double nn,  double bsr);
  void delta_t(double epoch);
  void delta_sst(double epoch, int NODE);
  void delta_el(double jd, double ii, double om, double ec,
                double ww, double ma, double nn, double bsr);
  void sxp4(double tsince);

  // vectors to mean SXP4 elements
  void rv2el(double* rr1, double* vv1);
  // vectors to equinocital elements
  void rv2eq(double* rr1, double* vv1);
  // equinocital elements to vectors
  void eq2rv(double* rr1, double* vv1);
  // mean elements to equinocital elements
  void el2eq(void);
  // equinocital elements to mean elements
  void eq2el(void);

  void convert(void);
  void print_rv(void);
  void print_el(void);

}; // end class Satellite
const int Satellite :: N ( 20);
const double Satellite :: xj3 ( -2.53881E-6);
const double Satellite :: ck2 ( 5.413079E-4);
const double Satellite :: ck4 ( 6.2098875E-7);
const double Satellite :: xke ( 0.0743669161331734132);    /* = (G*M)^(1/2)*(er/min)^(3/2) where G =
                               Newton's grav const, M = earth mass */
const double Satellite :: xkmper ( 6378.135);   // 6378.1363 equatorial earth radius, km

const double Satellite :: qoms2t ( 1.880279159015270643865e-9);
const double Satellite :: a3ovk2 ( -1.*xj3/ck2);
const double Satellite :: s ( 1.0122292801892716);


  // equinocital arrays
  double  G[23], H[23], dGh[23],dGk[23], dGalf[23], dGbet[23],
          Q[23][23], dQ[23][23], nK[23][23], ndK[23][23],
          K[23][23], dK[23][23], V[23][23];

////////////// math / vector functions /////////////////////////////////////////


  double acose(double x)
  {
     double rval;

     if( x >= 1.) rval = 0.;
     else if( x <= -1.) rval = pi;
     else rval = acos( x);
     return( rval);
  }
  //   v1 . v2
  double dot(double* v1, double* v2)
  {
    double sum = 0;
    for (int i = 0; i < 3; i++)
       sum += v1[i] * v2[i];
    return sum;
  }
  //   ||v||
  double norm(double* v)
  {
    return sqrt(dot(v, v));
  }
  //  a * v = av
  void smult(double a, double* v, double* av)
  {
    for (int i = 0; i < 3; i++)
      av[i] = a * v[i];
  }
  //   v1 + v2 = s
  void vadd(double* v1, double* v2, double* s)
  {
    for (int i = 0; i < 3; i++)
      s[i] = v1[i] + v2[i];
  }
  //   s = v1 + a * v2   (used for subtraction)
  inline void vmadd(double* v1, double* v2, double* s, double a)
  {
    for (int i = 0; i < 3; i++)
      s[i] = v1[i] + a * v2[i];
  }
  //   v1 x v2 = b
  void cross(double v1[3], double v2[3], double b[3])
  {
    b[0] = v1[1] * v2[2] - v1[2] * v2[1];
    b[1] = v1[2] * v2[0] - v1[0] * v2[2];
    b[2] = v1[0] * v2[1] - v1[1] * v2[0];
  }
  //   u = v / ||v||
  void unitv(double v[], double u[])
  {
    double no = norm(v);
    for (int i = 0; i < 3; i++)
      u[i] = v[i] / no;
  }
  //   mod(var, div) = remainder of var / div
  inline double mod(double var, double div)
  {
    return((var < 0) ? div + fmod(var, div) : fmod(var, div));
  }
  double fmod2p( const double x)
  {
     double twopi = 2 * pi;
     double rval = fmod( x, twopi);

     if(rval < 0.)
        rval += twopi;
     return(rval);
  }
  inline int fact(int x)
  {
     if (x<2)
     return 1;
     return fact(x-1)*x; // recursive call
  }
  inline char *s_in(char *prompt, char *buffer)
  {
    printf("%s", prompt);
    if (feof(stdin)) exit(0);
    gets(buffer);
    return buffer;
  }

////////////// member function SGP4 ////////////////////////////////////////////

  // four argument constructor => vectors passed as arguments
  // mean orbital elements calculated from vectors
  Satellite :: Satellite(double tm, double* rr1, double* vv1, double bsr) :
                            tle(tm),   rr2(rr1),    vv2(vv1), bstar(bsr)
  {
    Date t1(tle);       // create date object
    jd = t1.jd;         // convert to julian date
    thetag = t1.thetag; // Greenwich hour angle in radians
    rv2el(rr2, vv2);    // mean elements and state vectors rr, vv @ tle
    rr[0] = rr2[0];
    rr[1] = rr2[1];
    rr[2] = rr2[2];
    vv[0] = vv2[0];
    vv[1] = vv2[1];
    vv[2] = vv2[2];
    el2eq();            // equinocital elements at epoch
    convert();          // fit equinocital elements
  }

  // eight argument constructor => elements passed as arguments
  Satellite :: Satellite(double tm,  double ii, double om,  double ec,
                         double ww,  double ma, double nn,  double bsr) :
                         tle(tm),    xincl(ii), xnodeo(om), eo(ec),
                         omegao(ww), xmo(ma),   xno(nn),    bstar(bsr)
  {
    Date t1(tle);       // create date object
    jd = t1.jd;         // convert to julian date
    thetag = t1.thetag; // Greenwich hour angle in radians
    xincl *= de2ra;     // inclination
    xnodeo *= de2ra;    // right ascension of ascending node
    omegao *= de2ra;    // argument of the perigee
    xmo *= de2ra;       // mean anomaly
    xno *= nocon;       // mean motion
    sxp4(0.0);          // state vectors: rr, vv @ tle
    el2eq();            // equinocital elements at epoch
    convert();          // fit equinocital elements
  }

  // predict position state vectors at epoch + delta_t
  // accepts a julian date or tle date as epoch
  void Satellite :: delta_t(double epoch)
  {
    Date t2(epoch);
    yy = t2.yy;
    doy = t2.doy;
    hr = t2.hr;
    mn = t2.mn;
    ss = t2.ss;
    double time = (t2.jd - jd) * 1440;   // time is minutes since epoch
    sxp4(time);
  }

  // predict position state vectors at epoch + delta_t
  // accepts a julian date or tle date as epoch
  void Satellite :: delta_sst(double epoch, int NODE)
  {
    Date t1(epoch);
    double time = (t1.jd - jd) * 1440;   // time is minutes since epoch
    jd = sst(time, NODE);
    Date t2(jd);
    tle = t2.tle;
    thetag = t2.thetag;
  }

  // change element(s) of an existing Satellite object
  // elements passed in tle format, degrees & revolutions/day
  void Satellite :: delta_el(double jo, double ii, double om, double ec,
                             double ww, double ma, double nn, double bsr)
  {
    jd = jo;                // julian date
    xincl = ii * de2ra;     // inclination
    xnodeo = om * de2ra;    // right ascension of ascending node
    eo = ec;                // eccentricity
    omegao = ww * de2ra;    // argument of the perigee
    xmo = ma * de2ra;       // mean anomaly
    xno = nn * nocon;       // mean motion
    bstar = bsr;            // atmospheric drag
    sxp4(0.0);              // state vectors: rr, vv @ tle
    el2eq();                // equinocital elements at epoch
    eq[0] += DEQ[0];
    eq[1] += DEQ[1];
    eq[2] += DEQ[2];
    eq[3] += DEQ[3];
    eq[4] += DEQ[4];
    eq[5] += DEQ[5];
  }

///////////////////// Initialize /////////////////////

  void Satellite :: init(void)
  {
     double Area, Cd, Refl, Ms;
     Vxx();
                                            /*  AO-13 Physical data */
     Ms = 84;                               /*  Mass       kg */
     Area = 0.55E-6;                        /*  Area       km^2 */
     Cd = 1.0;                              /*  Drag coeff (typ. 1.0 - 2.0) */
     Kdr = 0.5*Cd*Area/Ms;                  /*  Drag factor */
     Refl = 1.5;                            /*  Refectivity constant */
     Ksrp = 4.51E-3*Refl*Area/Ms;           /*  km / s^2 */
     Ksrp = Ksrp*3600/xkmper;               /*  er / min^2 */

  }


//   gaussian elimination
void rref(int NN, double rv[6][7], double b[])    // scaled partial pivoting
{
  int i, j, s, ROW;
  double ss[NN], bin, mult;

////////////////// calculate permanent scale factors //////////////////////////

  for(i = 0; i < NN; i++)
  {
    ss[i] = abs(rv[i][0]);
    for(j = 1; j < NN; j++)
      if(ss[i] < abs(rv[i][j]))
      {
        ss[i] = abs(rv[i][j]);
      }
  }  // end for i

///////////////////// swap rows according to scale ////////////////////////////

  for(j = 0; j < NN-1; j++)
  {
    ROW = j;
    for(i = j + 1; i < NN; i++)
    {
      if(abs(rv[ROW][j] / ss[ROW]) < abs(rv[i][j] / ss[i]))
        ROW = i;
    }
    if(ROW != j)
    {
      for(s = j; s < NN+1; s++)    // swap rows
      {
        bin = rv[j][s];
        rv[j][s] = rv[ROW][s];
        rv[ROW][s] = bin;
      }
      bin = ss[j];                 // swap scales
      ss[j] = ss[ROW];
      ss[ROW] = bin;
    }  // end if
  }

///////////////////// forward elimination /////////////////////////////////////

  for(j = 0; j < NN; j++)
  {
    mult = rv[j][j];
    if(mult != 0.)
    {
      for(s = j; s < NN+1; s++)
      {
        rv[j][s] = rv[j][s] / mult;
      }
    }
    for(i = j + 1; i < NN; i++)
    {
      mult = rv[i][j];
      for(s = j + 1; s < NN+1; s++)
      {
        rv[i][s] = rv[i][s] - mult * rv[j][s];
      }
      rv[i][j] = 0;
    }  // end for i
  }  // end for j

///////////// test for singular matrix /////////////////////////////////////////

  bin = 1;
  for(i = 0; i < NN; i++)
  {
    bin *= rv[i][i];
  }
  if(bin == 0)
  {
    printf("Singular matrix");
    exit(0);
  }

////////////////// back sustitution ////////////////////////////////////////////

  b[NN-1] = rv[NN-1][NN] / rv[NN-1][NN-1];
  for(i = NN-1; i >= 0; i--)
  {
    bin = 0;
    for(s = i + 1; s < NN; s++)
      bin = bin + rv[i][s] * b[s];
    b[i] = (rv[i][NN] - bin) / rv[i][i];
  }
}  // end rref



  void Satellite :: convert(void)
  {
     printf("converting");
     int i, j, s, NN = 50, PASS = 0;          // counters
     double s_sum, sum, rms, rev, RREF;
     // small change amounts
     double delta, kdelta;
     double rv[6][7];               // normal matrix for rref
     double rdata[NN][4];
     double mdata[3*NN][7];         // define matrix for multiple regression
     double dat[3*NN][3];           // store satellite position
     double b[6];                   // output deltas, sat velocity
     double eq0, eq1, eq2, eq3, eq4, eq5;
     double bsum, brms, KDR;
     double jdo, tleo, step = 1.;

     // intitialize
     init();                          // SST initialization, Kdr
     KDR = 1;                         // direction change flags
     kdelta = Kdr*0.2;                // initial drag delta

     // store equinocital values
     eq0 = eq[0];
     eq1 = eq[1];
     eq2 = eq[2];
     eq3 = eq[3];
     eq4 = eq[4];
     eq5 = eq[5];

     // store mean element values
     DEL[0] = xincl;
     DEL[1] = xnodeo;
     DEL[2] = eo;
     DEL[3] = omegao;
     DEL[4] = xmo;
     DEL[5] = xno;

     // save equinocital state for reference
     DEQ[0] = eq[0];
     DEQ[1] = eq[1];
     DEQ[2] = eq[2];
     DEQ[3] = eq[3];
     DEQ[4] = eq[4];
     DEQ[5] = eq[5];

     // SXP4 epoch and predicted positions
     tleo = tle;
     jdo = jd;
     jdo += 10.;
     rev = .06 * nocon / xno;             // 16.66 points per orbit, 3 orbits

     for(i = 0; i < NN; i++)
     {
        jdo += rev;
        delta_t(jdo);                    // SXP4 pseudo obs  el -> rv
        rdata[i][0] = jdo;
        rdata[i][1] = rr[0];
        rdata[i][2] = rr[1];
        rdata[i][3] = rr[2];
     }
     jdo = jd;

     // find REF
     RREF  = 0.0;
     brms = 0.0;
     bsum = 0.0;
     for(i = 0; i < NN; i++)              // advance thru SXP4 points
     {
        // relocate satellite at new time
        delta_sst(rdata[i][0], 0);        // eq -> rv
        if(REF > RREF) RREF = REF;

        // find the rms
        bsum = (rr[0] - rdata[i][1]);
        brms += bsum*bsum;
        bsum = (rr[1] - rdata[i][2]);
        brms += bsum*bsum;
        bsum = (rr[2] - rdata[i][3]);
        brms += bsum*bsum;
     }
     sum = sqrt(brms/NN);                  // fit

     loop:

     // restore state
     eq[0] = eq0;
     eq[1] = eq1;
     eq[2] = eq2;
     eq[3] = eq3;
     eq[4] = eq4;
     eq[5] = eq5;
     jd = jdo;

////////////////////////////////////////////////////////////////////////////////

     // full equinocital differential correction

     // unperturbed sst propagation to SXP4 points
     for(i = 0; i < NN; i++)
     {
        delta_sst(rdata[i][0], 0);         // relocate satellite at new time
                                           // one pass thru SXP4 points
        // find the deltas and load into output matrix, dat
        mdata[3*i  ][6] = rdata[i][1] - rr[0];   // store delta_x
        mdata[3*i+1][6] = rdata[i][2] - rr[1];   // store delta_y
        mdata[3*i+2][6] = rdata[i][3] - rr[2];   // store delta_z

        // store these positions
        dat[3*i  ][0] = rr[0];        // store sat_x
        dat[3*i+1][1] = rr[1];        // store sat_y
        dat[3*i+2][2] = rr[2];        // store sat_z
     }

     delta = .001;                      // element change
     // perturbed sst propagation to SXP4 points
     for(j = 0; j < 6; j++)                // step thru equinocital perturbations
     {
        // restore state
        eq[0] = eq0;
        eq[1] = eq1;
        eq[2] = eq2;
        eq[3] = eq3;
        eq[4] = eq4;
        eq[5] = eq5;
        jd = jdo;

        eq[j] += delta;                    // delta element
        for(i = 0; i < NN; i++)            // advance thru NN points
        {
           delta_sst(rdata[i][0], 0);      // relocate satellite at new time

           // find the deltas and load into output matrix, mdata
           mdata[3*i  ][j] = (rr[0] - dat[3*i  ][0]) / delta;
           mdata[3*i+1][j] = (rr[1] - dat[3*i+1][1]) / delta;
           mdata[3*i+2][j] = (rr[2] - dat[3*i+2][2]) / delta;
        }
     }

     // multiple regression
     for(j = 0; j < 6; j++)
     {
        for(s = 0; s < 7; s++)
        {
           rv[j][s] = 0;
           for(i = 0; i < 3*NN; i++)
           {
              rv[j][s] = rv[j][s] + mdata[i][j] * mdata[i][s];
           }
        }
     }

     rref(6, rv, b);

     // restore state
     eq[0] = eq0;
     eq[1] = eq1;
     eq[2] = eq2;
     eq[3] = eq3;
     eq[4] = eq4;
     eq[5] = eq5;
     jd = jdo;
     // update trial state
     eq[0] += step*b[0];
     eq[1] += step*b[1];
     eq[2] += step*b[2];
     eq[3] += step*b[3];
     eq[4] += step*b[4];
     eq[5] += step*b[5];

     // test to see if trial state is better fit
     // find rms
     brms = 0.0;
     bsum = 0.0;
     for(i = 0; i < NN; i++)              // advance thru SXP4 points
     {
        delta_sst(rdata[i][0], 0);        // relocate satellite at new time

        // find the rms
        bsum = (rr[0] - rdata[i][1]);
        brms += bsum*bsum;
        bsum = (rr[1] - rdata[i][2]);
        brms += bsum*bsum;
        bsum = (rr[2] - rdata[i][3]);
        brms += bsum*bsum;

     }
     rms = sqrt(brms/NN);                 // fit

     if(sum > rms+1.E-11)                 // better fit
     {
        // update values
        eq0 += step*b[0];
        eq1 += step*b[1];
        eq2 += step*b[2];
        eq3 += step*b[3];
        eq4 += step*b[4];
        eq5 += step*b[5];

        sum = rms;

        printf(".");
        /*
        // display computed deltas for all 6 components
        printf("\ndQ   %12.9f\n", step*b[0]);
        printf("     %12.9f\n", step*b[1]);
        printf("     %12.9f\n", step*b[2]);
        printf("     %12.9f\n", step*b[3]);
        printf("     %12.9f\n", step*b[4]);
        printf("     %12.9f\n", step*b[5]);

        printf("rms  %12.9f\n", rms);
        */
        PASS++;
        if(PASS > 5)
        {
           PASS = 0;
           step *= 2.;
        }
        goto loop;
     }
     else if(step > 0.5)
     {
        step *= 0.5;
        goto loop;
     }
     s_sum = sum;
////////////////////////////////////////////////////////////////////////////////

     // Kdr

     if(RREF > 0.0)
     {

        k_loop:                             // drag loop

        Kdr += KDR*kdelta;
        if(Kdr < 1.E-10) Kdr = 1.E-10;
        if(Kdr > 1.E-7) Kdr = 1.E-7;
        if(fabs(kdelta) < 5.E-12) goto k_out;

        brms = 0.0;
        bsum = 0.0;
        for(i = 0; i < NN; i++)              // advance thru SXP4 points
        {
           delta_sst(rdata[i][0], 0);        // relocate satellite at new time

           // find the rms
           bsum = (rr[0] - rdata[i][1]);
           brms += bsum*bsum;
           bsum = (rr[1] - rdata[i][2]);
           brms += bsum*bsum;
           bsum = (rr[2] - rdata[i][3]);
           brms += bsum*bsum;

        }
        rms = sqrt(brms/NN);
        // restore state
        eq[0] = eq0;
        eq[1] = eq1;
        eq[2] = eq2;
        eq[3] = eq3;
        eq[4] = eq4;
        eq[5] = eq5;
        jd = jdo;

        if(rms > sum)
        {
           KDR *= -1;
           kdelta = kdelta*.5;              // zoom on Kdr
           goto k_loop;
        }
        if(rms < sum)                       // better drag coefficient
        {
           sum = rms;
           printf(",");
           /*
           printf("\nKdr   %g", Kdr);
           printf("\nrms  %12.9f\n", rms);
           */
        }
     }
     k_out:
     kdelta = Kdr*0.2;
     if(s_sum - sum > 1.e-9 ) goto loop;

////////////////////////////////////////////////////////////////////////////////

     // restore state
     eq[0] = eq0;
     eq[1] = eq1;
     eq[2] = eq2;
     eq[3] = eq3;
     eq[4] = eq4;
     eq[5] = eq5;
     jd = jdo;
     tle = tleo;
     if(REF == 0.0) bstar = 0.0;

     // compute el2eq deltas  (el2eq:eq[0] + DEQ[0] = eq[0])
     DEQ[0] = eq[0] - DEQ[0];
     DEQ[1] = eq[1] - DEQ[1];
     DEQ[2] = eq[2] - DEQ[2];
     DEQ[3] = eq[3] - DEQ[3];
     DEQ[4] = eq[4] - DEQ[4];
     DEQ[5] = eq[5] - DEQ[5];

     // compute eq2el deltas  (eq2el:xincl + DEL[0] = xincl)
     eq2el();
     DEL[0] -= xincl;
     DEL[1] -= xnodeo;
     DEL[2] -= eo;
     DEL[3] -= omegao;
     DEL[4] -= xmo;
     DEL[5] -= xno;

     printf("\n");
////////  end CONVERT  /////////////////////////////////////////////////////////

  }


  // vectors to SGP4 mean elements
void Satellite :: rvel(double* rr2, double* vv2)
{
    //double twopi = 2 * 3.14159265358979323846;
    int i;

    /* classical osculating orbit elements calculated from vectors rr2, vv2  */
    double xinck, xnodek, ek, mk, wk, xn, rk, uk, aodp, pl, rdotk, rfdotk, temp;

    double h[3], n[3], vec[3], vk[3];
    smult(1. / xke, vv2, vk);
    cross(rr2, vk, h);
    pl = dot(h, h);
    double vz[] = {0, 0, 1};
    double vy[3], t;
    cross(vz, h, n);
    if(n[0] == 0. && n[1] == 0.) n[0] = 1.;
    unitv(n, n);
    rk = norm(rr2);
    rdotk = dot(rr2, vv2) / rk;
    rfdotk = norm(h) * xke / rk;
    temp = dot(rr2, n) / rk;
    uk = acose(temp);
    if(rr2[2] < 0.) uk = twopi - uk;
    cross(vk, h, vz);
    smult(-1. / rk, rr2, vy);
    vadd(vz, vy, vec);
    ek = norm(vec);
    if( ek >= 1.) return;         // open orbit
    xnodek = atan2( n[1], n[0]);
    if( xnodek < 0.)
       xnodek += twopi;
    temp = sqrt( h[0] * h[0] + h[1] * h[1]);
    xinck =  atan2( temp, h[2]);
    temp = dot(vec, n) / ek;
    wk = acose( temp);
    if(vec[2] < 0.)
      wk = fmod2p(twopi - wk);
    aodp = pl / (1. - ek*ek);
    xn = xke * pow(aodp, -1.5);

    double cosio, sinio, sin2u, cos2u, temp1, temp2,
     rdot, rfdot, theta2, betal, x3thm1, x1mth2, x7thm1,
     esine, ecose, elsq, cosepw, sinepw, axn, ayn,
     cosu, sinu, capu, /*a3ovk2,*/ xlcof, aycof, aynl, xll,
     xl, a0, a1, a2, d0, d1, beta, beta2, r, u;

/*
In the first loop the osculating elements rk, uk, xnodek, xinck, rdotk,
and rfdotk are used as anchors to find the corresponding final SGP4
mean elements r, u, xnodeo, xincl, rdot, and rfdot.  Several other final
mean values based on these are also found: betal, cosio, sinio, theta2,
cos2u, sin2u, x3thm1, x7thm1, x1mth2.  In addition, the osculating values
initially held by aodp, pl, and xn are replaced by intermediate
(not osculating and not mean) values used by SGP4.  The loop converges
on the value of pl in about four iterations.
*/

    /*  seed value for first loop */
    xincl = xinck;
    u = uk;

    for (i = 0; i < 99; ++i)
    {
      a2 = pl;
      betal = sqrt(pl / aodp);
      temp1 = ck2  / pl;
      temp2 = temp1 / pl;
      cosio = cos(xincl);
      sinio = sin(xincl);
      sin2u = sin(2.*u);
      cos2u = cos(2.*u);
      theta2 = cosio * cosio;
      x3thm1 = 3. * theta2 - 1.;
      x1mth2 = 1. - theta2;
      x7thm1 = 7. * theta2 - 1.;
      r = (rk - .5 * temp1 * x1mth2 * cos2u)
         / (1. - 1.5 * temp2 * betal * x3thm1);
      u = uk + .25 * temp2 * x7thm1 * sin2u;
      xnodeo = xnodek - 1.5 * temp2 * cosio * sin2u;
      xincl = xinck - 1.5 * temp2 * cosio * sinio * cos2u;
      rdot = rdotk + xn * temp1 * x1mth2 * sin2u;
      rfdot = rfdotk - xn * temp1 * (x1mth2 * cos2u + 1.5 * x3thm1);
      temp = r * rfdot / xke;
      pl = temp * temp;

      // vis-viva equation
      temp = 2. / r - (rdot*rdot + rfdot*rfdot) / (xke*xke);
      aodp = 1. / temp;

      xn = xke * pow(aodp, -1.5);
      if(fabs(a2 - pl) < 1.e-13) break;
    }

/*
The next values are calculated from constants and a combination of mean
and intermediate quantities from the first loop.  These values all remain
fixed and are used in the second loop.
*/


    // preliminary values for the second loop
    ecose = 1. - r / aodp;
    esine = r * rdot / (xke * sqrt(aodp));   /* needed for Kepler's eqn.  */
    elsq = 1. - pl / aodp;               /* intermediate eccentricity squared */
    //a3ovk2 = -xj3 / ck2;
    xlcof = .125 * a3ovk2 * sinio * (3. + 5. * cosio)
          / (1. + cosio);
    aycof = .25 * a3ovk2 * sinio;
    temp1 = esine / (1. + sqrt(1. - elsq));
    cosu = cos(u);
    sinu = sin(u);

/*
The second loop normally converges in about six iterations to the final
mean value for the eccentricity, eo.  The mean perigee, omegao, is also
determined.  Cosepw and sinepw are found to twelve decimal places and
are used to calculate an intermediate value for the eccentric anomaly,
temp2.  Temp2 is then used in Kepler's equation to find an intermediate
value for the true longitude, capu.
*/
    /*  seed values for loop  */
    eo = sqrt(elsq);
    omegao = wk;
    axn = eo * cos(omegao);

    for (i = 0; i < 99; ++i)
    {
       a2 = eo;
       beta = 1. - eo*eo;
       temp = 1. / (aodp * beta);
       aynl = temp * aycof;
       ayn = eo * sin(omegao) + aynl;
       cosepw = r * cosu / aodp + axn - ayn * temp1;
       sinepw = r * sinu / aodp + ayn + axn * temp1;
       axn = cosepw * ecose + sinepw * esine;
       ayn = sinepw * ecose - cosepw * esine;
       omegao = fmod2p(atan2(ayn - aynl, axn));
       eo = .9*eo + .1*(axn / cos(omegao));
       if(eo > .999) eo = .999;
       if(fabs(a2 - eo) < 1.e-9) break;
    }

    temp2 = atan2(sinepw, cosepw);
    capu = temp2 - esine;             /* Kepler's equation */
    xll = temp * xlcof * axn;

    /* xll adjusts the intermediate true longitude,  */
    /* capu, to the mean true longitude, xl          */
    xl = capu - xll;

    xmo = fmod2p(xl - omegao);        /* mean anomaly */

/*
The third loop usually converges after three iterations to the
mean semi-major axis, a1, which is then used to find the mean motion, xno.
*/

    a0 = aodp;
    a1 = a0;
    beta2 = sqrt(beta);
    temp = 1.5 * ck2 * x3thm1 / (beta * beta2);
    for (i = 0; i < 99; ++i)
    {
       a2 = a1;
       d0 = temp / (a0*a0);
       a0 = aodp * (1. - d0);
       d1 = temp / (a1*a1);
       a1 = a0 / (1. - d1 / 3. - d1*d1 - 134. * d1*d1*d1 / 81.);
       if(fabs(a2 - a1) < 1.e-13) break;
    }
    xno = xke * pow(a1 , -1.5);

} /* end rvel  */

// vectors to SDP4 mean elements
void Satellite :: rv2el(double* rr, double* vv)
{
   double ik, ok, ek, wk, mk, nk;
   double iz, oz, ez, wz, mz, nz;
   double rr1[3], vv1[3];

   rvel(rr, vv);              // SGP4 elements

   if(xno / nocon < 6.4)
   {
      rr1[0] = rr[0];
      rr1[1] = rr[1];
      rr1[2] = rr[2];
      vv1[0] = vv[0];
      vv1[1] = vv[1];
      vv1[2] = vv[2];

      ik = xincl,
      ok = xnodeo,
      ek = eo,
      wk = omegao,
      mk = xmo,
      nk = xno;

      sdp4(0.0);              // SDP4 propagation to rr, vv
      rvel(rr, vv);           // SDP4 elements first correction

      xincl  = 2*ik - xincl;
      xnodeo = 2*ok - xnodeo;
      eo     = 2*ek - eo;
      omegao = 2*wk - omegao;
      xmo    = 2*mk - xmo;
      xno    = 2*nk - xno;

      // Loop ---->
      for(int i = 0; i < 2; i++)
      {
        // These elements are pretty close.  Save these
        // as the "pretty close" reference elements
        iz = xincl,
        oz = xnodeo,
        ez = eo,
        wz = omegao,
        mz = xmo,
        nz = xno;

        sdp4(0.0);              // SDP4 propagation to rr, vv
        rvel(rr, vv);           // SDP4 elements second correction

        // Correct the "pretty close" reference elements
        xincl  = iz - (xincl - ik);
        xnodeo = oz - (xnodeo - ok);
        eo     = ez - (eo - ek);
        omegao = wz - (omegao - wk);
        xmo    = mz - (xmo - mk);
        xno    = nz - (xno - nk);

      } //  <---- end Loop

      if(xincl < 0.) xincl  *= -1.;
      xincl = fmod2p(xincl);
      xnodeo = fmod2p(xnodeo);
      omegao = fmod2p(omegao);
      xmo    = fmod2p(xmo);

      rr[0] = rr1[0];
      rr[1] = rr1[1];
      rr[2] = rr1[2];
      vv[0] = vv1[0];
      vv[1] = vv1[1];
      vv[2] = vv1[2];
   }
} /* end rv2el  */

void Satellite :: sgp4(double tsince)
{
   // adapted from projectpluto.com

   // global tle variables
   // jd, xndt2o, xndd6o, bstar, c2;
   // xincl, xnodeo, eo, omegao, xmo, xno;

   double c1, c4, xnodcf, t2cof, aodp, cosio, sinio,
          omgdot, xmdot, xnodot, xnodp, c5, d2, d3,
          d4, delmo, eta, omgcof, sinmo, t3cof, t4cof,
          t5cof, xmcof, eosq, betao, betao2;
   int simple_flag;
   double MINIMAL_E = 1.e-4;
   double ECC_EPS = 1.e-6;     /* Too low for computing further drops. */

   double
         coef, coef1, eeta, etasq, perige, pinv, pinvsq,
         psisq, qoms24, s4, temp1, temp2, temp3,
         theta4, tsi, tsi_squared, x1mth2, x1m5th, xhdot1;

   const double a1 = pow(xke/ xno,two_thirds);
   double del1, ao, delo, x3thm1, tval, theta2;

   /* Recover original mean motion (xnodp) and   */
   /* semimajor axis (aodp) from input elements. */
   cosio = cos( xincl);
   theta2 =  cosio* cosio;
   x3thm1 = 3. *  theta2 - 1.;
   eosq =  eo* eo;
   betao2 = 1- eosq;
   betao = sqrt( betao2);
   tval = 1.5 * ck2 * x3thm1 / ( betao *  betao2);
   del1 = tval / (a1 * a1);
   ao = a1*(1-del1*(0.5*two_thirds+del1*(1.+134./81.*del1)));
   delo = tval / (ao * ao);
   xnodp =  xno/(1+delo);
   aodp = ao/(1-delo);

   x3thm1 = 3* theta2-1;
   /* For perigee below 156 km, the values */
   /* of s and qoms2t are altered.         */
   s4 = s;
   qoms24 = qoms2t;
   perige = ( aodp*(1- eo)-1.)*xkmper;
   if(perige < 156)
      {
      double temp_val, temp_val_squared;

      if(perige <= 98)
          s4 = 20;
      else
          s4 = perige-78.;
      temp_val = (120. -  s4) / xkmper;
      temp_val_squared = temp_val * temp_val;
      qoms24 = temp_val_squared * temp_val_squared;
       s4 =  s4/xkmper+1.;
      }  /* End of if(perige <= 156) */

   pinv = 1. / ( aodp* betao2);
   pinvsq = pinv * pinv;
   tsi = 1. / ( aodp -  s4);
   eta =  aodp* eo* tsi;
   etasq =  eta* eta;
   eeta =  eo* eta;
   psisq = fabs(1-etasq);
   tsi_squared =  tsi *  tsi;
   coef = qoms24 * tsi_squared * tsi_squared;
   coef1 =  coef / pow(psisq,3.5);
   c2 =  coef1 *  xnodp * ( aodp*(1+1.5*etasq+eeta*
   (4+etasq))+0.75*ck2* tsi/psisq*x3thm1*(8+3*etasq*(8+etasq)));
   c1 =  bstar*c2;
   sinio = sin( xincl);
   x1mth2 = 1- theta2;
   c4 = 2* xnodp* coef1* aodp* betao2*
        ( eta*(2+0.5*etasq)+ eo*(0.5+2*etasq)-2*ck2* tsi/
        ( aodp*psisq)*(-3*x3thm1*(1-2*eeta+etasq*
        (1.5-0.5*eeta))+0.75*x1mth2*(2*etasq-eeta*(1+etasq))*
        cos(2* omegao)));
   theta4 =  theta2* theta2;
   temp1 = 3*ck2*pinvsq* xnodp;
   temp2 = temp1*ck2*pinvsq;
   temp3 = 1.25*ck4*pinvsq*pinvsq* xnodp;
   xmdot =  xnodp+0.5*temp1* betao*
               x3thm1+0.0625*temp2* betao*
                    (13-78* theta2+137*theta4);
   x1m5th = 1-5* theta2;
   omgdot = -0.5*temp1*x1m5th+0.0625*temp2*
                     (7-114* theta2+395*theta4)+
                temp3*(3-36* theta2+49*theta4);
   xhdot1 = -temp1* cosio;
   xnodot = xhdot1+(0.5*temp2*(4-19* theta2)+
            2*temp3*(3-7* theta2))* cosio;
   xnodcf = 3.5* betao2*xhdot1*c1;
   t2cof = 1.5*c1;

   eeta = eo*eta;

   /* For perigee less than 220 kilometers, the "simple" flag is set */
   /* and the equations are truncated to linear variation in sqrt a  */
   /* and quadratic variation in mean anomaly.  Also, the c3 term,   */
   /* the delta omega term, and the delta m term are dropped.        */
   simple_flag = (aodp*(1- eo) < (220./xkmper+1.));
   if( !simple_flag)
      {
      const double c1sq = c1*c1;
      double temp;

      simple_flag = 0;
      delmo = 1. + eta * cos( xmo);
      delmo *= delmo * delmo;
      d2 = 4*aodp*tsi*c1sq;
      temp = d2* tsi*c1/3;
      d3 = (17*aodp+ s4)*temp;
      d4 = 0.5*temp* aodp* tsi*(221*aodp+31* s4)*c1;
      t3cof = d2+2*c1sq;
      t4cof = 0.25*(3*d3+c1*(12*d2+10*c1sq));
      t5cof = 0.2*(3*d4+12*c1*d3+6*d2*d2+15*c1sq*(2*d2+c1sq));
      sinmo = sin( xmo);
      if(  eo < MINIMAL_E)
         omgcof = xmcof = 0.;
      else
         {
         const double c3 =
              coef * tsi * a3ovk2 * xnodp * sinio /  eo;

         xmcof = -two_thirds * coef *  bstar / eeta;
         omgcof =  bstar*c3*cos( omegao);
         }
      } /* End of if (isFlagClear(SIMPLE_FLAG)) */
   etasq =  eta *  eta;
   c5 = 2* coef1* aodp * betao2*(1+2.75*(etasq+eeta)+eeta*etasq);

  double
        a, e, omega, omgadf,
        temp, tempa, tempe, templ, tsq,
        xl, xmdf, xmp, xnoddf, xnode;

  /* Update for secular gravity and atmospheric drag. */
  xmdf =  xmo+ xmdot*tsince;
  omgadf =  omegao+ omgdot*tsince;
  xnoddf =  xnodeo+ xnodot*tsince;
  omega = omgadf;
  xmp = xmdf;
  tsq = tsince*tsince;
  xnode = xnoddf+xnodcf*tsq;
  tempa = 1-c1*tsince;
  tempe =  bstar*c4*tsince;
  templ = t2cof*tsq;
  if( !simple_flag)
    {
      const double delomg = omgcof*tsince;
      double delm = 1. +  eta * cos(xmdf);
      double tcube, tfour;

      delm = xmcof * (delm * delm * delm - delmo);
      temp = delomg+delm;
      xmp = xmdf+temp;
      omega = omgadf-temp;
      tcube = tsq*tsince;
      tfour = tsince*tcube;
      tempa = tempa-d2*tsq-d3*tcube-d4*tfour;
      tempe = tempe+ bstar*c5*(sin(xmp)-sinmo);
      templ = templ+t3cof*tcube+tfour*(t4cof+tsince*t5cof);
    }; /* End of if (isFlagClear(SIMPLE_FLAG)) */

  a =  aodp*tempa*tempa;
  e =  eo-tempe;

  /* A highly arbitrary lower limit on e,  of 1e-6: */
  if( e < ECC_EPS)
     e = ECC_EPS;
  xl = xmp+omega+xnode+xnodp*templ;

  /* Long period periodics */
  double axn = e*cos(omega);
  temp = 1/(a*(1.-e*e));
  double xlcof = .125 * a3ovk2 * sinio * (3+5*cosio)/ (1. + cosio);
  double aycof = 0.25 * a3ovk2 * sinio;
  double xll = temp*xlcof*axn;
  double aynl = temp*aycof;
  double xlt = xl+xll;
  double ayn = e*sin(omega)+aynl;
  double elsq = axn*axn+ayn*ayn;
  double capu = fmod( xlt - xnode, twopi);
  double chicken_factor_on_eccentricity = 1.e-6;
  double epw = capu;
  double ecosE, esinE, pl, r;
  double betal;
  double u, sinu, cosu, sin2u, cos2u;
  double rk, uk, xnodek, xinck;
  double sinuk, cosuk, sinik, cosik, sinnok, cosnok, xmx, xmy;
  double sinEPW, cosEPW;
  double ux, uy, uz;
  int i, rval = 0;

  /* Dundee changes:  items dependent on cosio get recomputed: */
  double cosio_squared = cosio * cosio;
  x3thm1 = 3.0 * cosio_squared - 1.0;
  x1mth2 = 1.0 - cosio_squared;
  double x7thm1 = 7.0 * cosio_squared - 1.0;

         /* Added 29 Mar 2003,  modified 26 Sep 2006:  extremely    */
         /* decayed satellites can end up "orbiting" within the     */
         /* earth.  Eventually,  the semimajor axis becomes zero,   */
         /* then negative.  In that case,  or if the orbit is near  */
         /* to parabolic,  we zero the posn/vel and quit.  If the   */
         /* object has a perigee or apogee indicating a crash,  we  */
         /* just flag it.  Revised 28 Oct 2006.                     */

  if( elsq > 1. - chicken_factor_on_eccentricity)
     rval = -1;
  if( rval)
  {
     for( i = 0; i < 3; i++)
     {
        rr[i] = 0.;
        vv[i] = 0.;
     }
     exit( 0);
  }
  if( a * (1. - e) < 1. && a * (1. + e) < 1.)   /* entirely within earth */
     rval = -3;                      /* remember, e can be negative */
  if( a * (1. - e) < 1. || a * (1. + e) < 1.)   /* perigee within earth */
     rval = -4;
  /* Solve Kepler's' Equation */
  for( i = 0; i < 10; i++)
  {
     double newton_raphson_epsilon = 1e-12;
     double f, fdot, delta_epw;
     int do_second_order_newton_raphson = 1;

     sinEPW = sin(epw);
     cosEPW = cos(epw);
     ecosE = axn * cosEPW + ayn * sinEPW;
     esinE = axn * sinEPW - ayn * cosEPW;
     f = capu - epw + esinE;
     if ( fabs(f) < newton_raphson_epsilon) break;
     fdot = 1. - ecosE;
     delta_epw = f / fdot;
     if( !i)
     {
        double max_newton_raphson = 1.25 * fabs(e);

        do_second_order_newton_raphson = 0;
        if( delta_epw > max_newton_raphson)
            delta_epw = max_newton_raphson;
        else if( delta_epw < -max_newton_raphson)
                 delta_epw = -max_newton_raphson;
        else
            do_second_order_newton_raphson = 1;
     }
     if( do_second_order_newton_raphson)
         delta_epw = f / (fdot + 0.5*esinE*delta_epw);
                            /* f/(fdot - 0.5*fdotdot * f / fdot) */
     epw += delta_epw;
  }

  /* Short period preliminary quantities */
  temp = 1-elsq;
  pl = a*temp;
  r = a*(1-ecosE);
  temp2 = a / r;
  betal = sqrt(temp);
  temp = esinE/(1+betal);
  cosu = temp2 * (cosEPW - axn + ayn * temp);
  sinu = temp2 * (sinEPW - ayn - axn * temp);
  u = atan2( sinu, cosu);
  sin2u = 2*sinu*cosu;
  cos2u = 2*cosu*cosu-1;
  temp1 = ck2 / pl;
  temp2 = temp1 / pl;

  /* Update for short periodics */
  rk = r*(1-1.5*temp2*betal*x3thm1)+0.5*temp1*x1mth2*cos2u;
  uk = u-0.25*temp2*x7thm1*sin2u;
  xnodek = xnode+1.5*temp2*cosio*sin2u;
  xinck = xincl+1.5*temp2*cosio*sinio*cos2u;

  /* Orientation vectors */
  sinuk = sin(uk);
  cosuk = cos(uk);
  sinik = sin(xinck);
  cosik = cos(xinck);
  sinnok = sin(xnodek);
  cosnok = cos(xnodek);
  xmx = -sinnok*cosik;
  xmy = cosnok*cosik;
  ux = xmx*sinuk+cosnok*cosuk;
  uy = xmy*sinuk+sinnok*cosuk;
  uz = sinik*sinuk;

  /* Global Position, er, and velocity, er/min */
  rr[0] = rk*ux;
  rr[1] = rk*uy;
  rr[2] = rk*uz;

  double rdot = xke*sqrt(a)*esinE/r;
  double rfdot = xke*sqrt(pl)/r;
  double xn = xke/(a * sqrt(a));
  double rdotk = rdot-xn*temp1*x1mth2*sin2u;
  double rfdotk = rfdot+xn*temp1*(x1mth2*cos2u+1.5*x3thm1);
  double vx = xmx*cosuk-cosnok*sinuk;
  double vy = xmy*cosuk-sinnok*sinuk;
  double vz = sinik*cosuk;

  vv[0] = rdotk*ux+rfdotk*vx;
  vv[1] = rdotk*uy+rfdotk*vy;
  vv[2] = rdotk*uz+rfdotk*vz;
}


void Satellite :: sdp4(double tsince)
{
   // global tle variables
   // jd, xndt2o, xndd6o, bstar;
   // xincl, xnodeo, eo, omegao, xmo, xno;

   double ao, xnodp, delo, a1, del1, r1, temp, x3thm1, tval;

   // define in SDP4 and COMMON
   double c1, c4, xnodcf, t2cof ;

   // init_t
   double coef, coef1, tsi, s4, eta;

   // deep_arg
   double
   /* Common between SGP4 and SDP4: */
   aodp, cosio, sinio, omgdot, xmdot, xnodot,
   /* Used by dpinit */
   eosq, betao, theta2, sing, cosg, betao2,

   /* Used by dpsec and dpper */
   xll, omgadf, xnode, em, xinc, xn, t,

   /* 'd####' secular coeffs for 12-hour, e>.5 orbits: */
   d2201, d2211, d3210, d3222, d4410, d4422, d5220, d5232, d5421, d5433,
   atime, del2, del3, e3, ee2, omegaq, pe, pgh, ph, pinc, pl,
   se2, se3, sgh2, sgh3, sgh4, sh2, sh3, si2, si3, sl2, sl3, shdq,
   sl4, sse, ssg, ssh, ssi, ssl, thgr, xfact, xgh2, xgh3, xgh4, xh2,
   xh3, xi2, xi3, xl2, xl3, xl4, xlamo, xli, xni, xnq, zcosgl,
   zcoshl, zcosil, zmol, zmos, zsingl, zsinhl, zsinil;

   int i, resonance_flag, synchronous_flag;
   double
         eeta, etasq, perige, pinv, pinvsq,
         psisq, qoms24, temp1, temp2, temp3,
         theta4, tsi_squared, x1mth2, x1m5th, xhdot1;

   /* Recover original mean motion (xnodp) and   */
   /* semimajor axis (aodp) from input elements. */

   a1 = pow(xke / xno, 2/3.);
   cosio = cos(xincl);
   theta2 = cosio * cosio;
   x3thm1 = 3. * theta2 - 1.;
   eosq = eo * eo;
   betao2 = 1. - eosq;
   betao = sqrt(betao2);
   del1 = 1.5 * ck2 * x3thm1 / (a1 * a1 * betao * betao2);
   ao = a1 * (1. - del1 * (1/3. + del1 * (1. + 134.
    / 81. * del1)));
   delo = 1.5 * ck2 * x3thm1 / (ao * ao * betao * betao2);
   xnodp = xno / (1. + delo);
   aodp = ao / (1. - delo);

   /* For perigee below 156 km, the values */
   /* of s and qoms2t are altered.         */
   s4 = s;
   qoms24 = qoms2t;
   perige = (aodp*(1-eo)-1.)*xkmper;
   if(perige < 156)
   {
      double temp_val, temp_val_squared;

      if(perige <= 98)
         s4 = 20;
      else
         s4 = perige-78.;
      temp_val = (120. - s4) / xkmper;
      temp_val_squared = temp_val * temp_val;
      qoms24 = temp_val_squared * temp_val_squared;
      s4 = s4/xkmper+1.;
   }  /* End of if(perige <= 156) */

   pinv = 1. / (aodp*betao2);
   pinvsq = pinv * pinv;
   tsi = 1. / (aodp - s4);
   eta = aodp*eo*tsi;
   etasq = eta*eta;
   eeta = eo*eta;
   psisq = fabs(1-etasq);
   tsi_squared = tsi * tsi;
   coef = qoms24 * tsi_squared * tsi_squared;
   coef1 = coef / pow(psisq,3.5);
   c2 = coef1 * xnodp * (aodp*(1+1.5*etasq+eeta*
   (4+etasq))+0.75*ck2*tsi/psisq*x3thm1*(8+3*etasq*(8+etasq)));
   c1 = bstar*c2;
   sinio = sin(xincl);
   x1mth2 = 1-theta2;
   c4 = 2*xnodp*coef1*aodp*betao2*
        (eta*(2+0.5*etasq)+eo*(0.5+2*etasq)-2*ck2*tsi/
        (aodp*psisq)*(-3*x3thm1*(1-2*eeta+etasq*
        (1.5-0.5*eeta))+0.75*x1mth2*(2*etasq-eeta*(1+etasq))*
        cos(2*omegao)));
   theta4 = theta2*theta2;
   temp1 = 3*ck2*pinvsq*xnodp;
   temp2 = temp1*ck2*pinvsq;
   temp3 = 1.25*ck4*pinvsq*pinvsq*xnodp;
   xmdot = xnodp+0.5*temp1*betao*
               x3thm1+0.0625*temp2*betao*
                    (13-78*theta2+137*theta4);
   x1m5th = 1-5*theta2;
   omgdot = -0.5*temp1*x1m5th+0.0625*temp2*
                     (7-114*theta2+395*theta4)+
                temp3*(3-36*theta2+49*theta4);
   xhdot1 = -temp1*cosio;
   xnodot = xhdot1+(0.5*temp2*(4-19*theta2)+
           2*temp3*(3-7*theta2))*cosio;
   xnodcf = 3.5*betao2*xhdot1*c1;
   t2cof = 1.5*c1;
   /* end sxpx_common_init */

   sing = sin(omegao);
   cosg = cos(omegao);

   /* initialize */
   double zns  = 1.19459E-5;
   double zes  = 0.01675;
   double znl  = 1.5835218E-4;
   double zel  = 0.05490;
   double thdt = 4.3752691E-3;
   double dpsec_integration_step = 720.;
   int dpsec_integration_order = 2;

   /* dpinit */
   double sinq = sin(xnodeo);
   double cosq = cos(xnodeo);
   double aqnv = 1/aodp;
   double c1ss   =  2.9864797E-6;
   double days_since_1900 = jd - 2415020.;
                                     /* 1900 Jan 0.5 = JD 2415020. */
   double zcosi =  0.91744867;
   double zsini =  0.39785416;
   double zsing = -0.98088458;
   double zcosg =  0.1945905;
   double bfact, cc = c1ss, se;
   double ze = zes, zn = zns;
   double sgh, sh, si;
   double zsinh = sinq, zcosh = cosq;
   double sl;
   int iteration;

   /* thetaG Reference:  The 1992 Astronomical Almanac, page B6. */
   double omega_E = 1.00273790934;
                   /* Earth rotations per sidereal day (non-constant) */
   double UT = fmod( jd + .5, 1.);
   double seconds_per_day = 86400.;
   double jd_2000 = 2451545.0;   /* 1.5 Jan 2000 = JD 2451545. */
   double t_cen, GMST;

   t_cen = (jd - UT - jd_2000) / 36525.;
   GMST = 24110.54841 + t_cen * (8640184.812866 + t_cen *
                           (0.093104 - t_cen * 6.2E-6));
   GMST = fmod( GMST + seconds_per_day * omega_E * UT, seconds_per_day);
   if( GMST < 0.)
      GMST += seconds_per_day;
   thgr = twopi * GMST / seconds_per_day;
   /* end thetaG */

   xnq = xnodp;
   omegaq = omegao;

   double xnodce = 4.5236020 - 9.2422029E-4 * days_since_1900;
   double stem = sin(xnodce);
   double ctem = cos(xnodce);
   double c_minus_gam = 0.228027132 * days_since_1900 - 1.1151842;
   double gam = 5.8351514 + 0.0019443680 * days_since_1900;
   double zx, zy;

   zcosil = 0.91375164 - 0.03568096 * ctem;
   zsinil = sqrt(1. - zcosil * zcosil);
   zsinhl = 0.089683511 * stem / zsinil;
   zcoshl = sqrt(1. - zsinhl*zsinhl);
   zmol = fmod2p( c_minus_gam);
   zx = 0.39785416 * stem / zsinil;
   zy = zcoshl * ctem + 0.91744867 * zsinhl * stem;
   zx = atan2( zx, zy) + gam - xnodce;
   zcosgl = cos( zx);
   zsingl = sin( zx);
   zmos = fmod2p( 6.2565837 + 0.017201977 * days_since_1900);

   /* Do solar terms */

   /* There was previously some convoluted logic here,  but it boils    */
   /* down to this:  we compute the solar terms,  then the lunar terms. */
   /* On a second pass,  we recompute the solar terms,  taking advantage */
   /* of the improved data that resulted from computing lunar terms.     */
   for( iteration = 0; iteration < 2; iteration++)
   {
      double c1l = 4.7968065E-7;
      double a1 = zcosg * zcosh + zsing * zcosi * zsinh;
      double a3 = -zsing * zcosh + zcosg * zcosi * zsinh;
      double a7 = -zcosg * zsinh + zsing * zcosi * zcosh;
      double a8 = zsing * zsini;
      double a9 = zsing * zsinh + zcosg * zcosi * zcosh;
      double a10 = zcosg * zsini;
      double a2 = cosio * a7 + sinio * a8;
      double a4 = cosio * a9 + sinio * a10;
      double a5 = -sinio * a7 + cosio * a8;
      double a6 = -sinio * a9 + cosio * a10;
      double x1 = a1 * cosg + a2 * sing;
      double x2 = a3 * cosg + a4 * sing;
      double x3 = -a1 * sing + a2 * cosg;
      double x4 = -a3 * sing + a4 * cosg;
      double x5 = a5 * sing;
      double x6 = a6 * sing;
      double x7 = a5 * cosg;
      double x8 = a6 * cosg;
      double z31 = 12 * x1 * x1 - 3 * x3 * x3;
      double z32 = 24 * x1 * x2 - 6 * x3 * x4;
      double z33 = 12 * x2 * x2 - 3 * x4 * x4;
      double z11 = -6 * a1 * a5 + eosq * (-24 * x1 * x7 - 6 * x3 * x5);
      double z12 = -6 * (a1 * a6 + a3 * a5) +  eosq *
                (-24 * (x2 * x7 + x1 * x8) - 6 * (x3 * x6 + x4 * x5));
      double z13 = -6 * a3 * a6 + eosq * (-24 * x2 * x8 - 6 * x4 * x6);
      double z21 = 6 * a2 * a5 + eosq * (24 * x1 * x5 - 6 * x3 * x7);
      double z22 = 6 * (a4 * a5 + a2 * a6) +  eosq *
                (24 * (x2 * x5 + x1 * x6) - 6 * (x4 * x7 + x3 * x8));
      double z23 = 6 * a4 * a6 + eosq * (24 * x2 * x6 - 6 * x4 * x8);
      double s3 = cc / xnq;
      double s2 = -0.5 * s3 / betao;
      double s4 = s3 * betao;
      double s1 = -15 * eo * s4;
      double s5 = x1 * x3 + x2 * x4;
      double s6 = x2 * x3 + x1 * x4;
      double s7 = x2 * x4 - x1 * x3;
      double z1 = 3 * (a1 * a1 + a2 * a2) + z31 * eosq;
      double z2 = 6 * (a1 * a3 + a2 * a4) + z32 * eosq;
      double z3 = 3 * (a3 * a3 + a4 * a4) + z33 * eosq;

      z1 = z1 + z1 + betao2 * z31;
      z2 = z2 + z2 + betao2 * z32;
      z3 = z3 + z3 + betao2 * z33;
      se = s1*zn*s5;
      si = s2*zn*(z11+z13);
      sl = -zn*s3*(z1+z3-14-6*eosq);
      sgh = s4*zn*(z31+z33-6);

      sh = 0;
      shdq = 0;
      if(xincl >= 0.052359877)    /* pi / 60 radians = 3 degrees */
      {
        sh = -zn * s2 * (z21 + z23);
        shdq = sh / sinio;
      }

      ee2 = 2*s1*s6;
      e3 = 2*s1*s7;
      xi2 = 2*s2*z12;
      xi3 = 2*s2*(z13-z11);
      xl2 = -2*s3*z2;
      xl3 = -2*s3*(z3-z1);
      xl4 = -2*s3*(-21-9*eosq)*ze;
      xgh2 = 2*s4*z32;
      xgh3 = 2*s4*(z33-z31);
      xgh4 = -18*s4*ze;
      xh2 = -2*s2*z22;
      xh3 = -2*s2*(z23-z21);

      if( !iteration)   /* we compute lunar terms only on the first pass: */
      {
         sse = se;
         ssi = si;
         ssl = sl;
         ssh = shdq;
         ssg = sgh - cosio * ssh;
         se2 = ee2;
         si2 = xi2;
         sl2 = xl2;
         sgh2 = xgh2;
         sh2 = xh2;
         se3 = e3;
         si3 = xi3;
         sl3 = xl3;
         sgh3 = xgh3;
         sh3 = xh3;
         sl4 = xl4;
         sgh4 = xgh4;
         zcosg = zcosgl;
         zsing = zsingl;
         zcosi = zcosil;
         zsini = zsinil;
         zcosh = zcoshl*cosq + zsinhl*sinq;
         zsinh = sinq*zcoshl - cosq*zsinhl;
         zn = znl;
         cc = c1l;
         ze = zel;
      }
   }
   sse += se;
   ssi += si;
   ssl += sl;
   ssg += sgh - cosio * shdq;
   ssh += shdq;

         /* "if mean motion is 1.893053 to 2.117652 revs/day, and ecc >= .5" */
   if( xnq >= 0.00826 && xnq <= 0.00924 && eo >= .5)
   {           /* start of 12-hour orbit, e >.5 section */
      double root22 = 1.7891679E-6;
      double root32 = 3.7393792E-7;
      double root44 = 7.3636953E-9;
      double root52 = 1.1428639E-7;
      double root54 = 2.1765803E-9;
      double g201 = -0.306 - (eo - 0.64) * 0.440;
      double eoc = eo * eosq;
      double sini2 = sinio*sinio;
      double f220 = 0.75*(1+2*cosio+theta2);
      double f221 = 1.5 * sini2;
      double f321 = 1.875 * sinio * (1 - 2 *\
               cosio - 3 * theta2);
      double f322 = -1.875*sinio*(1+2*
               cosio-3*theta2);
      double f441 = 35 * sini2 * f220;
      double f442 = 39.3750 * sini2 * sini2;
      double f522 = 9.84375*sinio*(sini2*(1-2*cosio-5*
                 theta2)+0.33333333*(-2+4*cosio+
                 6*theta2));
      double f523 = sinio*(4.92187512*sini2*(-2-4*
                 cosio+10*theta2)+6.56250012
                 *(1+2*cosio-3*theta2));
      double f542 = 29.53125*sinio*(2-8*
                 cosio+theta2*
                 (-12+8*cosio+10*theta2));
      double f543 = 29.53125*sinio*(-2-8*cosio+
                 theta2*(12+8*cosio-10*
                 theta2));
      double g410, g422, g520, g521, g532, g533;
      double g211, g310, g322;
      double temp, temp1;

      resonance_flag = 1;       /* it _is_ resonant... */
      synchronous_flag = 0;     /* but it's not synchronous */
             /* Geopotential resonance initialization for 12 hour orbits: */
      if (eo <= 0.65)
      {
         g211 = 3.616-13.247*eo+16.290*eosq;
         g310 = -19.302+117.390*eo-228.419*
                    eosq+156.591*eoc;
         g322 = -18.9068+109.7927*eo-214.6334*
                    eosq+146.5816*eoc;
         g410 = -41.122+242.694*eo-471.094*
                    eosq+313.953*eoc;
         g422 = -146.407+841.880*eo-1629.014*
                    eosq+1083.435*eoc;
         g520 = -532.114+3017.977*eo-5740.032*   /* NOTE: was 5740 */
                    eosq+3708.276*eoc;
      }
      else
      {
         g211 = -72.099+331.819*eo-508.738*
                    eosq+266.724*eoc;
         g310 = -346.844+1582.851*eo-2415.925*
                    eosq+1246.113*eoc;
         g322 = -342.585+1554.908*eo-2366.899*
                    eosq+1215.972*eoc;
         g410 = -1052.797+4758.686*eo-7193.992*
                    eosq+3651.957*eoc;
         g422 = -3581.69+16178.11*eo-24462.77*
                    eosq+ 12422.52*eoc;
         if (eo <= 0.715)
            g520 = 1464.74-4664.75*eo+3763.64*eosq;
         else
            g520 = -5149.66+29936.92*eo-54087.36*
                         eosq+31324.56*eoc;
      } /* End if (eo <= 0.65) */

      if (eo < 0.7)
      {
         g533 = -919.2277+4988.61*eo-9064.77*
                    eosq+5542.21*eoc;
         g521 = -822.71072+4568.6173*eo-8491.4146*
                    eosq+5337.524*eoc;
         g532 = -853.666+4690.25*eo-8624.77*
                    eosq+ 5341.4*eoc;
      }
      else
      {
         g533 = -37995.78+161616.52*eo-229838.2*
                    eosq+109377.94*eoc;
         g521 = -51752.104+218913.95*eo-309468.16*
                    eosq+146349.42*eoc;
         g532 = -40023.88+170470.89*eo-242699.48*
                    eosq+115605.82*eoc;
      } /* End if (eo <= 0.7) */

      temp1 = 3 * xnq * xnq * aqnv * aqnv;
      temp = temp1*root22;
      d2201 = temp * f220 * g201;
      d2211 = temp * f221 * g211;
      temp1 *= aqnv;
      temp = temp1*root32;
      d3210 = temp * f321 * g310;
      d3222 = temp * f322 * g322;
      temp1 *= aqnv;
      temp = 2*temp1*root44;
      d4410 = temp * f441 * g410;
      d4422 = temp * f442 * g422;
      temp1 *= aqnv;
      temp = temp1*root52;
      d5220 = temp * f522 * g520;
      d5232 = temp * f523 * g532;
      temp = 2*temp1*root54;
      d5421 = temp * f542 * g521;
      d5433 = temp * f543 * g533;
      xlamo = xmo+xnodeo+xnodeo-thgr-thgr;
      bfact = xmdot + xnodot+
                   xnodot - thdt - thdt;
      bfact += ssl + ssh + ssh;
   }        /* end of 12-hour orbit, e >.5 section */
   else if( (xnq < 0.0052359877) && (xnq > 0.0034906585))
   {                        /* "if mean motion is .8 to 1.2 revs/day" */
      double q22    =  1.7891679E-6;
      double q31    =  2.1460748E-6;
      double q33    =  2.2123015E-7;
      double cosio_plus_1 = 1. + cosio;
      double g200 = 1+eosq*(-2.5+0.8125*eosq);
      double g300 = 1+eosq*(-6+6.60937*eosq);
      double f311 = 0.9375*sinio*sinio*
             (1+3*cosio)-0.75*cosio_plus_1;
      double g310 = 1+2*eosq;
      double f220 = 0.75 * cosio_plus_1 * cosio_plus_1;
      double f330 = 2.5 * f220 * cosio_plus_1;

      resonance_flag = synchronous_flag = 1;
      /* Synchronous resonance terms initialization */
      del1 = 3*xnq*xnq*aqnv*aqnv;
      del2 = 2*del1*f220*g200*q22;
      del3 = 3*del1*f330*g300*q33*aqnv;
      del1 = del1*f311*g310*q31*aqnv;
      xlamo = xmo+xnodeo+omegao-thgr;
      bfact = xmdot + omgdot + xnodot - thdt;
      bfact = bfact+ssl+ssg+ssh;
   } /* End of geosych case */
   else              /* it's neither a high-e 12-hr orbit nor a geosynch: */
      resonance_flag = synchronous_flag = 0;

   if( resonance_flag)
   {
      xfact = bfact-xnq;

      /* Initialize integrator */
      xli = xlamo;
      xni = xnq;
      atime = 0;
   }
   /* End case dpinit: */
   
   /* SDP4 */
  double
      a, tempa, tsince_squared,
      xl, xnoddf;

  /* Update for secular gravity and atmospheric drag */
  omgadf = omegao + omgdot * tsince;
  xnoddf = xnodeo + xnodot * tsince;
  tsince_squared = tsince*tsince;
  xnode = xnoddf + xnodcf * tsince_squared;
  xn = xnodp;

  /* Update for deep-space secular effects */
  xll = xmo + xmdot * tsince;
  t = tsince;

  /* dpsec */
  int final_integration_step = 0;

  xll += ssl*t;
  omgadf += ssg*t;
  xnode += ssh*t;
  em = eo+sse*t;
  xinc = xincl+ssi*t;
  if( !resonance_flag ) goto end_dpsec;

            /* If we're closer to t=0 than to the currently-stored data
               from the previous call to this function,  then we're
               better off "restarting",  going back to the initial data.
               The Dundee code rigs things up to _always_ take 720-minute
               steps from epoch to end time,  except for the final step.
               Easiest way to arrange similar behavior in this code is
               just to always do a restart  */

            /* Epoch restart */
  atime = 0;
  xni = xnq;
  xli = xlamo;


  while( !final_integration_step)
  {
     double xldot, xlpow = 1., delt_factor;
     double delt = t - atime;
     double derivs[20];
     int s = 0;

     if( delt > dpsec_integration_step)
        delt = dpsec_integration_step;
     else if( delt < -dpsec_integration_step)
        delt = -dpsec_integration_step;
     else
        final_integration_step = 1;

     /* compute_dpsec_derivs */

     double sin_li = sin( xli);
     double cos_li = cos( xli);
     double sin_2li = 2. * sin_li * cos_li;
     double cos_2li = 2. * cos_li * cos_li - 1.;

                /* Dot terms calculated,  using a lot of trig add/subtract */
                /* identities to reduce the computational load... at the   */
                /* cost of making the code somewhat hard to follow:        */
     if( synchronous_flag )
     {
        double c_fasx2 =  0.99139134268488593;
        double s_fasx2 =  0.13093206501640101;
        double c_2fasx4 =  0.87051638752972937;
        double s_2fasx4 = -0.49213943048915526;
        double c_3fasx6 =  0.43258117585763334;
        double s_3fasx6 =  0.90159499016666422;
        double sin_3li = sin_2li * cos_li + cos_2li * sin_li;
        double cos_3li = cos_2li * cos_li - sin_2li * sin_li;
        double term1a = del1 * (sin_li  * c_fasx2  - cos_li  * s_fasx2);
        double term2a = del2 * (sin_2li * c_2fasx4 - cos_2li * s_2fasx4);
        double term3a = del3 * (sin_3li * c_3fasx6 - cos_3li * s_3fasx6);
        double term1b = del1 * (cos_li  * c_fasx2  + sin_li  * s_fasx2);
        double term2b = 2. * del2 * (cos_2li * c_2fasx4 + sin_2li * s_2fasx4);
        double term3b = 3. * del3 * (cos_3li * c_3fasx6 + sin_3li * s_3fasx6);

        for( i = 0; i < dpsec_integration_order; i += 2)
        {
           derivs[s] = term1a + term2a + term3a;
           s++;
           derivs[s] = term1b + term2b + term3b;
           s++;
           if( i + 2 < dpsec_integration_order)
           {
              term1a = -term1a;
              term2a *= -4.;
              term3a *= -9.;
              term1b = -term1b;
              term2b *= -4.;
              term3b *= -9.;
           }
        }
     }        /* end of geosynch case */
     else
     {        /* orbit is a 12-hour resonant one: */
        double c_g22 =  0.87051638752972937;
        double s_g22 = -0.49213943048915526;
        double c_g32 =  0.57972190187001149;
        double s_g32 =  0.81481440616389245;
        double c_g44 = -0.22866241528815548;
        double s_g44 =  0.97350577801807991;
        double c_g52 =  0.49684831179884198;
        double s_g52 =  0.86783740128127729;
        double c_g54 = -0.29695209575316894;
        double s_g54 = -0.95489237761529999;
        double xomi = omegaq + omgdot * atime;
        double sin_omi = sin( xomi), cos_omi = cos( xomi);
        double sin_li_m_omi = sin_li * cos_omi - sin_omi * cos_li;
        double sin_li_p_omi = sin_li * cos_omi + sin_omi * cos_li;
        double cos_li_m_omi = cos_li * cos_omi + sin_omi * sin_li;
        double cos_li_p_omi = cos_li * cos_omi - sin_omi * sin_li;
        double sin_2omi = 2. * sin_omi * cos_omi;
        double cos_2omi = 2. * cos_omi * cos_omi - 1.;
        double sin_2li_m_omi = sin_2li * cos_omi - sin_omi * cos_2li;
        double sin_2li_p_omi = sin_2li * cos_omi + sin_omi * cos_2li;
        double cos_2li_m_omi = cos_2li * cos_omi + sin_omi * sin_2li;
        double cos_2li_p_omi = cos_2li * cos_omi - sin_omi * sin_2li;
        double sin_2li_p_2omi = sin_2li * cos_2omi + sin_2omi * cos_2li;
        double cos_2li_p_2omi = cos_2li * cos_2omi - sin_2omi * sin_2li;
        double sin_2omi_p_li = sin_li * cos_2omi + sin_2omi * cos_li;
        double cos_2omi_p_li = cos_li * cos_2omi - sin_2omi * sin_li;
        double term1a =
               d2201 * (sin_2omi_p_li*c_g22 - cos_2omi_p_li*s_g22)
             + d2211 * (sin_li * c_g22 - cos_li * s_g22)
             + d3210 * (sin_li_p_omi*c_g32 - cos_li_p_omi*s_g32)
             + d3222 * (sin_li_m_omi*c_g32 - cos_li_m_omi*s_g32)
             + d5220 * (sin_li_p_omi*c_g52 - cos_li_p_omi*s_g52)
             + d5232 * (sin_li_m_omi*c_g52 - cos_li_m_omi*s_g52);
        double term2a =
               d4410 * (sin_2li_p_2omi*c_g44 - cos_2li_p_2omi*s_g44)
             + d4422 * (sin_2li * c_g44 - cos_2li * s_g44)
             + d5421 * (sin_2li_p_omi*c_g54 - cos_2li_p_omi*s_g54)
             + d5433 * (sin_2li_m_omi*c_g54 - cos_2li_m_omi*s_g54);
        double term1b =
              (d2201 * (cos_2omi_p_li*c_g22 + sin_2omi_p_li*s_g22)
             + d2211 * (cos_li * c_g22 + sin_li * s_g22)
             + d3210 * (cos_li_p_omi*c_g32 + sin_li_p_omi*s_g32)
             + d3222 * (cos_li_m_omi*c_g32 + sin_li_m_omi*s_g32)
             + d5220 * (cos_li_p_omi*c_g52 + sin_li_p_omi*s_g52)
             + d5232 * (cos_li_m_omi*c_g52 + sin_li_m_omi*s_g52));
        double term2b = 2. *
              (d4410 * (cos_2li_p_2omi*c_g44 + sin_2li_p_2omi*s_g44)
             + d4422 * (cos_2li * c_g44 + sin_2li * s_g44)
             + d5421 * (cos_2li_p_omi*c_g54 + sin_2li_p_omi*s_g54)
             + d5433 * (cos_2li_m_omi*c_g54 + sin_2li_m_omi*s_g54));

        for( i = 0; i < dpsec_integration_order; i += 2)
        {
           derivs[s] = term1a + term2a;
           s++;
           derivs[s] = term1b + term2b;
           s++;
           if( i + 2 < dpsec_integration_order)
           {
              term1a = -term1a;
              term2a *= -4.;
              term1b = -term1b;
              term2b *= -4.;
           }
        }
     } /* End of 12-hr resonant case */
     /* end compute_dpsec_derivs */

     xldot = xni+xfact;

     xli += delt * xldot;
     xni += delt * derivs[0];
     delt_factor = delt;
     for( i = 2; i <= dpsec_integration_order; i++)
     {
        xlpow *= xldot;
        derivs[i - 1] *= xlpow;
        delt_factor *= delt / (double)i;
        xli += delt_factor * derivs[i - 2];
        xni += delt_factor * derivs[i - 1];
     }
     atime += delt;
  }

  xn = xni;

  temp = -xnode + thgr + t * thdt;

  xll = xli + temp
        + (synchronous_flag ? -omgadf : temp);
  end_dpsec:
  /*End case dpsec: */


  tempa = 1-c1*tsince;
  a = pow(xke/xn,two_thirds)*tempa*tempa;
  em -= bstar*c4*tsince;

  /* Update for deep-space periodic effects */
  xll += xnodp * t2cof * tsince_squared;

  /* dpper */
  double sinis, cosis;

        /* If the time didn't change by more than 30 minutes,      */
        /* there's no good reason to recompute the perturbations;  */
        /* they don't change enough over so short a time span.     */
        /* However,  the Dundee code _always_ recomputes,  so if   */
        /* we're attempting to replicate its results,  we've gotta */
        /* recompute everything,  too.                             */
  double zf, zm, sinzf, ses, sis, sil, sel, sll, sls;
  double f2, f3, sghl, sghs, shs, sh1;

        /* Update solar perturbations for time T: */
  zm =  zmos+zns* t;
  zf = zm+2*zes*sin(zm);
  sinzf = sin(zf);
  f2 = 0.5*sinzf*sinzf-0.25;
  f3 = -0.5*sinzf*cos(zf);
  ses =  se2*f2+ se3*f3;
  sis =  si2*f2+ si3*f3;
  sls =  sl2*f2+ sl3*f3+ sl4*sinzf;
  sghs =  sgh2*f2+ sgh3*f3+ sgh4*sinzf;
  shs =  sh2*f2+ sh3*f3;

        /* Update lunar perturbations for time T: */
  zm =  zmol+znl* t;
  zf = zm+2*zel*sin(zm);
  sinzf = sin(zf);
  f2 = 0.5*sinzf*sinzf-0.25;
  f3 = -0.5*sinzf*cos(zf);
  sel =  ee2*f2+ e3*f3;
  sil =  xi2*f2+ xi3*f3;
  sll =  xl2*f2+ xl3*f3+ xl4*sinzf;
  sghl =  xgh2*f2+ xgh3*f3+ xgh4*sinzf;
  sh1 =  xh2*f2+ xh3*f3;

       /* Sum the solar and lunar contributions: */
  pe = ses+sel;
  pinc = sis+sil;
  pl = sls+sll;
  pgh = sghs+sghl;
  ph = shs+sh1;

  xinc +=  pinc;
  sinis = sin(  xinc);
  cosis = cos(  xinc);

        /* Add solar/lunar perturbation correction to eccentricity: */
  em +=  pe;
  xll +=  pl;
  omgadf +=  pgh;
  if( xincl >= 0.2)
  {
     /* Apply periodics directly */
     double sinis = sin( xinc);
     double cosis = cos( xinc);
     double temp_val =  ph / sinis;

     omgadf -= cosis * temp_val;
     xnode += temp_val;
  }
  else
  {
     /* Apply periodics with Lyddane modification */
     double sinok = sin( xnode);
     double cosok = cos( xnode);
     double alfdp =  ph * cosok
                  + ( pinc * cosis + sinis) * sinok;
     double betdp = -  ph * sinok
                  + ( pinc * cosis + sinis) * cosok;
     double dls, delta_xnode;

     delta_xnode = atan2(alfdp,betdp) -  xnode;

      /* This is a patch to Lyddane modification suggested */
      /* by Rob Matson, streamlined very slightly by BJG, to */
      /* keep 'delta_xnode' between +/- 180 degrees: */

     if( delta_xnode < - pi)
        delta_xnode += twopi;
     else if( delta_xnode > pi)
        delta_xnode -= twopi;

     dls = - xnode * sinis *  pinc;

     omgadf += dls - cosis * delta_xnode;

     xnode += delta_xnode;
  } /* End case dpper: */

  xl = xll + omgadf + xnode;

               /* Dundee change:  Reset cosio,  sinio for new xinc: */
  cosio = cos( xinc);
  sinio = sin( xinc);

  /* Long period periodics */
  // double a3ovk2 = -1.*xj3/ck2;
  double axn = em*cos(omgadf);
  temp = 1/(a*(1.-em*em));
  double xlcof = .125 * a3ovk2 * sinio * (3+5*cosio)/ (1. + cosio);
  double aycof = 0.25 * a3ovk2 * sinio;
  xll = temp*xlcof*axn;       /* different xll from above use */
  double aynl = temp*aycof;
  double xlt = xl+xll;
  double ayn = em*sin(omgadf)+aynl;
  double elsq = axn*axn+ayn*ayn;
  double capu = fmod( xlt - xnode, twopi);
  double epw = capu;
  double ecosE, esinE, r;
  double betal;
  double u, sinu, cosu, sin2u, cos2u;
  double rk, uk, xnodek, xinck;
  double sinuk, cosuk, sinik, cosik, sinnok, cosnok, xmx, xmy;
  double sinEPW, cosEPW;
  double ux, uy, uz;
  int rval = 0;

  /* Dundee changes:  items dependent on cosio get recomputed: */
  double cosio_squared = cosio * cosio;
  x3thm1 = 3.0 * cosio_squared - 1.0;
  x1mth2 = 1.0 - cosio_squared;
  double x7thm1 = 7.0 * cosio_squared - 1.0;

  /* Solve Kepler's' Equation */
  for( i = 0; i < 10; i++)
  {
     double newton_raphson_epsilon = 1e-12;
     double f, fdot, delta_epw;
     int do_second_order_newton_raphson = 1;

     sinEPW = sin(epw);
     cosEPW = cos(epw);
     ecosE = axn * cosEPW + ayn * sinEPW;
     esinE = axn * sinEPW - ayn * cosEPW;
     f = capu - epw + esinE;
     if ( fabs(f) < newton_raphson_epsilon) break;
     fdot = 1. - ecosE;
     delta_epw = f / fdot;
     if( !i)
     {
        double max_newton_raphson = 1.25 * fabs(em);

        do_second_order_newton_raphson = 0;
        if( delta_epw > max_newton_raphson)
            delta_epw = max_newton_raphson;
        else if( delta_epw < -max_newton_raphson)
                 delta_epw = -max_newton_raphson;
        else
            do_second_order_newton_raphson = 1;
     }
     if( do_second_order_newton_raphson)
         delta_epw = f / (fdot + 0.5*esinE*delta_epw);
                            /* f/(fdot - 0.5*fdotdot * f / fdot) */
     epw += delta_epw;
  }

  /* Short period preliminary quantities */
  temp = 1-elsq;
  pl = a*temp;
  r = a*(1-ecosE);
  temp2 = a / r;
  betal = sqrt(temp);
  temp = esinE/(1+betal);
  cosu = temp2 * (cosEPW - axn + ayn * temp);
  sinu = temp2 * (sinEPW - ayn - axn * temp);
  u = atan2( sinu, cosu);
  sin2u = 2*sinu*cosu;
  cos2u = 2*cosu*cosu-1;
  temp1 = ck2 / pl;
  temp2 = temp1 / pl;

  /* Update for short periodics */
  rk = r*(1-1.5*temp2*betal*x3thm1)+0.5*temp1*x1mth2*cos2u;
  uk = u-0.25*temp2*x7thm1*sin2u;
  xnodek = xnode+1.5*temp2*cosio*sin2u;
  xinck = xinc+1.5*temp2*cosio*sinio*cos2u;

  /* Orientation vectors */
  sinuk = sin(uk);
  cosuk = cos(uk);
  sinik = sin(xinck);
  cosik = cos(xinck);
  sinnok = sin(xnodek);
  cosnok = cos(xnodek);
  xmx = -sinnok*cosik;
  xmy = cosnok*cosik;
  ux = xmx*sinuk+cosnok*cosuk;
  uy = xmy*sinuk+sinnok*cosuk;
  uz = sinik*sinuk;

  /* Global Position, er, and velocity, er/min */
  rr[0] = rk*ux;
  rr[1] = rk*uy;
  rr[2] = rk*uz;

  double rdot = xke*sqrt(a)*esinE/r;
  double rfdot = xke*sqrt(pl)/r;
  xn = xke/(a * sqrt(a));
  double rdotk = rdot-xn*temp1*x1mth2*sin2u;
  double rfdotk = rfdot+xn*temp1*(x1mth2*cos2u+1.5*x3thm1);
  double vx = xmx*cosuk-cosnok*sinuk;
  double vy = xmy*cosuk-sinnok*sinuk;
  double vz = sinik*cosuk;

  vv[0] = rdotk*ux+rfdotk*vx;
  vv[1] = rdotk*uy+rfdotk*vy;
  vv[2] = rdotk*uz+rfdotk*vz;

}   /* end sdp4 */

void Satellite :: sxp4(double tsince)
{
  double a1, r1, temp, del1, ao, delo, xnodp;

  /* Period > 225 minutes is deep space */
  a1 = pow(xke / xno, 2/3.);
  r1 = cos(xincl);
  temp = ck2 * 1.5 * (r1*r1 * 3.0 - 1.0) * pow( 1.0 - eo*eo, -1.5);
  del1 = temp / (a1*a1);
  ao = a1 * (1.0 - del1 * (1./3. + del1 * (del1 * 1.654320987654321 + 1.0)));
  delo = temp / (ao*ao);
  xnodp = xno / (delo + 1.0);

  /* Select a deep-space/near-earth ephemeris */
  /* If the object makes less than 6.4 revolutions around the earth... */
  if (twopi / (xnodp * 1440.) >= (1. / 6.4))
     sdp4(tsince);    /* yes,  it should be a deep-space (SDP4) ephemeris */
  else
     sgp4(tsince);

}

/* Does a checksum modulo 10 on the given line.
   Digits = their value, '-' = 1, all other chars = 0.
   Used to write a tle*/
void ccksum(char *line)
{
    int cksum, i;

    cksum = 0;
    for (i = 0; i < 68; i++) {
        char ich;
        ich = line[i];
        if (ich > '0' && ich <= '9') cksum += ich - '0';
        if (ich == '-') cksum += 1;
    }
    line[68] = (cksum % 10) + '0';
    line[69] = '\0';
}

// write TLE to screen
void Satellite :: print_el()
{
  char line1[80], line2[80];
  char ec_string[9];
  char xns_string[12];
  char bstar_string[13];
  char bstar_fract[7];
  char bstar_exp[3];

  sprintf(ec_string, "%.7f", eo);
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

  double xns = 2160 * bstar * xno/nocon * c2;

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

  sprintf(line1, "1 %05dU %-8s %014.8f %10s  00000-0 %6s%2s 0    00",
    ssn, desig, tle, xns_string, bstar_fract, bstar_exp);
  ccksum(line1);
  sprintf(line2, "2 %05d %8.4lf %8.4lf %.7s %8.4lf %8.4lf %11.8lf    00",
    ssn, xincl/de2ra, xnodeo/de2ra, ec_string, omegao/de2ra, xmo/de2ra, xno/nocon);
  ccksum(line2);

  printf("%s\n", line1);
  printf("%s\n", line2);
}

// display state vectors
void Satellite :: print_rv()
{
  int i;
  printf("\n");
  for(i = 0; i < 3; i++)
    printf("%12.9f ", rr[i]);
  printf("\n");
  for(i = 0; i < 3; i++)
    printf("%12.9f ", vv[i]);
  printf("\n");
}

///////////// SST //////////////////////////////////////////////////////////////

  // mean elements to equinocital elements
void Satellite :: el2eq(void)
{
/*
   eq[0],   // semimajor axis
   eq[1],   // eccentricity vector component
   eq[2],   // eccentricity vector component
   eq[3],   // ascending node vector component
   eq[4],   // ascending node vector component
   eq[5],   // mean longitude
   I,   // retrograde factor
   F,   // eccentric longitude
   L;   // true longitude

   xincl,     // inclination
   xnodeo,    // right ascension of ascending node
   eo,        // eccentricity
   omegao,    // argument of the perigee
   xmo,       // mean anomaly
   xno,       // mean motion, radians/min
*/
   int i;
   double temp, E, theta;
   I = xincl < pi ? 1. : -1;
   eq[0] = pow(xke/xno, 2./3.);
   eq[1] = eo*sin(omegao + I*xnodeo);
   eq[2] = eo*cos(omegao + I*xnodeo);
   temp = pow(tan(0.5*xincl), I);
   eq[3] = temp*sin(xnodeo);
   eq[4] = temp*cos(xnodeo);
   eq[5] = xmo + omegao + I*xnodeo;
   E = xmo;
   for(i = 0; i < 10; i++)
      E = E + (xmo + eo*sin(E) - E) / (1 - eo*cos(E));
   temp = (cos(E) - eo) / (1 - eo*cos(E));
   theta = acos(temp);
   if(E > pi)
      theta = 2*pi - theta;
   F = E + omegao + I*xnodeo;
   L = theta + omegao + I*xnodeo;
}

  // equinocital elements to mean elements
void Satellite :: eq2el(void)
{
   double xi, temp, temp1, temp2;
   xno    = xke / pow(eq[0], 1.5);
   eo     = sqrt(eq[1]*eq[1] + eq[2]*eq[2]);
   temp1  = eq[1] / eo;
   temp2  = eq[2] / eo;
   xi     = atan2(temp1, temp2);
   temp   = sqrt(eq[3]*eq[3] + eq[4]*eq[4]);
   xincl  = pi*(0.5*(1-I)) + 2*I*atan(temp);
   temp1  = eq[3] / temp;
   temp2  = eq[4] / temp;
   xnodeo = fmod2p(atan2(temp1, temp2));
   omegao = fmod2p(xi - I*xnodeo);
   xmo    = fmod2p(eq[5] - xi);
}

void Satellite :: rv2eq(double* rr1, double* vv1)
{
   double e[3], f[3], g[3], hh[3], w[3], tempv[3],
          temp, b, X, Y, sinF, cosF, sinL, cosL;

   smult(1. / xke, vv1, tempv);
   cross(rr1, tempv, hh);
   temp = hh[2] / norm(hh);
   I =  acose(temp);
   if(I > pi/2) I = -1.;
   else I = 1.;

   eq[0] = 1 / (2 / norm(rr1) - dot(vv1, vv1) / (xke*xke));
   cross(rr1, vv1, tempv);
   temp = 1 / norm(tempv);
   smult(temp, tempv, w);
   temp = (1 + I * w[2]);
   eq[3] =  w[0] / temp;
   eq[4] = -w[1] / temp;

   temp = 1 / (1 + eq[3]*eq[3] + eq[4]*eq[4]);
   f[0] =  temp * (1 - eq[3]*eq[3] + eq[4]*eq[4]);
   f[1] =  temp * (2*eq[3]*eq[4]);
   f[2] =  temp * (-2*I*eq[3]);
   g[0] =  temp * (2*I*eq[3]*eq[4]);
   g[1] =  temp * (I*(1 + eq[3]*eq[3] - eq[4]*eq[4]));
   g[2] =  temp * (2*eq[4]);
   X = dot(rr1, f);
   Y = dot(rr1, g);

   // w and hh are temporary vectors here
   cross(rr1, vv1, w);
   cross(vv1, w, hh);
   temp = -1 / norm(rr1);
   smult(temp, rr1, w);
   temp = 1 / (xke*xke);
   smult(temp, hh, tempv);
   vadd(w, tempv, e);
   eq[1] = dot(e, g);
   eq[2] = dot(e, f);

   temp = sqrt(1 - eq[1]*eq[1] - eq[2]*eq[2]);
   b = 1 / (1 + temp);
   temp = eq[0] * temp;
   sinF = eq[1] + ((1 - eq[1]*eq[1] * b)*Y - eq[1]*eq[2]*b*X) / (eq[0]*temp);
   cosF = eq[2] + ((1 - eq[2]*eq[2] * b)*X - eq[1]*eq[2]*b*Y) / (eq[0]*temp);
   F = atan2(sinF, cosF);
   eq[5] = F + eq[1]*cosF - eq[2]*sinF;

   temp = 1 - eq[1]*sinF - eq[2]*cosF;
   sinL = ((1 - eq[2]*eq[2] * b)*sinF + eq[1]*eq[2]*b*cosF - eq[1]) / temp;
   cosL = ((1 - eq[1]*eq[1] * b)*cosF + eq[1]*eq[2]*b*sinF - eq[2]) / temp;
   L = atan2(sinL, cosL);
}

void Satellite :: eq2rv(double* rr1, double* vv1)
{
   // rr1 and vv1 are output values
   // input elements are global

   double f[3], g[3], w[3], tempv1[3], tempv2[3],
      temp, tempb, b, ro, X, Xdot, Y, Ydot, sinF, cosF, sinL, cosL;
   int i;

   temp = 1 / (1 + eq[3]*eq[3] + eq[4]*eq[4]);
   f[0] =  temp * (1 - eq[3]*eq[3] + eq[4]*eq[4]);
   f[1] =  temp * (2*eq[3]*eq[4]);
   f[2] =  temp * (-2*I*eq[3]);
   g[0] =  temp * (2*I*eq[3]*eq[4]);
   g[1] =  temp * (I*(1 + eq[3]*eq[3] - eq[4]*eq[4]));
   g[2] =  temp * (2*eq[4]);
   w[0] =  temp * (2*eq[3]);
   w[1] =  temp * (-2*eq[4]);
   w[2] =  temp * (I*(1 - eq[3]*eq[3] - eq[4]*eq[4]));

   tempb = sqrt(1 - eq[1]*eq[1] - eq[2]*eq[2]);
   b = 1 / (1 + tempb);
   nn = xke / pow(eq[0], 1.5);

   // Keplers equation
   F = eq[5];
   for(i = 0; i < 10; i++)
   {
      sinF = sin(F);
      cosF = cos(F);
      temp = (1 - eq[1]*sinF - eq[2]*cosF);
      F = F - (F + eq[1]*cosF - eq[2]*sinF - eq[5]) / temp;
   }

   sinL = ((1 - eq[2]*eq[2]*b)*sinF + eq[1]*eq[2]*b*cosF - eq[1]) / temp;
   cosL = ((1 - eq[1]*eq[1]*b)*cosF + eq[1]*eq[2]*b*sinF - eq[2]) / temp;
   ro = eq[0] * temp;
   X = ro * cosL;
   Y = ro * sinL;
   Xdot = -nn*eq[0]*(eq[1] + sinL) / tempb;
   Ydot =  nn*eq[0]*(eq[2] + cosL) / tempb;

   smult(X, f, tempv1);
   smult(Y, g, tempv2);
   vadd(tempv1, tempv2, rr1);
   smult(Xdot, f, tempv1);
   smult(Ydot, g, tempv2);
   vadd(tempv1, tempv2, vv1);
}

/////////////////////  Matrices and derivatives   /////////////////////

void Satellite :: Vxx(void)
{
   int n, s;
   V[0][0] = 1.;
   for(s = 0; s < N; s++)
   {
      V[s+1][s+1] = V[s][s]/(2*s+2);
   }
   for(s = 0; s < N; s++)
   {
      for(n = s; n <= N; n += 2)
      {
         V[n+2][s] = -V[n][s]*(n-s+1) / (n+s+2);
      }
   }
}

void Satellite :: GdG(double eq[6])
{
   double temp1, temp2, temp3, temp4;
   int s;

   temp1 = eq[1]*alf;  // h
   temp2 = eq[2]*alf;  // k
   temp3 = eq[1]*bet;  // h
   temp4 = eq[2]*bet;  // k

   G[0]     = 1.;
   H[0]     = 0.;
   dGh[0]   = 0.;
   dGk[0]   = 0.;
   dGalf[0] = 0.;
   dGbet[0] = 0.;
   for(s = 1; s <= N; s++)
   {
      H[s]     = (temp1 - temp4)*G[s-1] + (temp2 + temp3)*H[s-1];
      G[s]     = (temp2 + temp3)*G[s-1] - (temp1 - temp4)*H[s-1];
      dGh[s]   = s*bet*G[s-1] - s*alf*H[s-1];
      dGk[s]   = s*alf*G[s-1] + s*bet*H[s-1];
      dGalf[s] = s*eq[2]*G[s-1] - s*eq[1]*H[s-1];
      dGbet[s] = s*eq[1]*G[s-1] + s*eq[2]*H[s-1];
   }
}

void Satellite :: QdQ(void)
{
   int n, s;
   Q[0][0] = 1.;
   for(s = 1; s <= N; s++)
   {
      Q[s][s]   = (2*s-1)*Q[s-1][s-1];
      Q[s][s-1] = (2*s-1)*gam*Q[s-1][s-1];
   }
   for(s = 0; s < N-1; s++)
   {
      for(n = s+2; n <= N; n++)
      {
         Q[n][s]   = ((2*n-1)*gam*Q[n-1][s] - (n+s-1)*Q[n-2][s]) / (n-s);
      }
   }
   for(s = 0; s < N-1; s++)
   {
      for(n = s+2; n <= N; n++)
      {
         dQ[n][s] = Q[n][s+1];
      }
   }
}

// negative superscript K(-n-1,s,0)
void Satellite :: nKdK(void)
{
   double temp, x = 1 / (B*B);
   int n, s;
   for(s = 1; s <= N; s++)
   {
      temp        = pow(.5*x,s-1);
      nK[s][s]    = 0.;
      nK[s][s-1]  = temp / B;
      ndK[s][s]   = 0.;
      ndK[s][s-1] = temp * (2.*s-1);
   }
   for(s = 0; s <= N-2; s++)
   {
      for(n = s+2; n <= N; n++)
      {
         nK[n][s]  = (n-1)*x*((2*n-3)*nK[n-1][s] - (n-2)*nK[n-2][s])
                   / ((n+s-1)*(n-s-1));
         ndK[n][s] = (n-1)*x*((2*n-3)*ndK[n-1][s] - (n-2)*ndK[n-2][s])
                   / ((n+s-1)*(n-s-1)) + 2*B*nK[n][s];
      }
   }
}

// non-negative superscript K(n,s,0)
void Satellite :: KdK(void)
{
   double x = B*B;
   int n, s;
   K[0][0]   =  1.;
   K[0][1]   = -1.;
   K[1][0]   =  0.;
   dK[0][0]  =  0.;
   dK[1][0]  =  0.;
   for(s = 2; s <= N; s++)
   {
      K[s-1][s]  = (1-2*s)*K[s-2][s-1]/s;
   }
   for(s = 1; s <= N; s++)
   {
      K[s][s]    = (2*s+1)*K[s-1][s]/(s+1);
      dK[s][s]   = 0.;
      dK[s-1][s] = 0.;
   }
   for(s = 1; s <= N-1; s++)
   {
      for(n = s+1; n <= N; n++)
      {
         K[n][s]  = (2*n+1)*K[n-1][s]/(n+1) - (n+s)*(n-s)*x*K[n-2][s]/(n*(n+1));
         dK[n][s] = (2*n+1)*dK[n-1][s]/(n+1) - (n+s)*(n-s)*x*dK[n-2][s]/(n*(n+1))
                    + 2*x*B*(n+s)*(n-s)*K[n-2][s]/(n*(n+1));
      }
   }
}

/////////////////////  Zonal   /////////////////////

/*   JGM3 earth model zonal harmonics, unnormalized.   */
double J[] = {
    0.,
    0.,
   -1.0826266905978e-003,
    2.5324353457544e-006,
    1.6193312050719e-006,
    2.2771610163688e-007,
   -5.3964849049834e-007,
    3.5136844210318e-007,
    2.0251871520885e-007,
    1.1936871324412e-007,
    2.4805686499981e-007,
   -2.4056521378881e-007,
    1.8191170311845e-007,
    2.0756773243261e-007,
   -1.1741738786363e-007,
    1.7627269667946e-008,
   -3.1194308417303e-008,
    1.0713059159610e-007,
    4.4216723715108e-008,
   -2.1973339642026e-008,
    1.2031461829791e-007 };

void Satellite :: zonal(double eq[6], double rate[6])
{

   // a, h, k, p, q, u
   double temp, temp1, temp2, temp3, tmp4;
   int n, s;

          A =  xke * sqrt(eq[0]);
          B =  sqrt(1 - eq[1]*eq[1] - eq[2]*eq[2]);
         CC =  1 + eq[3]*eq[3] + eq[4]*eq[4];
        alf = -2*I*eq[3] / CC;
        bet =  2*eq[4] / CC;
        gam =  (2 - CC)*I / CC;

      GdG(eq);
      QdQ();
      nKdK();

      temp1 = xke*xke / eq[0];
      temp2 = 1 / eq[0];
      temp3 = 1 / (B*B*B);

      dUa   = 0.;
      dUh   = 0.;
      dUk   = 0.;
      dUalf = 0.;
      dUbet = 0.;
      dUgam = 0.;

      for(s = 0; s < N-1; s++)
      {
         for(n = s+2; n <= N; n += 2)
         {
            tmp4   = -J[n]*(2-delos(s))*pow(temp2,n)*V[n][s];
            dUa   += tmp4*(n+1)*nK[n][s]*Q[n][s]*G[s];
            dUh   += tmp4*Q[n][s]*(nK[n][s]*dGh[s]+eq[1]*temp3*G[s]*ndK[n][s]);
            dUk   += tmp4*Q[n][s]*(nK[n][s]*dGk[s]+eq[2]*temp3*G[s]*ndK[n][s]);
            dUalf += tmp4*nK[n][s]*Q[n][s]*dGalf[s];
            dUbet += tmp4*nK[n][s]*Q[n][s]*dGbet[s];
            dUgam += tmp4*nK[n][s]*dQ[n][s]*G[s];
         }
      }
      dUa   *=  temp1 / eq[0];
      dUh   *= -temp1;
      dUk   *= -temp1;
      dUalf *= -temp1;
      dUbet *= -temp1;
      dUgam *= -temp1;

      temp1 =  alf*dUgam - gam*dUalf;
      temp2 =  bet*dUgam - gam*dUbet;
      temp3 =  1 / (A*B);
       temp =  eq[3]*temp1 - I*eq[4]*temp2;

    rate[0] +=  0.0;
    rate[1] +=  B*dUk/A + eq[2]*temp*temp3;
    rate[2] += -B*dUh/A - eq[1]*temp*temp3;
    rate[3] += -0.5*CC*temp2*temp3;
    rate[4] += -0.5*I*CC*temp1*temp3;
    rate[5] += -2*eq[0]*dUa/A + B*(eq[1]*dUh
            +  eq[2]*dUk)/(A*(1 + B)) + temp*temp3;

}

/////////////////////  end Zonal  //////////////////


///////////////////// Sun Moon /////////////////////

void Satellite :: sun(double ep)
{
   double e, m, lm;

   double t = (ep - 2451545.0) / 36525.0;

   e  = .4090928 - .013 * t;
   m  = 6.23999888 + 628.30193263678 * t;
   lm = 4.9382346 + m + .033413359 * sin(m) + 3.49E-4 * sin(2*m) + .024386 * t;

   Zsun[0] = cos(lm);
   Zsun[1] = sin(lm) * cos(e);
   Zsun[2] = sin(lm) * sin(e);
}

void Satellite :: moon(double ep)
{
   double m, lo, l, f, d, lm, b;

   double t = (ep - 2451545.0) / 36525.0;

   m  = 6.23999888 + 628.30193263678 * t;
   lo = 3.810336 + 8399.7091055 * t - .02438574 * t;
   l  = 2.3555476 + 8328.6914252 * t;
   f  = 1.627917986 + 8433.4661791 * t;
   d  = 5.19846789 + 7771.3771439 * t;
   lm = lo + .109762 * sin(l) + .003728 * sin(2*l) - .022233 * sin(l - 2*d)
        + .01149 * sin(2*d) - .003238 * sin(m) - .001997 * sin(2*f)
        - .0010278 * sin(2*l - 2*d) - .0009987 * sin(l + m - 2*d)
        + .0009308 * sin(l + 2*d) - .0007999 * sin(m - 2*d)
        + .0007175 * sin(l - m) - .000606 * sin(d) - .0005333 * sin(l + m)
        - .0002666 * sin(2*f - 2*d);
   b  = .0897874937 * sin(f + lm - lo + .0019974 * sin(2*f) + .0026228 * sin(m))
        - .00255 * sin(f - 2*d) + .0002133 * sin(l + f - 2*d)
        - .0001503 * sin(f - l - 2*d) - .0001211 * sin(f - 2*l)
        - .0001115 * sin(f + m - 2*d) + .0001018 * sin(f - l)
        + .0000533 * sin(f - m - 2*d);

   Zmoon[0] = cos(lm) * cos(b);
   Zmoon[1] = sin(lm) * cos(b);
   Zmoon[2] = sin(b);

}

void Satellite :: sun_moon(double ep, double eq[6], double rate[6])
{
   double
   MUsun  = 1841.334741711,
   MUmoon = 6.8024503292767e-5,
   Rsun   = 23454.79363478,
   Rmoon  = 60.268389059042;

   // a, h, k, p, q, u
   double f[3], g[3], w[3], temp, temp1, temp2, temp3, tmp4;
   temp = 1 / (1 + eq[3]*eq[3] + eq[4]*eq[4]);
   f[0] =  temp * (1 - eq[3]*eq[3] + eq[4]*eq[4]);
   f[1] =  temp * (2*eq[3]*eq[4]);
   f[2] =  temp * (-2*I*eq[3]);
   g[0] =  temp * (2*I*eq[3]*eq[4]);
   g[1] =  temp * (I*(1 + eq[3]*eq[3] - eq[4]*eq[4]));
   g[2] =  temp * (2*eq[4]);
   w[0] =  temp * (2*eq[3]);
   w[1] =  temp * (-2*eq[4]);
   w[2] =  temp * (I*(1 - eq[3]*eq[3] - eq[4]*eq[4]));

   // Sun
   sun(ep);

   alf = dot(Zsun, f);
   bet = dot(Zsun, g);
   gam = dot(Zsun, w);

   GdG(eq);
   QdQ();
   KdK();

   temp1 = MUsun / Rsun;
   temp2 = eq[0] / Rsun;
   temp3 = 1 / (B*B*B);

   dUa   = 0.;
   dUh   = 0.;
   dUk   = 0.;
   dUalf = 0.;
   dUbet = 0.;
   dUgam = 0.;

   int s, n;
   for(n = 2; n <= 8; n++)
   {
      for(s = n; s >= 0; s = s-2)
      {
         tmp4   = (2-delos(s))*pow(temp2,n)*V[n][s];
         dUa   += tmp4*n*K[n][s]*Q[n][s]*G[s];
         dUh   += tmp4*Q[n][s]*(K[n][s]*dGh[s] + eq[1]*temp3*G[s]*dK[n][s]);
         dUk   += tmp4*Q[n][s]*(K[n][s]*dGk[s] + eq[2]*temp3*G[s]*dK[n][s]);
         dUalf += tmp4*K[n][s]*Q[n][s]*dGalf[s];
         dUbet += tmp4*K[n][s]*Q[n][s]*dGbet[s];
         dUgam += tmp4*K[n][s]*dQ[n][s]*G[s];
      }
   }
   dUa   *=  temp1 / eq[0];
   dUh   *=  temp1;
   dUk   *=  temp1;
   dUalf *=  temp1;
   dUbet *=  temp1;
   dUgam *=  temp1;

   temp1 =  alf*dUgam - gam*dUalf;
   temp2 =  bet*dUgam - gam*dUbet;
   temp3 =  1 / (A*B);
   temp  =  eq[3]*temp1 - I*eq[4]*temp2;

   rate[0] +=  0.;
   rate[1] +=  B*dUk/A - eq[2]*temp*temp3;
   rate[2] += -B*dUh/A + eq[1]*temp*temp3;
   rate[3] += -0.5*CC*temp2*temp3;
   rate[4] += -0.5*I*CC*temp1*temp3;
   rate[5] += -2*eq[0]*dUa/A + B*(eq[1]*dUh
           +  eq[2]*dUk)/(A*(1 + B)) + temp*temp3;


   // Moon
   moon(ep);

   alf = dot(Zmoon, f);
   bet = dot(Zmoon, g);
   gam = dot(Zmoon, w);

   GdG(eq);
   QdQ();
   KdK();

   temp1 = MUmoon / Rmoon;
   temp2 = eq[0] / Rmoon;
   temp3 = 1 / (B*B*B);

   dUa   = 0.;
   dUh   = 0.;
   dUk   = 0.;
   dUalf = 0.;
   dUbet = 0.;
   dUgam = 0.;

   for(n = 2; n <= 8; n++)
   {
      for(s = n; s >= 0; s = s-2)
      {
         tmp4   = (2-delos(s))*pow(temp2,n)*V[n][s];
         dUa   += tmp4*n*K[n][s]*Q[n][s]*G[s];
         dUh   += tmp4*Q[n][s]*(K[n][s]*dGh[s] + eq[1]*temp3*G[s]*dK[n][s]);
         dUk   += tmp4*Q[n][s]*(K[n][s]*dGk[s] + eq[2]*temp3*G[s]*dK[n][s]);
         dUalf += tmp4*K[n][s]*Q[n][s]*dGalf[s];
         dUbet += tmp4*K[n][s]*Q[n][s]*dGbet[s];
         dUgam += tmp4*K[n][s]*dQ[n][s]*G[s];
      }
   }
   dUa   *=  temp1 / eq[0];
   dUh   *=  temp1;
   dUk   *=  temp1;
   dUalf *=  temp1;
   dUbet *=  temp1;
   dUgam *=  temp1;

   temp1 =  alf*dUgam - gam*dUalf;
   temp2 =  bet*dUgam - gam*dUbet;
   temp3 =  1 / (A*B);
   temp  =  eq[3]*temp1 - I*eq[4]*temp2;

   rate[0] +=  0.;
   rate[1] +=  B*dUk/A - eq[2]*temp*temp3;
   rate[2] += -B*dUh/A + eq[1]*temp*temp3;
   rate[3] += -0.5*CC*temp2*temp3;
   rate[4] += -0.5*I*CC*temp1*temp3;
   rate[5] += -2*eq[0]*dUa/A + B*(eq[1]*dUh
           +  eq[2]*dUk)/(A*(1 + B)) + temp*temp3;
}

/////////////////// end Sun Moon ///////////////////

////////////////////////  DRAG  ////////////////////

/* Modified Harris-Priester density by Montenbruck and Gill */
double Satellite :: rho(double h)
{
    int i, ih;
    double h_min, h_max, d_min, d_max, c_psi2;
    double ra_sun, dec_sun, c_dec;
    double uv[3];

    if(h > 2000.)
    {
       REF = 0.0;
       return 0.0;
    }

    // index search
    for(i = 0; i < 61; i++)
    {
        if(h > dens[i][0] && h < dens[i+1][0])
        {
            ih = i;
            break;
        }
    }

    // exponential density interpolation
    h_min = (dens[ih][0] - dens[ih+1][0])
            / log(dens[ih+1][1] / dens[ih][1]);
    h_max = (dens[ih][0] - dens[ih+1][0])
            / log(dens[ih+1][2] / dens[ih][2]);
    d_min = dens[ih][1] * exp((dens[ih][0] - h) / h_min);
    d_max = dens[ih][2] * exp((dens[ih][0] - h) / h_max);
    REF = (d_min + d_max) / 2.;   // average density for Kdr

    // Sun right ascension, declination
    ra_sun  = atan2(Zsun[1], Zsun[0]);
    dec_sun = atan2(Zsun[2], sqrt(Zsun[0]*Zsun[0] + Zsun[1]*Zsun[1]));

    // unit vector to apex of diurnal bulge
    c_dec = cos(dec_sun);
    uv[0] = c_dec * cos(ra_sun + 0.523599);
    uv[1] = c_dec * sin(ra_sun + 0.523599);
    uv[2] = sin(dec_sun);

    // cosine of half angle between satellite position vector
    // and apex of diurnal bulge
    c_psi2 = 0.5 + 0.5 * dot(rr, uv) / norm(rr);

    return (d_min + (d_max - d_min) * pow(c_psi2 ,3));
}

void Satellite :: drag(double ep, double eq[6])
{

   DG[0] = 0.0;
   DG[1] = 0.0;
   DG[2] = 0.0;
   DG[3] = 0.0;
   DG[4] = 0.0;
   DG[5] = 0.0;

   if((eq[0]*(1. - eo) - 1.)*xkmper > 2000.) return;  // perige > 2000 km

   double f[3], g[3], w[3], tempv1[3], tempv2[3],
      rr1[3], vv1[3], q[3], vrel[3], vmag, qmag, t,
      temp, ro, roa, X, Xdot, Y, Ydot, tempAB, hgt,
      sinL, cosL, tempv3[3], temp1, tempD, roa2;

   double dav[3], dhv[3], dkv[3], dpv[3], dqv[3], duv[3];

   double
      wrot  = .004375269,  // rotational velocity er/min
      twopi = pi * 2.,
      pio18 = pi / 18.,
      xkeo2 = 1 / (xke*xke);     // 1/mu

    A =  xke * sqrt(eq[0]);
    B =  sqrt(1 - eq[1]*eq[1] - eq[2]*eq[2]);
   CC =  1 + eq[3]*eq[3] + eq[4]*eq[4];

   tempAB = 1 / (A*B);
   tempD  = 1 / (twopi * B);
   temp   = 1 / CC;
   f[0] =  temp * (1 - eq[3]*eq[3] + eq[4]*eq[4]);
   f[1] =  temp * (2*eq[3]*eq[4]);
   f[2] =  temp * (-2*I*eq[3]);
   g[0] =  temp * (2*I*eq[3]*eq[4]);
   g[1] =  temp * (I*(1 + eq[3]*eq[3] - eq[4]*eq[4]));
   g[2] =  temp * (2*eq[4]);
   w[0] =  temp * (2*eq[3]);
   w[1] =  temp * (-2*eq[4]);
   w[2] =  temp * (I*(2 - CC));

   nn = xke / pow(eq[0], 1.5);       // er / min
   t  = ep + pi / (1440.*nn);        // add time for half orbit, days
   sun(t);                           // return global Zsun

   L += pi / 36.;
   // drag rates integral by midpoint rule
   for(int i = 0; i < 36; i++)    // 10 degree stepsize
   {

      sinL = sin(L);
      cosL = cos(L);
      roa  = B*B / (1 + eq[1]*sinL + eq[2]*cosL);
      ro   = eq[0] * roa;
      roa2 = roa*roa;
      X = ro * cosL;
      Y = ro * sinL;
      Xdot = -nn*eq[0]*(eq[1] + sinL) / B;
      Ydot =  nn*eq[0]*(eq[2] + cosL) / B;

      // rr1 and vv1
      smult(X, f, tempv1);
      smult(Y, g, tempv2);
      vadd(tempv1, tempv2, rr1);
      smult(Xdot, f, tempv1);
      smult(Ydot, g, tempv2);
      vadd(tempv1, tempv2, vv1);

      // vv + rr x w
      vrel[0] = vv1[0] + wrot * rr1[1];
      vrel[1] = vv1[1] - wrot * rr1[0];
      vrel[2] = vv1[2];

      vmag = norm(vrel);
      temp = rr1[2] / ro;
       hgt = (ro - (1 - temp * temp / 298.257)) * xkmper;
      qmag = -Kdr * rho(hgt) * vmag;

      // drag acceleration vector
      q[0] = qmag * vrel[0];
      q[1] = qmag * vrel[1];
      q[2] = qmag * vrel[2];

      // partials wrt velocity
      // dav
      temp = 2./A;
      smult(temp, vv1, dav);

      // dhv
      temp = (2.* Xdot * Y - X * Ydot) * xkeo2;
      smult(temp, f, tempv2);
      temp = -X * Xdot * xkeo2;
      smult(temp, g, tempv3);
      vadd(tempv2, tempv3, tempv1);
      temp1 = (I*eq[4]*Y - eq[3]*X) * tempAB;
      temp = eq[2] * temp1;
      smult(temp, w, tempv2);
      vadd(tempv1, tempv2, dhv);

      // dkv
      temp = (2.* X * Ydot - Xdot * Y) * xkeo2;
      smult(temp, g, tempv2);
      temp = -Y * Ydot * xkeo2;
      smult(temp, f, tempv3);
      vadd(tempv2, tempv3, tempv1);
      temp = -eq[1] * temp1;
      smult(temp, w, tempv2);
      vadd(tempv1, tempv2, dkv);

      // dpv
      temp = 0.5 * CC * Y * tempAB;
      smult(temp, w, dpv);

      // dqv
      temp = 0.5 * I * CC * X * tempAB;
      smult(temp, w, dqv);

      // duv
      temp = eq[2] / (1 + B);
      smult(temp, dhv, tempv1);
      temp = -eq[1] / (1 + B);
      smult(temp, dkv, tempv2);
      vadd(tempv1, tempv2, tempv3);
      temp = -2./A;
      smult(temp, rr1, tempv2);
      vadd(tempv2, tempv3, tempv1);
      temp = temp1 * B;
      smult(temp, w, tempv2);
      vadd(tempv1, tempv2, duv);

      temp   = roa2 * pio18 * tempD;
      DG[0] += temp * dot(dav, q);
      DG[1] += temp * dot(dhv, q);
      DG[2] += temp * dot(dkv, q);
      DG[3] += temp * dot(dpv, q);
      DG[4] += temp * dot(dqv, q);
      DG[5] += temp * dot(duv, q);

      L += pio18;

   } // end for
   L += pi / 36.;    // restore original L
}

/////////////////////  end DRAG  ///////////////////

//////////////////////  SRP  ///////////////////////

void Satellite :: srp(double ep, double eq[6], double rate[6])
{

   double
   MS = 84,        // kg
   AR = .55E-6,    // km^2
   CR = 1.5;       // reflectivity

   Ksrp = 4.51E-3*CR*AR*1.5*eq[0]/MS;     // km / s^2
   Ksrp = Ksrp*3600/xkmper;               // er / min^2

   //Ksrp = Ksrp*1.5*eq[0];

   // a, h, k, p, q, u
   double f[3], g[3], w[3], temp, tempAB, tempKP;
   temp = 1 / (1 + eq[3]*eq[3] + eq[4]*eq[4]);
   f[0] =  temp * (1 - eq[3]*eq[3] + eq[4]*eq[4]);
   f[1] =  temp * (2*eq[3]*eq[4]);
   f[2] =  temp * (-2*I*eq[3]);
   g[0] =  temp * (2*I*eq[3]*eq[4]);
   g[1] =  temp * (I*(1 + eq[3]*eq[3] - eq[4]*eq[4]));
   g[2] =  temp * (2*eq[4]);
   w[0] =  temp * (2*eq[3]);
   w[1] =  temp * (-2*eq[4]);
   w[2] =  temp * (I*(1 - eq[3]*eq[3] - eq[4]*eq[4]));

   // Sun
   sun(ep);

   alf = dot(Zsun, f);
   bet = dot(Zsun, g);
   gam = dot(Zsun, w);

   tempAB = 1 / (A*B);
   tempKP = eq[2]*eq[3] - I*eq[1]*eq[4];

   rate[0] +=  0.;
   rate[1] +=  Ksrp*(B*alf/A - eq[2]*gam*tempAB*tempKP);
   rate[2] += -Ksrp*(B*bet/A + eq[1]*gam*tempAB*tempKP);
   rate[3] +=  0.5*Ksrp*CC*eq[1]*gam*tempAB;
   rate[4] +=  0.5*Ksrp*CC*eq[2]*gam*tempAB*I;
   rate[5] += -(2+B)*Ksrp*(eq[2]*alf+eq[1]*bet)/(A*(1+B))-Ksrp*gam*tempAB*tempKP;

}


////////////////  end solar radiation pressure  //////////////////

void Satellite :: rates(double ep, double eq[6], double rate[6])
{

   // global drag acceleration vector
   rate[0] = DG[0];
   rate[1] = DG[1];
   rate[2] = DG[2];
   rate[3] = DG[3];
   rate[4] = DG[4];
   rate[5] = DG[5] + xke * pow(eq[0], -1.5);

   zonal(eq, rate);
   sun_moon(ep, eq, rate);
   srp(ep, eq, rate);

}

void Satellite :: rk4(double ep, double *eq, double dt)
{

   // single step
   int j;
   double eq1[6], rate[6], k1[6], k2[6], k3[6], k4[6];
   double dt1;

   // k1
   rates(ep, eq, rate);
   for(j = 0; j < 6; j++)
   {
      k1[j] = rate[j];
      eq1[j] = eq[j] + 0.5*dt*k1[j];
   }

   // k2
   dt1 = ep + 0.5*dt/1440.;
   rates(dt1, eq1, rate);
   for(j = 0; j < 6; j++)
   {
      k2[j] = rate[j];
      eq1[j] = eq[j] + 0.5*dt*k2[j];
   }

   // k3
   //dt1 = ep + 0.5*dt/1440.;
   rates(dt1, eq1, rate);
   for(j = 0; j < 6; j++)
   {
      k3[j] = rate[j];
      eq1[j] = eq[j] + dt*k3[j];
   }

   // k4
   dt1 = ep + dt/1440.;
   rates(dt1, eq1, rate);
   for(j = 0; j < 6; j++)
   {
      k4[j] = rate[j];
      eq[j] += dt*(k1[j] + 2*k2[j] + 2*k3[j] + k4[j])/6;
   }
}

double Satellite :: sst(double target, int NODE)
{
   double dt, dt1, jdo;

   // epochs are days, stepsize is minutes
   // 5 orbit stepsize
   dt  = 7200. * nocon / xno;   // stepsize = 5 revolutions (in minutes)
   if (target < 0.) dt *= -1;   // if step backwards
   dt1 = 0.0;
   jdo = jd;

   while(1)                     // loop until target time
   {

      // near target
      if(fabs(dt) > fabs(target - dt1))
      {
         // final step goes exactly to fixed target
         dt = (target - dt1);
         rk4(jdo, eq, dt);      // last step
         jdo += dt / 1440.;
         // if NODE, advance to ascending node
         if(NODE) goto node;
         else
         {
            eq2rv(rr, vv);
            return(jdo);
         }
      }  // end target

      dt1 += dt;                // elapsed minutes counter
      drag(jdo, eq);
      rk4(jdo, eq, dt);
      jdo += dt / 1440.;        // advance epoch one step (in days)
   }

   node:
   dt = 180. * nocon / xno;   // always go forward to node
   while (1)  // loop until target time
   {
      eq2rv(rr, vv);
      // find ascending node (can take more than one step)
      // If ascending part
      if (vv[2] > 0)
      {
         double tmp = -rr[2] / vv[2];
         if(fabs(tmp) < fabs(dt)) dt = tmp;

         // Close enough to equator?
         if (fabs(dt) < 0.0001)
         {
            return(jdo);
         }
      }  // approach the ascending node
      // solve by rk4
      rk4(jdo, eq, dt);
      // Update jdo, epoch
      jdo += dt / 1440.;
   }
}


