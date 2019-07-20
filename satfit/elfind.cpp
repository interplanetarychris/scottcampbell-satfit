// elfind.cpp
// by Scott Campbell campbel7@the-i.net

/******************************************************************************
*                                                                             *
*    Elfind finds and reads observations in the IOD format and computes a     *
*    set of orbital elements matching the observations.  The observations     *
*    must be made on a single pass from a single location.                    *
*                                                                             *
*    The default file for input is named "unid.txt".  However, any file       *
*    name may be used as a command line argument to elfind. i.e.              *
*         elfind input.txt                                                    *
*    would search for observations in the file "input.txt" and append the     *
*    computed TLE at the end of the same file.  If two command line arguments *
*    are found, elfind interprets the first as the input file and the second  *
*    as the output file. If elfind is executed with no command argument i.e.  *
*         elfind                                                              *
*    the file "unid.txt" is searched for properly formatted observations      *
*    and the computed TLE is appended to the end of the file "unid.txt".      *
*                                                                             *
*    A properly formatted observation is recorded as a single line in one of  *
*    the four IOD formats using right ascension and declination (1, 2, 3, 7)  *
*    Two observations gives the program enough information to find the        *
*    orbital elements of the circular earth orbit through those observations. *
*    Three observations finds all parameters but the atmoshperic drag.  A     *
*    guess is made for the drag element.                                      *
*                                                                             *
*******************************************************************************/

#include <cstdlib>  // for atof
#include <ctype.h>  // for isdigit
#include "date.h"   // for Date, Locate
using namespace std;

///////////////// PHYSICAL CONSTANTS //////////////////////////////////////////

/* dimensions & gravity of earth, World Geodetic System 1972 values */

const double
de2ra = .0174532925199433,
tu = .0093380977083333,   // time unit in days ~ 13 minutes
pi = 3.14159265358979323846,
xj3 = -2.53881E-6,
ck2 = 5.413079E-4,
xke = 7.43669161E-2;   /* = (G*M)^(1/2)*(er/min)^(3/2) where G =
                               Newton's grav const, M = earth mass */

///////////// declare global variables ////////////////////////////////////////

double jd, eo, xincl, xnodeo, omegao, xmo, xno;
double la, lo, hh;
int nobs, ssn;
char* fl;
char buf[6];
char desig[] = "0000000";
double rr2[3], vv2[3];

//////////// declare functions ////////////////////////////////////////////////

void get_obs(char *file_in, char iod_line[][81]);
void read_obs(char iod_line[][81],
            double odata[][3],
            double ll[][3],
            double rd[][3]);
void read_site(char line[]);
void so2rv(double od[][3], double rd[][3], double ll[][3]);
void so3rv(double od[][3], double rd[][3], double ll[][3]);
void write_tle(char *file_out);
inline char *s_in(char *prompt, char *buffer);

/////////////////// MAIN //////////////////////////////////////////////////////
int main(int argc, char *argv[])
{
  int a, b, c;
  char iod_line[50][81];
  char line[3][81];
  char buf[12];
  double odata[3][3];
  char file_in[80], file_out[80];

  if(argc == 2)
  {
    sprintf(file_in, argv[1]);
    sprintf(file_out, argv[1]);
  }
  else if(argc == 3)
  {
    sprintf(file_in, argv[1]);
    sprintf(file_out, argv[2]);
  }
  else          // default input/output files
  {
    sprintf(file_in, "unid.txt");
    sprintf(file_out, "unid.txt");
  }

  get_obs(file_in, iod_line);        // passes nobs globally

  // designate row numbers of obs to be used
  if(nobs == 2)
  {
    a = 1;
    b = 2;
    goto jump;
  }

  top:
  s_in("\nEnter Q or the row numbers of 2 or 3 observations: ", buf);
  if ((*buf & 0x5f) == 'Q') return 0;
  nobs = sscanf(buf, "%d %d %d", &a, &b, &c);

  jump:
  line[0] = iod_line[a - 1];
  line[1] = iod_line[b - 1];
  if(nobs == 3) line[2] = iod_line[c - 1];

  double ll[nobs][3];       // line of sight vectors; observer -> satellite
  double rd[nobs][3];       // topocentric vectors to observer positions

  read_obs(line, odata, ll, rd);

  if(nobs == 2) so2rv(odata, rd, ll);
  else so3rv(odata, rd, ll);

  jd = odata[1][0];
  write_tle(file_out);
  goto top;

}

////////////////// end main ///////////////////////////////////////////////////

inline char *s_in(char *prompt, char *buffer)
{
  printf("%s", prompt);
  if (feof(stdin)) exit(0);
  gets(buffer);
  return buffer;
}

inline int SGN(double var)
{
  return (var < 0) ? -1 : 1;
}
//   mod(var, div) = remainder of var / div
inline double mod(double var, double div)
{
  return((var < 0) ? div + fmod(var, div) : fmod(var, div));
}
//   v1 . v2
inline double dot(double v1[], double v2[])
{
  double sum = 0;
  for (int i = 0; i < 3; i++)
     sum += v1[i] * v2[i];
  return sum;
}
//   ||v||
inline double norm(double v[])
{
  return sqrt(dot(v, v));
}
//  av = a * v
inline void smult(double a, double v[], double av[])
{
  for (int i = 0; i < 3; i++)
    av[i] = a * v[i];
}
//   s = v1 + v2
inline void vadd(double v1[], double v2[], double s[])
{
  for (int i = 0; i < 3; i++)
    s[i] = v1[i] + v2[i];
}
//   s = v1 + a * v2   (used for subtraction)
inline void vmadd(double v1[], double v2[], double s[], double a)
{
  for (int i = 0; i < 3; i++)
    s[i] = v1[i] + a * v2[i];
}
//   u = v / ||v||
inline void unitv(double v[], double u[])
{
  double no = norm(v);
  for (int i = 0; i < 3; i++)
    u[i] = v[i] / no;
}
//   b = v1 x v2
inline void cross(double v1[3], double v2[3], double b[3])
{
  b[0] = v1[1] * v2[2] - v1[2] * v2[1];
  b[1] = v1[2] * v2[0] - v1[0] * v2[2];
  b[2] = v1[0] * v2[1] - v1[1] * v2[0];
}
double fmod2p( const double x)
{
   double twopi = 2 * pi;
   double rval = fmod( x, twopi);

   if(rval < 0.)
      rval += twopi;
   return(rval);
}

  // vectors to mean elements
void rv2el(double* rr2, double* vv2)
{
    int i;

    /* classical osculating orbit elements calculated from vectors rr2, vv2  */
    double xinck, xnodek, ek, mk, wk, xn, rk, uk, aodp, pl, rdotk, rfdotk, temp;
    double twopi = 2*pi;

    double h[3], n[3], vec[3], vk[3];
    smult(1. / xke, vv2, vk);
    cross(rr2, vk, h);
    pl = dot(h, h);
    double vz[] = {0, 0, 1};
    double vy[3], t;
    cross(vz, h, n);
    rk = norm(rr2);
    rdotk = dot(rr2, vv2) / rk;
    rfdotk = norm(h) * xke / rk;
    temp = dot(rr2, n) / rk / norm(n);
    if(fabs(temp) > 1) temp = SGN(temp);
    uk = acos(temp);
    if(rr2[2] < 0) uk = twopi - uk;
    cross(vk, h, vz);
    smult(-1 / rk, rr2, vy);
    vadd(vz, vy, vec);
    ek = norm(vec);
    temp = h[2] / norm(h);
    if(fabs(temp) > 1) temp = SGN(temp);
    xinck =  acos(temp);
    temp = n[0] / norm(n);
    if(fabs(temp) > 1) temp = SGN(temp);
    xnodek = acos(temp);
    if(n[1] < 0)
      xnodek = fmod2p(twopi - xnodek);
    temp = dot(vec, n) / ek / norm(n);
    if(fabs(temp) > 1) temp = SGN(temp);
    wk = acos(temp);
    if(vec[2] < 0)
      wk = fmod2p(twopi - wk);
    aodp = pl / (1 - ek*ek);
    xn = xke * pow(aodp, -1.5);

    double cosio, sinio, sin2u, cos2u, temp1, temp2,
     rdot, rfdot, theta2, betal, x3thm1, x1mth2, x7thm1,
     esine, ecose, elsq, cosepw, sinepw, axn, ayn,
     cosu, sinu, capu, a3ovk2, xlcof, aycof, aynl, xll,
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
      sin2u = sin(2*u);
      cos2u = cos(2*u);
      theta2 = cosio * cosio;
      x3thm1 = 3 * theta2 - 1.;
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
      temp = 2 / r - (rdot*rdot + rfdot*rfdot) / (xke*xke);
      aodp = 1 / temp;

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
    a3ovk2 = -xj3 / ck2;
    xlcof = .125 * a3ovk2 * sinio * (3. + 5. * cosio)
          / (1. + cosio);
    aycof = .25 * a3ovk2 * sinio;
    temp1 = esine / (1 + sqrt(1 - elsq));
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
       beta = 1 - eo*eo;
       temp = 1. / (aodp * beta);
       aynl = temp * aycof;
       ayn = eo * sin(omegao) + aynl;
       cosepw = r * cosu / aodp + axn - ayn * temp1;
       sinepw = r * sinu / aodp + ayn + axn * temp1;
       axn = cosepw * ecose + sinepw * esine;
       ayn = sinepw * ecose - cosepw * esine;
       omegao = fmod2p(atan2(ayn - aynl, axn));
       eo = axn / cos(omegao);
       if(fabs(a2 - eo) < 1.e-13) break;
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
       a0 = aodp * (1 - d0);
       d1 = temp / (a1*a1);
       a1 = a0 / (1 - d1 / 3 - d1*d1 - 134 * d1*d1*d1 / 81);
       if(fabs(a2 - a1) < 1.e-13) break;
    }
    xno = pow(a1 , -1.5) / (tu * twopi);

} /* end rv2el  */

/* find topocentric vector, rr, to satellite given the line of sight unit
vector, ll, from the observer to the satellite, the topocentric vector, rd,
to the observer, and the length of the vector, r. */
void so2r(double r, double rd[], double ll[], double rr[])
{
  double ang1, ang2, nrd, rho;

  nrd = norm(rd);
  ang1 = acos(dot(rd, ll) / nrd);
  if(ang1 < .001) rho = r - nrd;
  else
  {
    ang2 = asin(nrd * sin(ang1) / r);
    rho = r * sin(ang1 - ang2) / sin(ang1);
  }
  vmadd(rd, ll, rr, rho);
}

// fit circular orbit to two observations
void so2rv(double od[][3], double rd[][3], double ll[][3])
{
  double rro = 0;
  double rr = 1.1;
  double rr1[3], rrx[3], rry[3], vv, theta;

  double delt = od[1][0] - od[0][0];

  // use the first two positions in the array
  theta = acos(dot(ll[0], ll[1]));
  vmadd(ll[1], ll[0], rrx, -1);
  cross(rd[1], rrx, rry);
  double sin_phi = norm(rry) / (norm(rd[1])*norm(rrx));
  sin_phi = sin_phi * sin_phi;   // squared

  // initial estimate for rr
  for(int i = 0; i < 3; i++)
  {
     rr = 1 + delt * sin_phi * 53.575 / (tan(theta / 2) * sqrt(rr));
     if(rr > 7.) rr = 7.;
  }

  while(abs(rr - rro) > 1.E-8)
  {
    so2r(rr, rd[0], ll[0], rr1);
    so2r(rr, rd[1], ll[1], rr2);
    theta = acos(dot(rr1, rr2) / (rr*rr));
    vv = theta * rr * tu / delt;
    rro = 1 / (vv*vv);
    rr = .01 * (99 * rr + rro);   // weighted average 100:1
    if(rr > 8.)
    {
       rr = 8.;
       break;
    }
    if(rr < 1.002)
    {
       rr = 1.002;
       break;
    }
  }
  cross(rr1, rr2, rrx);           // placeholder for intermediate calculation
  cross(rrx, rr2, rry);           // placeholder for intermediate calculation
  unitv(rry, rrx);                // placeholder for intermediate calculation
  smult(vv, rrx, vv2);

  smult(xke, vv2, vv2);
  rv2el(rr2, vv2);

}

//   gaussian elimination
void rref(double m[][7], double b[])    // scaled partial pivoting
{
  int i, j, k, ROW;
  double s[6], bin, mult;

////////////////// calculate scale factors ////////////////////////////////////

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

///////////////////// forward elimination /////////////////////////////////////

    for(i = j + 1; i < 6; i++)
    {
      mult = m[i][j] / m[j][j];
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
    printf("Singular matrix\n");
    s_in("[exit]", buf);
    exit(0);
  }

////////////////// back sustitution //////////////////////////////////////////

  b[5] = m[5][6] / m[5][5];
  for(i = 4; i >= 0; i--)
  {
    bin = 0;
    for(k = i + 1; k < 6; k++)
      bin = bin + m[i][k] * b[k];
    b[i] = (m[i][6] - bin) / m[i][i];
  }
}  // end rref

// F and G series
void f8g(double rr1[], double vv1[], double delt, double& fr, double& gr)
{
  double r, u, p, q, p2, p4, q2, u2;
  double f[9], g[9];
  r = norm(rr1);
  u = 1 / (r*r*r);
  p = dot(rr1, vv1) / (r*r);
  q = dot(vv1, vv1) / (r*r) - u;
  p2 = p * p;
  p4 = p2 * p2;
  q2 = q * q;
  u2 = u * u;
  f[0] = 1;
  f[1] = 0;
  f[2] = -u / 2;
  f[3] = p * u / 2;
  f[4] = u * (u - 3 * (5 * p2 - q)) / 24;
  f[5] = -p * u * (u - 7 * p2 + 3 * q) / 8;
  f[6] = -u * (u2 - 6 * (35 * p2 - 4 * q) * u
         + 45 * (21 * p4 - 14 * p2 * q + q2)) / 720;
  f[7] = p * u * (u2 - 2 * (25 * p2 - 7 * q) * u
         + 5 * (33 * p4 - 30 * p2 * q + 5 * q2)) / 80;
  f[8] = u * (u2*u - 9 * (245 * p2 - 13 * q) * u2
         + 27 * (1925 * p4 - 910 * p2 * q + 41 * q2) * u
         - 315 * (429 * p4*p2 - 495 * p4 * q
         + 135 * p2 * q2 - 5 * q2*q)) / 40320;
  g[0] = 0;
  g[1] = 1;
  g[2] = 0;
  g[3] = -u / 6;
  g[4] = p * u / 4;
  g[5] = u * (u - 9 * (5 * p2 - q)) / 120;
  g[6] = -p * u * (u - 2 * (7 * p2 - 3 * q)) / 24;
  g[7] = -u * (u2 - 18 * (35 * p2 - 3 * q) * u
         + 225 * (21 * p4 - 14 * p2 * q + q2)) / 5040;
  g[8] = p * u * (u2 - 4 * (25 * p2 - 6 * q) * u
         + 15 * (33 * p4 - 30 * p2 * q + 5 * q2)) / 320;
  fr = 1 + delt*delt * (f[2] + delt * (f[3] + delt * (f[4]
       + delt * (f[5] + delt * (f[6] + delt * (f[7] + delt * f[8]))))));
  gr = delt * (1 + delt*delt * (g[3] + delt * (g[4] + delt * (g[5]
       + delt * (g[6] + delt * (g[7] + delt * g[8]))))));
}

// fit orbit to three observations
void so3rv(double od[][3], double rd[][3], double ll[][3])
{
  double delt1, delt3, fr, gr, f1, f3, g1, g3;
  double rr1[3], vv1[3], rvx[3], zz[6];
  //  initial orbit
  so2rv(od, rd, ll);
  //  refine initial orbit
  rr1 = rr2;
  vv1 = vv2;
  delt1 = (od[0][0] - od[1][0]) / tu;
  delt3 = (od[2][0] - od[1][0]) / tu;
  rr2[0] = 0; rr2[1] = 0; rr2[2] = 0;

  do
  {
    f8g(rr1, vv1, delt1, fr, gr);
    f1 = fr;
    g1 = gr;
    f8g(rr1, vv1, delt3, fr, gr);
    f3 = fr;
    g3 = gr;
    rr2 = rr1;
    vv2 = vv1;

    double mat[6][7] =
     {{f1 * ll[0][2],  0, -f1 * ll[0][0],  g1 * ll[0][2], 0,
      -g1 * ll[0][0], ll[0][2] * rd[0][0] - ll[0][0] * rd[0][2]},
      {0,             -f1 * ll[0][2],      f1 * ll[0][1], 0,
      -g1 * ll[0][2], g1 * ll[0][1],ll[0][1] * rd[0][2] - ll[0][2] * rd[0][1]},
      {ll[1][2],       0,     -ll[1][0],    0,       0,        0,
       ll[1][2] * rd[1][0] - ll[1][0] * rd[1][2]},
      {0,             -ll[1][2],           ll[1][1],      0,
       0,              0,    ll[1][1] * rd[1][2] - ll[1][2] * rd[1][1]},
      {f3 * ll[2][2],  0,                 -f3 * ll[2][0], g3 * ll[2][2],
       0,   -g3 * ll[2][0],   ll[2][2] * rd[2][0] - ll[2][0] * rd[2][2]},
      {0,             -f3 * ll[2][2],      f3 * ll[2][1], 0,   -g3 * ll[2][2],
       g3 * ll[2][1],   ll[2][1] * rd[2][2] - ll[2][2] * rd[2][1]}};
    rref(mat, zz);
    rr1[0] = zz[0]; rr1[1] = zz[1]; rr1[2] = zz[2];
    vv1[0] = zz[3]; vv1[1] = zz[4]; vv1[2] = zz[5];

    // averaging
    vmadd(rr1, rr2, rvx, 3);                  //  rr1 = (3*rr2 + rr1) / 4
    smult(.25, rvx, rr1);
    vmadd(vv1, vv2, rvx, 3);                  //  vv1 = (3*rr2 + rr1) / 4
    smult(.25, rvx, vv1);

    vmadd(rr1, rr2, rvx, -1);                 //  rx = rr1 - rr2

  }while(norm(rvx) > 1.E-12);

  smult(xke, vv2, vv2);
  rv2el(rr2, vv2);

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

// write TLE to output file and to screen
void write_tle(char *file_out)
{
  Date t2(jd);           // global jd computed in main
  double tle = t2.tle;   // epoch day of elements
  char ec_string[9];
  char line1[70];
  char line2[70];

  sprintf(ec_string, "%.7f", eo);
  ec_string[0] = ec_string[2];
  ec_string[1] = ec_string[3];
  ec_string[2] = ec_string[4];
  ec_string[3] = ec_string[5];
  ec_string[4] = ec_string[6];
  ec_string[5] = ec_string[7];
  ec_string[6] = ec_string[8];

  sprintf(line1, "1 %05dU %-8s %014.8f 0.00000073  00000-0  50000-4 0    00",
           ssn, desig, tle);
  ccksum(line1);
  sprintf(line2, "2 %05d %8.4lf %8.4lf %.7s %8.4lf %8.4lf %11.8lf    00",
           ssn,
           xincl / de2ra,
           xnodeo / de2ra,
           ec_string,
           omegao / de2ra,
           xmo / de2ra,
           xno);
  ccksum(line2);
  FILE *fp;
  fp = fopen(file_out, "a");
  fprintf(fp, "\n%s\n", line1);
  fprintf(fp, "%s\n", line2);
  fclose(fp);
  printf("\n%s\n", line1);
  printf("%s\n", line2);
}

// load IOD lines from input file_in
void get_obs(char *file_in, char iod_line[][81])
{
  FILE *fp;
  if ((fp = fopen(file_in, "r")) == NULL)
  {
    printf("can't open %s\n", file_in);
    s_in("[exit]", buf);
    exit(0);
  }
  nobs = 0;

  while(fgets(iod_line[nobs], 80, fp))     // find observation entries
  {
    if(strlen(iod_line[nobs]) > 60 && strlen(iod_line[nobs]) != 70)
    {
      printf("(%2d) %s",nobs + 1, iod_line[nobs]);
      nobs++;
    }
  }
  fclose(fp);
  if(nobs < 2)
  {
    printf("Not enough observations in file %s \n", file_in);
    s_in("[exit]", buf);
    exit(0);
  }
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

// decodes the iod_line data
void read_obs(char iod_line[][81],
            double odata[][3],
            double ll[][3],
            double rd[][3])
{
  int year, month, day, hour, min, sign, format, age, epoch, obscode;
  double sec, ra, mm, ss, dc, dd, ds;
  double t, a, b, c, dcx;
  char sn[1], sbuf[5], dbuf[2];

  for(int i = 0; i < nobs; i++)
  {
    sscanf(iod_line[i], "%5d", &ssn);
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
      if(fabs(dc) > 1.4) dc = c / fabs(c) * acos(sqrt(a*a + b*b));
    }

    // precession from J2000
    t = (t1.jd - 2451545) / 36525;
    csi = (2306.2181 + .30188 * t + .017998 *t*t) * t * de2ra / 3600;
    zet = (2306.2181 + 1.09468 * t + .018203 *t*t) * t * de2ra / 3600;
    the = (2004.3109 - .42665 * t - .041833 *t*t) * t * de2ra / 3600;
    a = cos(dc) * sin(ra + csi);
    b = cos(the) * cos(dc) * cos(ra + csi)
        - sin(the) * sin(dc);
    c = sin(the) * cos(dc) * cos(ra + csi)
        + cos(the) * sin(dc);
    ra = atan(a / b);          // ra - zet
    if(b < 0) ra += pi;
    ra += zet;                 // right ascension, radians
    dc = asin(c);
    if(fabs(dc) > 1.4) dc = c / fabs(c) * acos(sqrt(a*a + b*b));

    // line-of-sight vectors
    ll[i][0] = cos(dc) * cos(ra);
    ll[i][1] = cos(dc) * sin(ra);
    ll[i][2] = sin(dc);

    odata[i][0] = t1.jd;                   // julian date
    odata[i][1] = ra;                      // ra radians (observed)
    odata[i][2] = dc;                      // dc radians (observed)

    read_site(iod_line[i]);
    Locate rre(t1.jd, la, lo, hh);
    rd[i] =  rre.rre;
  } // end for
}

