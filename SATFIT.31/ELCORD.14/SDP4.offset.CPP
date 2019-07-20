/* SDP4.C */


#include "elcor.h"   /* the global header */

static double c2, a3ovk2 = -xj3 / ck2, two_thirds = 2./3.;

void sgp4(double tsince)
{
   // global tle variables
   // jd, xndt2o, xndd6o, bstar;
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
  xl = xmp+omega+xnode+ xnodp*templ;
  if( tempa < 0.)       /* force negative a,  to indicate error condition */
     a = -a;

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

  if( a < 0.)
     rval = -2;
  if( elsq > 1. - chicken_factor_on_eccentricity)
     rval = -1;
  if( rval)
     exit( 0);
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

  double rdot = xke*sqrt(a)*esinE/r;
  double rfdot = xke*sqrt(pl)/r;
  double xn = xke/(a * sqrt(a));
  double rdotk = rdot-xn*temp1*x1mth2*sin2u;
  double rfdotk = rfdot+xn*temp1*(x1mth2*cos2u+1.5*x3thm1);
  double vx = xmx*cosuk-cosnok*sinuk;
  double vy = xmy*cosuk-sinnok*sinuk;
  double vz = sinik*cosuk;

   /* position */

   sat.x = rk * ux;
   sat.y = rk * uy;
   sat.z = rk * uz;

   satdot.x = rdotk * ux + rfdotk * vx;
   satdot.y = rdotk * uy + rfdotk * vy;
   satdot.z = rdotk * uz + rfdotk * vz;

}


void sdp4(double jd, double tsince)
{
  // global tle variables
  // jd, xndt2o, xndd6o, bstar;
  // xincl, xnodeo, eo, omegao, xmo, xno;

  double ao, xnodp, delo, a1, del1, r1, temp, x3thm1, tval;

  // define in SDP4 and COMMON
  double c1, c2, c4, xnodcf, t2cof ;

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
  se2, se3, sgh2, sgh3, sgh4, sh2, sh3, si2, si3, sl2, sl3,
  sl4, sse, ssg, ssh, ssi, ssl, thgr, xfact, xgh2, xgh3, xgh4, xh2,
  xh3, xi2, xi3, xl2, xl3, xl4, xlamo, xli, xni, xnq, xqncl, zcosgl,
  zcoshl, zcosil, zmol, zmos, zsingl, zsinhl, zsinil;

  int i, resonance_flag, synchronous_flag;

  /* sxpx_common_init( params, tle, &init, deep_arg); */
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
  /* end sxpall_common_init */

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
        x3thm1+0.0625*temp2*betao*(13-78*theta2+137*theta4);
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
  xqncl = xincl;
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
     if (xqncl < pi / 60.)      /* pi / 60 radians = 3 degrees */
         sh = 0;
     else
         sh = -zn*s2*(z21+z23);
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
        ssh = (sinio ? sh / sinio : 0.);
        ssg = sgh-cosio*ssh;
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
        zcosh = zcoshl*cosq+zsinhl*sinq;
        zsinh = sinq*zcoshl-cosq*zsinhl;
        zn = znl;
        cc = c1l;
        ze = zel;
     }
  } /* end iteration */

  sse += se;
  ssi += si;
  ssl += sl;
  ssg += sgh;
  if(sinio)
  {
     ssg -= sh * cosio / sinio;
     ssh += sh / sinio;
  }

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
     double f321 = 1.875 * sinio * (1 - 2 * cosio - 3 * theta2);
     double f322 = -1.875*sinio*(1+2*cosio-3*theta2);
     double f441 = 35 * sini2 * f220;
     double f442 = 39.3750 * sini2 * sini2;
     double f522 = 9.84375*sinio*(sini2*(1-2*cosio-5*
                 theta2)+0.33333333*(-2+4*cosio+
                 6*theta2));
     double f523 = sinio*(4.92187512*sini2*(-2-4*
                 cosio+10*theta2)+6.56250012
                 *(1+2*cosio-3*theta2));
     double f542 = 29.53125*sinio*(2-8*
                 cosio+theta2*(-12+8*cosio+10*theta2));
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
     bfact = xmdot + xnodot + xnodot - thdt - thdt;
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
     int k = 0;

     if( delt > dpsec_integration_step)
        delt = dpsec_integration_step;
     else if( delt < -dpsec_integration_step)
        delt = -dpsec_integration_step;
     else
        final_integration_step = 1;

     /* compute_dpsec_derivs  */

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
           derivs[k] = term1a + term2a + term3a;
           k++;
           derivs[k] = term1b + term2b + term3b;
           k++;
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
           derivs[k] = term1a + term2a;
           k++;
           derivs[k] = term1b + term2b;
           k++;
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

  /* Deep_dpper( tle, deep_arg); */
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

  xinc += pinc;
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
     double betdp = -ph * sinok
                  + ( pinc * cosis + sinis) * cosok;
     double dls, delta_xnode;

     delta_xnode = atan2(alfdp,betdp) - xnode;

      /* This is a patch to Lyddane modification suggested */
      /* by Rob Matson, streamlined very slightly by BJG, to */
      /* keep 'delta_xnode' between +/- 180 degrees: */

     if( delta_xnode < -pi)
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

  double rdot = xke*sqrt(a)*esinE/r;
  double rfdot = xke*sqrt(pl)/r;
  xn = xke/(a * sqrt(a));
  double rdotk = rdot-xn*temp1*x1mth2*sin2u;
  double rfdotk = rfdot+xn*temp1*(x1mth2*cos2u+1.5*x3thm1);
  double vx = xmx*cosuk-cosnok*sinuk;
  double vy = xmy*cosuk-sinnok*sinuk;
  double vz = sinik*cosuk;

   /* position / velocity */
   sat.x  =  rk * ux;          /* er  */
   sat.y  =  rk * uy;
   sat.z  =  rk * uz;
   satdot.x  =  (rdotk * ux + rfdotk * vx);        /* er / min  */
   satdot.y  =  (rdotk * uy + rfdotk * vy);
   satdot.z  =  (rdotk * uz + rfdotk * vz);

}   /* end sdp4 */

void sxp4(double jd, double tsince)
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
  if (2 * 3.14159265358979323846 / (xnodp * 1440.) >= (1. / 6.4))
     sdp4(jd, tsince);    /* yes,  it should be a deep-space (SDP4) ephemeris */
  else
     sgp4(tsince);

}

