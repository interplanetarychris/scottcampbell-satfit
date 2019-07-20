// date.h
// by Scott Campbell campbel7@the-i.net

/****************************** Class Date ************************************
*                                                                             *
*                         Three Object Constructors                           *
*                                                                             *
*  Date t1;  Creates a date object, t1,  initialized to computer time at the  *
*            instant of the creation of the object.                           *
*                                                                             *
*  Date t1(time);   Creates a date object, t1, initialized to the time,       *
*                   which can be in either Julian date or TLE format.         *
*                                                                             *
*  Date t1(year, month, day, hour, min, sec);                                 *
*            Creates a date object, t1, initialized to the calendar date and  *
*            time passed by the six calendar variables                        *
*                                                                             *
*                                                                             *
*                         Three Member Functions                              *
*                                                                             *
*  t1.now();   Re-initializes an existing date object, t1, to the computer    *
*              time at the instant this command is executed by the computer.  *
*                                                                             *
*  t1.input();   Screen input of calendar date variables with current values  *
*                presented as default. Just press ENTER to accept current     *
*                value.  Useful to examine current values at run time.        *
*                                                                             *
*  t1.print();   Screen output of calendar date variables.                    *
*                                                                             *
*                                                                             *
*                         Ten Output Values                                   *
*                                                                             *
*  t1.thetag;  Sidereal time                                                  *
*  t1.jd;      Julian date                                                    *
*  t1.mjd;     Modified Julian Date                                           *
*  t1.tle;     Date in TLE format, only for years  2000 < year < 2100         *
*  t1.doy;     day-of-year                                                    *
*  t1.yy;      year                                                           *
*  t1.mm;      month                                                          *
*  t1.dd;      day of month, Greenwich                                        *
*  t1.hr;      hour, Greenwich                                                *
*  t1.mn;      minute                                                         *
*  t1.ss;      seconds                                                        *
*                                                                             *
*  all output values are immediately available after any initialization       *
*******************************************************************************/


/***************************** Class Locate ***********************************
*                                                                             *
*                        Two Object Constructors                              *
*                                                                             *
*  Locate rr(jd);   Creates a location object, rr, initialized to the default *
*                   coordinates, home(). Default values may be edited in the  *
*                   source code for home() at the very bottom of this         *
*                   program listing.                                          *
*                                                                             *
*  Locate rr(jd, latitude, longitude, height);                                *
*            Creates a location object, rr, initialized to the location at    *
*            lat(degrees, N > 0), long(degrees, E > 0), height(meters MSL),   *
*            at Julian date, jd.                                              *
*                                                                             *
*                         Two Member Functions                                *
*                                                                             *
*  rr.date(jd);   Moves existing location vector to the new location as a     *
*                 function of the Julian date.                                *
*                                                                             *
*  rr.input();   Screen input of location parameters with current values      *
*                presented as default. Just press ENTER to accept current     *
*                value.  Useful to examine current values at run time.        *
*                                                                             *
*                                                                             *
*                          One Output Value                                   *
*                                                                             *
*  rr.rre;  Topocentric location vector to observers position adjusted for    *
*           precession to epoch 2000                                          *
*                                                                             *
*                                                                             *
*   Note that the units of the input/output values are different from the     *
*   units as used internally by the class.                                    *
*   Input/Output Units:                 Internal Units:                       *
*     latitude,  degrees                  latitude,  radians                  *
*     longitude, degrees                  longitude, degrees                  *
*     height,    meters                   height,    earth radii              *
*                                                                             *
*******************************************************************************/

#include <cstdio>   // for printf, stdin
#include <cstdlib>  // for atof
#include <cmath>    // for fmod in sidereal()
#include <time.h>   // for now()
#include <string.h> // for strlen

using namespace std;

////////////// define classes /////////////////////////////////////////////////

class Date
{

private:

  void julian();           // calculate julian date
  void sidereal();         // calculate Greenwich sidereal hour angle
  void jcal();             // calendar date from julian date
  void calc_tle();         // calculate date in tle format and day-of-year
  void tle2cal();          // calendar date and doy from TLE date

  char *s_in(char *prompt, char *buffer)    // input utility
  {
    printf("%s", prompt);
    if (feof(stdin)) exit(0);
    fgets(buffer,sizeof(buffer),stdin);
    return buffer;
  }

  // returns remainder of division  var / div
  double mod(double var, double div)
  {
    return((var < 0) ? div + fmod(var, div) : fmod(var, div));
  }

public:
  double yy, mm, dd, hr, mn, ss;  // calendar date variables
  double thetag;                  // sidereal time in degrees
  double jd, mjd;                 // Julian date, modified julian date
  double tle;                     // date in TLE format
  int doy;                        // day of the year

  Date();                         // no argument constructor = computer time
  Date(double jdi);               // single argument constructor
  Date(double year,
       double month,
       double day,
       double hour,
       double min,
       double sec );              // six argument constructor
  void now();                     // initialize to GMT
  void input();                   // input calendar date values
  void print();                   // output calendar date values

}; // end class Date  /////////////////////////////////////////////////////////

  Date :: Date()                   // no argument constructor = computer time
  {
    now();
    julian();
    sidereal();
    calc_tle();
  }

  Date :: Date(double jdi) :  jd(jdi)     // single argument constructor
  {
    if(jd < 2400000)              // this date is in tle format
    {
      tle = jd;
      tle2cal();
      julian();
      sidereal();
    }
    else                          // this date is julian
    {
      jcal();
      sidereal();
      calc_tle();
      mjd = jd - 2400000.5;
    }
  }

  // six argument constructor
  Date :: Date(double year,
       double month,
       double day,
       double hour,
       double min,
       double sec ) :
       yy(year),
       mm(month),
       dd(day),
       hr(hour),
       mn(min),
       ss(sec)
  {
    julian();
    sidereal();
    jcal();
    calc_tle();
  }

  void Date :: julian()            // calculate julian date
  {
    double yr = yy;
    double mon = mm;
    if( mm < 3 )
    {
      yr = yy - 1;
      mon = mm + 12;
    }
    double c = (long)(yr / 100);
    jd = (long)(365.25 * (yr + 4716)) + (long)(30.6001 * (mon + 1))
         + dd + 2 - c + (long)(c / 4) - 1524.5;
    jd = jd + (hr + (mn + ss / 60) / 60) / 24;
    mjd = jd - 2400000.5;
  }

  void Date :: sidereal()          // calculate Greenwich sidereal hour angle
  {
    double t = (jd - 2451545) / 36525;
    thetag = 280.46061837 + 360.98564736629 * (jd - 2451545)
             + .000387933 * (t * t) - (t * t * t) / 38710000;
    thetag = mod(thetag, 360);     // degrees
  }

  void Date :: jcal()              // calendar date from julian date
  {
    long a, aa, b, c, d, e, z;
    double f;
    z = (long)(jd + .5);
    f =  jd + .5 - z;
    if(z < 2299161)
      a = z;
    else
    {
      aa = (long)((z - 1867216.25) / 36524.25);
      a = z + 1 + aa - aa / 4;
    }
    b = a + 1524;
    c = (long)((b - 122.1) / 365.25);
    d = (long)(365.25 * c);
    e = (long)((b - d) / 30.6001);
    if(e < 13.5)
      mm = e - 1;
    else
      mm = e - 13;
    if(mm > 2.5)
      yy = c - 4716;
    else
      yy = c - 4715;
    dd = b - d - (long)(30.6001 * e) + f;
    hr = dd - (long)dd;
    dd -= hr;
    hr *= 24;
    mn = hr - (long)hr;
    hr -= mn;
    mn *= 60;
    ss = mn - (long)mn;
    mn -= ss;
    ss *= 60;
    if(ss >= 59.9999)
    {
       mn += 1;
       ss -= ss;
    }
    if(mn >= 60)
    {
       hr += 1;
       mn -= mn;
    }
  }

  void Date :: calc_tle()      // calculate date in tle format and day-of-year
  {
    int k = ( ((int)yy%4 == 0 && (int)yy%100 != 0) || (int)yy%400 == 0) ? 1 : 2;
    doy = (int)(275 * mm / 9) - k * (int)(mm / 12 + .75) + (int)dd - 30;
    tle = doy + hr / 24 + mn / 1440 + ss / 86400;
    if(yy >= 2000) tle += 1000 * (yy - 2000);
    else tle += 1000 * (yy - 1900);
  }

  void Date :: tle2cal()      //  calendar date and doy from TLE date
  {
    yy = (int)(tle * .001);
    double dx = tle - (int)tle;             // place holder for decimal day
    doy = (int)((int)tle - yy * 1000);
    yy += 2000;                             // update 01/01/2100
    int k = ( ((int)yy%4 == 0 && (int)yy%100 != 0) || (int)yy%400 == 0) ? 1 : 2;
    mm = (int)(9 * (k + doy) / 275 + .98) + 1;
    if(doy < 32)
      mm = 1;
    dd = doy - (int)(275 * mm / 9) + k * (int)(mm / 12 + .75) + 30;
    hr = (int)(dx * 24);
    dx = dx * 24 - hr;                     // decimal hours
    mn = (int)(dx * 60);
    dx = dx * 60 - mn;                     // decimal minutes
    ss = dx * 60;                          // seconds and decimal seconds
  }

  void Date :: now()          // initialize to GMT
  {
    time_t rawtime;
    tm * ptm;
    time ( &rawtime );
    ptm = gmtime ( &rawtime );
    yy = ptm->tm_year + 1900;
    mm = ptm->tm_mon + 1;
    dd = ptm->tm_mday;
    hr = ptm->tm_hour;
    mn = ptm->tm_min;
    ss = ptm->tm_sec;
    julian();
    sidereal();
    calc_tle();
  }

  void Date :: input()        // input calendar date values
  {
    char buf[20];
    printf("Year   [%.0f", yy);
    if (strlen(s_in((char *)"]: ", buf)))
      yy = atof(buf);
    printf("Month  [%.0f", mm);
    if (strlen(s_in((char *)"]: ", buf)))
      mm = atof(buf);
    printf("Day    [%.0f", dd);
    if (strlen(s_in((char *)"]: ", buf)))
      dd = atof(buf);
    printf("Hour   [%.0f", hr);
    if (strlen(s_in((char *)"]: ", buf)))
      hr = atof(buf);
    printf("Minute [%.0f", mn);
    if (strlen(s_in((char *)"]: ", buf)))
      mn = atof(buf);
    printf("Second [%.0f", ss);
    if (strlen(s_in((char *)"]: ", buf)))
      ss = atof(buf);
    julian();
    sidereal();
    calc_tle();
  }

  void Date :: print()        // output calendar date values
  {
  printf("Year   %.0f\n", yy);
  printf("Month  %.0f\n", mm);
  printf("Day    %.0f\n", dd);
  printf("Hour   %.0f\n", hr);
  printf("Minute %.0f\n", mn);
  printf("Second %.2f\n", ss);
  }

//////////////////// Class Locate /////////////////////////////////////////////
class Locate
{
private:

  static const double de2ra,
                     xkmper;   // equatorial earth radius, km

  double la,   // latitude, degrees
         lo,   // east longitude (west is negative), degrees
         hh,   // height MSL, meters
         jd;   // Julian date

  void rgeo();                 // vector to topocentric position
  void convert();
  void home();

  char *s_in(char *prompt, char *buffer)
  {
    printf("%s", prompt);
    if (feof(stdin)) exit(0);
    fgets(buffer,sizeof(buffer),stdin);
    return buffer;
  }

  // returns remainder of division  var / div
  double mod(double var, double div)
  {
    return((var < 0) ? div + fmod(var, div) : fmod(var, div));
  }

public:
  double rre[3];

// single argument constructor takes julian date and uses default location
  Locate(double jdi);
// 4 - argument constructor to input date and location parameters
 Locate(double jdi, double lat, double lon, double hite);
// input new location variables on existing object
  void input();
  // location as a function of date
  void date(double jdi);

};  // end class Locate ///////////////////////////////////////////////////////
const double Locate::de2ra(.0174532925199433);
const double Locate::xkmper(6378.1363);   // equatorial earth radius, km

// single argument constructor takes julian date and uses default location
  Locate :: Locate(double jdi) :  jd(jdi)
  {
     home();                 // default location parameters
     rgeo();
  }

// 4 - argument constructor to input date and location parameters
 Locate :: Locate(double jdi, double lat, double lon, double hite) :
                      jd(jdi),    la(lat),    lo(lon),    hh(hite)
 {
    convert();
    rgeo();
 }

  // input new location variables on existing object
  void Locate :: input()
  {                         // i.e. change location of the object
    char buf[20];
    printf("\n\n");
    printf("Latitude ( + N ) [%6.2f", la / de2ra);
    if (strlen(s_in((char *)"]: ", buf)))
    {
      la = atof(buf);
      la *= de2ra;              // convert to radians
    }
    printf("E longitude      [%6.2f", lo);
    if (strlen(s_in((char *)"]: ", buf)))
    {
      lo = atof(buf);
      lo = mod(lo, 360);       // convert w long(-) to e long(+)
    }
    // swap printf if units in feet
//    printf("Altitude MSL(ft) [%6.2f", hh * xkmper / .0003048);
    printf("Altitude MSL(m)  [%6.2f", hh * xkmper * 1000);
    if (strlen(s_in((char *)"]: ", buf)))
    {
      hh = atof(buf);
//      hh *= .3048 ;             // feet to meters <- comment out for meters
      hh *= .001;               // meters to kilometers
      hh /= xkmper;             // kilometers to earth radii
    }
    rgeo();
  }

  // location as a function of date
  void Locate :: date(double jdi)
  {
    jd = jdi;
    rgeo();
  }

  void Locate :: rgeo()        // vector to topocentric position
  {
    double ff, f, gc, gs, t;
    double theta0;

    Date t1(jd);               // calculate sidereal time
    theta0 = mod(t1.thetag + lo, 360) * de2ra ;
    ff = 1 / 298.2572;
    f = sqrt(1 - ff * (2 - ff) * sin(la)*sin(la));
    gc = 1 / f + hh;
    gs = (1 - ff)*(1 - ff) / f + hh;

    rre[0] = gc * cos(theta0) * cos(la);
    rre[1] = gc * sin(theta0) * cos(la);
    rre[2] = gs * sin(la);

  } // end rgeo

void Locate :: convert()
{
//   hh *= .3048 ;             // feet to meters <- comment out for meters
   hh *= .001;                 // meters to kilometers
   hh /= xkmper;               // kilometers to earth radii
   la *= de2ra;                // convert to radians
   lo = mod(lo, 360);          // convert w long(-) to e long(+)
}

void Locate :: home()        // default location values, Beeville TX
{
   la = 28.4861;             // north latitude is positive
   lo = 262.1806;            // west longitude is positive
   hh = 107;                 // elevation MSL in meters (see convert)
   convert();
}