
#include "elcor.h"   /* global header */

/* Read file of elements and observations,
   set globals nobs and n accordingly */

static int error;
static long getlong(char *buf, int start, int stop);
static float getfloat(char *buf, int start, int stop, int numint);

static double f = 3.35278e-3;   /* flattening of earth */

void getobs(char *filename) /* Input file name (argv[1]) */
{
    FILE *fp;

    if ((fp = fopen(filename, "r")) == NULL) {
        printf("can't open %s\n", filename);
        s_in("[exit]", buf);
        exit(2);
    }

    /* load elements */
    switch (getel(fp)) {
    case -1:
        printf("Incorrect format for elements\n");
        exit(0);
    case 0:
        /* printf("Elements loaded\n"); */
        printf("\n");
        break;
    case 1:
        printf("caution - deep space\n");
    }

    /* Read and process the observations */
    nobs = 0;
    max_obs_day = 0.0;

    /* Loop until all obs processed */
    while (1) {
        char inp_str[81];
        int site, year, month, day, hour, minute;
        float second;
        int type, obs_epoch;
        long ncatin;
        char desig[9];
        int rahr;
        float ramin, rasec;
        float raminfr;
        float decdeg, decmin, decsec, dec;
        float wgt;
        int i;
        double sinlat, coslat, a, b, c, s, h, t, otime, lhaa;
        double ra, dc, temp1, temp2;
        double obs_day_time;

        if (fgets(inp_str, 80, fp) == NULL) break;

        wgt = 1.0;   /* default if field is blank or alpha */

        /* Reset the "error" flag */
        error = 0;

/*0         1         2         3         4         5         6     */
/*01234567890123456789012345678901234567890123456789012345678901234 */
/*25017 97 064A   5918 G 19990806205132250 17 35 1726600+347700 56 S+025 04 */

        /* Get the site number */
        site = getlong(inp_str, 16, 19);

        if (error != 0) continue;
        if (site == 0) continue;

        /* Codes 1,2,3,7 decoded below */
        type = inp_str[44] - '0';
        if (type < 1 || (type > 3 && type != 7)) continue;

        /* Code 4 is epoch 1950, 5 is epoch 2000 */
        if (inp_str[45] == '4')
            obs_epoch = 4;
        else if (inp_str[45] == '5')
            obs_epoch = 5;
        else
            continue;

        /* Get the catalog number */
        ncatin = getlong(inp_str, 0, 4);
        if (ncatin == 0) continue;

        /* Catalog number must match the current elset! */
        if (ncatin != sat_ncat) continue;

        /* Get the designation */
        memcpy(desig, &inp_str[6], 8);
        desig[8] = '\0';

/*0         1         2         3         4         5         6         7 */
/*01234567890123456789012345678901234567890123456789012345678901234567890123 */
/*25017 97 064A   5919 G 19990806205132250 17 35 1726600+347700 56 S+025 04 */
/*24680 96 072A   5919 G 19990806213419500 56 35 1418600+157800 56 S+060 10 */
/*22251 92 083A   5918 P 19991007192000800 57 35 1708196+104900 56 S+060 10 */
/*24680 96 072A   9999 G 20010522205854380 76 25 1322879-005970 28 S */
/*24680 96 072A   2453 B 20010519094216470 57 75 0951105+242870 57 S+040 05 */

        /* Get the year, month, day, hour, minute, second */
        year = getlong(inp_str, 23, 26);
        month = getlong(inp_str, 27, 28);
        day = getlong(inp_str, 29, 30);
        hour = getlong(inp_str, 31, 32);
        minute = getlong(inp_str, 33, 34);
        second = getfloat(inp_str, 35, 39, 2);
        if (error != 0) continue;

        /* Handle the two different types of RA */
        switch (type) {

            case 1:
                /* Format 1: RA/DEC = HHMMSSs+DDMMSS */
            case 7:
                /* Format 7: RA/DEC = HHMMSSs+DDdddd */
                rahr = getlong(inp_str, 47, 48);
                ramin = getlong(inp_str, 49, 50);
                rasec = getfloat(inp_str, 51, 53, 2);
                ramin = ramin + rasec / 60.0;
                break;

            case 2:
                /* Format 2: RA/DEC = HHMMmmm+DDMMmm */
            case 3:
                /* Format 3: RA/DEC = HHMMmmm+DDdddd */
                rahr = getlong(inp_str, 47, 48);
                ramin = getfloat(inp_str, 49, 53, 2);
                break;

        }

        /* Handle the three different types of Dec */
        switch (type) {

            case 1:
                /* Format 1: RA/DEC = HHMMSSs+DDMMSS */
                decdeg = getlong(inp_str, 55, 56);
                decmin = getlong(inp_str, 57, 58);
                decsec = getlong(inp_str, 59, 60);
                sscanf(inp_str + 59, "%2s", &buf);
                decsec = atof(buf) * pow(10., 2 - strlen(buf));
                dec = decdeg + (decmin + decsec / 60.0) / 60.0;
                break;

            case 2:
                /* Format 2: RA/DEC = HHMMmmm+DDMMmm */
                decdeg = getlong(inp_str, 55, 56);
                sscanf(inp_str + 57, "%4s", &buf);
                decmin = atof(buf) / pow(10., strlen(buf) - 2);
                dec = decdeg + decmin / 60.0;
                break;

            case 3:
                /* Format 3: RA/DEC = HHMMmmm+DDdddd */
            case 7:
                /* Format 7: RA/DEC = HHMMSSs+DDdddd */
                dec = getfloat(inp_str, 55, 60, 2);
                break;
        }

        if (error != 0) continue;

        /* Check for negative sign for dec */
        if (inp_str[54] == '-')
            dec = -dec;

        /* printf("%s", inp_str); */

        /* Lookup this site */
        for (i = 0; i < num_sites; i++) {
           if (sitenum[i] == site) goto found;
        }
        printf("obs file site %d not in sites file\n", site);
        continue;

found:

        /* printf("%d %s %d %d %d %d:%d:%.1f %dHr %.2fMn %.2fdec %d %ld %s\n",
           site, siteabbr[i], year, month, day, hour, minute, second,
           rahr, ramin, dec, obs_epoch, ncatin, desig);   */

        if (nobs >= MAXOBS) {
            printf("Can't handle more than %d observations\n", MAXOBS);
            exit(2);
        }

        /* Create the "label" */
        sprintf(label[nobs], "(%2d) %s %2d %2d %02d:%02d:%02d %2dH %4.1fM%5.1f",
           nobs + 1, siteabbr[i], month, day, hour, minute, (int)second,
           rahr, ramin, dec);

        temp1 = xlat[i] * de2ra;    /* latitude */
        sinlat = sin(temp1);        /* sin latitude */
        coslat = cos(temp1);        /* cos latitude */
        h = xhgt[i] * 1.568e-7;     /* meters -> earth radii */

        /* Preliminary quantities involving observer's position with respect
        to geocenter (formulas from Baker & Makemson) */
        c = 1. / sqrt(1. - (f + f - f * f) * sinlat * sinlat);
        s = c * (1. - f) * (1. - f);
        /* distance from polar axis to observer */
        temp1 = (c + h) * coslat;
        /* north distance from equatorial plane to observer */
        temp2 = (s + h) * sinlat;

        /* Compute and store the tobs of the observation, then obtain local
        hour angle of Aries. */
        otime = (double)julday(year, month, day) +
                (hour + (minute + second / 60.0) / 60.0) / 24.0;
        tobs[nobs] = otime;
        /* Note long west is negative */
        lhaa = thetag(otime) + xlong[i] * de2ra;

        /* location of observer at time of observation */
        place[nobs].x = temp1 * cos(lhaa);
        place[nobs].y = temp1 * sin(lhaa);
        place[nobs].z = temp2;

        /* observed position of satellite */
        dc = dec * de2ra;      /* declination */
        ra = (rahr + ramin / 60.0) / 24. * twopi;   /* Right Asc. */

        /* Precess observation from obs_epoch to mean equator
        at epoch of elements */
        double  csi = .0055878713278878,
                zet = .0055888307019922,
                the = .0048580335354883;

        if(obs_epoch == 4)   // precess from B1950 to J2000
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
        t = (otime - 1545) / 36525;
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

        /* line-of-sight vectors */
        obs[nobs].x = cos(dc) * cos(ra);
        obs[nobs].y = cos(dc) * sin(ra);
        obs[nobs].z = sin(dc);

        /* store square of weight */
        W[nobs] = wgt * wgt;

        if (max_obs_day < otime)
           max_obs_day = otime;

        /* Increment count of obs */
        nobs++;
    }

    fclose(fp);
    printf("%d observations found\n", nobs);
    n = nobs * 3;   /* each obs yields an x, y, and z */
}

static long getlong(char *buf, int start, int stop)
{
    long value;
    int length, i;

    value = 0;
    length = stop - start + 1;

    for (i = 0; i < length; i++) {
        int ich;

        ich = buf[start + i];
        if (ich < '0' || ich > '9') {
            error = 1;
            return(-1); /* error */
        }
        value = value * 10 + ich - '0';
    }
    return(value);
}

static float getfloat(char *buf, int start, int stop, int numint)
{
    float value, factor;
    int length, i;

    value = 0.0;
    factor = 0.1;
    length = stop - start + 1;

    for (i = 0; i < length; i++) {
        int ich;

        ich = buf[start + i];

        /* Ignore trailing blanks */
        if (i >= numint && ich == ' ') break;

        if (ich < '0' || ich > '9') {
            error = 1;
            return(-1.0); /* error */
        }

        if (i < numint)
            value = value * 10.0 + ich - '0';
        else {
            value = value + ((float)(ich - '0') * factor);
            factor *= 0.1; /* decrease scale factor */
        }
    }
    return(value);
}
