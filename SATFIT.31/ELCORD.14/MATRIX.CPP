/* MATRIX.C
by Paul S. Hirose 1990 Sep 30
mathematical routines for ELCOR program
*/

#include "elcor.h"   /* global header */

/*########################### DATA ###########################*/

static double dscale[NUMEL];   /* scaling factor for derivatives */

/*###################### LOCAL FUNCTIONS ######################*/

/* dot product */
static double dot(double *fltp1, double *fltp2);
static void invert();  /* matrix inversion */
static void pvecs(struct vector *vecp);   /* print vectors */
static void sub(double *dp1, double *dp2, int u, double scale);     /* subtract matrix rows */
static void swap(double *dp1, double *dp2);    /* swap matrix rows */
static void xyzres(struct vector *res);  /* compute & print residuals */

/*########################### CODE ###########################*/

/* Compute & store (in row r of a p x n matrix) the weighted partial
derivatives of each observed quantity with respect to the orbital element
corresponding to tomod[r].  We will store in ATW, which is a p x (n+1) matrix,
in such fashion that column 0 of ATW is unused. */

void deriv(int r)
{
   double delta, orig;
   double *ATWp, *fltp, *bp;
   unsigned u, u3;
   struct elemen *elemp;

   ATWp = ATW + r * (n + 1) + 1;      /* point to active row */
   elemp = elarr + tomod[r];   /* element we're working on */

   /* request a delta value, add it to the element */
   orig = *elemp->elp;      /* save value of element */
   /*** no request!!! printf("delta ");
   el_in(99, elemp, &delel[r], &iflag); ***/
   delta = delel[r];
   *elemp->elp = orig + delta;

   /* place (old residuals) - (new residuals) in active row */
   xyzres((struct vector *) ATWp);      /* new residuals */
   fltp = ATWp;         /* point to active row */
   bp = (double *)b;      /* point to old residuals */
   for (u = 0; u < n; u++) {
      *fltp = *bp++ - *fltp;   /* old resid - new resid */
      fltp++;
   }
   /*** pvecs(ATWp); /* print delta residuals */

   /* form scaled derivatives */
   bp = W;      /* point to array of weights */
   u3 = 3;
   for (u = 0; u < n; u++) {
      *ATWp *= *bp * 100.;
      ++ATWp;
      if (--u3 == 0) {
         ++bp;      /* point to next weight */
         u3 = 3;
      }
   }
   dscale[r] = 100. * delta;   /* scaling factor */
   *elemp->elp = orig;   /* restore old value of element */
}

/* Returns dot product of fltp1[] and fltp2[], each of which is n long. */
static double dot(double *fltp1, double *fltp2)
{
   double sum;
   unsigned u;

   sum = 0.;
   for (u = 0; u < n; u++)
      sum += *fltp1++ * *fltp2++;
   return sum;
}

/* Invert a p x p matrix.  On entry the matrix to invert must be in the left
half of PX2P.  On exit, its inverse is in the right half. */
static void invert(void)
{
   double temp, pivot, *coljp, *rowip, *rowjp, *diagp;
   unsigned i, j, twop;

   if (p == 1) {      /*  invert 1 x 1 matrix */
      PX2P[1] = 1. / *PX2P;
      return;
   }

   twop = p << 1;      /* 2 * p */

   /* Put an identity matrix into right half of PX2P.  For this section,
   consider PX2P[0][0] bottom right corner. */
   coljp = PX2P + p;      /* &PX2P[p-1][p-1] */
   for (i = p; i--; ) {   /* do rows p-1 to 0 */
      for (j = p; j--; ) {   /* do columns p-1 to 0 */
         if (j == i)
            *coljp = 1.;
         else
            *coljp = 0.;
         ++coljp;   /* next column */
      }
      coljp += p;      /* next row */
   }
   /* Convert left half of PX2P to lower triangular.  Revert to
   conventional notation of PX2P[0][0] as top left. */
   rowjp = PX2P + twop * p;      /* &PX2P[p][0] */
   diagp = rowjp + p;         /* &PX2P[p][p] */
   for (j = p; --j; ) {   /* col p-1 to 1 */
      rowjp -= twop;      /* &PX2P[j][0] */
      diagp -= twop + 1;   /* &PX2P[j][j] */
      /* find biggest absolute value in column j row 0 to j */
      pivot = 0.;      /* biggest number found so far */
      coljp = diagp;      /* &PX2P[j][j] */
      rowip = rowjp;
      for (i = j + 1; i--; ) {   /* row j to 0 */
         temp = fabs(*coljp);
         if (temp > pivot) {
            rowip = coljp - j;   /* &PX2P[i][0] */
            pivot = temp;
         } coljp -= twop;   /* move up one row */
      }
      /* Swap the row with row j, if it isn't already there */
      if (rowip != rowjp)
         swap(rowjp, rowip);

      /* zero all entries in col j row 0 to j - 1 */
      rowip = rowjp;      /* &PX2P[j][0] */
      coljp = diagp;      /* &PX2P[j][j] */
      pivot = *diagp;      /* PX2P[j][j] */
      for (i = j; i--; ) {   /* row j-1 to 0 */
         rowip -= twop;      /* &PX2P[i][0] */
         coljp -= twop;      /* &PX2P[i][j] */
         temp = *coljp / pivot;   /* scaling factor */
         sub(rowip, rowjp, twop, temp);
      }
   }

   /* Left half of PX2P is lower triangular.  Use row operations to
   convert it to a diagonal matrix.  To save time, we will only perform
   the row additions on the right half of the matrix. */

   rowjp = PX2P + p;   /* &PX2P[0][p] ( the pivot row) */
   diagp = PX2P;      /* &PX2P[0][0] */
   for (j = 0; j != p - 1; ++j) {   /* column 0 to p-2 */
      pivot = *diagp;
      rowip = rowjp;      /* &PX2P[j][p] */
      coljp = diagp;      /* &PX2P[j][j] */
      for (i = j + 1; i != p; ++i) {      /* row j+1 to p-1 */
         rowip += twop;         /* &PX2P[i][p] */
         coljp += twop;         /* &PX2P[i][j] */
         temp = *coljp / pivot;
         sub(rowip, rowjp, p, temp);
      }
      diagp += twop + 1;
      rowjp += twop;
   }

   /* Left half of PX2P is diagonal.  Divide each row in the right half
   of PX2P by its diagonal entry in the left half. */

   diagp = PX2P;      /* &PX2P[0][0] */
   rowip = PX2P + p;   /* &PX2P[0][p] */
   for (i = 0; i != p; ++i) {   /* row 0 to p-1 */
      temp = 1. / *diagp;
      coljp = rowip;      /* &PX2P[i][p] */
      for (j = p; j--; ) {   /* do p columns */
         *coljp *= temp;
         ++coljp;   /* next column */
      }
      diagp += twop + 1;   /* next diagonal entry */
      rowip += twop;
   }
}

/* Perform least squares adjustment of the orbital elements */
void lsqr(void)
{
   double temp, *PX2Pp;
   double *ATWp, *fltp, *Wp;
   unsigned col, row, u, u3;

   /* Form the left half of PX2P */
   for (col = 0; col != p; ++col) {
      /* form AT */
      ATWp = ATW + (n + 1) * col + 1;
      fltp = AT;
      Wp = W;
      u3 = 0;
      for (u = 0; u < n; u++) {
         *fltp++ = *ATWp++ / *Wp;   /* remove weight */
         if (++u3 == 3) {
            u3 = 0;
            Wp++;
         }
      }

      /* form a column in PX2P */
      ATWp = ATW + 1;
      PX2Pp = PX2P + col;
      for (row = 0; row < p; row++) {
         *PX2Pp = dot(ATWp, AT);
         ATWp += n + 1;      /* down one row */
         PX2Pp += p + p;      /* down one row */
      }
   }

   invert();   /* invert PX2P */

   /* compute (right half of PX2P) * ATW */
   for (col = 0; col != n; col++) {
      /* form a column */
      PX2Pp = PX2P + p;      /* row 0 col p */
      fltp = ATW + col;      /* destination column */
      for (row = 0; row != p; row++) {
         /* form a dot product */
         ATWp = ATW + col + 1;   /* top of column */
         temp = 0.;
         for (u = 0; u < p; u++) {   /* p times */
            /* accumulate product in temp */
            temp += *PX2Pp++ * *ATWp;
            ATWp += n + 1;   /* down a row */
         }
         *fltp = temp;      /* store dot product */
         PX2Pp += p;      /* next row */
         fltp += n + 1;      /* down a row */
      }
   }

   /* make adjustments to the elements */
   ATWp = ATW;
   for (u = 0; u != p; u++) {
      *elarr[tomod[u]].elp += dot(ATWp, (double *)b) * dscale[u];
      ATWp += n + 1;      /* next row */
   }
}

/* Print the magnitudes of an array of vectors of length nobs.  Each is
   labeled with the corresponding observation label. */
static void pvecs(struct vector *vecp)
{
   int u, line;

   line = 0;
   for (u = 0; u < nobs; u++) {
      printf("%-25s %7.3f\n", label[u],
        sqrt(vecp->x * vecp->x + vecp->y * vecp->y +
             vecp->z * vecp->z) / de2ra);
      vecp++;
   }
}

/* Store and print residuals and RMS residuals */
void prnres(void)
{
   double temp;
   double *floatp, *Wp;
   unsigned u, u3;

   // printf("\nResiduals:\n");
   xyzres(b);   /* compute & store residuals */
   if(out) pvecs(b);  /* print residuals in degrees */

   /* compute the weighted RMS of residuals */
   oldrms = newrms;
   newrms = 0.;
   floatp = (double *)b;   /* array of residuals */
   Wp = W;         /* array of weights */
   u3 = 0;
   for (u = 0; u < n; u++) {
      temp = *floatp++ / de2ra;
      newrms += temp * temp;  /*   * *Wp; */
      if (++u3 == 3) {
         u3 = 0;
         Wp++;
      }
   }
   newrms = sqrt(newrms / nobs);
}

/* Matrix row subtraction.  dp1[] = dp1[] - dp2[] * scale */
static void sub(double *dp1, double *dp2, int u, double scale)
{
   do {
      *dp1 -= *dp2++ * scale;
      ++dp1;
   } while (--u);
}

/* Swap 2 arrays of double of length 2*p. */
static void swap(double *dp1, double *dp2)
{
   double temp;
   unsigned u = p << 1;

   do {
      temp = *dp1;
      *dp1++ = *dp2;
      *dp2++ = temp;
   } while (--u);
}

/* Computes residuals (observed position - predicted position) for all
observations, stores residuals in res[]. */
static void xyzres(struct vector *res)
{
   double epoch, jdi;
   double x, y, z, r;
   int u;

   xmo = qmo - omegao;   /* recover mean anomaly */

   /* epoch of the elset (MJD) */
   epoch = julday(epoch_year, 1, 0) + epoch_day;
   jdi = epoch + 2450000;

   for (u = 0; u < nobs; u++) {
      sxp4(jdi - 0.5, (tobs[u] - epoch) * xmnpda);   /* predict position of sat */

      /* topocentric coordinates (sat xyz - observer xyz) */
      x = sat.x - place[u].x;
      y = sat.y - place[u].y;
      z = sat.z - place[u].z;

      r = sqrt(x*x + y*y + z*z);   /* range */

      /* convert to unit vector, subtract from observed vector */
      res->x = obs[u].x - x / r;
      res->y = obs[u].y - y / r;
      res->z = obs[u].z - z / r;

      res++;
   }
}
