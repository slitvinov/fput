#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <math.h>
#include <stdio.h>

enum { n = 32 };
#define PI (3.141592653589793)

static int func(double t, const double x[], double f[], void *params) {
  int i;
  double xm;
  double x0;
  double xp;
  double dp;
  double dm;
  double a = *(double *)params;

  (void)(t);
  for (i = 0; i < n; i++)
    f[i] = x[i + n];

  for (i = 0; i < n; i++) {
    xm = i - 1 >= 0 ? x[i - 1] : 0;
    x0 = x[i];
    xp = i + 1 < n ? x[i + 1] : 0;
    dp = xp - x0;
    dm = x0 - xm;
    f[i + n] = xp + xm - 2 * x0 + a * (dp * dp - dm * dm);
  }

  return GSL_SUCCESS;
}

int main(void) {

  double a;
  double a0;
  double a1;
  double E;
  double si;
  double t;
  double T;
  double w;
  double x[2 * n];
  int i;
  int k;
  int l;
  int s;
  int ns;

  a = 1.0 / 4;
  for (i = 1; i <= n; i++) {
    x[i - 1] = sin(i * PI / n);
    x[i + n - 1] = 0;
  }

  t = 0;
  w = 2 * sqrt(a) * sin(PI / (2 * (n + 1)));
  T = 2 * PI / w;
  ns = 10000;
  gsl_odeiv2_system sys = {func, NULL, 2 * n, &a};
  gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_y_new(
      &sys, gsl_odeiv2_step_rk4, T / ns, 1e-8, 1e-8);
  for (;;) {
    printf("%.16e ", t);
    for (k = 1; k <= n; k++) {
      a0 = 0;
      a1 = 0;
      for (l = 1; l <= n; l++) {
        si = sin(k * l * PI / n);
        a0 += x[l - 1] * si;
        a1 += x[l + n - 1] * si;
      }
      si = sin(PI * k / (2 * n));
      E = a1 * a1 / 2 + 8 * a * a0 * a0 * si * si;
      printf(" %.16e", E);
    }
    printf("\n");
    s = gsl_odeiv2_driver_apply_fixed_step(d, &t, T / ns, ns, x);
    if (s != GSL_SUCCESS) {
      fprintf(stderr, "error: driver returned %d\n", s);
      exit(1);
    }
    if (t > 30000 * sqrt(1.0 / 8))
      break;
  }

  gsl_odeiv2_driver_free(d);
}

/*

cc main.c `gsl-config --libs`

*/
