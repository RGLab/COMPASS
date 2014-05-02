#include <math.h>

double digamma(double x) {
  double r, f1, f2, f3, f4;
  r = 0;

  while (x < 5) {
    r -= 1 / x;
    x += 1;
  }

  f1 = 1 / (x * x);
  f2 = f1 * f1;
  f3 = f2 * f1;
  f4 = f2 * f2;

  static const double c1 = -1 / 12.0;
  static const double c2 = 1 / 120.0;
  static const double c3 = -1 / 252.0;
  static const double c4 = 1 / 240.0;
  static const double c5 = -1 / 132.0;
  static const double c6 = 691 / 32760.0;
  static const double c7 = -1 / 12.0;
  static const double c8 = 3617 / 8160.0;

  double A = c1 + c2 * f1 + c3 * f2 + c4 * f3;
  double B = c5 * f1 + c6 * f2 + c7 * f3 + c8 * f4;
  return r + log(x) - 0.5 / x + f1 * A + f4 * B;
}
