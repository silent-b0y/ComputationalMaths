#include <iostream>
#include "Forsythe.h"

int main() {
  Float* X = new Float[6]{ -1.0, -0.9, -0.8, -0.7, -0.6, -0.5 };
  Float* Y = new Float[6]{ 0.5440, -0.4121, -0.9894, -0.6570, 0.2794, 0.9589 };
  Float* B = new Float[5];
  Float* C = new Float[5];
  Float* D = new Float[5];
  Spline(6, X, Y, B, C, D);
  Float err = 0.000001;
  Float a = -1.0;
  Float b = -0.5;
  Float c = 0.0;
  Float f_c = 0.0;
  Float f_a = 0.0;
  while (ABS(b - a) > err) {
    c = (a + b) / 2.0;
    f_c = SEval(6, c, X, Y, B, C, D) - 1.8 * c * c;
    f_a = SEval(6, a, X, Y, B, C, D) - 1.8 * a * a;
    if (SIGN(f_a) == SIGN(f_c)) {
      a = c;
    }
    else {
      b = c;
    }
    std::cout << std::fixed << "c = " << c << "   b - a = " << b - a << "\n";
  }
  std::cout << "Result = " << c;
  return 0;
}
