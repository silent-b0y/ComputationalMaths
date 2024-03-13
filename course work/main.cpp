#include <iostream>
#include <cmath>
#include <iomanip>
#include "Forsythe.h"

Float f(Float x) {
  return std::log(1 + x) / (1 + x * x);
}

Float fZero(Float x) {
  return std::exp(x) - 2 * (x - 1) * (x - 1);
}

rkf parameters;
unsigned char* work = new unsigned char[6 * 2 * sizeof(Float) + sizeof(struct rkf_inside)];
Float C;
Float L;
Float L1;
Float L3;
Float R;
Float R1;
Float R2;
Float R3;
Float E1;
Float E2;
Float Umod[11];
Float Uexp[] = { -1.0, 7.777, 12.017, 10.701, 5.407, -0.843, -5.159, -6.015, -3.668, 0.283, 3.829 };

void Frkf(Float t, Float* X, Float* dX) {
  dX[0] = 1 / L1 * (E1 - E2 - X[2] + X[1] * R2 - X[0] * (R1 + R2));
  dX[1] = 1 / L3 * (E2 + X[2] + X[0] * R2 - X[1] * (R2 + R3));
  dX[2] = 1 / C * (X[0] - X[1]);
}

Float F(Float p) {
  C = p;
  Float sum = 0;
  Float X[3];
  X[0] = E1 / R1;
  X[1] = 0;
  X[2] = -E2;
  parameters.Y = X;
  parameters.t = 0.0;
  parameters.flag = 1;
  Umod[0] = X[2];
  std::cout << std::fixed << "t = " << 0.0 << " | " << Umod[0] << " | " << Umod[0] << "\n";
  for (int i = 1; i < 11; i++) {
    parameters.tout = parameters.t + 0.0001;
    rkf45(&parameters);
    std::cout << "t = "  <<  i / 10.0 << " | " << std::setw(9) << std::right << Uexp[i] << " | " << X[2] << "\n";
    Umod[i] = X[2];
    parameters.flag = 1;
  }
  for (int i = 0; i < 11; i++) {
    sum += std::pow(Uexp[i] - Umod[i], 2);
  }
  std::cout << "\n";
  return sum;
}


int main() {
  Float A[3][3] = { { 16.0, -18.0, 24.0 },
                  { -18.0, 49.0, -42.0 },
                  { 24.0, -42.0, 46.0 } };
  Float B[3] = { 304.0,
                218.0,
                166.0 };
  Float* cond = new Float();
  int* ipvt = new int();
  Decomp(3, *A, cond, ipvt);
  Solve(3, *A, B, ipvt);
  R = B[0];
  R2 = B[1];
  E2 = B[2];
  R1 = R;
  R3 = R;
  Float* errest = new Float();
  int* nofun = new int();
  Float* flag = new Float();
  L = 0.1469517 * Quanc8(f, 0.0, 1.0, 0.000001, 0.0, errest, nofun, flag);
  E1 = 18.75217 * Zeroin(fZero, 0.0, 1.7, 0.000001);
  L1 = L;
  L3 = L;
  parameters.f = Frkf;
  parameters.neqn = 3;
  parameters.re = 0.00000000001;
  parameters.ae = 0.00000000001;
  parameters.work = work;
  std:: cout << "C = " << FMin(F, (Float)0.0000005, (Float)0.000002, (Float)0.000000000001) * 1000000 << "\n\n";

  E1 += E1 * 0.01;
  E2 += E2 * 0.01;
  std::cout << "E1 = " << E1 << " E2 = " << E2 << "\n\n";
  Float X[3];
  X[0] = E1 / R1;
  X[1] = 0;
  X[2] = -E2;
  parameters.Y = X;
  parameters.t = 0.0;
  parameters.flag = 1;
  Umod[0] = X[2];
  std::cout << std::fixed << "t = " << 0.0 << " | " << Umod[0] << " | " << Umod[0] << "\n";
  for (int i = 1; i < 11; i++) {
    parameters.tout = parameters.t + 0.0001;
    rkf45(&parameters);
    std::cout << "t = " << i / 10.0 << " | " << std::setw(9) << std::right << Uexp[i] << " | " << X[2] << "\n";
    Umod[i] = X[2];
    parameters.flag = 1;
  }
  std::cout << "\n\n";

  C += C * 0.01;
  std::cout << "C = " << C * 1000000 << "\n\n";
  X[0] = E1 / R1;
  X[1] = 0;
  X[2] = -E2;
  parameters.Y = X;
  parameters.t = 0.0;
  parameters.flag = 1;
  Umod[0] = X[2];
  std::cout << std::fixed << "t = " << 0.0 << " | " << Umod[0] << " | " << Umod[0] << "\n";
  for (int i = 1; i < 11; i++) {
    parameters.tout = parameters.t + 0.0001;
    rkf45(&parameters);
    std::cout << "t = " << i / 10.0 << " | " << std::setw(9) << std::right << Uexp[i] << " | " << X[2] << "\n";
    Umod[i] = X[2];
    parameters.flag = 1;
  }
  std::cout << "\n";
  delete errest;
  delete nofun;
  delete flag;
  delete cond;
  delete ipvt;
  delete[] work;
  return 0;
}
