#include <iostream>
#include "Forsythe.h"

Float* matrixA(int N, Float Xn) {
  Float* matrix = new Float[N * N];
  int index = 0;
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      if (i == 0) {
        matrix[index] = 1;
      }
      else {
        if (j != N - 1) {
          matrix[index] = std::pow(j + 1, i);
        }
        else {
          matrix[index] = std::pow(Xn, i);
        }
      }
      index += 1;
    }
  }
  return matrix;
}

Float* matrixB(int N) {
  Float* matrix = new Float[N];
  for (int i = 0; i < N; i++) {
    matrix[i] = std::pow(2, i + 1) + std::cos(i + 1);
  }
  return matrix;
}

int main() {
  Float* A;
  Float* B;
  Float* cond = new Float();
  int* ipvt = new int[5];
  Float Xn[4] = { 1.1, 1.01, 1.001, 1.0001 };
  for (int N = 5; N < 10; N += 2) {
    std::cout << "N = " << N << "\n";
    for (int i = 0; i < 4; i++) {
      A = matrixA(N, Xn[i]);
      B = matrixB(N);
      Decomp(N, A, cond, ipvt);
      Solve(N, A, B, ipvt);
      for (int j = 0; j < N; j++) {
        std::cout << std::fixed << B[j] << "  ";
      }
      std::cout << "| Xn = " << Xn[i] << "\n" << "cond = " << *cond << "\n";
      delete[] A;
      delete[] B;
    }
  }
  delete cond;
  delete[] ipvt;
  return 0;
}
