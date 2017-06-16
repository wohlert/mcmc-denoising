#include <cmath>
#include <random>
#include <fstream>
#include <iostream>
//#include "omp.h"

using namespace std;

template<size_t H, size_t W>
unsigned int neighbourhood(const unsigned int (&X)[H][W],
                           const unsigned int i, const unsigned int j,
                           const unsigned int colour)
{
  unsigned int same_colours = 0;

  // Left, right, up, down
  if (i-1 > -1)  { if (X[i-1][j] == colour) same_colours++; }
  if (i+1 <  H)  { if (X[i+1][j] == colour) same_colours++; }
  if (j+1 <  W)  { if (X[i][j+1] == colour) same_colours++; }
  if (j-1 > -1)  { if (X[i][j-1] == colour) same_colours++; }

  // Diagonal
  if (i-1 > -1 && j-1 > -1)  { if (X[i-1][j-1] == colour) same_colours++; }
  if (i+1 <  H && j-1 > -1)  { if (X[i+1][j-1] == colour) same_colours++; }
  if (i+1 <  H && j+1 <  W)  { if (X[i+1][j+1] == colour) same_colours++; }
  if (i-1 > -1 && j+1 <  W)  { if (X[i-1][j+1] == colour) same_colours++; }

  // Extremal
  if (i-2 > -1)  { if (X[i-2][j] == colour) same_colours++; }
  if (i+2 <  H)  { if (X[i+2][j] == colour) same_colours++; }
  if (j+2 <  W)  { if (X[i][j+2] == colour) same_colours++; }
  if (j-2 > -1)  { if (X[i][j-2] == colour) same_colours++; }

  return same_colours;
}

template<size_t H, size_t W>
void metropolisHastings(unsigned int (&X)[H][W], const float (&Y)[H][W],
                        const unsigned int iterations,
                        const float noise, const float beta)
{
  unsigned int inverse, k;
  float h, hInverse;
  float d, dInverse;
  float p, u;

  // Use Mersenne Twister random number generator
  random_device device;
  mt19937 generator(device());
  uniform_real_distribution<double> U(0, 1);

  for (k = 0; k < iterations; k++) {
    if (k % 5 == 0) {
      cout << "iteration: " << k << endl;
    }

    for (size_t i = 0; i < H; i++) {
      for (size_t j = 0; j < W; j++) {
        inverse  = 1 - X[i][j];

        h        = -1.0/(2*noise) * pow((Y[i][j] - X[i][j]), 2);
        hInverse = -1.0/(2*noise) * pow((Y[i][j] - inverse), 2);

        d        = beta * neighbourhood(X, i, j, X[i][j]) + h;
        dInverse = beta * neighbourhood(X, i, j, inverse) + hInverse;

        u = U(generator);
        p = exp(fmin(0, dInverse - d));

        if (u < p) {
          X[i][j] = inverse;
        }
      }
    }
  }
  cout << "iteration: " << k << endl;
}

int main(int argc, char const *argv[]) {
  unsigned int iterations = 1;
  float noise = 0.1;
  float beta  = 1.0;

  if (argc > 1) iterations = atoi(argv[1]);
  if (argc > 2) noise = atof(argv[2]);
  if (argc > 3) beta  = atof(argv[3]);

  // Load in files
  ifstream inputFile("horse.txt", ifstream::in);
  ofstream outputFile("denoised.txt", ofstream::out);

  // Instantiate Mersenne Twister
  random_device device;
  mt19937 generator(device());
  uniform_int_distribution<unsigned int> U(0, 1);

  // Iniitalise X and Y
  const size_t width  = 400;
  const size_t height = 328;

  const size_t window_width  = width/2;
  const size_t window_height = height/2;


  float        Y[height][width];
  unsigned int X[window_height][window_width][4];
  float        Z[window_height][window_width][4];

  for (size_t i = 0; i < height; i++) {
    for (size_t j = 0; j < width; j++) {
      inputFile >> Y[i][j];
    }
  }

  size_t k = 0;
  int k1, k2;
  k1 = k2 = 1;

  for (size_t i = 0; i < height; i++) {
    for (size_t j = 0; j < width; j++) {
      cout << k << endl;
      if (i == window_height && k1 == 1) {
        k++;
        k1 = 0;
      }
      if (j == window_width && k2 == 1) {
        k++;
        k2 = 0;
      }
      X[i % window_height][j % window_width][k] = U(generator);
      Z[i % window_height][j % window_width][k] = Y[i][j];
    }
  }

  /*
  // Perform computation
  metropolisHastings(X, Y, iterations, noise, beta);

  // Write output to file
  if (outputFile.is_open()) {
    for (size_t i = 0; i < height; i++) {
      for (size_t j = 0; j < width; j++) {
        outputFile << X[i][j] << " ";
      }
    }
  }

  */
  inputFile.close();
  outputFile.close();

  return 0;
}
