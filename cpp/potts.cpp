#include <cmath>
#include <random>
#include <fstream>
#include <iostream>

using namespace std;

float phi(float u) {
  //return log(1 + pow(u, 2));
  return u;
}

template<size_t H, size_t W>
unsigned int neighbourhood(const double (&X)[H][W],
                           const unsigned int i, const unsigned int j,
                           const unsigned int colour)
{
  double colour_sum = 0;
  double w1 = 1;
  double w2 = 0.5;
  double w3 = 0.25;
  double delta = 100;

  // Left, right, up, down
  if (i-1 > -1)  { colour_sum += w1 * phi(X[i-1][j] - colour)/delta; }
  if (i+1 <  H)  { colour_sum += w1 * phi(X[i+1][j] - colour)/delta; }
  if (j+1 <  W)  { colour_sum += w1 * phi(X[i][j+1] - colour)/delta; }
  if (j-1 > -1)  { colour_sum += w1 * phi(X[i][j-1] - colour)/delta; }

  // Diagonal
  if (i-1 > -1 && j-1 > -1)  { colour_sum += w2 * phi(X[i-1][j-1] - colour)/delta; }
  if (i+1 <  H && j-1 > -1)  { colour_sum += w2 * phi(X[i+1][j-1] - colour)/delta; }
  if (i+1 <  H && j+1 <  W)  { colour_sum += w2 * phi(X[i+1][j+1] - colour)/delta; }
  if (i-1 > -1 && j+1 <  W)  { colour_sum += w2 * phi(X[i-1][j+1] - colour)/delta; }

  // Extremal
  if (i-2 > -1) { colour_sum += w3 * phi(X[i-2][j] - colour)/delta; }
  if (i+2 <  H) { colour_sum += w3 * phi(X[i+2][j] - colour)/delta; }
  if (j+2 <  W) { colour_sum += w3 * phi(X[i][j+2] - colour)/delta; }
  if (j-2 > -1) { colour_sum += w3 * phi(X[i][j-2] - colour)/delta; }

  // First derivative
  return colour_sum;
}

template<size_t H, size_t W>
void metropolisHastings(double (&X)[H][W], const double (&Y)[H][W],
                        const unsigned int iterations,
                        const double noise, const double beta)
{
  unsigned int inverse, k;
  double h, hInverse;
  double d, dInverse;
  double p, u;

  // Use Mersenne Twister random number generator
  random_device device;
  mt19937 generator(device());
  uniform_real_distribution<double> U(0, 1);
  normal_distribution<double>       N(0, noise);

  for (k = 0; k < iterations; k++) {
    if (k % 5 == 0) {
      cout << "iteration: " << k << endl;
    }

    for (size_t i = 0; i < H; i++) {
      for (size_t j = 0; j < W; j++) {
        double normal = N(generator);
        inverse  = (X[i][j] + normal);

        h        = -1.0/(2*noise) * pow((Y[i][j] - X[i][j]), 2);
        hInverse = -1.0/(2*noise) * pow((Y[i][j] - inverse), 2);

        d        = beta * neighbourhood(X, i, j, X[i][j]) + h;
        dInverse = beta * neighbourhood(X, i, j, inverse) + hInverse;

        u = U(generator);
        p = exp(fmin(0, dInverse - d));

        if (u < p) {
          X[i][j] = neighbourhood(X, i, j, inverse);
        }
      }
    }
  }
  cout << "iteration: " << k << endl;
}

int main(int argc, char const *argv[]) {
  unsigned int iterations = 1;
  double noise = 0.1;
  double beta  = 1.0;

  if (argc > 1) iterations = atoi(argv[1]);
  if (argc > 2) noise = atof(argv[2]);
  if (argc > 3) beta  = atof(argv[3]);

  // Load in files
  ifstream inputFile("data/noisy.txt", ifstream::in);
  ofstream outputFile("denoised.txt", ofstream::out);

  // Instantiate Mersenne Twister
  random_device device;
  mt19937 generator(device());
  uniform_real_distribution<double> U(0, 1);

  // Iniitalise X and Y
  const size_t height = 400;
  const size_t width  = 328;

  double X[height][width];
  double Y[height][width];

  if (inputFile.is_open()) {
    for (size_t i = 0; i < height; i++) {
      for (size_t j = 0; j < width; j++) {
        X[i][j] = U(generator);
        inputFile >> Y[i][j];
      }
    }
  }

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

  inputFile.close();
  outputFile.close();

  return 0;
}
