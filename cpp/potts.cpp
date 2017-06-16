#include <cmath>
#include <random>
#include <fstream>
#include <iostream>

using namespace std;

const size_t height = 400;
const size_t width  = 328;

float phi(float u) {
  return log(1 + u*u);
}

template<size_t H, size_t W>
void saveImage(const double (&Z)[H][W], const unsigned int k)
{
  string filename = "history/history_" + std::to_string(k) + ".txt";
  ofstream saveFile(filename, ofstream::out);

  if (saveFile.is_open()) {
    for (size_t i = 0; i < H; i++) {
      for (size_t j = 0; j < W; j++) {
        saveFile << Z[i][j] << " ";
      }
    }
  }

  saveFile.close();
}

template<size_t H, size_t W>
double neighbourhood(const double (&X)[H][W], const double colour,
                     const int i, const int j)
{
  double colour_sum = 0;
  double w1    = 1;
  double w2    = 0.5;
  double w3    = 0.25;
  double delta = 50;
  double diff;

  // Left, right, up, down
  //if ((W/2 - i) > -1) {
    //diff = (X[i-1][j] - colour)/delta;

    if (X[i-1][j] != colour) {
      diff = 1;
    } else {
      diff = 0;
    }
    colour_sum += diff;
  //}

  //if (i+W < H/2) {
    //diff = (X[i+1][j] - colour)/delta;

    if (X[i+1][j] != colour) {
      diff = 1;
    } else {
      diff = 0;
    }

    colour_sum += diff;
  //}

  //if (j+1 < W/2) {
    //diff = (X[i][j+1] - colour)/delta;

    if (X[i][j+1] != colour) {
      diff = 1;
    } else {
      diff = 0;
    }

    colour_sum += diff;
  //}
  //if ((H/2 - i) > -1) {
    //diff = (X[i][j-1] - colour)/delta;

    if (X[i][j-1] != colour) {
      diff = 1;
    } else {
      diff = 0;
    }

    colour_sum += diff;
  //}

  colour_sum = w1 * colour_sum;

  /*
  // Diagonal
  if (i-1 > -1 && j-1 > -1)  { colour_sum += phi(X[i-1][j-1] - colour)/delta; }
  if (i+1 <  H && j-1 > -1)  { colour_sum += phi(X[i+1][j-1] - colour)/delta; }
  if (i+1 <  H && j+1 <  W)  { colour_sum += phi(X[i+1][j+1] - colour)/delta; }
  if (i-1 > -1 && j+1 <  W)  { colour_sum += phi(X[i-1][j+1] - colour)/delta; }

  // Extremal
  if (i-2 > -1) { colour_sum += phi(X[i-2][j] - colour)/delta; }
  if (i+2 <  H) { colour_sum += phi(X[i+2][j] - colour)/delta; }
  if (j+2 <  W) { colour_sum += phi(X[i][j+2] - colour)/delta; }
  if (j-2 > -1) { colour_sum += phi(X[i][j-2] - colour)/delta; }
  */

  // First derivative
  return colour_sum;
}

template<size_t H, size_t W>
void metropolisHastings(double (&X)[H][W], const double (&Y)[H][W],
                        const unsigned int iterations,
                        const double noise, const double beta)
{
  unsigned int k;
  double inverse;
  double h, hInverse;
  double d = 0;
  double dInverse = 0;
  double p, u;

  // Use Mersenne Twister random number generator
  random_device device;
  mt19937 generator(device());
  uniform_real_distribution<double> U(0, 1);
  normal_distribution<double>       N(0, noise);

  double Z[height][width];

  for (k = 0; k < iterations; k++) {
    if (k % 5 == 0) {
      cout << "iteration: " << k << endl;
    }

    for (size_t i = 0; i < H; i++) {
      for (size_t j = 0; j < W; j++) {

        inverse  = (X[i][j] + N(generator));

        h        = -1.0/(2*noise) * pow((Y[i][j] - X[i][j]), 2);
        hInverse = -1.0/(2*noise) * pow((Y[i][j] - inverse), 2);

        d        = beta * neighbourhood(X, X[i][j], i, j) + h;
        dInverse = beta * neighbourhood(X, inverse, i, j) + hInverse;

        u = U(generator);
        p = exp(fmin(0, dInverse - d));

        if (u < p) {
          X[i][j] = Y[i][j];
        }
        Z[i][j] = neighbourhood(Y, Y[i][j], i, j);
      }
    }
    saveImage(Z, k);
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
  ifstream inputFile("data/original.txt", ifstream::in);

  // Instantiate Mersenne Twister
  random_device device;
  mt19937 generator(device());
  uniform_real_distribution<double> U(0, 1);

  // Inititalise X and Y
  double X[height][width];
  double Y[height][width];

  for (size_t i = 0; i < height; i++) {
    for (size_t j = 0; j < width; j++) {
      X[i][j] = U(generator);
    }
  }

  if (inputFile.is_open()) {
    for (size_t i = 0; i < height; i++) {
      for (size_t j = 0; j < width; j++) {
        inputFile >> Y[i][j];
      }
    }
  }

  // Perform computation
  metropolisHastings(X, Y, iterations, noise, beta);

  ofstream outputFile("denoised.txt", ofstream::out);
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
