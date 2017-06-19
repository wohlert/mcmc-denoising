#include <cmath>
#include <random>
#include <fstream>
#include <iostream>

using namespace std;

const size_t height = 300;
const size_t width  = 465;

/**
 * arraystd
 *
 * Calculates the standard deviation
 * of an array.
 */
template<size_t H, size_t W>
double arraystd(const double (&X)[H][W])
{
  double mean = 0;
  double variance = 0;

  for (size_t i = 0; i < H; i++) {
    for (size_t j = 0; j < W; j++) {
      mean += X[i][j];
    }
  }

  mean /= (H*W);

  for (size_t i = 0; i < H; i++) {
    for (size_t j = 0; j < W; j++) {
      variance += pow((X[i][j] - mean), 2);
    }
  }

  variance /= (H*W);

  return sqrt(variance);
}

/**
 * normpdf
 *
 * Evaluates the pdf of a normal distribution
 * with mean=0 and variance=variance
 */
double normpdf(double x, double variance)
{
  return exp(-pow(x, 2)/(2*variance))/sqrt(2*M_PI*variance);
}

/**
 * discretise
 */
double discretise(double x, size_t bins)
{
  double colour = round(x*bins)/bins;
  if (colour < (1/bins)) {
    colour = 1/bins;
  }
  if (colour > (1 - (1/bins))) {
    colour = (1 - (1/bins));
  }

  return colour;
}

template<size_t H, size_t W>
double neighbourhood(const double (&X)[H][W], const double colour,
                     const int i, const int j,
                     const double noise, const double beta)
{
  double colour_sum = 0;
  unsigned int indicator;


  if (i-1 > -1) {
    indicator = (int) (colour == X[i-1][j]);
    colour_sum += indicator;
  }

  if (i+1 < H) {
    indicator = (int) (colour == X[i+1][j]);
    colour_sum += indicator;
  }

  if (j+1 < W) {
    indicator = (int) (colour == X[i][j+1]);
    colour_sum += indicator;
  }

  if (j-1 > -1) {
    indicator = (int) (colour == X[i][j-1]);
    colour_sum += indicator;
  }

  if (i-1 > -1 && i+1 < H && j+1 < W && j-1 > -1) {
    colour_sum = normpdf(X[i][j] - colour, noise) * exp(beta * colour_sum);
  }

  return colour_sum;
}

// Noise should be the std of the image
template<size_t H, size_t W>
void metropolisHastings(double (&X)[H][W], const double (&Y)[H][W],
                        const unsigned int iterations, size_t bins,
                        const double noise, const double beta)
{
  unsigned int k;
  double candidate;
  double p;

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

        // Choose random candidate from X
        candidate  = (double)((rand()+1) % bins)/bins;

        //cout << candidate << endl;

        p = fmin(1, neighbourhood(X, candidate, i, j, noise, beta)/
                    neighbourhood(X, X[i][j],   i, j, noise, beta));

        if (U(generator) < p) {
          X[i][j] = candidate;
        }
      }
    }
  }
  cout << "iteration: " << k << endl;
}

int main(int argc, char const *argv[]) {
  size_t iterations = 1;
  size_t bins = 10;
  double beta  = 100;

  if (argc > 1) iterations = atoi(argv[1]);
  if (argc > 2) beta = atof(argv[2]);

  // Load in files
  ifstream inputFile("data/noisy.txt", ifstream::in);

  // Inititalise X and Y
  double X[height][width];
  double Y[height][width];

  if (inputFile.is_open()) {
    for (size_t i = 0; i < height; i++) {
      for (size_t j = 0; j < width; j++) {
        inputFile >> Y[i][j];
      }
    }

    for (size_t i = 0; i < height; i++) {
      for (size_t j = 0; j < width; j++) {
        X[i][j] = discretise(Y[i][j], bins);
      }
    }
  }

  // Perform computation
  metropolisHastings(X, Y, iterations, bins, arraystd(Y), beta);

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
