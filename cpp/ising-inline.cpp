#include <fstream>
#include <iostream>
#include <memory>
#include <cmath>
#include <random>

#include "ising.hpp"

float Ising::initialise(float value)
{
  std::random_device device;
  std::mt19937 generator(device());
  std::uniform_real_distribution<double> U(0, 1);

  return U(generator);
}

Ising::Ising(std::string input, std::string output, size_t h, size_t w, float b)
{
  beta = b;
  outputFile = output;

  std::random_device device;
  std::mt19937 generator(device());
  std::uniform_real_distribution<double> U(0, 1);

  std::ifstream inputFile(input, std::ifstream::in);

  height = h;
  width  = w;

  X.resize(height);
  Y.resize(height);

  // Initialise X and Y
  for (size_t i = 0; i < height; i++) {
    X[i].resize(width);
    Y[i].resize(width);

    for (size_t j = 0; j < width; j++) {
      X[i][j] = initialise(0.0);
      inputFile >> Y[i][j];
    }
  }

  inputFile.close();
}

Ising::~Ising()
{

}

/**
 * Uses Metropolis-Hastings to solve the image.
 */
void Ising::solve(const unsigned int iterations)
{
  unsigned int inverse, k;
  float h, hInverse;
  float d, dInverse;
  float p;

  std::random_device device;
  std::mt19937 generator(device());
  std::uniform_real_distribution<double> U(0, 1);

  for (k = 0; k < iterations; k++) {
    if (k % 5 == 0) {
      std::cout << "iteration: " << k << std::endl;
    }

    for (size_t i = 0; i < height; i++) {
      for (size_t j = 0; j < width; j++) {

        inverse  = round(1 - X[i][j]);

        h        = -1.0/(2*noise) * pow((Y[i][j] - X[i][j]), 2);
        hInverse = -1.0/(2*noise) * pow((Y[i][j] - inverse), 2);

        Point coords = { i, j };

        d        = beta * neighbourhood(coords, X[i][j]) + h;
        dInverse = beta * neighbourhood(coords, inverse) + hInverse;

        p = exp(fmin(0, dInverse - d));

        if (U(generator) < p) {
          X[i][j] = inverse;
        }
      }
    }
  }
  std::cout << "iteration: " << k << std::endl;


  std::ofstream output(outputFile, std::ofstream::out);

  if (output.is_open()) {
    for (size_t i = 0; i < height; i++) {
      for (size_t j = 0; j < width; j++) {
        output << X[i][j] << " ";
      }
    }
  }

  output.close();
}

float Ising::neighbourhood(const Point coords, const float colour)
{
  unsigned int same_colours = 0.0;
  unsigned int i = coords.x;
  unsigned int j = coords.y;

  // Left, right, up, down
  if (i-1 > -1)     { if (X[i-1][j] == colour) same_colours++; }
  if (i+1 < height) { if (X[i+1][j] == colour) same_colours++; }
  if (j+1 < width)  { if (X[i][j+1] == colour) same_colours++; }
  if (j-1 > -1)     { if (X[i][j-1] == colour) same_colours++; }

  // Diagonal
  if (i-1 > -1 && j-1 > -1)         { if (X[i-1][j-1] == colour) same_colours++; }
  if (i+1 < height && j-1 > -1)     { if (X[i+1][j-1] == colour) same_colours++; }
  if (i+1 < height && j+1 < width)  { if (X[i+1][j+1] == colour) same_colours++; }
  if (i-1 > -1 && j+1 <  width)     { if (X[i-1][j+1] == colour) same_colours++; }

  // Extremal
  if (i-2 > -1)     { if (X[i-2][j] == colour) same_colours++; }
  if (i+2 < height) { if (X[i+2][j] == colour) same_colours++; }
  if (j+2 < width)  { if (X[i][j+2] == colour) same_colours++; }
  if (j-2 > -1)     { if (X[i][j-2] == colour) same_colours++; }

  return same_colours;
}

int main(int argc, char const *argv[]) {
  std::unique_ptr<Ising> ising(new Ising("data/noisy.txt", "denoised.txt", 328, 400));

  ising->solve(50);

  return 0;
}
