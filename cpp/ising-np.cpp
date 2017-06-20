#include <fstream>
#include <iostream>
#include <memory>
#include <cmath>
#include <random>

#include "ising.hpp"

Ising::Ising(std::vector<std::vector<float>> y, float b, float sigma)
{
  beta = b;
  noise = sigma;
  Y = y;

  height = Y.size();
  width  = Y[0].size();
}

Ising::~Ising()
{

}

/**
 * Uses Metropolis-Hastings to solve the image.
 */
std::vector<std::vector<float>> Ising::solve(std::vector<std::vector<float>> x, const unsigned int iterations)
{
  X = x;

  unsigned int inverse;
  float h, hInverse;
  float d, dInverse;
  float p;

  for (size_t i = 0; i < height; i++) {
    for (size_t j = 0; j < width; j++) {
      X[i][j] = U(generator);
    }
  }

  for (size_t k = 0; k < iterations; k++) {
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

  return X;
}

unsigned int Ising::neighbourhood(const Point coords, const float colour)
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
