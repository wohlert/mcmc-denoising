#include <memory>
#include <cmath>
//#include <random>
#include <memory>
#include <iostream>
#include <cstdlib>
#include <stdexcept>

#include "denoising.hpp"

using namespace denoising;

float normpdf(float x, float variance)
{
  return exp(-pow(x, 2)/(2*variance))/sqrt(2*M_PI*variance);
}

float discretise(float x, int bins)
{
  float colour = round(x*bins)/bins;
  if (colour < (1/bins)) {
    colour = 1/bins;
  }
  if (colour > (1 - (1/bins))) {
    colour = (1 - (1/bins));
  }

  return colour;
}

Potts::Potts(const std::vector<std::vector<float>> y, float b, float sigma, int nBins)
{
  beta = b;
  noise = sigma;
  bins = nBins;

  Y = y;

  height = (int) Y.size();
  width  = (int) Y[0].size();

  X.resize(height);

  for (int i = 0; i < height; i++) {
    X[i].resize(width);
  }
}

Potts::~Potts()
{

}

/**
 * Uses Metropolis-Hastings to solve the image.
 */
std::vector<std::vector<float>> Potts::solve(const unsigned int iterations)
{
  for (int i = 0; i < height; i++) {
    for (int j = 0; j < width; j++) {
      X[i][j] = discretise(Y[i][j], bins);
    }
  }

  float candidate;
  float p = 0;

  for (size_t k = 0; k < iterations; k++) {
    for (int i = 0; i < height; i++) {
      for (int j = 0; j < width; j++) {

        // Choose random candidate from X
        candidate = (float)((rand()+1) % bins)/bins;

        Point coords = { i, j };
        p = fmin(1, neighbourhood(coords, candidate)/neighbourhood(coords, X[i][j]));

        if (((float) rand() / (RAND_MAX)) < p) {
          X[i][j] = candidate;
        }
      }
    }
  }

  return X;
}

float Potts::neighbourhood(const Point coords, const float colour)
{
  float colourSum = 0;
  unsigned int indicator;

  int i = coords.x;
  int j = coords.y;

  try {
    indicator = (int) (colour == X.at(i-1).at(j));
    colourSum += indicator;
  } catch (const std::exception& e) {}

  try {
    indicator = (int) (colour == X.at(i+1).at(j));
    colourSum += indicator;
  } catch (const std::exception& e) {}

  try {
    indicator = (int) (colour == X.at(i).at(j+1));
    colourSum += indicator;
  } catch (const std::exception& e) {}

  try {
    indicator = (int) (colour == X.at(i).at(j-1));
    colourSum += indicator;
  } catch (const std::exception& e) {}

  try {
    colourSum = normpdf(X.at(i).at(j) - colour, noise) * exp(beta * colourSum);
  } catch (const std::exception& e) {}

  return colourSum;
}
