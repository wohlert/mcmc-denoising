#include <cmath>
//#include <random>
#include <cstdlib>
#include <stdexcept>
#include <memory>

#include "denoising.hpp"

using namespace denoising;

float normpdf(float x, float variance)
{
  return -pow(x, 2)/(2*variance);
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
std::vector<std::vector<float>> Potts::metropolisHastings(const unsigned int iterations)
{
  history.clear();

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
        p = exp(fmin(0, neighbourhood(coords, candidate) - neighbourhood(coords, X[i][j])))/sqrt(2*M_PI*noise);
        //p = hook->update(p);

        if (((float) rand() / (RAND_MAX)) < p) {
          X[i][j] = candidate;
        }
      }
    }
    //hook->iterate(k, X);
    history.push_back(X);
  }

  return X;
}

std::vector<std::vector<float>> Potts::metropolisHastings(const unsigned int iterations, std::vector<std::vector<float>> x)
{
  history.clear();

  X = x;

  float candidate;
  float p = 0;

  for (size_t k = 0; k < iterations; k++) {
    for (int i = 0; i < height; i++) {
      for (int j = 0; j < width; j++) {

        // Choose random candidate from X
        candidate = (float)((rand()+1) % bins)/bins;

        Point coords = { i, j };
        p = exp(fmin(0, neighbourhood(coords, candidate) - neighbourhood(coords, X[i][j])))/sqrt(2*M_PI*noise);
        //p = hook->update(p);

        if (((float) rand() / (RAND_MAX)) < p) {
          X[i][j] = candidate;
        }
      }
    }
    //hook->iterate(k, X);
    history.push_back(X);
  }

  return X;
}

std::vector<std::vector<float>> Potts::MAP(const unsigned int iterations, const float tInit, const float diffusion)
{
  history.clear();
  float t = tInit;

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
        p = exp(fmin(0, neighbourhood(coords, candidate) - neighbourhood(coords, X[i][j]))/t)/sqrt(2*M_PI*noise);

        if (((float) rand() / (RAND_MAX)) < p) {
          X[i][j] = candidate;
        }
      }
    }
    history.push_back(X);
    t = tInit * pow(diffusion, k);
  }

  return X;
}

std::vector<std::vector<float>> Potts::MAP(const unsigned int iterations, std::vector<std::vector<float>> x, const float tInit, const float diffusion)
{
  history.clear();
  float t = tInit;

  X = x;

  float candidate;
  float p = 0;

  for (size_t k = 0; k < iterations; k++) {
    for (int i = 0; i < height; i++) {
      for (int j = 0; j < width; j++) {

        // Choose random candidate from X
        candidate = (float)((rand()+1) % bins)/bins;

        Point coords = { i, j };
        p = exp(fmin(0, neighbourhood(coords, candidate) - neighbourhood(coords, X[i][j]))/t)/sqrt(2*M_PI*noise);

        if (((float) rand() / (RAND_MAX)) < p) {
          X[i][j] = candidate;
        }
      }
    }
    history.push_back(X);
    t = tInit * pow(diffusion, k);
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
    colourSum = normpdf(X.at(i).at(j) - colour, noise) + beta * colourSum;
  } catch (const std::exception& e) {}

  return colourSum;
}

/*
int main(int argc, char const *argv[]) {
  std::vector<std::vector<float>> Y;
  int n = 100;

  Y.resize(n);
  for (size_t i = 0; i < n; i++) {
    Y[i].resize(n);
    for (size_t j = 0; j < n; j++) {
      Y[i][j] = i + j;
    }
  }

  std::unique_ptr<Potts> potts(new Potts(Y, 10, 0.1, 10));
  potts->MAP(100, 4);

  return 0;
}
*/
