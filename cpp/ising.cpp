#include <memory>
#include <cmath>
#include <stdexcept>
//#include <random>
#include <cstdlib>

#include "denoising.hpp"

using namespace denoising;

Ising::Ising(const std::vector<std::vector<float>> y, float b, float sigma)
{
  beta = b;
  noise = sigma;
  Y = y;

  height = (int) Y.size();
  width  = (int) Y[0].size();

  X.resize(height);

  for (int i = 0; i < height; i++) {
    X[i].resize(width);
  }
}

Ising::~Ising()
{

}

/**
 * Uses Metropolis-Hastings to solve the image.
 */
std::vector<std::vector<int>> Ising::metropolisHastings(const unsigned int iterations)
{
  unsigned int inverse;
  float h, hInverse;
  float d, dInverse;
  float p;

  //std::random_device device;
  //std::mt19937 generator(device());
  //std::uniform_real_distribution<float> U(0, 1);

  for (int i = 0; i < height; i++) {
    for (int j = 0; j < width; j++) {
      X[i][j] = ((float) rand() / (RAND_MAX));
    }
  }

  for (size_t k = 0; k < iterations; k++) {
    for (int i = 0; i < height; i++) {
      for (int j = 0; j < width; j++) {

        inverse  = round(1 - X[i][j]);

        h        = -1.0/(2*noise) * pow((Y[i][j] - X[i][j]), 2);
        hInverse = -1.0/(2*noise) * pow((Y[i][j] - inverse), 2);

        Point coords = { i, j };

        d        = beta * neighbourhood(coords, X[i][j]) + h;
        dInverse = beta * neighbourhood(coords, inverse) + hInverse;

        p = exp(fmin(0, dInverse - d));

        if (((float) rand() / (RAND_MAX)) < p) {
          X[i][j] = inverse;
        }
      }
    }
  }

  return X;
}

int Ising::neighbourhood(const Point coords, const float colour)
{
  int same_colours = 0.0;
  int i = coords.x;
  int j = coords.y;

  // Left, right, up, down
  try { if (X.at(i-1).at(j) == colour) same_colours++; } catch (const std::out_of_range& oor) {}
  try { if (X.at(i+1).at(j) == colour) same_colours++; } catch (const std::out_of_range& oor) {}
  try { if (X.at(i).at(j-1) == colour) same_colours++; } catch (const std::out_of_range& oor) {}
  try { if (X.at(i).at(j+1) == colour) same_colours++; } catch (const std::out_of_range& oor) {}

  // Diagonal
  try { if (X.at(i-1).at(j-1) == colour) same_colours++; } catch (const std::out_of_range& oor) {}
  try { if (X.at(i+1).at(j-1) == colour) same_colours++; } catch (const std::out_of_range& oor) {}
  try { if (X.at(i+1).at(j+1) == colour) same_colours++; } catch (const std::out_of_range& oor) {}
  try { if (X.at(i-1).at(j+1) == colour) same_colours++; } catch (const std::out_of_range& oor) {}

  // Extremal
  try { if (X.at(i-2).at(j) == colour) same_colours++; } catch (const std::out_of_range& oor) {}
  try { if (X.at(i+2).at(j) == colour) same_colours++; } catch (const std::out_of_range& oor) {}
  try { if (X.at(i).at(j-2) == colour) same_colours++; } catch (const std::out_of_range& oor) {}
  try { if (X.at(i).at(j+2) == colour) same_colours++; } catch (const std::out_of_range& oor) {}

  return same_colours;
}
