#include <vector>

struct Point {
  size_t x;
  size_t y;
};

class Ising {
protected:
  float beta;
  float noise;

  std::vector<std::vector<unsigned int>> X;
  std::vector<std::vector<float>> Y;

  size_t height, width;

public:
  Ising(std::vector<std::vector<float>> y, float b=1.0, float sigma=0.1);
  ~Ising();

  unsigned int neighbourhood(const Point coords, const float colour);
  std::vector<std::vector<float>> solve(std::vector<std::vector<unsigned int>> x, const unsigned int iterations);
};
