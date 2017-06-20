#include <vector>

struct Point {
  size_t x;
  size_t y;
};

class Ising {
protected:
  float beta;
  float noise = 0.1;

  std::string outputFile;

  std::vector<std::vector<unsigned int>> X;
  std::vector<std::vector<float>> Y;

  size_t height, width;

  float initialise(float value);

public:
  Ising(std::string input, std::string output, size_t h, size_t w, float b=1.0);
  ~Ising();

  float neighbourhood(const Point coords, const float colour);
  void solve(const unsigned int iterations);
};
