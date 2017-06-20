#include <vector>

namespace denoising {
  struct Point {
    int x;
    int y;
  };

  class Ising {
  protected:
    float beta;
    float noise;

    std::vector<std::vector<int>> X;
    std::vector<std::vector<float>> Y;

    int height, width;

    int neighbourhood(const Point coords, const float colour);

  public:
    Ising(const std::vector<std::vector<float>> y, float b=1.0, float sigma=0.1);
    ~Ising();
    std::vector<std::vector<int>> solve(const unsigned int iterations);
  };

  class Potts {
  protected:
    float beta;
    float noise;
    int bins;

    std::vector<std::vector<float>> X;
    std::vector<std::vector<float>> Y;

    int height, width;

    float neighbourhood(const Point coords, const float colour);

  public:
    Potts(const std::vector<std::vector<float>> y, float b=1.0, float sigma=0.1, int nBins=10);
    ~Potts();
    std::vector<std::vector<float>> solve(const unsigned int iterations);
  };
}
