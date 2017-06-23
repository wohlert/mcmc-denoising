#include <vector>

namespace denoising {
  struct Point {
    int x;
    int y;
  };

  /*
  class BaseHook {
  public:
    BaseHook() {};
    ~BaseHook();
    std::vector<std::vector<std::vector<float>>> history;

    inline float update(float value) {
      return value;
    }
    void iterate(unsigned int iteration, const std::vector<std::vector<float>> X) {}
  };

  class MAPHook : public BaseHook {
  private:
    float t, init, diffusion;

  public:
    MAPHook(float tInit, float diff=0.995) : BaseHook() {
      init = tInit;
      t = tInit;
      diffusion = diff;
    }
    ~MAPHook();
    inline float update(float value) {
      return value * exp(-t);
    }
    inline void iterate(unsigned int iteration, std::vector<std::vector<float>> X) {
      t = init * pow(diffusion, iteration);
    }
  };

  class SaveHook : public BaseHook {
  public:

    SaveHook() : BaseHook() {}
    ~SaveHook();
    inline void iterate(unsigned int iteration, std::vector<std::vector<float>> X) {
      history.push_back(X);
    }
  };
  */

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
    std::vector<std::vector<int>> metropolisHastings(const unsigned int iterations);
    std::vector<std::vector<int>> gibbs(const unsigned int iterations);
  };

  class Potts {
  protected:
    float beta;
    float noise;
    int bins;

    std::vector<std::vector<std::vector<float>>> history;

    std::vector<std::vector<float>> X;
    std::vector<std::vector<float>> Y;

    int height, width;

    float neighbourhood(const Point coords, const float colour);

  public:
    Potts(const std::vector<std::vector<float>> y, float b=1.0, float sigma=0.1, int nBins=10);
    ~Potts();
    std::vector<std::vector<float>> metropolisHastings(const unsigned int iterations);
    std::vector<std::vector<float>> metropolisHastings(const unsigned int iterations, std::vector<std::vector<float>> x);
    std::vector<std::vector<float>> MAP(const unsigned int iterations, const float tInit=4, const float diffusion=0.995);
    std::vector<std::vector<float>> MAP(const unsigned int iterations, std::vector<std::vector<float>> x, const float tInit=4, const float diffusion=0.995);
    inline std::vector<std::vector<std::vector<float>>> getHistory() { return history; }
  };
}
