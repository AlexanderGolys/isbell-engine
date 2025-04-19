#include <cassert>
#include <iostream>
#include "../common/specific.hpp"
#include "../geometry/pde_dicrete.hpp"

using namespace glm;

void fftInverseTest()
{
	auto fn = DiscreteRealFunction(max(15-abs(X_R*X_R-16), -1), vec2(0, 7), 16);
  auto fn_fft = fn.fft();
  auto fn_back = fn_fft.ifft().re();

  assert(nearlyEqual(fn.getDomain(), fn_back.getDomain()));
  assert(fn.samples() == fn_back.samples());
  assert(fn.samples() == fn_fft.samples());
  assert(fn.samples() == 16);
  assert(nearlyEqual(fn.sampling_step(), fn_back.sampling_step()));
  assert((fn-fn_back).L2_norm() < 0.01);
  std::cout << "OK: FFT inverse test passed" << std::endl;
}


void fftRealDomainSymmetryDFTTest()
{
  auto fn = DiscreteRealFunction(X_R*X_R-1, vec2(0, 1), 16);
  auto fn_fft = fn.fft();

  for (int i = 1; i < fn_fft.samples(); ++i)
  {
    auto c = fn_fft[i];
    auto c_conj = fn_fft[fn_fft.samples()-i].conj();
    assert(c.nearlyEqual(c_conj));
  }

  std::cout << "OK: FFT real symmetry test passed" << std::endl;
}

void paddingTest()
{
  auto fn = DiscreteRealFunction(X_R*X_R-1, vec2(0, 1), 50);
  auto padded = fn.two_sided_zero_padding(64);

  assert(padded.samples() == 64);
  for (int i = 0; i < 7; ++i)
  {
    assert(nearlyEqual(padded[i], 0.f));
    assert(nearlyEqual(padded[63-i], 0.f));
  }
  for (int i = 0; i < 50; ++i)
  {
    assert(nearlyEqual(padded[i+7], fn[i]));
  }
  std::cout << "OK: zero padding" << std::endl;
}



void gaborTest()
{
  auto fn = DiscreteRealFunction(X_R*X_R-1, vec2(0, 1), 100);
  auto gab = DiscreteGaborTransform(10).transform(fn);

  std::cout << "OK: Gabor transform" << std::endl;

}


  int main()
  {
    fftInverseTest();
    fftRealDomainSymmetryDFTTest();
    paddingTest();
    gaborTest();
    return 0;
  }
