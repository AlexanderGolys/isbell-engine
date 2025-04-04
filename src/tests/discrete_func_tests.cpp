#include "../geometry/pde_dicrete.hpp"
#include "../common/specific.hpp"
#include <cassert>
#include <iostream>

using namespace glm;

void fftInverseTest()
{
	auto fn = DiscreteRealFunction(max(1-abs(3*X_R-3), 0), vec2(0, 9), 200);
  auto fn_fft = fn.fft();
  auto fn_back = fn_fft.ifft().re();

  assert(nearlyEqual(fn.getDomain(), fn_back.getDomain()));
  assert(fn.samples() == fn_back.samples());
  assert(fn.samples() == fn_fft.samples());
  assert(fn.samples() == 200);
  assert(nearlyEqual(fn.sampling_step(), fn_back.sampling_step()));
  assert((fn-fn_back).L2_norm() < 0.01);
  std::cout << "OK: FFT inverse test passed" << std::endl;
}

void fftShiftsTest()
{
  auto fn = DiscreteRealFunction(max(1-abs(3*X_R-3), 0), vec2(0, 9), 200);
  auto fn_fft = fn.fft();
  auto fn_back = fn_fft.ifft().re();

  assert(nearlyEqual(fn.getDomain(), fn_back.getDomain()));
  assert(fn.samples() == fn_back.samples());
  assert(fn.samples() == fn_fft.samples());
  assert(fn.samples() == 200);
  assert(nearlyEqual(fn.sampling_step(), fn_back.sampling_step()));
  assert((fn-fn_back).L2_norm() < 0.01);
  std::cout << "OK: FFT inverse test passed" << std::endl;
}




  int main(void)
  {
    fftInverseTest();
    return 0;
  }
