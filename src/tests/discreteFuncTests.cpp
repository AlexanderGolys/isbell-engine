#include <cassert>
#include <iostream>
#include "../engine/specific.hpp"
#include "../geometry/pdeDiscrete.hpp"

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

  printf("OK: FFT inverse test passed");
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

  printf("OK: FFT real symmetry test passed");
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
    assert(nearlyEqual(padded[i+7], fn[i]));

  printf("OK: zero padding");
}



void gaborTest()
{
  auto fn = DiscreteRealFunction(X_R*X_R-1, vec2(0, 1), 100);
  auto gab = DiscreteGaborTransform(10).transform(fn);

  printf("OK: Gabor transform");

}


void quaternionTest() {
  auto i = Quaternion::i();
  auto j = Quaternion::j();
  auto one = Quaternion::one();
  auto q = Quaternion(.5, .5, .5, .5);
  auto ii = i*i;
  assert(nearlyEqual(q.norm(), 1.f));
  assert(nearlyEqual(i*i, j*j));
  assert(nearlyEqual(i*i, -one));
  assert(nearlyEqual(i*j, -j*i));

  printf("OK: Quaternions arithmetic");

}

  int main()
  {
    fftInverseTest();
    fftRealDomainSymmetryDFTTest();
    paddingTest();
    gaborTest();
    quaternionTest();
    return 0;
  }
