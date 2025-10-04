#pragma once
#include "unittests.hpp"
#include "../utils/func.hpp"
#include "../utils/integralTransforms.hpp"
#include "../utils/logging.hpp"

using namespace glm;

inline bool fftInverseTest()
{
    bool passed = true;
    DiscreteRealFunction fn = DiscreteRealFunction(max(15-abs(X_R*X_R-16), -1), vec2(0, 7), 16);
    DiscreteComplexFunction fn_fft = fn.fft();
    DiscreteRealFunction fn_back = fn_fft.ifft().re();

    passed &= assertNearlyEqual_UT(fn.getDomain(), fn_back.getDomain());
    passed &= assertEqual_UT(fn.samples(), fn_back.samples());
    passed &= assertEqual_UT(fn.samples(), fn_fft.samples());
    passed &= assertEqual_UT(fn.samples(), 16);
    passed &= assertNearlyEqual_UT(fn.sampling_step(), fn_back.sampling_step());
    passed &= assertLess_UT((fn-fn_back).L2_norm(), 0.01);

    return passed;

}


inline bool fftRealDomainSymmetryDFTTest()
{

  bool passed = true;
  auto fn = DiscreteRealFunction(X_R*X_R-1, vec2(0, 1), 16);
  auto fn_fft = fn.fft();

  for (int i = 1; i < fn_fft.samples(); ++i)
  {
    auto c = fn_fft[i];
    auto c_conj = fn_fft[fn_fft.samples()-i].conj();
    passed &= assertNearlyEqual_UT(c_conj, c);
  }

  return passed;

}

inline bool paddingTest()
{

  bool passed = true;
  auto fn = DiscreteRealFunction(X_R*X_R-1, vec2(0, 1), 50);
  auto padded = fn.two_sided_zero_padding(64);

  passed &= assertEqual_UT(padded.samples(), 64);

  for (int i = 0; i < 7; ++i)
  {
    passed &= assertNearlyEqual_UT(padded[i], 0.f);
    passed &= assertNearlyEqual_UT(padded[63-i], 0.f);
  }

  for (int i = 0; i < 50; ++i)
    passed &= assertNearlyEqual_UT(padded[i+7], fn[i]);

  return passed;
}


inline bool gaborTest()
{

  auto fn = DiscreteRealFunction(X_R*X_R-1, vec2(0, 1), 100);
  auto gab = DiscreteGaborTransform(10).transform(fn);

  return true;

}


inline bool quaternionTest() {

  bool passed = true;
  auto i = Quaternion::i();
  auto j = Quaternion::j();
  auto one = Quaternion::one();
  auto q = Quaternion(.5, .5, .5, .5);
  auto ii = i*i;
  passed &= assertNearlyEqual_UT(q.norm(), 1.f);
  passed &= assertNearlyEqual_UT(i*i, j*j);
  passed &= assertNearlyEqual_UT(i*i, -one);
  passed &= assertNearlyEqual_UT(i*j, -j*i);

  return passed;
}

inline UnitTestResult discreteFuncTests__all()
{
	UnitTestResult result;
	result.runTest(fftInverseTest);
	result.runTest(fftRealDomainSymmetryDFTTest);
	result.runTest(paddingTest);
	result.runTest(gaborTest);
	result.runTest(quaternionTest);

	return result;

}
