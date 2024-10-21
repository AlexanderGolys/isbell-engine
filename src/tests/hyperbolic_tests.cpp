
#include "src/common/hyperbolic.hpp"
#include <cassert>
#include <iostream>

using namespace std;

bool isAdditionH(Mob m) {
  return abs(m.a - ONE) < 1e-6 && abs(m.d - ONE) < 1e-6 && abs(m.c) < 1e-6;
}

bool isMultiplicationH(Mob m) {
  return abs(m.b) < 1e-6 && abs(m.d - ONE) < 1e-6 && abs(m.c) < 1e-6;
}

bool isInAutH(Mob m) {
  return abs(m.a.im()) < 0.0001 && abs(m.b.im()) < 0.0001 && abs(m.c.im()) < 0.0001 && abs(m.d.im()) < 0.0001;
}
bool isInAutD(Mob m) {
  return abs(m.a.conj() - m.d) < 0.0001 && abs(m.b.conj() - m.c) < 0.0001;
}



void cayleyTransformTest()
  {
     auto M = CayleyTransform;
     assert(M.a == ONE);
     Complex z = Complex(1, 3);
     assert(M.mobius(z) == (z - I)/(z + I));
     cout << "Cayley transform test passed" << endl;
  }

bool nearlyEqual(Complex a, Complex b) {
  return abs(a-b)<1e-6;
}

void MeromorphismsTests()
{
    auto exp_f = EXP;
    auto sqrt_f = SQRT;

    assert( nearlyEqual(sqrt_f(ONE*4), ONE*2) );
    assert(nearlyEqual(exp_f.compose(sqrt_f)(9), exp(1)*exp(1)*exp(1)) );
    assert(nearlyEqual(sqrt_f.compose(exp_f)(2), exp(1)));
    assert(nearlyEqual(sqrt_f.inv(I), -1) );
    assert(nearlyEqual(Id(I*8 + 2), I*8 + 2) );
    assert(nearlyEqual(Biholomorphism::linear(I, I)(I*3), I - 3) );
    cout << "Composition, inverse tests passed" << endl;

}


void FuchsianGroupTest()
{
  auto G1 = Ga(1, 4);
  auto G2 = Gm(1, 2);

  auto addMatricesD = G1.generateElementsD(13);
  auto multMatricesD = G2.generateElementsD(5);
  auto addMatricesH = G1.generateElementsH(13);
  auto multMatricesH = G2.generateElementsH(5);

  assert(addMatricesD.size() == 13);
  assert(multMatricesD.size() == 5);
  assert(addMatricesH.size() == 13);
  assert(multMatricesH.size() == 5);

  for (auto M : addMatricesD)
  {
    assert(isInAutD(M));
  }

  for (auto M : multMatricesD)
  {
    assert(isInAutD(M));
  }

  for (auto M : addMatricesH)
  {
    assert(isInAutH(M));
    assert(isAdditionH(M));
  }

  for (auto M : multMatricesH)
  {
    assert(isInAutH(M));
    assert(isMultiplicationH(M));
  }

  cout << "Additive and multiplicative Fuchsian groups seem to work." << endl;


}


  int main(void)
  {
    cayleyTransformTest();
    MeromorphismsTests();
    FuchsianGroupTest();
    return 0;
  }