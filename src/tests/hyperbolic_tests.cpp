
#include "src/common/hyperbolic.hpp"
#include <cassert>
#include <iostream>

using namespace glm;

float angle(vec3 a, vec3 b)
{
  return acos(dot(a, b)/(length(a)*length(b)));
}


void rotationsTest()
  {
    float t1 = 0;
    float t2 = 1;

    vec3 n2 = vec3(0, 1, 1);
    vec3 n3 = vec3(1, -4, 1)/100.f;
    vec3 n4 = vec3(0, -2, -2);

    assert(nearlyEqual(rotationMatrix3(0), rotationMatrix3(TAU)));
    assert(nearlyEqual(rotationMatrix3(n2, 0), rotationMatrix3(n2, TAU)));

    mat3 m11 = rotationMatrix3(t1);
    mat3 m12 = rotationMatrix3(n2, t1);
    mat3 m13 = rotationMatrix3(n3, t1);
    mat3 m14 = rotationMatrix3(n4, t1);
    mat3 m21 = rotationMatrix3(t2);
    mat3 m22 = rotationMatrix3(n2, t2);
    mat3 m23 = rotationMatrix3(n3, t2);
    mat3 m24 = rotationMatrix3(n4, t2);

    float ang = angle(e1, m22*e1);
    assert(nearlyEqual(ang, t2));

    std::vector<mat3> all = {m11, m12, m13, m14, m21, m22, m23, m24};

    for (auto m: all)
    {
      assert(nearlyEqual(abs(det(m)), 1));
    }

    std::cout << "Rotation matrix tests passed" << std::endl;
  }


  int main(void)
  {
    rotationsTest();
    return 0;
  }