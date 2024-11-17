#include "rigid.hpp"
#include <glm/glm.hpp>

using namespace glm;
using std::vector;

glm::vec3 RigidBody::calculateCenterOfMass(float t) {
    auto res = (_masses*_positions)(t);
    return vec3(res[0], res[1], res[2]);
}

float RigidBody::kineticEnergy(float t) {
    return float(BigMatrix(_angularVelocity(t)).transpose()*_I(t)*_angularVelocity(t)) + 0.5*float((_masses*_V*_V)(t));
}

float RigidBody::potentialEnergy(float t) const {
    return _totalEnergy - kineticEnergy(t);
}

void RigidBody::setTotalEnergy(float totalEnergy) {
    _totalEnergy = totalEnergy;
}

void RigidBody::addForceField(const ForceField &field) {
    throw std::logic_error("Not implemented");
}

glm::vec3 RigidBody::accumulateForces(glm::vec3 pos) {
    throw std::logic_error("Not implemented");
}

void RigidBody::update(float t) {
    throw std::logic_error("Not implemented");
}

vec3 RigidBody::position(float t, int i) {
    return cast_vec3(_positions(t)[i]);
}

glm::vec3 RigidBody::positionRelCm(float t, int i) {
    return cast_vec3(_positions(t)[i]);
}

void RigidBody::changeFrame(function<function<glm::mat3(const IndexedTriangle &)>(float)> M_t) {
    // pos = [M_t](float time, TRG tr, HOM(float, R3) relativeTo) {return relativeTo(t) + M_t(time).inv() * pos(t, tr);};
    throw std::runtime_error("changeFrame not done yet");
}

glm::vec3 RigidBody::linearVelocity(float t, int i) {
    return _dPos(t)[i] + glm::cross(angularVelocity(t), cast_vec3((_positions(t)[i]));
}

glm::vec3 RigidBody::linearMomentum(float t, int i) {
    return (linearVelocity(t, i) + glm::cross(angularVelocity(t), position(t, i)))*_masses[i];
}


glm::vec3 RigidBody::totalLinearMomentum(float t) {
    return cast_vec3(_masses*_V(t)) + glm::cross(BigMatrix(_masses(t))*angularVelocity(t), positions(t));
}
//
// glm::vec3 RigidBody::angularMomentum(float t) {
//     return InertiaTensor(t)*angularVelocity(t);
// }
//
// glm::vec3 RigidBody::angularVelocity(float t) {
//     return _angularVelocity(t);
// }
//
// void RigidBody::calculateInertiaTensor() {
//     return _I =
// }
//
// void RigidBody::calculateMomenta() {
//     _L = _I*angularVelocity
//
// }
