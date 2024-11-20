#include "rigid.hpp"
#include <glm/glm.hpp>
#include <glm/detail/_vectorize.hpp>

using namespace glm;
using std::vector;

glm::vec3 RigidBody::calculateCenterOfMass() {
    auto m =  BigMatrix(_masses).transpose()*_positions;
    return vec3(m);
}

float RigidBody::kineticEnergy() {
	if (!_angularVelocity.has_value()) throw std::logic_error("No angular velocity");
    return float(_angularVelocity*BigMatrix(_I.or_else(mat3(0)))*_angularVelocity) + 0.5*float((BigMatrix(_masses)*_V).transpose()*_V);
}

float RigidBody::potentialEnergy() const {
    return _totalEnergy - kineticEnergy();
}

void RigidBody::setTotalEnergy(float totalEnergy) {
    _totalEnergy = totalEnergy;
}

void RigidBody::addForceField(const ForceField &field) {
    throw std::logic_error("Not implemented");
}

RigidCollisionForce::RigidCollisionForce(const glm::vec3 &normal, const float magnitude, const std::weak_ptr<RigidBody> &body1,
        const std::weak_ptr<RigidBody> &body2): RigidForce(ForceType::COLLISION),
                                                normal(normal),
                                                magnitude(magnitude),
                                                body1(body1),
                                                body2(body2) {}

void RigidCollisionForce::movementEvent(int notifierBodyIndex) {
    if (notifierBodyIndex == 0 && !body2.expired() && !notifierBody2.expired()) {
        notifierBody2.lock()->operator();
        return;
    }
    if (notifierBodyIndex == 1 && !body1.expired() && !notifierBody1.expired()) {
        notifierBody1.lock()->operator();
        return;
    }
    noLongerValid();
    this->~RigidCollisionForce();
}

function<BigMatrix(float)> dMdt(const function<BigMatrix(float)> &M, float delta) {
    return [M, delta](float t) {return (M(t) - M(t-delta))/delta;};
}

RigidBody::RigidBody(const BigMatrix &positions, const vec69 &masses, const BigMatrix &linearVelocities, glm::vec3 centerOfWorld, glm::vec3 angularVelocity, float time):
	_positions(positions),
	time(time) {
			_masses       = masses;
			_dPos         = HM_NO;
			_ddPos        = HM_NO;
			_angularVelocity = angularVelocity;
			_centerOfMass = calculateCenterOfMass();
			eps           = 0.01;
			_V			  =	linearVelocities;
			_P            = calculateMomentum();
			_L            = calculateAngularMomentum();
			_I            = calculateInertiaTensor();
			_potentialEnergy = calculatePotentialEnergy();
			_kineticEnergy = calculateKineticEnergy();
			_totalEnergy  = _potentialEnergy.has_value() && _kineticEnergy.has_value() ? _potentialEnergy.or_else(0.f) + _kineticEnergy.or_else(0.f) : HM_NO;

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


vec3 RigidBody::totalLinearMomentum(float t) {
    return cast_vec3(_masses*_V(t)) + glm::cross(BigMatrix(_masses(t))*angularVelocity(t), _positions(t));
}

mat3 RigidBody::calculateInertiaTensor() { throw std::logic_error("Not implemented"); }
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
