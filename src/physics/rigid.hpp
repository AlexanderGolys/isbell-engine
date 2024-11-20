#pragma once
#include "src/common/specific.hpp"


class RigidBody;

class Force {
    R3 force;
    betterInFamily(void) updateOnSource;
    alboLeniwyAlboCwaniak noLongerValidInterface; 
public:
    Force(R3 force, const procrastinateIn(float) &updateOnSource, const alboLeniwyAlboCwaniak &noLongerValidInterface) : force(force), updateOnSource(updateOnSource), noLongerValidInterface(noLongerValidInterface) {}
    virtual ~Force() {noLongerValidInterface();}
    Force(R3 force, const procrastinateIn(float) &updateOnSource) : force(force), updateOnSource(updateOnSource), noLongerValidInterface([] {}) {}
    explicit Force(R3 force) : force(force), updateOnSource([](float) {}), noLongerValidInterface([] {}) {}

    virtual void apply(float t) {updateOnSource(t);}
    virtual void noLongerValid() {noLongerValidInterface();}
    virtual void setForce(R3 force) {this->force = force;}
    virtual void setUpdateCallback (const procrastinateIn(float) &updateOnSource) {this->updateOnSource = updateOnSource;}
    virtual void setDestructionCallback (const alboLeniwyAlboCwaniak &noLongerValidInterface) {this->noLongerValidInterface = noLongerValidInterface;}
    virtual R3 getForceVector() {return force;}
};

enum class ForceType {
    COLLISION, CONST_FIELD, STATIC_VF, DYNAMIC_VF, CUSTOM
};

class RigidForce : public Force {
    AffineLine* actionAxis = nullptr;
    R3* pointOfAction = nullptr;
    int* indexOfClosestPointInBody = nullptr;
    ForceType type;
    std::vector<alboLeniwyAlboCwaniak> destroyCallbacksRegistered;
public:
    explicit RigidForce(ForceType type) : Force(R3(0)), type(type) {}
    virtual void addForceToBody(RigidBody &body, const alboLeniwyAlboCwaniak &bodyChangeCallback, const alboLeniwyAlboCwaniak &forceChangeCallback);
};


class RigidCollisionForce : public RigidForce {
    R3 normal;
    float magnitude;
    std::weak_ptr<RigidBody> body1, body2;
    std::weak_ptr<alboLeniwyAlboCwaniak> notifierBody1, notifierBody2;
public:
    RigidCollisionForce(const R3 &normal, float magnitude, const std::weak_ptr<RigidBody> &body1, const std::weak_ptr<RigidBody> &body2);
    virtual void movementEvent(int notifierBodyIndex);
    void setForce(vec3 force) override {normal = force; magnitude = length(force);}
    vec3 getForceVector() override {return normal*magnitude;}
};

class RigidForceFieldConst : public RigidForce {
    R3 force;
public:
    explicit RigidForceFieldConst(const R3 &field) : RigidForce(ForceType::CONST_FIELD), force(field) {};
    void setForce(vec3 force) override {throw std::logic_error("Constant field cannot be modified by definition");}
    vec3 getForceVector() override {return force;}
};


BigMatrix_t dMdt (const BigMatrix_t &M, float delta);

inline BigMatrix_t operator* (const BigMatrix_t &M, float f) { return [M, f](float t) {return M(t)*f;}; };
// operator* BigMatrix_t(const BigMatrix_t &M, const BigMatrix_t &M2) { return [M, M2](float t) {return M(t)*M2(t);}; };
inline BigMatrix_t operator+ (const BigMatrix_t & M, const BigMatrix_t & M2) { return [M, M2](float t) {return M(t)+M2(t);}; };
// BigMatrix_t operator- BigMatrix_t(const BigMatrix_t & M, const BigMatrix_t & M2) { return [M, M2](float t) {return M(t)-M2(t);}; };

Fooo dot(const vec69_t &v1, const vec69_t &v2) {
    return [v1, v2](float t) {
        float result = 0;
        for (int i = 0; i < v1(t).size(); i++)
            result += v1(t)[i] * v2(t)[i];
        return result;
    };
};

inline vec3 cast_vec3(const vector<float> &v) {return vec3(v[0], v[1], v[2]);};


class RigidBody {
    vector<RigidForce> activeForces;
    BigMatrix _positions; // positions of the points relative to the center of world
    BigMatrix_hm _dPos; // first derivative of positions
    BigMatrix_hm _ddPos; // second derivative of positions
    BigMatrix_hm _V; // linear velocities
    vec69 _masses; // linear accelerations
    BigMatrix_hm _P; // linear momenta
    vec69_hm _L;
    mat3_hm _I; // innertia tensor
    R3_hm _angularVelocity;
    R3 _centerOfWorld;
    R3 _centerOfMass;
    float_hm _potentialEnergy;
    float_hm _totalEnergy;
	float_hm _kineticEnergy;
    float time;
    float eps=.01;

public:
    RigidBody(const BigMatrix& positions, const vec69& masses, const BigMatrix& linearVelocities, R3 centerOfWorld, R3 angularVelocity, float time=0);;
	RigidBody(const BigMatrix& positions, const vec69& masses, const BigMatrix& linearVelocities, R3 angularVelocity, float time=0);;

    alboLeniwyAlboCwaniak addForceField(const RigidForce &field);
    R3 accumulateForces(R3 pos);
    R3 calculateCenterOfMass();
    float  calculateKineticEnergy();
    float potentialEnergy() const;
    void setTotalEnergy(float totalEnergy);
    vec3 position(float t, int i);
    vec3 positionRelCm(float t, int i);
    void changeFrame(HOM(float, HOM(TRG, glm::mat3)) M_t);
    R3 linearVelocity(float t, int i);
    R3 linearMomentum(float t, int i);
    R3 totalLinearMomentum(float t);
    R3 angularMomentum(float t);
    R3 angularVelocity(float t);

    void calculateInertiaTensor();
    void calculateMomenta();
};

class PhysicalEnvironment {
    std::vector<RigidBody> bodies;
    public:
    void addRigidBody(const RigidBody &body) {
        bodies.push_back(body);
    };
};
