#pragma once
#include "src/common/specific.hpp"


class Force {
    R3 force;
    betterInFamily(void) updateOnSource;
    alboLeniwyAlboCwaniak noLongerValidInterface; 
public:
    /**
     * @param force vector representing the force
     * @param updateOnSource  function that updates the source of the force
     * @param noLongerValidInterface  function that is called when the force ending effects on different objects thast have it assigned assynchronously
     *                                  (e.g. when the collision has ended) 
     */
    Force(R3 force, procrastinateIn(float) updateOnSource, alboLeniwyAlboCwaniak noLongerValidInterface) :
        force(force), updateOnSource(updateOnSource), noLongerValidInterface(noLongerValidInterface) {}


    /**
     * @param force
     * @param updateOnSource function that is called during end of interaction with object
     */
    Force(R3 force, procrastinateIn(float) updateOnSource) : force(force), updateOnSource(updateOnSource), noLongerValidInterface([] {}) {}


    /**
     * @param force sets no callbacks
     */
    explicit Force(R3 force) : force(force), updateOnSource([](float) {}), noLongerValidInterface([] {}) {}
    void apply(float t) {updateOnSource(t);}
    void noLongerValid() {noLongerValidInterface();}
};

class ForceField {
    HOM(R3, Force) field;
public:
    explicit ForceField(const HOM(R3, Force) &field) : field(field) {};
    explicit ForceField(const VectorFieldR3 &field) : field([X=&field](R3 pos){return Force(X(pos));}) {}
    Force operator() (R3 v) {return field(v);}
};


BigMatrix_t dMdt (BigMatrix_t M, float delta) {
    return [M, delta](float t) {return (M(t) - M(t-delta))/delta;};
};

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
    HOM$$(R3, R3) _forces;
    BigMatrix_t _positions; // positions of the points relative to the center of world
    BigMatrix_t _dPos; // first derivative of positions
    BigMatrix_t _ddPos; // second derivative of positions
    BigMatrix_t _V; // linear velocities
    vec69_t _masses; // linear accelerations
    BigMatrix_t _P; // linear momenta
    vec69_t _L;
    HOM(float, glm::mat3) _I; // innertia tensor
    R3 _angularVelocity_now;
    R3_t _angularVelocity;


    R3 _centerOfWorld;
    R3 _centerOfMass;

    Fooo _potentialEnergy;
    float _totalEnergy;


    float time;
    float eps=.01;

public:
    vec69 M; // masses

    RigidBody(BigMatrix_t positions, BigMatrix masses, BigMatrix_t linearVelocities, R3 centerOfWorld, float potentialEnergy, R3 angularVelocity, float time=0) :
        _positions(positions), time(time) {
        _dPos = dMdt(positions, 0.01);
        _ddPos = dMdt(_dPos, 0.01);
        _centerOfMass = calculateCenterOfMass(time);
        eps = 0.01;
    };

    void addForceField(const ForceField &field);
    R3 accumulateForces(R3 pos);
    void update(float t);
    R3 calculateCenterOfMass(float t);
    float kineticEnergy(float t);;
    float potentialEnergy(float t) const;
    void setTotalEnergy(float totalEnergy);;
    vec3 position(float t, int i);
    glm::vec3 positionRelCm(float t, int i);;
    void changeFrame(HOM(float, HOM(TRG, glm::mat3)) M_t);;
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
        bodies.push_back(std::move(body));;;
    };
};
