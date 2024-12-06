#pragma once
#include "solidMeshes.hpp"

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

//Fooo dot(const vec69_t &v1, const vec69_t &v2) {
//    return [v1, v2](float t) {
//        float result = 0;
//        for (int i = 0; i < v1(t).size(); i++)
//            result += v1(t)[i] * v2(t)[i];
//        return result;
//    };
//};

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


class RigidBodyTriangulated2D {
	std::shared_ptr<WeakSuperMesh> mesh;
	vec3 centerOfMass = vec3(0);
	mat3 I = mat3(0);
	vec3 angularVelocity;
	vec3 linearVelocityCM;
	vec3 angularAcceleration;
	vec3 linearAccelerationCM;
public:
	RigidBodyTriangulated2D(std::shared_ptr<WeakSuperMesh> mesh, vec3 angularVelocity,  vec3 linearVelocity,  vec3 angularAcceleration,  vec3 linearAcceleration);

	void update(float dt);

	void setAngularVelocity(const vec3 &omega) {angularVelocity = omega;}
	void setLinearVelocity(const vec3 &v) {linearVelocityCM = v;}
	void setAngularAcceleration(const vec3 &alpha) {angularAcceleration = alpha;}
	void setLinearAcceleration(const vec3 &a) {linearAccelerationCM = a;}
	void addAngularVelocity(const vec3 &omega) {angularVelocity += omega;}
	void addLinearVelocity(const vec3 &v) {linearVelocityCM += v;}
	void addAngularAcceleration(const vec3 &alpha) {angularAcceleration += alpha;}
	void addLinearAcceleration(const vec3 &a) {linearAccelerationCM += a;}

	void approximateInertiaTensor(vec3 p);
	void approximateInertiaTensorCM() {approximateInertiaTensor(centerOfMass);}
	void calculateCenterOfMass();

	vec3 getCm() const {return centerOfMass;}
	mat3 getI() const {return I;}
};

class RollingBody {
	SmoothParametricCurve boundary; // polar b(phi)
	SmoothParametricCurve floor; // arclen f(s)
	vec3 gravity;
	vec3 centerOfMass;
	float angularVelocity = 0;
	float distanceTravelled=0;
	float polarParam=0;
	float I_center;

	PolyGroupID boundaryID;
	PolyGroupID floorID;
	PolyGroupID centerID;

public:
	std::shared_ptr<WeakSuperMesh> mesh;
	std::shared_ptr<WeakSuperMesh> floormesh;
	std::shared_ptr<WeakSuperMesh> centermesh;

	RollingBody(std::shared_ptr<WeakSuperMesh> mesh, SmoothParametricCurve boundary, SmoothParametricCurve floor, vec3 gravity);
	RollingBody(SmoothParametricCurve boundary, SmoothParametricCurve floor, vec3 polarConeCenter, vec3 gravity, int n, int m, float pipe_r);

	float rollingPathLen(float phi);
	void roll(float dt);
	float calculate_angularAcc();
	void update_omega(float dt);
	float I_zz(vec3 p);
	float I_0();
	vec3 contactPoint() {return floor(distanceTravelled);}
	vec3 r_contact() {return centerOfMass-contactPoint();}
	void shift(vec3 v);
	void rotate(float angle);
	void step(float dt);
	float rotationAngle(float d1, float d2);

	void shift(vec3 v, WeakSuperMesh &mesh, WeakSuperMesh &centermesh);
	void rotate(float angle, WeakSuperMesh &mesh, WeakSuperMesh &centermesh);
	void step(float dt, WeakSuperMesh &mesh, WeakSuperMesh &centermesh);
	void roll(float dt, WeakSuperMesh &mesh, WeakSuperMesh &centermesh);


	vec3 getCM() {return centerOfMass;}
	float getAngularVelocity() {return angularVelocity;}
	std::shared_ptr<WeakSuperMesh> getMesh() {return std::make_shared<WeakSuperMesh>(*mesh);}
	std::shared_ptr<WeakSuperMesh> getCenterMesh() {return std::make_shared<WeakSuperMesh>(*centermesh);}
	std::shared_ptr<WeakSuperMesh> getFloorMesh() {return std::make_shared<WeakSuperMesh>(*floormesh);}

};


class RigidBody3D {
	vec3 centerOfMass = vec3(0);
	mat3 I = mat3(0);
	vec3 angularVelocity;
	vec3 linearVelocityCM;
	vec3 angularAcceleration;
	vec3 linearAccelerationCM;
	float mass;
	mat3 R;

public:
	RigidBody3D(vec3 cm, mat3 I, vec3 angularVelocity,
		vec3 linearVelocity,  vec3 angularAcceleration,  vec3 linearAcceleration, float mass, mat3 R);

	void rotateAroundCM(mat3 M, WeakSuperMesh &mesh);
	void rotate(mat3 M, vec3 p, WeakSuperMesh &mesh);
	vec3 angularMomentum(vec3 rotationCenter);
	float  kineticEnergy(vec3 rotationCenter);
	mat3 I_p(vec3 p);
	mat3 spin();
	mat3 dRdt() {return spin()*R;}
	mat3 I_rotated(const mat3 &R) {return R*I*transpose(R);}
	vec3 alpha();
	void freeMotion(float dt, WeakSuperMesh &mesh);
	void shift(vec3 v, WeakSuperMesh &mesh);
	vec3 angularMomentum() {return I*angularVelocity;}

	vec3 getCM() const {return centerOfMass;}
};

std::pair<RigidBody3D, WeakSuperMesh> boxRigid(vec3 size, vec3 center, float mass, vec3 velocity, vec3 angularVelocity, vec3 angularAcceleration, vec3 linearAcceleration, const mat3 &R);
