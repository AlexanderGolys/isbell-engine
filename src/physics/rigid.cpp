#include "rigid.hpp"


using namespace glm;
using std::vector, std::string, std::shared_ptr, std::unique_ptr, std::pair, std::make_unique, std::make_shared;

RigidBodyTriangulated2D::RigidBodyTriangulated2D(std::shared_ptr<WeakSuperMesh> mesh, vec3 angularVelocity,  vec3 linearVelocity,  vec3 angularAcceleration,  vec3 linearAcceleration):
mesh(mesh), angularVelocity(angularVelocity), linearVelocityCM(linearVelocity), angularAcceleration(angularAcceleration), linearAccelerationCM(linearAcceleration) {
	calculateCenterOfMass();
	approximateInertiaTensorCM();
}

void RigidBodyTriangulated2D::update(float dt) {
	angularVelocity += angularAcceleration*dt;
	linearVelocityCM += linearAccelerationCM*dt;
	mesh->deformWithAmbientMap(SpaceAutomorphism::translation(linearVelocityCM*dt));
	centerOfMass += linearVelocityCM*dt;
	if (norm(angularVelocity) > 0)
		mesh->deformWithAmbientMap(SpaceAutomorphism::rotation(normalise(angularVelocity), norm(angularVelocity)*dt, centerOfMass));
}

void RigidBodyTriangulated2D::approximateInertiaTensor(vec3 p) {
	mat3 result     = mat3(0);
	float totalArea = 0;
	for (auto& id: mesh->getPolyGroupIDs())
		for (auto& tr: mesh->getTriangles(id)){
			vec3 center = tr.center();
			vec3 r      = center - p;
			result += mat3(r.y*r.y + r.z*r.z, -r.x*r.y, -r.x*r.z,
						   -r.x*r.y, r.x*r.x + r.z*r.z, -r.y*r.z,
						   -r.x*r.z, -r.y*r.z, r.x*r.x + r.y*r.y)*tr.area();
			totalArea += tr.area();
		}
	I = result/totalArea;
}

void RigidBodyTriangulated2D::calculateCenterOfMass() {
	vec3 sum  = vec3(0);
	float totalArea = .01;
	for (auto& id: mesh->getPolyGroupIDs())
		for (auto& tr: mesh->getTriangles(id)){
			vec3 center = tr.center();
			sum += center*tr.area();
			totalArea += tr.area();
		}
	centerOfMass = sum/totalArea;
}




RollingBody::RollingBody(SmoothParametricCurve boundary, SmoothParametricCurve floor, vec3 polarConeCenter, vec3 gravity, int n, int m, float pipe_r)
		: mesh(make_shared<WeakSuperMesh>()), boundary(boundary), floor(floor), gravity(gravity),
boundaryID(randomID()), floorID(randomID()), centerID(randomID()){


	mesh->addUniformSurface(polarCone(boundary, polarConeCenter), n, m, boundaryID);
	centerOfMass = mesh->centerOfMass();
	floormesh = make_shared<WeakSuperMesh>(floor.pipe(pipe_r, false), n, m, floorID);
	centermesh = make_shared<WeakSuperMesh>(icosphere(.02, 2, centerOfMass, centerID));
	I_center = I_0();
}

float RollingBody::rollingPathLen(float phi) {
	return boundary.length(polarParam, polarParam + phi, 10);
}

void RollingBody::roll(float dt) {
	float angle = angularVelocity*dt;
	float pathLen = rollingPathLen(angle);
	distanceTravelled += pathLen;
	polarParam += angle;
	vec3 contact = contactPoint();
	vec3 bd = boundary(polarParam);
	shift(contact - bd);
	float rotAngle = rotationAngle(polarParam-angle, polarParam);
	rotate(rotAngle);
}

float RollingBody::rotationAngle(float d1, float d2)
{
//	vec3 t2 = boundary.tangent(d2);
	float a = angle(boundary.tangent(d1), boundary.tangent(d2));
	if (a < -100 || a > 260)
		return 0;
	return a;
}



float RollingBody::calculate_angularAcc() {
	return norm(cross(gravity, r_contact())/I_zz(contactPoint()));
}

void RollingBody::update_omega(float dt) {
	float inc = calculate_angularAcc()*dt;
	if (inc < -260 || inc > 260)
		return;
	angularVelocity += calculate_angularAcc()*dt;
}

float RollingBody::I_zz(vec3 p) {
	return I_center + norm2(p-centerOfMass);
}

float RollingBody::I_0(){
	float result = 0;
	float totalArea = 0;
	for (auto& tr: mesh->getTriangles(boundaryID)){
		vec3 center = tr.center();
		vec3 r      = center - centerOfMass;
		result += tr.area()*(r.x*r.x + r.y*r.y);
		totalArea += tr.area();
	}
	return result/totalArea;
}

void RollingBody::shift(vec3 v) {
	centerOfMass += v;
	mesh->shift(v, boundaryID);
	boundary = boundary.shift(v);
}

void RollingBody::rotate(float angle) {
	SpaceAutomorphism rot = SpaceAutomorphism::rotation(angle).applyWithShift(contactPoint());
	centerOfMass = rot(centerOfMass);
	mesh->deformWithAmbientMap(boundaryID, rot);
	boundary = boundary.precompose(rot);
}

void RollingBody::step(float dt) {
	update_omega(dt);
	roll(dt);
}

void RollingBody::shift(vec3 v, WeakSuperMesh &mesh, WeakSuperMesh &centermesh) {
	centerOfMass += v;
	mesh.shift(v, boundaryID);
	centermesh.shift(v, centerID);
	boundary = boundary.shift(v);
}
void RollingBody::rotate(float angle, WeakSuperMesh &mesh, WeakSuperMesh &centermesh) {
	SpaceAutomorphism rot = SpaceAutomorphism::rotation(angle).applyWithShift(contactPoint());
	centerOfMass = rot(centerOfMass);
	mesh.deformWithAmbientMap(boundaryID, rot);
	centermesh.deformWithAmbientMap(centerID, rot);
	boundary = boundary.precompose(rot);
}
void RollingBody::step(float dt, WeakSuperMesh &mesh, WeakSuperMesh &centermesh) {
	update_omega(dt);
	roll(dt, mesh, centermesh);
}

void RollingBody::roll(float dt, WeakSuperMesh &mesh, WeakSuperMesh &centermesh) {
	if (dt < 1e-6)
		return;
	float angle = angularVelocity*dt;
	float pathLen = rollingPathLen(angle);
	distanceTravelled += pathLen;
	polarParam += angle;
	vec3 contact = contactPoint();
	vec3 bd = boundary(polarParam);
	shift(contact - bd, mesh, centermesh);
	float rotAngle = rotationAngle(polarParam-angle, polarParam);
	rotate(rotAngle, mesh, centermesh);
}

RigidBody3D::RigidBody3D( vec3 cm, mat3 I, vec3 angularVelocity, vec3 linearVelocity, vec3 angularAcceleration, vec3 linearAcceleration, float mass, mat3 R)
: centerOfMass(cm), I(I), angularVelocity(angularVelocity), linearVelocityCM(linearVelocity),
angularAcceleration(angularAcceleration), linearAccelerationCM(linearAcceleration), mass(mass), R(mat3(1)) {}

void RigidBody3D::rotateAroundCM(mat3 M, WeakSuperMesh &mesh) {
	mesh.deformWithAmbientMap(SpaceAutomorphism::linear(M).applyWithShift(centerOfMass));

}

vec3 RigidBody3D::angularMomentum(vec3 rotationCenter) {
	return I_p(rotationCenter)*angularVelocity;
}

float RigidBody3D::kineticEnergy(vec3 rotationCenter) {
	return 0.5f * dot(angularMomentum(rotationCenter), angularVelocity);
}

mat3 RigidBody3D::I_p(vec3 p) {
	vec3 d = p-centerOfMass;
	return I + mass*mat3(d.y*d.y + d.z*d.z, -d.x*d.y, -d.x*d.z, -d.x*d.y, d.x*d.x + d.z*d.z, -d.y*d.z, -d.x*d.z, -d.y*d.z, d.x*d.x + d.y*d.y);
}


mat3 RigidBody3D::spin() {
	return spinTensor(angularVelocity);
}

vec3 RigidBody3D::alpha() {
	mat3 I1 = I_rotated(R);
	return -inverse(I1) * cross(angularVelocity, I1*angularVelocity);
}

void RigidBody3D::freeMotion(float dt, WeakSuperMesh &mesh) {
	mat3 R0 = R;
	mat3 I_n = R*I*transpose(R);
	R += spinTensor(inverse(I_n)*angularMomentum(centerOfMass))*R*dt;
	shift(linearVelocityCM*dt, mesh);

	linearVelocityCM += linearAccelerationCM*dt;
	rotateAroundCM(R*inverse(R0), mesh);
//	linearVelocityCM += linearAccelerationCM*dt;
}

void RigidBody3D::shift(vec3 v, WeakSuperMesh &mesh) {
	mesh.shift(v);
	centerOfMass += v;
}

std::pair<RigidBody3D, WeakSuperMesh> boxRigid(vec3 size, vec3 center, float mass, vec3 velocity, vec3 angularVelocity, vec3 angularAcceleration, vec3 linearAcceleration,
	const mat3 &R) {
	mat3 I =size.x*size.y*size.z*mass* mat3(1.f/12.f*(size.y*size.y + size.z*size.z), 0, 0,
				  0, 1.f/12.f*(size.x*size.x + size.z*size.z), 0,
				  0, 0, 1.f/12.f*(size.x*size.x + size.y*size.y));
	return {RigidBody3D(center, I, angularVelocity, velocity, angularAcceleration, linearAcceleration, mass, R),
			box(size, center, randomID())};
}
