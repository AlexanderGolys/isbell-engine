#include "rigid.hpp"


using namespace glm;
using std::vector;

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
