#include "scene3D.hpp"

Camera::Camera(vec3 position, vec3 lookAt, vec3 upVector, float fov, float aspectRatio, float nearClip, float farClip)
: position(position), lookAt(lookAt), upVector(upVector), fov(fov), aspectRatio(aspectRatio), nearClip(nearClip), farClip(farClip) {}

mat4 Camera::projectionMatrix() const {
	return glm::perspective(fov, aspectRatio, nearClip, farClip);
}

mat4 Camera::viewMatrix() const {
	return glm::lookAt(position, lookAt, upVector);
}

mat4 Camera::vp() const {
	return projectionMatrix() * viewMatrix();
}

mat4 Camera::mvp(const mat4& modelTransform) const {
	return vp() * modelTransform;
}
