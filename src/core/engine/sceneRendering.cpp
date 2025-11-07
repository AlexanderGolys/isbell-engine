#include "sceneRendering.hpp"

#include "SDFObjects.hpp"
#include "SDFObjects.hpp"
#include "SDFObjects.hpp"
#include "SDFObjects.hpp"
#include "SDFObjects.hpp"
#include "SDFObjects.hpp"


CameraSettings::CameraSettings(float fov_x, float aspectRatio, float clippingRangeMin, float clippingRangeMax)
: fov_x(fov_x), aspectRatio(aspectRatio), clippingRangeMin(clippingRangeMin), clippingRangeMax(clippingRangeMax) {
	projMatrix = perspective(fov_x, aspectRatio, clippingRangeMin, clippingRangeMax);
}

Camera::Camera(const CameraSettings& settings, vec3 position, vec3 lookAtPos, vec3 upVector)
: position(position), lookAtPos(lookAtPos), upVector(upVector), settings(settings), _v(0), _p(0), _vp(0) {
	_v = lookAt(position, lookAtPos, upVector);
	_p = settings.projMatrix;
	_vp = _p * _v;
}

mat4 Camera::v() const {
	return _v;
}

mat4 Camera::p() const {
	return _p;
}

void Camera::setPosition(vec3 newPosition) {
	position = newPosition;
	_v = lookAt(position, lookAtPos, upVector);
	_vp = _p * _v;
}

void Camera::setLookAt(vec3 newLookAt) {
	lookAtPos = newLookAt;
	_v = lookAt(position, lookAtPos, upVector);
	_vp = _p * _v;
}

void Camera::setUpVector(vec3 newUp) {
	upVector = newUp;
	_v = lookAt(position, lookAtPos, upVector);
	_vp = _p * _v;
}

void Camera::setTransform(vec3 newPosition, vec3 newLookAt, vec3 newUp) {
	position = newPosition;
	lookAtPos = newLookAt;
	upVector = newUp;
	_v = lookAt(position, lookAtPos, upVector);
	_vp = _p * _v;
}

mat4 Camera::vp() const {
	return _vp;
}


mat4 Camera::mvp(const mat4& modelTransform) const {
	return vp() * modelTransform;
}

