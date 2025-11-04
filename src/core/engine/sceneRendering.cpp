#include "sceneRendering.hpp"


Camera::Camera() {
	this->lookAtFunc = make_shared<SmoothParametricCurve>(SmoothParametricCurve::constCurve(vec3(0)));
	this->up = [](float t) {
		return vec3(0, 0, 1);
	};
	this->fov_x = 45.0f;
	this->aspectRatio = 16.0f / 9.0f;
	this->clippingRangeMin = 0.1f;
	this->clippingRangeMax = 100.0f;
	this->trajectory = make_shared<SmoothParametricCurve>(SmoothParametricCurve::constCurve(vec3(2.0f, 3.0f, 1.0f)));
	this->moving = false;
	this->projectionMatrix = perspective(fov_x, aspectRatio, clippingRangeMin, clippingRangeMax);
}

Camera::Camera(vec3 position, vec3 lookAtPos, vec3 upVector, float fov_x, float aspectRatio, float clippingRangeMin, float clippingRangeMax) {
	this->lookAtFunc = make_shared<SmoothParametricCurve>(SmoothParametricCurve::constCurve(lookAtPos));
	this->up = [upVector](float t) {
		return upVector;
	};
	this->fov_x = fov_x;
	this->aspectRatio = aspectRatio;
	this->clippingRangeMin = clippingRangeMin;
	this->clippingRangeMax = clippingRangeMax;
	this->trajectory = make_shared<SmoothParametricCurve>(SmoothParametricCurve::constCurve(position));
	this->moving = false;
	this->projectionMatrix = perspective(fov_x, aspectRatio, clippingRangeMin, clippingRangeMax);
}

Camera::Camera(float radius, float speed, float height, vec3 lookAtPos, vec3 upVector, float fov_x, float aspectRatio, float clippingRangeMin, float clippingRangeMax)
: Camera(make_shared<SmoothParametricCurve>([speed, radius, height](float t) {
	return vec3(radius * cos(speed * t), radius * sin(speed * t), height);
}), lookAtPos, upVector, fov_x, aspectRatio, clippingRangeMin, clippingRangeMax) {}

Camera::Camera(const shared_ptr<SmoothParametricCurve>& trajectory, vec3 lookAtPos, vec3 upVector, float fov_x, float aspectRatio, float clippingRangeMin, float clippingRangeMax) {
	this->lookAtFunc = make_shared<SmoothParametricCurve>(SmoothParametricCurve::constCurve(lookAtPos));
	this->up = [upVector](float t) {
		return upVector;
	};
	this->fov_x = fov_x;
	this->aspectRatio = aspectRatio;
	this->clippingRangeMin = clippingRangeMin;
	this->clippingRangeMax = clippingRangeMax;
	this->trajectory = trajectory;
	this->moving = true;
	this->projectionMatrix = perspective(fov_x, aspectRatio, clippingRangeMin, clippingRangeMax);
}

Camera::Camera(const shared_ptr<SmoothParametricCurve>& trajectory, const shared_ptr<SmoothParametricCurve>& lookAtPos, vec3 upVector, float fov_x, float aspectRatio,
			   float clippingRangeMin, float clippingRangeMax) {
	this->lookAtFunc = lookAtPos;
	this->up = [upVector](float t) {
		return upVector;
	};
	this->fov_x = fov_x;
	this->aspectRatio = aspectRatio;
	this->clippingRangeMin = clippingRangeMin;
	this->clippingRangeMax = clippingRangeMax;
	this->trajectory = trajectory;
	this->moving = true;
	this->projectionMatrix = perspective(fov_x, aspectRatio, clippingRangeMin, clippingRangeMax);
}

Camera::Camera(const shared_ptr<SmoothParametricCurve>& trajectory, const shared_ptr<SmoothParametricCurve>& lookAtPos, const std::function<vec3(float)>& upVector, float fov_x,
			   float aspectRatio, float clippingRangeMin, float clippingRangeMax) {
	this->lookAtFunc = lookAtPos;
	this->up = upVector;
	this->fov_x = fov_x;
	this->aspectRatio = aspectRatio;
	this->clippingRangeMin = clippingRangeMin;
	this->clippingRangeMax = clippingRangeMax;
	this->trajectory = trajectory;
	this->moving = true;
	this->projectionMatrix = perspective(fov_x, aspectRatio, clippingRangeMin, clippingRangeMax);
}

vec3 Camera::position(float t) const {
	return trajectory->operator()(t);
}

vec3 Camera::lookAtPoint(float t) const {
	return lookAtFunc->operator()(t);
}

vec3 Camera::upVector(float t) const {
	return up(t);
}


mat4 Camera::viewMatrix(float t) const {
	return lookAt(position(t), lookAtPoint(t), upVector(t));
}

mat4 Camera::vp(float t) const {
	return projectionMatrix * viewMatrix(t);
}

mat4 Camera::mvp(float t, const mat4& modelTransform) const {
	return vp(t) * modelTransform;
}

