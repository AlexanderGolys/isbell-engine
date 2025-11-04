#pragma once
#include "smoothParametric.hpp"

class Camera {
public:
	float fov_x;
	float aspectRatio;
	float clippingRangeMin;
	float clippingRangeMax;
	bool moving;
	shared_ptr<SmoothParametricCurve> trajectory;
	shared_ptr<SmoothParametricCurve> lookAtFunc;
	std::function<vec3(float)> up;
	mat4 projectionMatrix;

	Camera();

	Camera(vec3 position, vec3 lookAtPos, vec3 upVector = vec3(0, 0, 1), float fov_x = PI / 4, float aspectRatio = 16 / 9.f, float clippingRangeMin = .01f,
		   float clippingRangeMax = 100.f);

	Camera(float radius, float speed, float height, vec3 lookAtPos = vec3(0), vec3 upVector = vec3(0, 0, 1), float fov_x = PI / 4, float aspectRatio = 16 / 9.f,
		   float clippingRangeMin = .01f, float clippingRangeMax = 100.f);

	Camera(const shared_ptr<SmoothParametricCurve>& trajectory, vec3 lookAtPos, vec3 upVector, float fov_x = PI / 4, float aspectRatio = 16 / 9.f, float clippingRangeMin = .01f,
		   float clippingRangeMax = 100.f);

	Camera(const shared_ptr<SmoothParametricCurve>& trajectory, const shared_ptr<SmoothParametricCurve>& lookAtPos, vec3 upVector, float fov_x = PI / 4,
		   float aspectRatio = 16 / 9.f, float clippingRangeMin = .01f, float clippingRangeMax = 100.f);

	Camera(const shared_ptr<SmoothParametricCurve>& trajectory, const shared_ptr<SmoothParametricCurve>& lookAtPos, const std::function<vec3(float)>& upVector,
		   float fov_x = PI / 4, float aspectRatio = 16 / 9.f, float clippingRangeMin = .01f, float clippingRangeMax = 100.f);

	vec3 position(float t) const;
	vec3 lookAtPoint(float t) const;
	vec3 upVector(float t) const;
	mat4 mvp(float t, const mat4& modelTransform) const;
	mat4 viewMatrix(float t) const;
	mat4 vp(float t) const;
};
