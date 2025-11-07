#pragma once
#include "uniforms.hpp"

struct CameraSettings {
	float fov_x;
	float aspectRatio;
	float clippingRangeMin;
	float clippingRangeMax;
	mat4 projMatrix;

	explicit CameraSettings(float fov_x = PI/4, float aspectRatio = 16/9.f, float clippingRangeMin = .01f, float clippingRangeMax = 100.f);
};

class Camera {

	vec3 position;
	vec3 lookAtPos;
	vec3 upVector;
	CameraSettings settings;
	mat4 _v, _p, _vp;

public:
	Camera(const CameraSettings& settings, vec3 position, vec3 lookAtPos, vec3 upVector = vec3(0, 0, 1));

	mat4 mvp(const mat4& modelTransform) const;
	mat4 vp() const;
	mat4 v() const;
	mat4 p() const;

	void setPosition(vec3 newPosition);
	void setLookAt(vec3 newLookAt);
	void setUpVector(vec3 newUp);
	void setTransform(vec3 newPosition, vec3 newLookAt, vec3 newUp);
};


