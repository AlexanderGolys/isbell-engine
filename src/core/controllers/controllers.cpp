#include "controllers.hpp"


CameraController::CameraController(sptr<Camera> camera, std::function<void(sptr<Camera>, float, float)> controlFunction):	controlFunction(controlFunction), camera(camera) {}

void CameraController::update(float t, float dt) {
	controlFunction(camera, t, dt);
}

CameraPositionController::CameraPositionController(sptr<Camera> camera, std::function<vec3(vec3, float, float)> controlFunction):	CameraController(camera, [f=controlFunction](sptr<Camera> cam, float t, float dt) {
	cam->set_position(f(cam->get_position(), t, dt));
}) {}
