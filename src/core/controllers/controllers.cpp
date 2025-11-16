#include "controllers.hpp"




CameraController::CameraController(sptr<Camera> camera, BIHOM(CameraTransform, TimeStep, CameraTransform) controlFunction)
:	Controller(camera, [f=controlFunction](sptr<Camera> cam, TimeStep timeStep) {
	cam->set_transform(f(cam->get_transform(), timeStep));
}) {}

ImpulseOnKeyComponent::ImpulseOnKeyComponent(const Distribution& impulseDistribution, KeyCode key, KeyCode keyReverse, float regenerationTime)
: impulseAccumulator(make_shared<DefferedResponseAccumulator>()),
	impulseDistribution(impulseDistribution), forwardTimer(regenerationTime), reverseTimer(regenerationTime), key(key), keyReverse(keyReverse) {}

void ImpulseOnKeyComponent::onKeyPressed(KeyCode keyCode, TimeStep timeStep) {
	if (keyCode == key and forwardTimer.available(timeStep.t))
		impulseAccumulator->addImpulse(impulseDistribution);
	else if (keyCode == keyReverse and reverseTimer.available(timeStep.t))
		impulseAccumulator->addImpulse(-impulseDistribution);
}

float ImpulseOnKeyComponent::step(float dt) const {
	return impulseAccumulator->step(dt);
}

ImpulseOnScrollComponent::ImpulseOnScrollComponent(const Distribution& impulseDistribution, float scrollSensitivity, bool xAxis):
impulseAccumulator(make_unique<DefferedResponseAccumulator>()),
impulseDistribution(impulseDistribution),
scrollSensitivity(scrollSensitivity),
xAxis(xAxis) {}

void ImpulseOnScrollComponent::onScroll(float xOffset, float yOffset, TimeStep timeStep) {
	if (xAxis and not nearlyEqual(xOffset, 0.f))
		impulseAccumulator->addImpulse(impulseDistribution * xOffset * scrollSensitivity);
	else if (not xAxis and not nearlyEqual(yOffset, 0.f))
		impulseAccumulator->addImpulse(impulseDistribution * yOffset * scrollSensitivity);
}

float ImpulseOnScrollComponent::step(float dt) const {
	return impulseAccumulator->step(dt);
}

CameraTransformFromImpulseOnKey::CameraTransformFromImpulseOnKey(sptr<Camera> camera, const BIHOM(CameraTransform, float, CameraTransform)& onImpulse, const Distribution& impulseDistribution, KeyCode key, KeyCode keyReverse, float regenerationTime)
: CameraController(camera, [this](CameraTransform transform, TimeStep timeStep) {
	  return onImpulseCall(transform, step(timeStep.dt));
  }),
	ImpulseOnKeyComponent(impulseDistribution, key, keyReverse, regenerationTime),
	onImpulse(onImpulse) {}

CameraTransform CameraTransformFromImpulseOnKey::onImpulseCall(CameraTransform transform, float impulse) const {
	return onImpulse(transform, impulse);
}

CameraTransformFromImpulseOnScroll::CameraTransformFromImpulseOnScroll(sptr<Camera> camera, const std::function<CameraTransform(CameraTransform, float)>& onImpulse,
	const Distribution& impulseDistribution, float scrollSensitivity, bool xAxis)
:	CameraController(camera, [this](CameraTransform transform, TimeStep timeStep) {
		  return onImpulseCall(transform, step(timeStep.dt)); }),
	ImpulseOnScrollComponent(impulseDistribution, scrollSensitivity, xAxis),
	onImpulse(onImpulse) {}

CameraTransform CameraTransformFromImpulseOnScroll::onImpulseCall(CameraTransform transform, float impulse) const {
	return onImpulse(transform, impulse);
}

CameraHorisontalRotationOnKey::CameraHorisontalRotationOnKey(sptr<Camera> camera, const Distribution& impulseDistribution, KeyCode key, KeyCode keyReverse, float regenerationTime)
: CameraTransformFromImpulseOnKey(camera,[](CameraTransform transform, float angle) {
	  mat3 M = rotationMatrix(transform.upVector, angle);
	  return CameraTransform(
		  M * (transform.position - transform.lookAtPos) + transform.lookAtPos,
		  transform.lookAtPos,
		  transform.upVector
	  );}, impulseDistribution, key, keyReverse, regenerationTime) {}

CameraHorisontalShiftOnKey::CameraHorisontalShiftOnKey(sptr<Camera> camera, const Distribution& impulseDistribution, KeyCode key, KeyCode keyReverse, float regenerationTime)
: CameraTransformFromImpulseOnKey(camera,[](CameraTransform transform, float angle) {
	vec3 shift = transform.getRightVector() * angle;
	  return CameraTransform(
		  transform.position + shift,
		  transform.lookAtPos + shift,
		  transform.upVector
	  );}, impulseDistribution, key, keyReverse, regenerationTime) {}

CameraVerticalRotationOnKey::CameraVerticalRotationOnKey(sptr<Camera> camera, const Distribution& impulseDistribution, KeyCode key, KeyCode keyReverse, float regenerationTime)
: CameraTransformFromImpulseOnKey(camera, [](CameraTransform transform, float angle) {
	 vec3 v = normalize(cross(transform.lookAtPos - transform.position, transform.upVector));
	 mat3 M = rotationMatrix(v, angle);
	 return CameraTransform(
		 M * (transform.position - transform.lookAtPos) + transform.lookAtPos,
		 transform.lookAtPos,
		 M * transform.upVector
	 );
 }, impulseDistribution, key, keyReverse, regenerationTime) {}

CameraVerticalShiftOnKey::CameraVerticalShiftOnKey(sptr<Camera> camera, const Distribution& impulseDistribution, KeyCode key, KeyCode keyReverse, float regenerationTime)
: CameraTransformFromImpulseOnKey(camera, [](CameraTransform transform, float angle) {
	 return CameraTransform(
		 transform.position + transform.upVector * angle,
		 transform.lookAtPos + transform.upVector * angle,
		 transform.upVector
	 );
 }, impulseDistribution, key, keyReverse, regenerationTime) {}
CameraZoomOnScroll::CameraZoomOnScroll(sptr<Camera> camera, const Distribution& impulseDistribution, float scrollSensitivity)
: CameraTransformFromImpulseOnScroll(camera, [](CameraTransform transform, float zoomAmount) {
	 vec3 dir = transform.lookAtPos - transform.position;
	 float dist = length(dir);
	 dir = dir / dist;
	 dist = std::max(0.001f, dist - zoomAmount);
	 return CameraTransform(
		 transform.lookAtPos - dir * dist,
		 transform.lookAtPos,
		 transform.upVector
	 );
 }, impulseDistribution, scrollSensitivity, false) {}

CameraRotateOnScroll::CameraRotateOnScroll(sptr<Camera> camera, const Distribution& impulseDistribution, float scrollSensitivity)
: CameraTransformFromImpulseOnScroll(camera, [](CameraTransform transform, float s) {
	 mat3 M = rotationMatrix(transform.getDirection(), s);
	 return CameraTransform(
		 transform.position,
		 transform.lookAtPos,
		 M * transform.upVector
	 );
 }, impulseDistribution, scrollSensitivity, true) {}
