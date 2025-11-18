#pragma once
#include "controllers.hpp"
#include "sceneRendering.hpp"


class CameraController : public Controller<Camera> {
public:
	CameraController(sptr<Camera> camera, BIHOM(CameraTransform, TimeStep, CameraTransform) controlFunction);
};

class CameraTransformFromImpulseOnKey : public CameraController, public ImpulseOnKeyComponent {
	BIHOM(CameraTransform, float, CameraTransform) onImpulse;
public:
	CameraTransformFromImpulseOnKey(sptr<Camera> camera, const BIHOM(CameraTransform, float, CameraTransform)& onImpulse, const Distribution& impulseDistribution, KeyCode key, KeyCode keyReverse, float regenerationTime=0);
	CameraTransform onImpulseCall(CameraTransform transform, float impulse) const;
};

class CameraTransformFromImpulseOnScroll : public CameraController, public ImpulseOnScrollComponent {
	BIHOM(CameraTransform, float, CameraTransform) onImpulse;
public:
	CameraTransformFromImpulseOnScroll(sptr<Camera> camera, const BIHOM(CameraTransform, float, CameraTransform)& onImpulse, const Distribution& impulseDistribution, float scrollSensitivity, bool xAxis=true);
	CameraTransform onImpulseCall(CameraTransform transform, float impulse) const;
};

class CameraHorisontalRotationOnKey : public CameraTransformFromImpulseOnKey {
public:
	CameraHorisontalRotationOnKey(sptr<Camera> camera, const Distribution& impulseDistribution, KeyCode key, KeyCode keyReverse, float regenerationTime=0);
};

class CameraHorisontalShiftOnKey : public CameraTransformFromImpulseOnKey {
public:
	CameraHorisontalShiftOnKey(sptr<Camera> camera, const Distribution& impulseDistribution, KeyCode key, KeyCode keyReverse, float regenerationTime=0);
};

class CameraVerticalRotationOnKey : public CameraTransformFromImpulseOnKey {
public:
	CameraVerticalRotationOnKey(sptr<Camera> camera, const Distribution& impulseDistribution, KeyCode key, KeyCode keyReverse, float regenerationTime=0);
};

class CameraVerticalShiftOnKey : public CameraTransformFromImpulseOnKey {
public:
	CameraVerticalShiftOnKey(sptr<Camera> camera, const Distribution& impulseDistribution, KeyCode key, KeyCode keyReverse, float regenerationTime=0);
};

class CameraZoomOnScroll : public CameraTransformFromImpulseOnScroll {
public:
	CameraZoomOnScroll(sptr<Camera> camera, const Distribution& impulseDistribution, float scrollSensitivity);
};

class CameraRotateOnScroll : public CameraTransformFromImpulseOnScroll {
public:
	CameraRotateOnScroll(sptr<Camera> camera, const Distribution& impulseDistribution, float scrollSensitivity);
};


