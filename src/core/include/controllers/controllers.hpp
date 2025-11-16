#pragma once
#include "accumulators.hpp"
#include "listeners.hpp"
#include "sceneRendering.hpp"

template <typename T>
class Controller : public UpdateComponent {
	sptr<T> object;
	BIHOM(sptr<T>, TimeStep, void) controlFunction;

public:
	Controller(sptr<T> object, BIHOM(sptr<T>, TimeStep, void) controlFunction) : object(object), controlFunction(controlFunction) {}
	void update(TimeStep timeStep) final { controlFunction(object, timeStep); }
	void init() final {}
};

class CameraController : public Controller<Camera> {
public:
	CameraController(sptr<Camera> camera, BIHOM(CameraTransform, TimeStep, CameraTransform) controlFunction);
};

class ImpulseOnKeyComponent : public KeyPressedListener {
	sptr<DefferedResponseAccumulator> impulseAccumulator;
	Distribution impulseDistribution;
	WaitTimer forwardTimer, reverseTimer;
	KeyCode key, keyReverse;

public:
	ImpulseOnKeyComponent(const Distribution& impulseDistribution, KeyCode key, KeyCode keyReverse, float regenerationTime=0);
	void onKeyPressed(KeyCode keyCode, TimeStep timeStep) final;
	float step(float dt) const;
};

class ImpulseOnScrollComponent : public ScrollListener {
	uptr<DefferedResponseAccumulator> impulseAccumulator;
	Distribution impulseDistribution;
	float scrollSensitivity;
	bool xAxis;

public:
	ImpulseOnScrollComponent(const Distribution& impulseDistribution, float scrollSensitivity, bool xAxis=true);

	void onScroll(float xOffset, float yOffset, TimeStep timeStep) final;
	float step(float dt) const;
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