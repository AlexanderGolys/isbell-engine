#pragma once
#include "accumulators.hpp"
#include "sceneRendering.hpp"

class Controller : public LayerComponent {
public:
	void init() final {}
	void setDuringRender() const final {}
};

class CameraController : public Controller {
	std::function<void(sptr<Camera>, float, float)> controlFunction;
	CONST_PROPERTY(sptr<Camera>, camera);

public:
	CameraController(sptr<Camera> camera, std::function<void(sptr<Camera>, float, float)> controlFunction);
	void update(float t, float dt) override;
};

class CameraPositionController : public CameraController {
public:
	CameraPositionController(sptr<Camera> camera, std::function<vec3(vec3, float, float)> controlFunction);
};

class CameraRotatorFromImpulse : public CameraController {
	sptr<DefferedResponseAccumulator> impulseAccumulator;
public:
	CameraRotatorFromImpulse(sptr<Camera> camera, sptr<DefferedResponseAccumulator> impulseAccumulator)
	: CameraController(camera, [&impulseAccumulator](sptr<Camera> cam, float t, float dt) {
		float angle = impulseAccumulator->step(dt);
		mat3 rot = rotationMatrix(cam->get_upVector(), angle);
		vec3 pos = rot*(cam->get_position() - cam->get_lookAtPos()) + cam->get_lookAtPos();
		cam->set_position(pos);
	}), impulseAccumulator(impulseAccumulator) {}
};

class DefferedImpulseFromKeypress : public EventListener {
	sptr<DefferedResponseAccumulator> impulseAccumulator;
	Distribution impulseDistribution;
	KeyCode key;
	float regenerationTime;
	float lastImpulseTime;
public:
	DefferedImpulseFromKeypress(sptr<DefferedResponseAccumulator> impulseAccumulator, const Distribution& impulseDistribution, KeyCode key, float regenerationTime = 2.f)
	: impulseAccumulator(impulseAccumulator), impulseDistribution(impulseDistribution), key(key), regenerationTime(regenerationTime), lastImpulseTime(-regenerationTime) {}

	void onEvent(const Event& event, float t, float dt) override {
		if (event.getEventType() == EventType::KeyPressed) {
			auto& e = static_cast<const KeyPressedEvent&>(event);
			if (e.get_keycode() == key) {
				if (t - lastImpulseTime < regenerationTime)
					return;
				lastImpulseTime = t;
				impulseAccumulator->addImpulse(impulseDistribution);
			}
		}
	}
};