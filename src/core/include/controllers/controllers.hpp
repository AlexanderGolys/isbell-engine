#pragma once
#include "accumulators.hpp"
#include "listeners.hpp"
#include "renderLayers.hpp"

template <typename Func, typename Obj>
concept controller_function =
	same_as<Func, BIHOM(sptr<Obj>&, TimeStep, void)> ||
	same_as<Func, HOM(TimeStep, Obj)> ||
	same_as<Func, HOM(TimeStep, sptr<Obj>)>;


template <typename T>
class Controller : public UpdateComponent {
	sptr<T> object;
	BIHOM(sptr<T>&, TimeStep, void) controlFunction;

public:
	Controller(sptr<T> object, BIHOM(sptr<T>&, TimeStep, void) controlFunction);
	Controller(sptr<T> object, HOM(TimeStep, T) setterFunction);
	Controller(sptr<T> object, HOM(TimeStep, sptr<T>) setterFunction);

	void update(TimeStep timeStep) final;
	void init() override {}
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





// ------------------- Implementation ------------------ //


template <typename T>
Controller<T>::Controller(sptr<T> object, std::function<void(sptr<T>&, TimeStep)> controlFunction): object(object), controlFunction(controlFunction) {}

template <typename T>
Controller<T>::Controller(sptr<T> object, std::function<T(TimeStep)> setterFunction): object(object) {
	controlFunction = [setterFunction](sptr<T> &obj, TimeStep timeStep) {
		obj = make_shared<T>(setterFunction(timeStep));
	};
}

template <typename T>
Controller<T>::Controller(sptr<T> object, std::function<sptr<T>(TimeStep)> setterFunction): object(object) {
	controlFunction = [setterFunction](sptr<T>& obj, TimeStep timeStep) {
		obj = setterFunction(timeStep);
	};
}

template <typename T>
void Controller<T>::update(TimeStep timeStep) { controlFunction(object, timeStep); }
