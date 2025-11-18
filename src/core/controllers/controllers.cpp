#include "controllers.hpp"






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
