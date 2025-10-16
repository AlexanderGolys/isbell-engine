#pragma once
#include "func.hpp"


class Layer {
public:
	virtual ~Layer() = default;

	virtual void update(float t, float dt);
	virtual void render();
	virtual void onEvent(int event);
};
