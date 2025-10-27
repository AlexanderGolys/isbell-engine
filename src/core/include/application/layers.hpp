#pragma once
#include "logging.hpp"
#include "events.hpp"


class LayerABC {
public:
	virtual ~LayerABC();

	virtual void onUpdate(float time, float dt) {}
	virtual void onEvent(Event& event) {}
	virtual void onRender() = 0;
	virtual void attach() {}
	virtual void detach() {}
};

template <typename L>
concept ImplementsLayer = std::is_base_of_v<LayerABC, L>;

