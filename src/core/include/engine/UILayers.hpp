#pragma once
#include "renderLayers.hpp"


class HoverListener : public EventListener {
	HOM(vec2, bool) hoverCheck;
	bool isHovered = false;
	bool isClicked = false;
	uint elementID;
public:
	HoverListener(uint elementID, const HOM(vec2, bool)& hoverCheck);
	bool listensToEventType(EventType type) const final;
	void onEvent(sptr<Event> event, TimeStep timeStep) final;
};

struct ElementIDGenerator {
	static uint currentID;
	static uint generateID();
};

class UILayer : public GenericMeshLayer {
	CONST_PROPERTY(uint, elementID);
	sptr<HOM(vec2, bool)> hoverCheck;
public:
	UILayer(sptr<ShaderProgram> shader, sptr<GeometricData> mesh, sptr<HOM(vec2, bool)> hoverCheck);
};