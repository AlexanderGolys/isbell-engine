#include "UILayers.hpp"


HoverListener::HoverListener(uint elementID, const std::function<bool(vec2)>& hoverCheck): hoverCheck(hoverCheck), elementID(elementID) {}

bool HoverListener::listensToEventType(EventType type) const {
	return type == EventType::MOUSE_MOVED or type == EventType::MOUSE_BUTTON_PRESSED or type == EventType::MOUSE_BUTTON_RELEASED;
}

void HoverListener::onEvent(sptr<Event> event, TimeStep timeStep) {
	if (event->getEventType() != EventType::MOUSE_MOVED) {
		auto e = static_pointer_cast<MouseMovedEvent>(event);
		bool currentlyHovered = hoverCheck(e->get_position());
		if (currentlyHovered && !isHovered) {
			isHovered = true;
			ElementHoveredEvent::emmit(elementID);
		} else if (!currentlyHovered && isHovered) {
			isHovered = false;
			ElementUnhoveredEvent::emmit(elementID);
		}
	}
	else if (event->getEventType() == EventType::MOUSE_BUTTON_PRESSED) {
		auto e = static_pointer_cast<MouseButtonPressedEvent>(event);
		vec2 mousePos = e->get_position();
		if (e->get_button() == GLFW_MOUSE_BUTTON_LEFT and hoverCheck(mousePos)) {
			isClicked = true;
			ElementClickedEvent::emmit(elementID);
		}
	}
	else if (event->getEventType() == EventType::MOUSE_BUTTON_RELEASED and isClicked) {
		auto e = static_pointer_cast<MouseButtonReleasedEvent>(event);
		if (e->get_button() == GLFW_MOUSE_BUTTON_LEFT) {
			vec2 mousePos = e->get_position();
			isClicked = false;
			bool inside = hoverCheck(mousePos);
			if (inside)
				ElementUnclickedInsideEvent::emmit(elementID);
			else
				ElementUnclickedOutsideEvent::emmit(elementID);
		}
	}
}

uint ElementIDGenerator::generateID() {
	return currentID++;
}

UILayer::UILayer(sptr<ShaderProgram> shader, sptr<GeometricData> mesh, sptr<HOM(vec2, bool)> hoverCheck)
: GenericMeshLayer(shader, mesh),
elementID(ElementIDGenerator::generateID()),
hoverCheck(hoverCheck)
{
	if (hoverCheck)
		addEventListener(make_shared<HoverListener>(elementID, *hoverCheck));
}

uint ElementIDGenerator::currentID = 0;
