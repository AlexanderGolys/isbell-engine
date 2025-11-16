#include "listeners.hpp"

#include "eventQueue.hpp"

bool KeyPressedListener::listensToEventType(EventType type) const {return type == EventType::KEY_PRESSED;}

void KeyPressedListener::onEvent(sptr<Event> event, TimeStep timeStep) {
	auto e = static_pointer_cast<KeyPressedEvent>(event);
	onKeyPressed(e->get_keycode(), timeStep);
}

bool ScrollListener::listensToEventType(EventType type) const {
	return type == EventType::MOUSE_SCROLLED;
}

void ScrollListener::onEvent(sptr<Event> event, TimeStep timeStep) {
	auto scrollEvent = static_pointer_cast<MouseScrolledEvent>(event);
	onScroll(scrollEvent->get_xOffset(), scrollEvent->get_yOffset(), timeStep);
}
