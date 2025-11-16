#include "event.hpp"

#include "eventQueue.hpp"

void WindowResizeEvent::emmit(uint width, uint height) {
	EventQueue::push(make_shared<WindowResizeEvent>(width, height));
}

void KeyPressedEvent::emmit(KeyCode keycode) {
	EventQueue::push(make_shared<KeyPressedEvent>(keycode));
}

void KeyReleasedEvent::emmit(KeyCode keycode) {
	EventQueue::push(make_shared<KeyReleasedEvent>(keycode));
}

void KeyRepeatEvent::emmit(KeyCode keycode) {
	EventQueue::push(make_shared<KeyRepeatEvent>(keycode));
}

void MouseButtonPressedEvent::emmit(MouseButtonCode button, vec2 position) {
	EventQueue::push(make_shared<MouseButtonPressedEvent>(button, position));
}

void MouseButtonReleasedEvent::emmit(MouseButtonCode button, vec2 position) {
	EventQueue::push(make_shared<MouseButtonReleasedEvent>(button, position));
}

void MouseButtonRepeatEvent::emmit(MouseButtonCode button, vec2 position) {
	EventQueue::push(make_shared<MouseButtonRepeatEvent>(button, position));
}

void MouseMovedEvent::emmit(vec2 position) {
	EventQueue::push(make_shared<MouseMovedEvent>(position));
}

void MouseScrolledEvent::emmit(float xOffset, float yOffset) {
	EventQueue::push(make_shared<MouseScrolledEvent>(xOffset, yOffset));
}

void KeyCombinationEvent::emmit(KeyCode key1, KeyCode key2) {
	EventQueue::push(make_shared<KeyCombinationEvent>(key1, key2));
}

void KeyCombinationReleasedEvent::emmit(KeyCode key1, KeyCode key2) {
	EventQueue::push(make_shared<KeyCombinationReleasedEvent>(key1, key2));
}

void ElementHoveredEvent::emmit(uint elementID) {
	EventQueue::push(make_shared<ElementHoveredEvent>(elementID));
}

void ElementUnhoveredEvent::emmit(uint elementID) {
	EventQueue::push(make_shared<ElementUnhoveredEvent>(elementID));
}

void ElementClickedEvent::emmit(uint elementID) {
	EventQueue::push(make_shared<ElementClickedEvent>(elementID));
}

void ElementUnclickedInsideEvent::emmit(uint elementID) {
	EventQueue::push(make_shared<ElementUnclickedInsideEvent>(elementID));
}

void ElementUnclickedOutsideEvent::emmit(uint elementID) {
	EventQueue::push(make_shared<ElementUnclickedOutsideEvent>(elementID));
}
