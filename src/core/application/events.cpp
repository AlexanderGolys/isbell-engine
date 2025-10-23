#include "events.hpp"


bool Event::isHandled() const {
	return handled;
}

void Event::setHandled() {
	handled = true;
}

WindowResizeEvent::WindowResizeEvent(int width, int height)
: width(width), height(height) {}

EventType WindowResizeEvent::getType() const {
	return EventType::EVENT_WINDOW_RESIZE;
}

string WindowResizeEvent::getName() const {
	return "WindowResizeEvent (" + std::to_string(width) + ", " + std::to_string(height) + ")";
}

WindowMoveEvent::WindowMoveEvent(int xPos, int yPos): xPos(xPos), yPos(yPos) {}

EventType WindowMoveEvent::getType() const {
	return EventType::EVENT_WINDOW_MOVE;
}
string WindowMoveEvent::getName() const {
	return "WindowMoveEvent (" + std::to_string(xPos) + ", " + std::to_string(yPos) + ")";
}

EventType WindowCloseEvent::getType() const {
	return EventType::EVENT_WINDOW_CLOSE;
}

string WindowCloseEvent::getName() const {
	return "WindowCloseEvent";
}

KeyEvent::KeyEvent(int keyCode)
: keyCode(static_cast<Key>(keyCode)) {}

KeyEvent::KeyEvent(Key keyCode)
: keyCode(keyCode) {}

KeyPressedEvent::KeyPressedEvent(int keyCode, int repeatCount)
: KeyEvent(keyCode), repeatCount(repeatCount) {}

KeyPressedEvent::KeyPressedEvent(Key keyCode, int repeatCount)
: KeyEvent(keyCode), repeatCount(repeatCount) {}

EventType KeyPressedEvent::getType() const {
	return EventType::EVENT_KEY_PRESSED;
}

string KeyPressedEvent::getName() const {
	return "KeyPressedEvent(" + keyToString(keyCode) + ", " + std::to_string(repeatCount) + ")";
}

EventType KeyReleasedEvent::getType() const {
	return EventType::EVENT_KEY_RELEASED;
}

string KeyReleasedEvent::getName() const {
	return "KeyReleasedEvent(" + keyToString(keyCode) + ")";
}

MouseButtonEvent::MouseButtonEvent(MouseButton button): button(button) {}

EventType MouseButtonPressedEvent::getType() const {
	return EventType::EVENT_MOUSE_BUTTON_PRESSED;
}

string MouseButtonPressedEvent::getName() const {
	return "MouseButtonPressed(" + mouseButtonToString(button) + ")";
}

EventType MouseButtonReleasedEvent::getType() const {
	return EventType::EVENT_MOUSE_BUTTON_RELEASED;
}
string MouseButtonReleasedEvent::getName() const {
	return "MouseButtonReleased(" + mouseButtonToString(button) + ")";
}
MouseMovedEvent::MouseMovedEvent(float xPos, float yPos): xPos(xPos), yPos(yPos) {}
EventType MouseMovedEvent::getType() const {
	return EventType::EVENT_MOUSE_MOVED;
}
string MouseMovedEvent::getName() const {
	return "MouseMovedEvent(" + std::to_string(xPos) + ", " + std::to_string(yPos) + ")";
}
