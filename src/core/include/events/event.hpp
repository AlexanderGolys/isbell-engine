#pragma once
#include "exceptions.hpp"
#include "keyCodes.hpp"


enum class EventType : uint {
	WindowClose = 0,
	WindowResize = 1,
	KeyPressed = 2,
	KeyReleased = 3,
	MouseButtonPressed = 4,
	MouseButtonReleased = 5,
	MouseMoved = 6,
	MouseScrolled = 7
};

class Event {
public:
	virtual ~Event() = default;
	virtual EventType getEventType() const = 0;
};

class WindowCloseEvent : public Event {
public:
	WindowCloseEvent() = default;
	EventType getEventType() const override {
		return EventType::WindowClose;
	}
};

class WindowResizeEvent : public Event {
	CONST_PROPERTY(uint, width);
	CONST_PROPERTY(uint, height);
public:
	WindowResizeEvent(uint width, uint height) : width(width), height(height) {}
	EventType getEventType() const override {
		return EventType::WindowResize;
	}
};

class KeyPressedEvent : public Event {
	CONST_PROPERTY(KeyCode, keycode);
public:
	explicit KeyPressedEvent(KeyCode keycode) : keycode(keycode) {}
	EventType getEventType() const override {
		return EventType::KeyPressed;
	}
};
class KeyReleasedEvent : public Event {
	CONST_PROPERTY(KeyCode, keycode);
public:
	explicit KeyReleasedEvent(KeyCode keycode) : keycode(keycode) {}
	EventType getEventType() const override {
		return EventType::KeyReleased;
	}
};
class MouseButtonPressedEvent : public Event {
	CONST_PROPERTY(MouseButton, button);
public:
	explicit MouseButtonPressedEvent(MouseButton button) : button(button) {}
	EventType getEventType() const override {
		return EventType::MouseButtonPressed;
	}
};

class MouseButtonReleasedEvent : public Event {
	CONST_PROPERTY(MouseButton, button);
public:
	explicit MouseButtonReleasedEvent(MouseButton button) : button(button) {}
	EventType getEventType() const override {
		return EventType::MouseButtonReleased;
	}
};

class MouseMovedEvent : public Event {
	CONST_PROPERTY(float, x);
	CONST_PROPERTY(float, y);
public:
	MouseMovedEvent(float x, float y) : x(x), y(y) {}
	EventType getEventType() const override {
		return EventType::MouseMoved;
	}
};

class MouseScrolledEvent : public Event {
	CONST_PROPERTY(float, xOffset);
	CONST_PROPERTY(float, yOffset);
public:
	MouseScrolledEvent(float xOffset, float yOffset) : xOffset(xOffset), yOffset(yOffset) {}
	EventType getEventType() const override {
		return EventType::MouseScrolled;
	}
};

