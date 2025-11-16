#pragma once
#include "keyCodes.hpp"
#include "logging.hpp"


enum class EventType : uint {
	WINDOW_CLOSE = 0,
	WINDOW_RESIZE = 1,

	KEY_PRESSED = 2,
	KEY_RELEASED = 3,
	KEY_REPEAT = 4,

	MOUSE_BUTTON_PRESSED = 5,
	MOUSE_BUTTON_RELEASED = 6,
	MOUSE_BUTTON_REPEAT = 7,
	MOUSE_MOVED = 8,
	MOUSE_SCROLLED = 9,

	KEY_COMBINATION_PRESSED = 10,
	KEY_COMBINATION_RELEASED = 11,

	ELEMENT_HOVERED = 12,
	ELEMENT_UNHOVERED = 13,
	ELEMENT_CLICKED = 14,
	ELEMENT_UNCLICKED_INSIDE = 15,
	ELEMENT_UNCLICKED_OUTSIDE = 16,
};

#define EVENT_TYPE(type) EventType getEventType() const override { return EventType::type; }


class Event {
public:
	virtual ~Event() = default;
	virtual EventType getEventType() const = 0;
};

class KeyEvent : public Event {
	CONST_PROPERTY(KeyCode, keycode);
public:
	explicit KeyEvent(KeyCode keycode) : keycode(keycode) {}
};

class MouseButtonEvent : public Event {
	CONST_PROPERTY(MouseButtonCode, button);
	CONST_PROPERTY(vec2, position);
public:
	explicit MouseButtonEvent(MouseButtonCode button, vec2 position) : button(button), position(position) {}
};


class WindowCloseEvent : public Event {
public:
	WindowCloseEvent() = default;
	EVENT_TYPE(WINDOW_CLOSE);
};


class WindowResizeEvent : public Event {
	CONST_PROPERTY(uint, width);
	CONST_PROPERTY(uint, height);
public:
	WindowResizeEvent(uint width, uint height) : width(width), height(height) {}
	EVENT_TYPE(WINDOW_RESIZE);
	static void emmit(uint width, uint height);
};


class KeyPressedEvent : public KeyEvent {
public:
	using KeyEvent::KeyEvent;
	EVENT_TYPE(KEY_PRESSED);

	static void emmit(KeyCode keycode);
};


class KeyReleasedEvent : public KeyEvent {
public:
	using KeyEvent::KeyEvent;
	EVENT_TYPE(KEY_RELEASED);

	static void emmit(KeyCode keycode);
};

class KeyRepeatEvent : public KeyEvent {
public:
	using KeyEvent::KeyEvent;
	EVENT_TYPE(KEY_REPEAT);

	static void emmit(KeyCode keycode);
};

class MouseButtonPressedEvent : public MouseButtonEvent {
public:
	using MouseButtonEvent::MouseButtonEvent;
	EVENT_TYPE(MOUSE_BUTTON_PRESSED);
	static void emmit(MouseButtonCode button, vec2 position);
};


class MouseButtonReleasedEvent : public MouseButtonEvent {
public:
	using MouseButtonEvent::MouseButtonEvent;
	EVENT_TYPE(MOUSE_BUTTON_RELEASED);

	static void emmit(MouseButtonCode button, vec2 position);
};

class MouseButtonRepeatEvent : public MouseButtonEvent {
public:
	using MouseButtonEvent::MouseButtonEvent;
	EVENT_TYPE(MOUSE_BUTTON_REPEAT);

	static void emmit(MouseButtonCode button, vec2 position);
};


class MouseMovedEvent : public Event {
	CONST_PROPERTY(vec2, position);
public:
	explicit MouseMovedEvent(vec2 position) : position(position) {}
	EVENT_TYPE(MOUSE_MOVED);

	static void emmit(vec2 position);
};


class MouseScrolledEvent : public Event {
	CONST_PROPERTY(float, xOffset);
	CONST_PROPERTY(float, yOffset);
public:
	MouseScrolledEvent(float xOffset, float yOffset) : xOffset(xOffset), yOffset(yOffset) {}
	EVENT_TYPE(MOUSE_SCROLLED);

	static void emmit(float xOffset, float yOffset);
};

class KeyCombinationEvent : public Event {
	CONST_PROPERTY(KeyCode, key1);
	CONST_PROPERTY(KeyCode, key2);
public:
	KeyCombinationEvent(KeyCode key1, KeyCode key2) : key1(key1), key2(key2) {}
	EVENT_TYPE(KEY_COMBINATION_PRESSED);

	static void emmit(KeyCode key1, KeyCode key2);
};

class KeyCombinationReleasedEvent : public Event {
	CONST_PROPERTY(KeyCode, key1);
	CONST_PROPERTY(KeyCode, key2);
public:
	KeyCombinationReleasedEvent(KeyCode key1, KeyCode key2) : key1(key1), key2(key2) {}
	EVENT_TYPE(KEY_COMBINATION_RELEASED);

	static void emmit(KeyCode key1, KeyCode key2);
};

class ElementHoveredEvent : public Event {
	CONST_PROPERTY(uint, elementID);
public:
	explicit ElementHoveredEvent(uint elementID) : elementID(elementID) {}
	EVENT_TYPE(ELEMENT_HOVERED);

	static void emmit(uint elementID);
};

class ElementUnhoveredEvent : public Event {
	CONST_PROPERTY(uint, elementID);
public:
	explicit ElementUnhoveredEvent(uint elementID) : elementID(elementID) {}
	EVENT_TYPE(ELEMENT_UNHOVERED);

	static void emmit(uint elementID);
};

class ElementClickedEvent : public Event {
	CONST_PROPERTY(uint, elementID);
public:
	explicit ElementClickedEvent(uint elementID) : elementID(elementID) {}
	EVENT_TYPE(ELEMENT_CLICKED);
	static void emmit(uint elementID);
};

class ElementUnclickedInsideEvent : public Event {
	CONST_PROPERTY(uint, elementID);
public:
	explicit ElementUnclickedInsideEvent(uint elementID) : elementID(elementID) {}
	EVENT_TYPE(ELEMENT_UNCLICKED_INSIDE);
	static void emmit(uint elementID);
};

class ElementUnclickedOutsideEvent : public Event {
	CONST_PROPERTY(uint, elementID);
public:
	explicit ElementUnclickedOutsideEvent(uint elementID) : elementID(elementID) {}
	EVENT_TYPE(ELEMENT_UNCLICKED_OUTSIDE);
	static void emmit(uint elementID);
};
