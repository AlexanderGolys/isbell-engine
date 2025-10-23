#pragma once
#include "keycodes.hpp"

enum class EventType {
    None = 0,
    EVENT_WINDOW_CLOSE,
	EVENT_WINDOW_RESIZE,
	EVENT_WINDOW_MOVE,
	EVENT_KEY_PRESSED,
	EVENT_KEY_RELEASED,
	EVENT_MOUSE_BUTTON_PRESSED,
	EVENT_MOUSE_BUTTON_RELEASED,
	EVENT_MOUSE_MOVED,
};

class Event {
    bool handled = false;
public:
    virtual ~Event() = default;
    virtual EventType getType() const = 0;
    virtual string getName() const = 0;
    bool isHandled() const;
    void setHandled();
};

class WindowResizeEvent : public Event {
public:
    int width, height;
    WindowResizeEvent(int width, int height);
    EventType getType() const override;
    string getName() const override;
};

class WindowMoveEvent : public Event {
public:
	int xPos, yPos;
	WindowMoveEvent(int xPos, int yPos);

	EventType getType() const override;
	string getName() const override;
};

class WindowCloseEvent : public Event {
public:
	WindowCloseEvent() = default;
	EventType getType() const override;
	string getName() const override;
};

class KeyEvent : public Event {
public:
    Key keyCode;
    explicit KeyEvent(int keyCode);
    explicit KeyEvent(Key keyCode);
};

class KeyPressedEvent : public KeyEvent {
public:
    int repeatCount;
	explicit KeyPressedEvent(int keyCode, int repeatCount=1);
	explicit KeyPressedEvent(Key keyCode, int repeatCount=1);

    EventType getType() const override;
    string getName() const override;
};

class KeyReleasedEvent : public KeyEvent {
public:
    using KeyEvent::KeyEvent;

    EventType getType() const override;
    string getName() const override;
};

class MouseButtonEvent : public Event {
public:
	MouseButton button;
	explicit MouseButtonEvent(MouseButton button);
};

class MouseButtonPressedEvent : public MouseButtonEvent {
public:
	using MouseButtonEvent::MouseButtonEvent;

	EventType getType() const override;
	string getName() const override;
};

class MouseButtonReleasedEvent : public MouseButtonEvent {
public:
	using MouseButtonEvent::MouseButtonEvent;

	EventType getType() const override;
	string getName() const override;
};

class MouseMovedEvent : public Event {
public:
	float xPos, yPos;
	MouseMovedEvent(float xPos, float yPos);

	EventType getType() const override;
	string getName() const override;
};

