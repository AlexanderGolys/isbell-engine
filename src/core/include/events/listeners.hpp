#pragma once
#include "clock.hpp"
#include "event.hpp"

class EventListener {
public:
	virtual ~EventListener() = default;
	virtual void onEvent(sptr<Event> event, TimeStep timeStep) = 0;
	virtual bool listensToEventType(EventType type) const = 0;
};

class KeyPressedListener : public EventListener {
public:
	virtual void onKeyPressed(KeyCode key, TimeStep timeStep) = 0;

	bool listensToEventType(EventType type) const final;
	void onEvent(sptr<Event> event, TimeStep timeStep) final;
};

class ScrollListener : public EventListener {
public:
	virtual void onScroll(float xOffset, float yOffset, TimeStep timeStep) = 0;

	bool listensToEventType(EventType type) const final;
	void onEvent(sptr<Event> event, TimeStep timeStep) final;
};

