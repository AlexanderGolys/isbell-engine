#pragma once
#include "event.hpp"
#include "exceptions.hpp"


class EventQueue {
	vector<Event> events;
	EventQueue() = default;
	static sptr<EventQueue> instance;
  public:
	EventQueue(const EventQueue&) = delete;
	EventQueue& operator=(const EventQueue&) = delete;
	void pushEvent(const Event& event);

	vector<Event>::iterator begin();
	vector<Event>::iterator end();
	void clearQueue();

	static sptr<EventQueue> getInstance();
	static void push(const Event& event);

};
