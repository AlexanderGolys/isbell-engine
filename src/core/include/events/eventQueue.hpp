#pragma once
#include "event.hpp"
#include "exceptions.hpp"


class EventQueue {
	vector<sptr<Event>> events;
	EventQueue() {}
	static sptr<EventQueue> instance;

  public:
	EventQueue(const EventQueue&) = delete;
	EventQueue& operator=(const EventQueue&) = delete;
	void pushEvent(sptr<Event> event);

	vector<sptr<Event>>::iterator begin();
	vector<sptr<Event>>::iterator end();
	void clearQueue();

	static sptr<EventQueue> getInstance();
	static void push(sptr<Event> event);
};
