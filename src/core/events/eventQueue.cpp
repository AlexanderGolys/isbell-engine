#include "eventQueue.hpp"

void EventQueue::pushEvent(sptr<Event> event) {
	events.push_back(event);
}

vector<sptr<Event>>::iterator EventQueue::begin() {
	return events.begin();
}

vector<sptr<Event>>::iterator EventQueue::end() {
	return events.end();
}

void EventQueue::clearQueue() {
	events.clear();
}



sptr<EventQueue> EventQueue::getInstance() {
	if (!instance)
		instance = sptr<EventQueue>(new EventQueue());
	return instance;
}

void EventQueue::push(sptr<Event> event) {
	getInstance()->pushEvent(event);
}

sptr<EventQueue> EventQueue::instance = nullptr;
