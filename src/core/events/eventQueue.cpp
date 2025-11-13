#include "eventQueue.hpp"

void EventQueue::pushEvent(const Event& event) {
	events.push_back(event);
}

vector<Event>::iterator EventQueue::begin() {
	return events.begin();
}

vector<Event>::iterator EventQueue::end() {
	return events.end();
}

void EventQueue::clearQueue() { events.clear(); }



sptr<EventQueue> EventQueue::getInstance() {
	if (!instance)
		instance = make_shared<EventQueue>();
	return instance;
}

void EventQueue::push(const Event& event) {
	getInstance()->pushEvent(event);
}

sptr<EventQueue> EventQueue::instance = nullptr;
