#pragma once
#include "events.hpp"
#include "window.hpp"
#include "logging.hpp"

class Layer {
public:
    virtual ~Layer();

    virtual void onUpdate(float time, float dt) {}
    virtual void onEvent(Event &event) {}
    virtual void onRender() = 0;
	virtual void attach() {}
	virtual void detach() {}

};


struct AppSettings {
	WindowSettings windowSettings;
};

class App {
    AppSettings settings;
	shared_ptr<Window> window;
	vector<unique_ptr<Layer>> layers;
	bool running = false;
public:
	explicit App(const AppSettings &settings = AppSettings());
	void run();
	void stop();
	shared_ptr<Window> getWindow() const;
	void pushLayer(unique_ptr<Layer> layer);
	void onEvent(Event &event) {
		for (auto &layer : layers)
			layer->onEvent(event);
		if (event.getType() == EventType::EVENT_WINDOW_CLOSE)
			stop();
	}
};
