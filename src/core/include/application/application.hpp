#pragma once
#include "layers.hpp"
#include "window.hpp"


struct AppSettings {
	WindowSettings windowSettings;
	bool measureFPS;
	float fpsUpdateInterval; // seconds
};

class App {
	AppSettings settings;
	unique_ptr<Window> window;
	vector<unique_ptr<LayerABC>> layers;
	bool running;
	float time, dt; // ms
	float fpsFlushAccumulator = 0;
	int frameCount = 0;
	float worstFrameTime = 0;

public:
	explicit App(const AppSettings& settings);
	void run();
	void stop();
	void fpsUpdate();
	void onEvent(Event& event);

	template <ImplementsLayer L, class... Args>
	void pushLayer(Args&&... args);
};


template <ImplementsLayer L, class... Args>
void App::pushLayer(Args&&... args) {
	auto layer = make_unique<L>(std::forward<Args>(args)...);
	layer->attach();
	layers.push_back(std::move(layer));
}
