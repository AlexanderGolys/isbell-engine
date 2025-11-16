#pragma once

#include "clock.hpp"
#include "colors.hpp"
#include "eventQueue.hpp"
#include "renderLayers.hpp"
#include "window.hpp"


struct RenderSettings {
	WindowSettings windowSettings;
	Color bgColor;
	bool alphaBlending = true;
	bool depthTest = true;
	bool timeUniform = true;
	bool measureFPSStats = true;
	float FPSAverageWindowSeconds = 2.f;
	float speed = 1.0f;
};


class Renderer {
	RenderSettings settings;
	unique_ptr<Window> window;
	vector<sptr<Layer>> layerStack;
	FPSClock fpsClock;
	sptr<EventQueue> eventQueue;

public:
	explicit Renderer(const RenderSettings& settings);
	virtual ~Renderer() = default;

	void addLayer(sptr<Layer> layer);
	virtual void initRendering();
	void update(TimeStep timeStep) const;
	void renderStep() const;
	void eventHandlingStep(TimeStep timeStep) const;
	virtual void mainLoop();
	void run();
};
