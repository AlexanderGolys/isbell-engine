#pragma once

#include <cstdio>
#include <cstdlib>
#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include "clock.hpp"
#include "indexedRendering.hpp"
#include "renderingUtils.hpp"
#include "renderLayers.hpp"
#include "window.hpp"


struct RenderSettings {
	WindowSettings windowSettings;
	vec4 bgColor;
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
	END(float) animSpeed;
	BIHOM(float, float, void) perFrameFunction;
	FPSClock fpsClock;

public:
	explicit Renderer(const RenderSettings& settings);
	virtual ~Renderer();

	void addLayer(sptr<Layer> layer);
	virtual void initRendering();
	void update();
	void renderStep() const;
	virtual int mainLoop();
};
