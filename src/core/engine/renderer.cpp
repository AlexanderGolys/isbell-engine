#include "renderer.hpp"
#include "glCommand.hpp"

Renderer::Renderer(const RenderSettings& settings): settings(settings), eventQueue(EventQueue::getInstance()) {
	logging::Logger::init();
	GLFWCommand::init();
	window = make_unique<Window>(settings.windowSettings);
	GLCommand::init();
	GLCommand::enableDepth();
	GLCommand::enableBlending();
	fpsClock = FPSClock(settings.measureFPSStats ? settings.FPSAverageWindowSeconds : -1.f);
}

void Renderer::addLayer(sptr<Layer> layer) {
	layerStack.push_back(layer);
}

void Renderer::initRendering() {
	for (const auto& layer : layerStack)
		layer->init();
	LOG("OpenGL renderer initialized");
}

void Renderer::update(float t, float delta) const {
	for (const auto& layer : layerStack)
		layer->updateStep(t, delta);
}

void Renderer::renderStep() const {
	GLCommand::clearScreen(settings.bgColor);
	for (const auto& layer : layerStack)
		layer->renderStep();
	window->renderFramebufferToScreen();
	glfwPollEvents();
}

void Renderer::eventHandlingStep(float t, float delta) const {
	for (const auto& event : eventQueue->begin(), eventQueue->end())
		for (const auto& layer : layerStack)
			layer->onEvent(event, t, delta);
	eventQueue->clearQueue();
}

void Renderer::mainLoop() {
	while (window->isOpen()) {
		auto [t, dt] = fpsClock.tick();
		update(t, dt);
		eventHandlingStep(t, dt);
		renderStep();
	}
}

void Renderer::run() {
	LOG("Initializing rendering...");
	LOG("--------------------------------");
	initRendering();
	LOG("Rendering initialised. Starting main loop...");
	LOG("--------------------------------");
	mainLoop();
}
