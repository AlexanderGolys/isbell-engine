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

void Renderer::update(TimeStep timeStep) const {
	for (const auto& layer : layerStack)
		layer->updateStep(timeStep);
}

void Renderer::renderStep() const {
	GLCommand::clearScreen(settings.bgColor);
	for (const auto& layer : layerStack)
		layer->renderStep();
	window->renderFramebufferToScreen();
	glfwPollEvents();
}

void Renderer::eventHandlingStep(TimeStep timeStep) const {
	for (auto event : *eventQueue)
		for (const auto& layer : layerStack)
			layer->onEvent(event, timeStep);
	eventQueue->clearQueue();
}

void Renderer::mainLoop() {
	while (window->isOpen()) {
		TimeStep t = fpsClock.tick();
		eventHandlingStep(t);
		update(t);
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
