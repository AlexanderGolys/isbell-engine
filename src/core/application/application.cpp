#include "application.hpp"
#include "GL/glew.h"


App::App(const AppSettings& settings)
: running(false), time(0), dt(0) {
	logging::Logger::init();
	glfwSetErrorCallback([](int error, const char* description) {
		LOG_ERROR("GLFW Error (" + std::to_string(error) + "): " + string(description));
	});
	glfwInit();
	this->settings = settings;
	this->window = make_unique<Window>(settings.windowSettings);

	glewExperimental = true;
	if (glewInit() != GLEW_OK)
		THROW(SystemError, "Failed to initialize GLEW");
}

void App::run() {
	running = true;
	time = glfwGetTime();

	while (running) {
		if (not window->isOpen()) {
			stop();
			break;
		}
		float new_time = glfwGetTime();
		dt = new_time - time;
		time = new_time;

		for (auto& layer : layers)
			layer->onUpdate(time, dt);

		for (auto& layer : layers)
			layer->onRender();

		if (settings.measureFPS)
			fpsUpdate();

		window->update();
	}
}

void App::stop() {
	running = false;
}


void App::fpsUpdate() {
	fpsFlushAccumulator += dt;
	if (fpsFlushAccumulator >= settings.fpsUpdateInterval) {
		int fps = frameCount / settings.fpsUpdateInterval;
		int worstFps = 1000 / worstFrameTime;
		LOG("FPS: " + to_string(fps) + " (worst drop: " + to_string(worstFps) + ")");
		fpsFlushAccumulator = 0;
		frameCount = 0;
		worstFrameTime = 0;
	}
	worstFrameTime = std::max(worstFrameTime, dt);
	frameCount++;
}


void App::onEvent(Event& event) {
	for (auto& layer : layers)
		layer->onEvent(event);
	if (event.getType() == EventType::EVENT_WINDOW_CLOSE)
		stop();
}
