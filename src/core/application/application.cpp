#include "application.hpp"
#include "GL/glew.h"


App::App(const AppSettings& settings) {
	logging::Logger::init();
	glfwSetErrorCallback([](int error, const char* description) {
		LOG_ERROR("GLFW Error (" + std::to_string(error) + "): " + string(description));
	});
	glfwInit();
	this->settings = settings;
	this->window = make_shared<Window>(settings.windowSettings);

	glewExperimental = true;
	if (glewInit() != GLEW_OK)
		THROW(SystemError, "Failed to initialize GLEW");

}

void App::run() {
	running = true;
	float last_time = glfwGetTime();
	float time = 0;
	float dt = 0;

	while (running) {
		if (not window->isOpen()) {
			stop();
			break;
		}
		time = glfwGetTime();
		dt = time - last_time;

		for (auto &layer : layers)
			layer->onUpdate(time, dt);

		for (auto &layer : layers)
			layer->onRender();

		last_time = time;
		window->update();
	}
}

void App::stop() {
	running = false;
}

shared_ptr<Window> App::getWindow() const {
	return window;
}
void App::pushLayer(unique_ptr<Layer> layer) {
	layer->attach();
	layers.push_back(std::move(layer));
}
