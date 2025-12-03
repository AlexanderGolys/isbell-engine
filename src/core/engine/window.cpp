#include "window.hpp"

#include "event.hpp"
#include "exceptions.hpp"
#include "glCommand.hpp"

Resolution::Resolution(TypicalResolution r) {
	switch (r) {
		case TypicalResolution::p480:
			width = 854;
			height = 480;
			break;
		case TypicalResolution::HD720:
			width = 1280;
			height = 720;
			break;
		case TypicalResolution::FHD:
			width = 1920;
			height = 1080;
			break;
		case TypicalResolution::HD2K:
			width = 2560;
			height = 1440;
			break;
		case TypicalResolution::UHD:
			width = 3840;
			height = 2160;
			break;
	}
	THROW(ValueError, "Unknown resolution.");
}

WindowSettings::WindowSettings(Resolution resolution, string_cr windowTitle, bool stickyKeys, bool stickyMouseButtons)
: resolution(resolution), windowTitle(windowTitle), stickyKeys(stickyKeys), stickyMouseButtons(stickyMouseButtons) {}

Window::Window(WindowSettings settings) : settings(settings) {
	window = GLFWCommand::createWindow(settings.resolution, settings.windowTitle.c_str());
	GLFWCommand::setWindowData(window, &this->settings);

	glfwSetInputMode(window, GLFW_STICKY_KEYS, settings.stickyKeys ? GLFW_TRUE : GLFW_FALSE);
	glfwSetInputMode(window, GLFW_STICKY_MOUSE_BUTTONS, settings.stickyMouseButtons ? GLFW_TRUE : GLFW_FALSE);

	glfwSetWindowSizeCallback(window, [](GLFWwindow* w, int width, int height){
		WindowResizeEvent::emmit(width, height);
		WindowSettings& winData = GLFWCommand::getWindowData(w);
		winData.width = width;
		winData.height = height;
	});

	glfwSetWindowCloseCallback(window, [](GLFWwindow* w){
		WindowCloseEvent::emmit();
		WindowSettings& winData = GLFWCommand::getWindowData(w);
		winData.open = false;
	});

	glfwSetKeyCallback(window, [](GLFWwindow* w, int key, int scancode, int action, int mods){
		if (action == GLFW_PRESS)
			KeyPressedEvent::emmit(key);
		else if (action == GLFW_RELEASE)
			KeyReleasedEvent::emmit(key);
		else if (action == GLFW_REPEAT)
			KeyRepeatEvent::emmit(key);

	});

	glfwSetMouseButtonCallback(window, [](GLFWwindow* w, int button, int action, int mods){
		vec2 mousePos = GLFWCommand::getCursorPosition(w);
		if (action == GLFW_PRESS) {
			MouseButtonPressedEvent::emmit(button, mousePos);
		}
		else if (action == GLFW_RELEASE){
			MouseButtonReleasedEvent::emmit(button, mousePos);
		}
		else if (action == GLFW_REPEAT)
			MouseButtonRepeatEvent::emmit(button, mousePos);

	});

	glfwSetCursorPosCallback(window, [](GLFWwindow* w, double xpos, double ypos){
		MouseMovedEvent::emmit(vec2(xpos, ypos));
	});

	glfwSetScrollCallback(window, [](GLFWwindow* w, double xoffset, double yoffset){
		MouseScrolledEvent::emmit(xoffset, yoffset);
	});
}

Window::~Window() {
	destroy();
}

void Window::destroy() const {
	GLFWCommand::destroyWindow(window);
	GLFWCommand::terminate();
}

void Window::renderFramebufferToScreen() const {
	GLFWCommand::swapFrameBuffers(window);
}

bool Window::isOpen() const {
	return GLFWCommand::isWindowOpen(window);
}

void Window::initViewport() const {
	GLCommand::initViewport(settings.width, settings.height);
}
