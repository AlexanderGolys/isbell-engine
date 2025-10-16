#pragma once
#include "window.hpp"
#include "layer.hpp"

class Application {
	unique_ptr<Window> window;
	vector<shared_ptr<Layer>> layers;
	bool running;
};
