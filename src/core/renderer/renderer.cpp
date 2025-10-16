#include "renderer.hpp"

RenderAPI Renderer::getAPI() {
	return currentAPI;
}

void Renderer::setAPI(RenderAPI api) {
	currentAPI = api;
}
