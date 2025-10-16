#pragma once
#include "exceptions.hpp"

enum RenderAPI {
	NO_RENDERING,
	OPENGL,
	VULKAN
};

class Renderer {
	static RenderAPI currentAPI;

public:
	static RenderAPI getAPI();
	static void setAPI(RenderAPI api);
};
