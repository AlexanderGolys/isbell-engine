#include "renderer.hpp"
#include <stdio.h>
#include <string>
#include <map>
#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <chrono>
#include <ranges>
#include <sstream>

#include <string.h>

#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include "configFiles.hpp"
#include "specific.hpp"

using namespace glm;

//
// MeshLayer::MeshLayer(const shared_ptr<ShaderProgram>& shader, const shared_ptr<MaterialPhong>& material, const shared_ptr<IndexedMesh>& mesh)
// : customStep([](float){}), elementBuffer(nullptr)
// {
// 	this->shader = shader;
// 	this->material = material;
// 	this->mesh = mesh;
// 	vao = make_shared<VertexArray>();
// }
//
//
//
//
//
// void MeshLayer::initMaterialTextures() const {
// 	material->initTextures();
// 	shader->initTextureSampler(material->texture_ambient.get());
// 	shader->initTextureSampler(material->texture_diffuse.get());
// 	shader->initTextureSampler(material->texture_specular.get());
// }
//
// void MeshLayer::bindTextures() const {
// 	material->bindTextures();
// 	shader->setTextureSampler(material->texture_ambient.get());
// 	shader->setTextureSampler(material->texture_diffuse.get());
// 	shader->setTextureSampler(material->texture_specular.get());
// }
//
//
//
// void MeshLayer::initAttributes() {
// 	shader->bind();
// 	vao->bind();
// 	attributes.push_back(make_shared<AttributeBuffer>("position", VEC3, 0));
// 	attributes.push_back(make_shared<AttributeBuffer>("normal", VEC3, 1));
// 	attributes.push_back(make_shared<AttributeBuffer>("uv", VEC2, 2));
// 	int nextInputNumber = 2;
// 	for (const string& attribute : mesh->getActiveExtraBuffers()) {
// 		attributes.push_back(make_shared<AttributeBuffer>(attribute, VEC4, ++nextInputNumber));
// 	}
//
// 	for (const auto& attribute : attributes) {
// 		vao->addAttributeBuffer(attribute);
//
// 		attribute->load(
// 			mesh->getBufferLocation(attribute->name),
// 			mesh->getBufferSize(attribute->name)
// 		);
// 	}
//
// 	elementBuffer = make_shared<ElementBuffer>();
// 	vao->addElementBuffer(elementBuffer);
// 	elementBuffer->load(
// 		mesh->bufferIndexLocation(),
// 		mesh->bufferIndexSize()
// 	);
// }
//
//
// void MeshLayer::addUniform(string uniformName, GLSLType uniformType, shared_ptr<std::function<void(float, shared_ptr<ShaderProgram>)>> setter) {
// 	this->uniforms[uniformName] = uniformType;
// 	this->uniformSetters[uniformName] = std::move(setter);
// 	this->shader->initUniforms({{uniformName, uniformType}});
// }
//
// void MeshLayer::addConstFloats(const std::map<string, float>& uniforms) {
// 	for (auto uniform : uniforms)
// 		addUniform(uniform.first, FLOAT, make_shared<std::function<void(float, shared_ptr<ShaderProgram>)>>([uniform](float t, const shared_ptr<ShaderProgram>& shader) {
// 			shader->setUniform(uniform.first, uniform.second);
// 		}));
// }
//
// void MeshLayer::addConstVec4(const string& name, vec4 value) {
// 	auto uniformSetter = [value, name](float t, const shared_ptr<ShaderProgram>& shader) {
// 		shader->setUniform(name, value);
// 	};
// 	addUniform(name, VEC4, make_shared<std::function<void(float, shared_ptr<ShaderProgram>)>>(uniformSetter));
// }
//
// void MeshLayer::addConstColor(const string& name, vec4 value) {
// 	addConstVec4(name, value);
// }
//
// void MeshLayer::addMaterialUniforms() {
// 	vec4 intenc = material->compressIntencities();
// 	auto materialSetter = [intenc](float t, const shared_ptr<ShaderProgram>& s) {
// 		s->setUniform("intencities", intenc);
// 	};
// 	addUniform("intencities", VEC4, make_shared<std::function<void(float, shared_ptr<ShaderProgram>)>>(materialSetter));
//
// 	auto textureSetterAmbient = [m=material](float t, const shared_ptr<ShaderProgram>& s) {
// 		m->texture_ambient->bind();
// 		s->setTextureSampler(m->texture_ambient.get());
// 	};
// 	addUniform("texture_ambient", SAMPLER2D, make_shared<std::function<void(float, shared_ptr<ShaderProgram>)>>(textureSetterAmbient));
//
// 	auto textureSetterDiff = [m=material](float t, const shared_ptr<ShaderProgram>& s) {
// 		m->texture_diffuse->bind();
// 		s->setTextureSampler(m->texture_diffuse.get());
// 	};
// 	addUniform("texture_diffuse", SAMPLER2D, make_shared<std::function<void(float, shared_ptr<ShaderProgram>)>>(textureSetterDiff));
//
// 	auto textureSetterSpec = [m=material](float t, const shared_ptr<ShaderProgram>& s) {
// 		m->texture_specular->bind();
// 		s->setTextureSampler(m->texture_specular.get());
// 	};
// 	addUniform("texture_specular", SAMPLER2D, make_shared<std::function<void(float, shared_ptr<ShaderProgram>)>>(textureSetterSpec));
// }
//
// void MeshLayer::setUniforms(float t) {
// 	for (const auto& key : uniforms | std::views::keys) {
// 		auto name = key;
// 		(*uniformSetters[name])(t, shader);
// 	}
// }
//
// void MeshLayer::addCameraUniforms(const shared_ptr<Camera>& camera) {
//
// 	auto MVPsetter = [camera, this](float t, const shared_ptr<ShaderProgram>& shader) {
// 		mat4 mvp = camera->mvp(t, mat4(mat3(1)));
// 		shader->setUniform("mvp", mvp);
// 	};
//
// 	auto positionSetter = [camera, this](float t, const shared_ptr<ShaderProgram>& shader) {
// 		vec3 camPos = camera->position(t);
// 		shader->setUniform("camPosition", camPos);
// 	};
// 	addUniform("mvp", MAT4, make_shared<std::function<void(float, shared_ptr<ShaderProgram>)>>(MVPsetter));
// 	addUniform("camPosition", VEC3, make_shared<std::function<void(float, shared_ptr<ShaderProgram>)>>(positionSetter));
// }
//
// void MeshLayer::addLightUniform(const Light& pointLight, int lightIndex) {
// 	string lightName = "light" + to_string(lightIndex);
// 	auto lightSetter = [pointLight, lightName](float t, const shared_ptr<ShaderProgram>& shader) {
// 		mat4 lightMat = pointLight.compressToMatrix();
// 		shader->setUniform(lightName, lightMat);
// 	};
// 	addUniform(lightName, MAT4, make_shared<std::function<void(float, shared_ptr<ShaderProgram>)>>(lightSetter));
// }
//
// void MeshLayer::addLightsUniforms(const vector<Light>& lights) {
// 	for (int i = 1; i <= lights.size(); i++)
// 		addLightUniform(lights[i - 1], i);
// }
//
//
// void MeshLayer::setScene(const shared_ptr<Camera>& cam, const vector<Light>& lights) {
// 	shader->bind();
// 	initAttributes();
// 	addMaterialUniforms();
// 	addCameraUniforms(cam);
// 	addLightsUniforms(lights);
// }
//
// void MeshLayer::addUniforms(const std::map<string, GLSLType>& uniforms, std::map<string, shared_ptr<std::function<void(float, shared_ptr<ShaderProgram>)>>> setters) {
// 	for (const auto& uniform : uniforms)
// 		addUniform(uniform.first, uniform.second, setters[uniform.first]);
// }
//
// void MeshLayer::addCustomAction(const std::function<void(float)>& action) {
// 	this->customStep = action;
// }
//
//
// void MeshLayer::renderStep(float t) {
// 	shader->bind();
// 	vao->bind();
// 	customStep(t);
// 	setUniforms(t);
// 	glDrawElements(GL_TRIANGLES, mesh->bufferIndexSize(), GL_UNSIGNED_INT, nullptr);
// }



Renderer::Renderer(const RenderSettings& settings)
: Renderer(settings.speed, settings.bgColor, settings.screenshotDirectory.string(), settings.screenshotFrequency) {
	this->settings = settings;
	this->animSpeed = [settings](float t) {
		return settings.speed;
	};
}

Renderer::~Renderer() {
	if (window != nullptr)
		window->destroy();
}

void Renderer::addLayer(std::shared_ptr<Layer> layer) {
	layerStack.push_back(std::move(layer));
}

void Renderer::initRendering() {
	for (const auto& layer : layerStack)
		layer->init();
}

void Renderer::update() {
	auto [t, dt] = fpsClock.tick();
	perFrameFunction(t, dt);
	for (const auto& layer : layerStack)
		layer->updateStep(t, dt);
}

void Renderer::renderStep() const {
	GLCommand::clearScreen(settings.bgColor);
	for (const auto& layer : layerStack)
		layer->renderStep();
	window->renderFramebufferToScreen();
	glfwPollEvents();
}

void Renderer::initMainWindow(int width, int height, const char* title) {
	this->window = make_unique<Window>(width, height, title);
	glewExperimental = true;
	THROW_IF(glewInit() != GLEW_OK, SystemError, "Failed to initialize GLEW");
	LOG("Window initialized with size: " + to_string(width) + "x" + to_string(height));
}

void Renderer::initMainWindow(Resolution resolution, const char* title) {
	initMainWindow(predefinedWidth(resolution), predefinedHeight(resolution), title);
}

void Renderer::addRenderingStep(shared_ptr<MeshLayer> renderingStep) {
	layerStack.push_back(std::move(renderingStep));
}

void Renderer::addMeshStep(const ShaderProgram& shader, const shared_ptr<IndexedMesh>& model, const shared_ptr<MaterialPhong>& material) {
	auto renderingStep = make_shared<MeshLayer>(make_shared<ShaderProgram>(shader), material, model);
	addRenderingStep(renderingStep);
}


void Renderer::setCamera(const shared_ptr<Camera>& camera) {
	this->camera = camera;
}

void Renderer::setLights(const vector<Light>& lights) {
	this->lights = lights;
}

void Renderer::setLightWithMesh(const Light& light, const shared_ptr<MaterialPhong>& material, const ShaderProgram& shader, float radius) {
	this->lights.push_back(light);
	addMeshStep(shader, make_shared<IndexedMesh>(icosphere(radius, 2, light.getPosition(), randomID())), material);
}

void Renderer::setLightsWithMesh(const vector<Light>& lights, const shared_ptr<MaterialPhong>& material, const ShaderProgram& shader, float radius) {
	for (const auto& light : lights)
		setLightWithMesh(light, material, shader, radius);
}

void Renderer::setLightWithMesh(const Light& light, float ambient, float diff, float spec, float shine, const ShaderProgram& shader, float radius) {
	setLightWithMesh(light, make_shared<MaterialPhong>(light.getColor(), light.getColor(), WHITE, ambient, diff, spec, shine), shader, radius);
}

void Renderer::setLightsWithMesh(const vector<Light>& lights_, float ambient, float diff, float spec, float shine, const ShaderProgram& shader, float radius) {
	for (const auto& light : lights_)
		setLightWithMesh(light, ambient, diff, spec, shine, shader, radius);
}

float Renderer::initFrame() {


	float real_dt = 0;
	while (real_dt * settings.maxFPS < 1.0f) {
		auto t = glfwGetTime();
		real_dt += t - last_time_capture;
		last_time_capture = t;
	}
	dt = real_dt * animSpeed(time);
	time += dt;
	return time;
}


float Renderer::lastDeltaTime() const {
	return dt;
}

void Renderer::addPerFrameUniforms(const std::map<string, GLSLType>& uniforms,
								   const std::map<string, shared_ptr<std::function<void(float, shared_ptr<ShaderProgram>)>>>& setters) const {
	for (const auto& renderingStep : layerStack)
		renderingStep->addUniforms(uniforms, setters);
}

void Renderer::addPerFrameUniform(const string& uniformName, GLSLType uniformType, const shared_ptr<std::function<void(float, shared_ptr<ShaderProgram>)>>& setter) const {
	for (const auto& renderingStep : layerStack)
		renderingStep->addUniform(uniformName, uniformType, setter);
}

void Renderer::addSurfaceFamilyDeformer(SurfaceParametricPencil& pencil, IndexedMesh& surface) {
	addCustomAction([&pencil, &surface](float t) {
		surface.adjustToNewSurface(pencil(t));
	});
}

void Renderer::initRendering() {
	window->initViewport();

	glShadeModel(GL_FLAT);

	if (settings.depthTest) {
		glEnable(GL_DEPTH_TEST);
		glDepthFunc(GL_LESS);
	}

	if (settings.alphaBlending) {
		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	}

	if (settings.timeUniform) {
		addTimeUniform();
	}

	for (const auto& renderingStep : layerStack)
		renderingStep->setScene(camera, lights);

	LOG("OpenGL renderer initialized: " + string((const char *)glGetString(GL_RENDERER)));
}

void Renderer::addConstUniform(const string& uniformName, GLSLType uniformType, shared_ptr<std::function<void(shared_ptr<ShaderProgram>)>> setter) const {
	std::function<void(float, shared_ptr<ShaderProgram>)> setterWrapper = [setter](float t, const shared_ptr<ShaderProgram>& shader) {
		(*setter)(shader);
	};
	for (const auto& renderingStep : layerStack)
		renderingStep->addUniform(uniformName, uniformType, make_shared<std::function<void(float, shared_ptr<ShaderProgram>)>>(setterWrapper));
}

void Renderer::addTimeUniform() const {
	addPerFrameUniform("time", FLOAT, make_shared<std::function<void(float, shared_ptr<ShaderProgram>)>>([](float t, const shared_ptr<ShaderProgram>& shader) {
		shader->setUniform("time", t);
	}));
}

void Renderer::addConstFloats(const std::map<string, float>& uniforms) const {
	for (const auto& renderingStep : layerStack)
		renderingStep->addConstFloats(uniforms);
}

void Renderer::addCustomAction(std::function<void(float)> action) {
	perFrameFunction = [a=perFrameFunction, n=action](float t, float delta) {
		a(t, delta);
		n(t);
	};
}

void Renderer::addCustomAction(std::function<void(float, float)> action) {
	perFrameFunction = [a=perFrameFunction, n=action](float t, float delta) {
		a(t, delta);
		n(t, delta);
	};
}

void Renderer::nonlinearSpeed(const std::function<float(float)>& speed) {
	animSpeed = speed;
}

void Renderer::addConstUniforms(const std::map<string, GLSLType>& uniforms, std::map<string, shared_ptr<std::function<void(shared_ptr<ShaderProgram>)>>> setters) const {
	for (const auto& uniform : uniforms)
		addConstUniform(uniform.first, uniform.second, setters[uniform.first]);
}

void Renderer::renderAllSteps() {
	for (const auto& renderingStep : layerStack)
		renderingStep->renderStep(time);
}

int Renderer::mainLoop() {
	initRendering();
	fpsClock.reset();
	while (window->isOpen()) {
		auto [time, dt] = fpsClock.tick();
		initFrame();
		perFrameFunction(time, dt);
		renderAllSteps();
		window->renderFramebufferToScreen();
	}
	return window->destroy();
}
