#include "glslUtils.hpp"
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

#include <stdlib.h>
#include <string.h>

#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include "configFiles.hpp"
#include "specific.hpp"

using namespace glm;


// void LOG(const string &message) {
// 	std::cout << message << std::endl;
// }
//
// void LOG_ERROR(string message) {
// 	std::cerr << message << std::endl;
// }

int predefinedWidth(Resolution res) {
	switch (res) {
	case FHD:
		return 1920;
	case UHD:
		return 3840;
	case HD2K:
		return 2560;
	}
	throw UnknownVariantError("Resolution not recognized", __FILE__, __LINE__);
}

int predefinedHeight(Resolution res) {
	switch (res) {
	case FHD:
		return 1080;
	case UHD:
		return 2160;
	case HD2K:
		return 1440;
	}
	throw UnknownVariantError("Resolution not recognized", __FILE__, __LINE__);
}


GLenum Shader::getTypeFromExtension(const string& extension) {
	string nodot = extension;
	if (nodot[0] == '.')
		nodot = nodot.substr(1);

	if (nodot == "vert")
		return GL_VERTEX_SHADER;
	if (nodot == "frag")
		return GL_FRAGMENT_SHADER;
	if (nodot == "geom" || nodot == "geo")
		return GL_GEOMETRY_SHADER;
	throw SystemError("Unknown shader getExtension: ." + nodot, __FILE__, __LINE__);
}


GLuint Shader::getID() const {
	return shaderID;
}

void Shader::compile() {
	shaderID = glCreateShader(shaderType);
	GLint Result = GL_FALSE;
	int InfoLogLength;

	char const* sourcePtr = code.c_str();

	glShaderSource(shaderID, 1, &sourcePtr, nullptr);
	glCompileShader(shaderID);
	glGetShaderiv(shaderID, GL_COMPILE_STATUS, &Result);
	glGetShaderiv(shaderID, GL_INFO_LOG_LENGTH, &InfoLogLength);
	if (InfoLogLength > 0) {
		vector<char> error(InfoLogLength + 1);
		glGetShaderInfoLog(shaderID, InfoLogLength, nullptr, &error[0]);
		throw SystemError(&error[0], __FILE__, __LINE__);
	}
}

string Shader::getCode() const {
	return code;
}

Shader::Shader(const string& code, GLenum shaderType)
: shaderType(shaderType), code(code) {
	shaderID = glCreateShader(shaderType);
}


Shader::Shader(CodeFileDescriptor&& file) {
	shaderType = getTypeFromExtension(file.extension());
	code = file.getCode();
	shaderID = glCreateShader(shaderType);
}

Shader::Shader(CodeFileDescriptor& file)
: Shader(file.getCode(), getTypeFromExtension(file.extension())) {}

Shader::Shader(const Shader& other)
: shaderType(other.shaderType), shaderID(other.shaderID), code(other.code) {}

Shader::Shader(Shader&& other) noexcept
: shaderType(other.shaderType), shaderID(other.shaderID), code(std::move(other.code)) {}

Shader& Shader::operator=(const Shader& other) {
	if (this == &other)
		return *this;
	shaderType = other.shaderType;
	shaderID = other.shaderID;
	code = other.code;
	return *this;
}

Shader& Shader::operator=(Shader&& other) noexcept {
	if (this == &other)
		return *this;
	shaderType = other.shaderType;
	shaderID = other.shaderID;
	code = std::move(other.code);
	return *this;
}


Window::Window(int width, int height, const char* title) {
	// glfwInit();
	glfwWindowHint(GLFW_SAMPLES, 4);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
	this->width = width;
	this->height = height;
	this->aspectRatio = (float)width / (float)height;
	this->window = glfwCreateWindow(width, height, title, nullptr, nullptr);
	if (!this->window) {
		// glfwTerminate(); // removed; handled by Renderer
		exit(2136);
	}
	glfwMakeContextCurrent(this->window);
	glfwGetFramebufferSize(window, &width, &height);
}

Window::Window(Resolution resolution, const char* title) {
	// glfwInit();
	glfwWindowHint(GLFW_SAMPLES, 4);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
	width = predefinedWidth(resolution);
	height = predefinedHeight(resolution);
	aspectRatio = (float)width / height;
	window = glfwCreateWindow(width, height, title, nullptr, nullptr);
	if (!window) {
		// glfwTerminate(); // removed; handled by Renderer
		throw SystemError("GLFW window creation failed", __FILE__, __LINE__);
	}
	glfwMakeContextCurrent(window);
}

Window::~Window() {
	destroy();
}

int Window::destroy() {
	glfwDestroyWindow(this->window);
	// glfwTerminate(); // removed; handled by Renderer
	return 0;
}

void Window::renderFramebufferToScreen() {
	glfwSwapBuffers(this->window);
	glfwPollEvents();
}

void Window::showCursor() {
	glfwSetInputMode(this->window, GLFW_CURSOR, GLFW_CURSOR_NORMAL);
}

void Window::disableCursor() {
	glfwSetInputMode(this->window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
}

void Window::hideCursorWithinWindow() {
	glfwSetInputMode(this->window, GLFW_CURSOR, GLFW_CURSOR_HIDDEN);
}

void Window::stickyKeys(bool sticky) {
	if (sticky)
		glfwSetInputMode(this->window, GLFW_STICKY_KEYS, GL_TRUE);
	else
		glfwSetInputMode(this->window, GLFW_STICKY_KEYS, GL_FALSE);
}

void Window::stickyMouseButtons(bool sticky) {
	if (sticky)
		glfwSetInputMode(this->window, GLFW_STICKY_MOUSE_BUTTONS, GL_TRUE);
	else
		glfwSetInputMode(this->window, GLFW_STICKY_MOUSE_BUTTONS, GL_FALSE);
}

void Window::setCallbacks(const GLFWkeyfun* keyCallback, const GLFWcharfun* charCallback, const GLFWmousebuttonfun* mouseButtonCallback, GLFWcursorposfun* cursorPosCallback,
						  GLFWcursorenterfun* cursorEnterCallback, GLFWscrollfun* scrollCallback, GLFWdropfun* dropCallback) {
	if (keyCallback != nullptr)
		glfwSetKeyCallback(this->window, *keyCallback);
	if (charCallback != nullptr)
		glfwSetCharCallback(this->window, *charCallback);
	if (mouseButtonCallback != nullptr)
		glfwSetMouseButtonCallback(this->window, *mouseButtonCallback);
	if (cursorPosCallback != nullptr)
		glfwSetCursorPosCallback(this->window, *cursorPosCallback);
	if (cursorEnterCallback != nullptr)
		glfwSetCursorEnterCallback(this->window, *cursorEnterCallback);
	if (scrollCallback != nullptr)
		glfwSetScrollCallback(this->window, *scrollCallback);
	if (dropCallback != nullptr)
		glfwSetDropCallback(this->window, *dropCallback);
}

bool Window::isOpen() {
	return !glfwWindowShouldClose(this->window);
}

void Window::initViewport() {
	glViewport(0, 0, this->width, this->height);
}


void ShaderProgram::linkShaders() {
	GLint Result = GL_FALSE;
	int InfoLogLength;

	vertexShader.compile();
	fragmentShader.compile();
	if (geometryShader)
		geometryShader->compile();


	GLuint VertexShaderID = vertexShader.getID();
	GLuint FragmentShaderID = fragmentShader.getID();
	GLuint GeometryShaderID = geometryShader ? geometryShader->getID() : 0;

	programID = glCreateProgram();

	glAttachShader(programID, VertexShaderID);
	glAttachShader(programID, FragmentShaderID);
	if (geometryShader)
		glAttachShader(programID, GeometryShaderID);

	glLinkProgram(programID);
	glGetProgramiv(programID, GL_LINK_STATUS, &Result);
	glGetProgramiv(programID, GL_INFO_LOG_LENGTH, &InfoLogLength);

	if (InfoLogLength > 0) {
		vector<char> ProgramErrorMessage(InfoLogLength + 1);
		glGetProgramInfoLog(programID, InfoLogLength, NULL, &ProgramErrorMessage[0]);
		throw SystemError(&ProgramErrorMessage[0], __FILE__, __LINE__);
	}

	glDetachShader(programID, VertexShaderID);
	if (geometryShader)
		glDetachShader(programID, GeometryShaderID);
	glDetachShader(programID, FragmentShaderID);

	glDeleteShader(VertexShaderID);
	if (geometryShader)
		glDeleteShader(GeometryShaderID);
	glDeleteShader(FragmentShaderID);

	shaderType = geometryShader ? GEOMETRY1 : CLASSIC;
}


ShaderProgram::ShaderProgram(const Shader& vertexShader, const Shader& fragmentShader, const Shader& geometryShader)
: vertexShader(vertexShader), fragmentShader(fragmentShader), geometryShader(geometryShader) {
	linkShaders();
}

ShaderProgram::ShaderProgram(const string& vertexPath, const string& fragPath)
: vertexShader(CodeFileDescriptor(vertexPath, false)), fragmentShader(CodeFileDescriptor(fragPath, false)), geometryShader(std::nullopt) {
	linkShaders();
}

ShaderProgram::ShaderProgram(const Shader& vertexShader, const Shader& fragmentShader)
: vertexShader(vertexShader), fragmentShader(fragmentShader), geometryShader(std::nullopt) {
	linkShaders();
}


ShaderProgram::~ShaderProgram() {
	glDeleteProgram(this->programID);
}

void ShaderProgram::use() {
	glUseProgram(this->programID);
}

void ShaderProgram::initUniforms(const unordered_map<string, GLSLType>& uniforms) {
	for (auto uni : uniforms)
		this->uniformTypes[uni.first] = uni.second;

	for (const auto& key : uniforms | std::views::keys) {
		GLuint uniformLocation = glGetUniformLocation(this->programID, key.c_str());
		this->uniformLocations[key] = uniformLocation;
	}
}

void ShaderProgram::initTextureSampler(const Texture* texture) {
	string name = texture->samplerName;
	uniformTypes[name] = SAMPLER2D;
	GLuint uniformLocation = glGetUniformLocation(this->programID, texture->samplerName);
	uniformLocations[name] = uniformLocation;
}

void ShaderProgram::setTextureSampler(const Texture* texture) const {
	string name = texture->samplerName;
	GLuint uniformLocation = this->uniformLocations.at(name);
	glUniform1i(uniformLocation, texture->abs_slot);
}

void ShaderProgram::setUniforms(const unordered_map<string, const GLfloat*>& uniformValues) {
	for (auto const& uniform : uniformValues)
		setUniform(uniform.first, uniform.second);
}

void ShaderProgram::setUniform(const string& uniformName, const GLfloat* uniformValue) {
	GLSLType uniformType = this->uniformTypes[uniformName];
	GLuint uniformLocation = this->uniformLocations[uniformName];

	switch (uniformType) {
	case FLOAT: glUniform1fv(uniformLocation, 1, uniformValue);
		break;
	case INT: glUniform1iv(uniformLocation, 1, (GLint*)uniformValue);
		break;
	case VEC2: glUniform2fv(uniformLocation, 1, uniformValue);
		break;
	case VEC3: glUniform3fv(uniformLocation, 1, uniformValue);
		break;
	case VEC4: glUniform4fv(uniformLocation, 1, uniformValue);
		break;
	case MAT2: glUniformMatrix2fv(uniformLocation, 1, GL_FALSE, uniformValue);
		break;
	case MAT3: glUniformMatrix3fv(uniformLocation, 1, GL_FALSE, uniformValue);
		break;
	case MAT4: glUniformMatrix4fv(uniformLocation, 1, GL_FALSE, uniformValue);
		break;

	default:
		throw std::invalid_argument("Uniform type not recognized");
	}
}

void ShaderProgram::setUniform(const string& uniformName, float uniformValue) {
	if (this->uniformTypes[uniformName] != FLOAT)
		throw std::invalid_argument("Uniform type must be FLOAT");
	GLuint uniformLocation = this->uniformLocations[uniformName];
	glUniform1f(uniformLocation, uniformValue);
}

void ShaderProgram::setUniform(const string& uniformName, int uniformValue) {
	if (this->uniformTypes[uniformName] != INT)
		throw std::invalid_argument("Uniform type must be INT");
	GLuint uniformLocation = this->uniformLocations[uniformName];
	glUniform1i(uniformLocation, uniformValue);
}

void ShaderProgram::setUniform(const string& uniformName, vec2 uniformValue) {
	if (this->uniformTypes[uniformName] != VEC2)
		throw std::invalid_argument("Uniform type must be VEC2");
	GLuint uniformLocation = this->uniformLocations[uniformName];
	glUniform2f(uniformLocation, uniformValue.x, uniformValue.y);
}

void ShaderProgram::setUniform(const string& uniformName, vec3 uniformValue) {
	if (this->uniformTypes[uniformName] != VEC3)
		throw std::invalid_argument("Uniform type must be VEC3");
	GLuint uniformLocation = this->uniformLocations[uniformName];
	glUniform3f(uniformLocation, uniformValue.x, uniformValue.y, uniformValue.z);
}

void ShaderProgram::setUniform(const string& uniformName, vec4 uniformValue) {
	if (this->uniformTypes[uniformName] != VEC4)
		throw std::invalid_argument("Uniform type must be VEC4");
	GLuint uniformLocation = this->uniformLocations[uniformName];
	glUniform4f(uniformLocation, uniformValue.x, uniformValue.y, uniformValue.z, uniformValue.w);
}

void ShaderProgram::setUniform(const string& uniformName, mat2 uniformValue) {
	if (this->uniformTypes[uniformName] != MAT2)
		throw std::invalid_argument("Uniform type must be MAT2");
	GLuint uniformLocation = this->uniformLocations[uniformName];
	glUniformMatrix2fv(uniformLocation, 1, GL_FALSE, &uniformValue[0][0]);
}

void ShaderProgram::setUniform(const string& uniformName, mat3 uniformValue) {
	if (this->uniformTypes[uniformName] != MAT3)
		throw std::invalid_argument("Uniform type must be MAT3");
	GLuint uniformLocation = this->uniformLocations[uniformName];
	glUniformMatrix3fv(uniformLocation, 1, GL_FALSE, &uniformValue[0][0]);
}

void ShaderProgram::setUniform(const string& uniformName, mat4 uniformValue) {
	if (this->uniformTypes[uniformName] != MAT4)
		throw std::invalid_argument("Uniform type must be MAT4");
	GLuint uniformLocation = this->uniformLocations[uniformName];
	glUniformMatrix4fv(uniformLocation, 1, GL_FALSE, &uniformValue[0][0]);
}

void ShaderProgram::setUniform(const string& uniformName, float x, float y) {
	if (this->uniformTypes[uniformName] != VEC2)
		throw std::invalid_argument("Uniform type must be VEC2");
	GLuint uniformLocation = this->uniformLocations[uniformName];
	glUniform2f(uniformLocation, x, y);
}

void ShaderProgram::setUniform(const string& uniformName, float x, float y, float z) {
	if (this->uniformTypes[uniformName] != VEC3)
		throw std::invalid_argument("Uniform type must be VEC3");
	GLuint uniformLocation = this->uniformLocations[uniformName];
	glUniform3f(uniformLocation, x, y, z);
}

void ShaderProgram::setUniform(const string& uniformName, float x, float y, float z, float w) {
	if (this->uniformTypes[uniformName] != VEC4)
		throw std::invalid_argument("Uniform type must be VEC4");
	GLuint uniformLocation = this->uniformLocations[uniformName];
	glUniform4f(uniformLocation, x, y, z, w);
}




RenderingStep::RenderingStep(const shared_ptr<ShaderProgram>& shader, const shared_ptr<MaterialPhong>& material, const shared_ptr<IndexedMesh>& mesh)
: customStep([](float){}), elementBuffer(nullptr)
{
	this->shader = shader;
	this->material = material;
	this->mesh = mesh;
	vao = make_shared<VertexArray>();
}


RenderingStep::~RenderingStep() {
	for (auto attribute : attributes)
		attribute.reset();
	shader.reset();
}



void RenderingStep::initMaterialTextures() const {
	material->initTextures();
	shader->initTextureSampler(material->texture_ambient.get());
	shader->initTextureSampler(material->texture_diffuse.get());
	shader->initTextureSampler(material->texture_specular.get());
}

void RenderingStep::bindTextures() const {
	material->bindTextures();
	shader->setTextureSampler(material->texture_ambient.get());
	shader->setTextureSampler(material->texture_diffuse.get());
	shader->setTextureSampler(material->texture_specular.get());
}



void RenderingStep::initAttributes() {
	shader->use();
	vao->bind();
	attributes.push_back(make_shared<AttributeBuffer>("position", VEC3, 0));
	attributes.push_back(make_shared<AttributeBuffer>("normal", VEC3, 1));
	attributes.push_back(make_shared<AttributeBuffer>("uv", VEC2, 2));
	int nextInputNumber = 2;
	for (const string& attribute : mesh->getActiveExtraBuffers()) {
		attributes.push_back(make_shared<AttributeBuffer>(attribute, VEC4, ++nextInputNumber));
	}

	for (const auto& attribute : attributes) {
		vao->addAttributeBuffer(*attribute);

		attribute->load(
			mesh->getBufferLocation(attribute->name),
			mesh->getBufferSize(attribute->name)
		);
	}

	elementBuffer = make_shared<ElementBuffer>();
	vao->addElementBuffer(*elementBuffer);
	elementBuffer->load(
		mesh->bufferIndexLocation(),
		mesh->bufferIndexSize()
	);
}


void RenderingStep::addUniform(string uniformName, GLSLType uniformType, shared_ptr<std::function<void(float, shared_ptr<ShaderProgram>)>> setter) {
	this->uniforms[uniformName] = uniformType;
	this->uniformSetters[uniformName] = std::move(setter);
	this->shader->initUniforms({{uniformName, uniformType}});
}

void RenderingStep::addConstFloats(const std::map<string, float>& uniforms) {
	for (auto uniform : uniforms)
		addUniform(uniform.first, FLOAT, make_shared<std::function<void(float, shared_ptr<ShaderProgram>)>>([uniform](float t, const shared_ptr<ShaderProgram>& shader) {
			shader->setUniform(uniform.first, uniform.second);
		}));
}

void RenderingStep::addConstVec4(const string& name, vec4 value) {
	auto uniformSetter = [value, name](float t, const shared_ptr<ShaderProgram>& shader) {
		shader->setUniform(name, value);
	};
	addUniform(name, VEC4, make_shared<std::function<void(float, shared_ptr<ShaderProgram>)>>(uniformSetter));
}

void RenderingStep::addConstColor(const string& name, vec4 value) {
	addConstVec4(name, value);
}

void RenderingStep::addMaterialUniforms() {
	vec4 intenc = material->compressIntencities();
	auto materialSetter = [intenc](float t, const shared_ptr<ShaderProgram>& s) {
		s->setUniform("intencities", intenc);
	};
	addUniform("intencities", VEC4, make_shared<std::function<void(float, shared_ptr<ShaderProgram>)>>(materialSetter));

	auto textureSetterAmbient = [m=material](float t, const shared_ptr<ShaderProgram>& s) {
		m->texture_ambient->bind();
		s->setTextureSampler(m->texture_ambient.get());
	};
	addUniform("texture_ambient", SAMPLER2D, make_shared<std::function<void(float, shared_ptr<ShaderProgram>)>>(textureSetterAmbient));

	auto textureSetterDiff = [m=material](float t, const shared_ptr<ShaderProgram>& s) {
		m->texture_diffuse->bind();
		s->setTextureSampler(m->texture_diffuse.get());
	};
	addUniform("texture_diffuse", SAMPLER2D, make_shared<std::function<void(float, shared_ptr<ShaderProgram>)>>(textureSetterDiff));

	auto textureSetterSpec = [m=material](float t, const shared_ptr<ShaderProgram>& s) {
		m->texture_specular->bind();
		s->setTextureSampler(m->texture_specular.get());
	};
	addUniform("texture_specular", SAMPLER2D, make_shared<std::function<void(float, shared_ptr<ShaderProgram>)>>(textureSetterSpec));
}

void RenderingStep::setUniforms(float t) {
	for (const auto& key : uniforms | std::views::keys) {
		auto name = key;
		(*uniformSetters[name])(t, shader);
	}
}

void RenderingStep::addCameraUniforms(const shared_ptr<Camera>& camera) {

	auto MVPsetter = [camera, this](float t, const shared_ptr<ShaderProgram>& shader) {
		mat4 mvp = camera->mvp(t, mat4(mat3(1)));
		shader->setUniform("mvp", mvp);
	};

	auto positionSetter = [camera, this](float t, const shared_ptr<ShaderProgram>& shader) {
		vec3 camPos = camera->position(t);
		shader->setUniform("camPosition", camPos);
	};
	addUniform("mvp", MAT4, make_shared<std::function<void(float, shared_ptr<ShaderProgram>)>>(MVPsetter));
	addUniform("camPosition", VEC3, make_shared<std::function<void(float, shared_ptr<ShaderProgram>)>>(positionSetter));
}

void RenderingStep::addLightUniform(const Light& pointLight, int lightIndex) {
	string lightName = "light" + to_string(lightIndex);
	auto lightSetter = [pointLight, lightName](float t, const shared_ptr<ShaderProgram>& shader) {
		mat4 lightMat = pointLight.compressToMatrix();
		shader->setUniform(lightName, lightMat);
	};
	addUniform(lightName, MAT4, make_shared<std::function<void(float, shared_ptr<ShaderProgram>)>>(lightSetter));
}

void RenderingStep::addLightsUniforms(const vector<Light>& lights) {
	for (int i = 1; i <= lights.size(); i++)
		addLightUniform(lights[i - 1], i);
}


void RenderingStep::init(const shared_ptr<Camera>& cam, const vector<Light>& lights) {
	shader->use();
	initAttributes();
	addMaterialUniforms();
	addCameraUniforms(cam);
	addLightsUniforms(lights);
}

void RenderingStep::addUniforms(const std::map<string, GLSLType>& uniforms, std::map<string, shared_ptr<std::function<void(float, shared_ptr<ShaderProgram>)>>> setters) {
	for (const auto& uniform : uniforms)
		addUniform(uniform.first, uniform.second, setters[uniform.first]);
}

void RenderingStep::addCustomAction(const std::function<void(float)>& action) {
	this->customStep = action;
}


void RenderingStep::renderStep(float t) {
	shader->use();
	vao->bind();
	customStep(t);
	setUniforms(t);
	glDrawElements(GL_TRIANGLES, mesh->bufferIndexSize(), GL_UNSIGNED_INT, nullptr);
}

RenderSettings::RenderSettings(vec4 bgColor, bool alphaBlending, bool depthTest, bool timeUniform, float speed, int maxFPS, bool takeScreenshots, Resolution resolution,
							   float screenshotFrequency, const string& windowTitle)
: bgColor(bgColor), alphaBlending(alphaBlending), depthTest(depthTest), timeUniform(timeUniform), speed(speed), maxFPS(maxFPS), resolution(resolution),
  windowTitle(windowTitle), takeScreenshots(takeScreenshots), screenshotFrequency(screenshotFrequency) {
	if (takeScreenshots)
		throw NotImplementedError("Screenshot functionality not implemented yet.", __FILE__, __LINE__);
	this->screenshotFrequency = -1;

	screenshotDirectory = ConfigFile().getScreenshotsDir();
}

Renderer::Renderer(float animSpeed, vec4 bgColor, const string& screenshotDirectory, float screenshotFrequency)
: settings(bgColor, true, true, true, animSpeed, 120, false, UNKNOWN, screenshotFrequency) {
	logging::Logger::init();

	glfwSetErrorCallback([](int error, const char* description) {
		LOG_ERROR("GLFW Error (" + to_string(error) + "): " + string(description));
	});
	THROW_IF(not glfwInit(), SystemError, "Failed to initialize GLFW");
	this->window = nullptr;
	this->camera = nullptr;
	this->time = 0;
	this->lights = vector<Light>();
	this->animSpeed = [animSpeed](float t) {
		return animSpeed;
	};
	this->perFrameFunction = [](float t, float delta) {};
}

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
	if (camera != nullptr)
		camera.reset();
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

void Renderer::initMainWindow() {
	initMainWindow(predefinedWidth(settings.resolution), predefinedHeight(settings.resolution), settings.windowTitle.c_str());
}

void Renderer::resetTimer() {
	last_time_capture = glfwGetTime();
	time = 0;
}

void Renderer::addRenderingStep(shared_ptr<RenderingStep> renderingStep) {
	renderingSteps.push_back(std::move(renderingStep));
}

void Renderer::addMeshStep(const ShaderProgram& shader, const shared_ptr<IndexedMesh>& model, const shared_ptr<MaterialPhong>& material) {
	auto renderingStep = make_shared<RenderingStep>(make_shared<ShaderProgram>(shader), material, model);
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
	glClearColor(settings.bgColor.x, settings.bgColor.y, settings.bgColor.z, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

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
	for (const auto& renderingStep : renderingSteps)
		renderingStep->addUniforms(uniforms, setters);
}

void Renderer::addPerFrameUniform(const string& uniformName, GLSLType uniformType, const shared_ptr<std::function<void(float, shared_ptr<ShaderProgram>)>>& setter) const {
	for (const auto& renderingStep : renderingSteps)
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

	for (const auto& renderingStep : renderingSteps)
		renderingStep->init(camera, lights);

	LOG("OpenGL renderer initialized: " + string((const char *)glGetString(GL_RENDERER)));
}

void Renderer::addConstUniform(const string& uniformName, GLSLType uniformType, shared_ptr<std::function<void(shared_ptr<ShaderProgram>)>> setter) const {
	std::function<void(float, shared_ptr<ShaderProgram>)> setterWrapper = [setter](float t, const shared_ptr<ShaderProgram>& shader) {
		(*setter)(shader);
	};
	for (const auto& renderingStep : renderingSteps)
		renderingStep->addUniform(uniformName, uniformType, make_shared<std::function<void(float, shared_ptr<ShaderProgram>)>>(setterWrapper));
}

void Renderer::addTimeUniform() const {
	addPerFrameUniform("time", FLOAT, make_shared<std::function<void(float, shared_ptr<ShaderProgram>)>>([](float t, const shared_ptr<ShaderProgram>& shader) {
		shader->setUniform("time", t);
	}));
}

void Renderer::addConstFloats(const std::map<string, float>& uniforms) const {
	for (const auto& renderingStep : renderingSteps)
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
	for (const auto& renderingStep : renderingSteps)
		renderingStep->renderStep(time);
}

int Renderer::mainLoop() {
	initRendering();
	resetTimer();
	while (window->isOpen()) {
		initFrame();
		perFrameFunction(time, dt);
		renderAllSteps();
		window->renderFramebufferToScreen();
	}
	return window->destroy();
}
