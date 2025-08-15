
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

#include "specific.hpp"

using namespace glm;


void setUniformTextureSampler(GLuint programID, Texture *texture, int textureSlot)
{
	GLuint samplerUniform = glGetUniformLocation(programID, texture->samplerName);
	glUniform1i(samplerUniform, textureSlot);
    texture->bind();
}

int predefinedWidth(Resolution res) {
	switch (res) {
		case FHD:
			return 1920;
		case UHD:
			return 3840;
		case HD2K:
			return 2560;
	}
	throw UnknownVariantError("Resolution not recognized");
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
	throw UnknownVariantError("Resolution not recognized");

}

size_t sizeOfGLSLType(GLSLType type)
{
    switch (type)
	{
	case FLOAT:
		return sizeof(float);
	case INT:
		return sizeof(int);
	case VEC2:
		return sizeof(vec2);
	case VEC3:
		return sizeof(vec3);
	case VEC4:
		return sizeof(vec4);
	case MAT2:
		return sizeof(mat2);
	case MAT3:
		return sizeof(mat3);
	case MAT4:
		return sizeof(mat4);
	case SAMPLER1D:
	case SAMPLER2D:
	case SAMPLER3D:
		return sizeof(GLuint);
	}
	throw UnknownVariantError("GLSL type not recognized");
}

int lengthOfGLSLType(GLSLType type)
{
   switch (type)
	{
	case FLOAT:
	case INT:
   	case SAMPLER1D:
	case SAMPLER2D:
	case SAMPLER3D:
		return 1;
	case VEC2:
		return 2;
	case VEC3:
		return 3;
	case VEC4:
	case MAT2:
		return 4;
	case MAT3:
		return 9;
	case MAT4:
		return 16;
	}
	std::cout << "Error: unknown GLSLType" << std::endl;
	return -1;
}

GLenum primitiveGLSLType(GLSLType type)
{
    switch (type)
    {
    case FLOAT:
    case VEC2:
    case VEC3:
    case VEC4:
    case MAT2:
    case MAT3:
    case MAT4:
        return GL_FLOAT;
    case INT:
    case SAMPLER1D:
    case SAMPLER2D:
    case SAMPLER3D:
        return GL_INT;

    }
    std::cout << "Error: unknown GLSLType" << std::endl;
    return -1;
}


GLenum Shader::getTypeFromExtension(const string &extension) {
	string nodot = extension;
	if (nodot[0] == '.')
		nodot = nodot.substr(1);

	if (nodot == "vert")
		return GL_VERTEX_SHADER;
	if (nodot == "frag")
		return GL_FRAGMENT_SHADER;
	if (nodot == "geom" || nodot == "geo")
		return GL_GEOMETRY_SHADER;
	throw SystemError("Unknown shader extension: ." + nodot);
}


GLuint Shader::getID() const {
	return shaderID;
}

void Shader::compile() {
	shaderID = glCreateShader(shaderType);
	GLint Result = GL_FALSE;
	int InfoLogLength;

	char const *sourcePtr = code.c_str();

	glShaderSource(shaderID, 1, &sourcePtr, nullptr);
	glCompileShader(shaderID);
	glGetShaderiv(shaderID, GL_COMPILE_STATUS, &Result);
	glGetShaderiv(shaderID, GL_INFO_LOG_LENGTH, &InfoLogLength);
	if (InfoLogLength > 0)
	{
		vector<char> error(InfoLogLength + 1);
		glGetShaderInfoLog(shaderID, InfoLogLength, nullptr, &error[0]);
		throw SystemError(&error[0]);
	}
}

string Shader::getCode() const { return code; }

Shader::Shader(const string &code, GLenum shaderType)
: shaderType(shaderType), code(code)
{
	shaderID = glCreateShader(shaderType);
}


Shader::Shader(CodeFileDescriptor &&file) {
	shaderType = getTypeFromExtension(file.extension());
	code = file.getCode();
	shaderID = glCreateShader(shaderType);
}

Shader::Shader(CodeFileDescriptor &file) : Shader(file.getCode(), getTypeFromExtension(file.extension())) {}

Shader::Shader(const Shader &other): shaderType(other.shaderType),
									 shaderID(other.shaderID),
									 code(other.code) {}

Shader::Shader(Shader &&other) noexcept: shaderType(other.shaderType),
										 shaderID(other.shaderID),
										 code(std::move(other.code)) {}

Shader & Shader::operator=(const Shader &other) {
	if (this == &other)
		return *this;
	shaderType = other.shaderType;
	shaderID   = other.shaderID;
	code       = other.code;
	return *this;
}
Shader & Shader::operator=(Shader &&other) noexcept {
	if (this == &other)
		return *this;
	shaderType = other.shaderType;
	shaderID   = other.shaderID;
	code       = std::move(other.code);
	return *this;
}

ConfigFile::ConfigFile(const string &path): CodeFileDescriptor("config.txt", path, false) {
	int i = 0;
	config = std::map<string, string>();
	while (readLine(i).contains(':')) {
		string key = readLine(i).substr(0, readLine(i).find(':'));
		string value = readLine(i).substr(readLine(i).find(':') + 1);
		std::erase(value, '\t');
		std::erase(value, ' ');
		std::erase(key, ':');
		config[key] = value;
		i++;
	}
}


ConfigFile::ConfigFile(): ConfigFile(DEFAULT_CONFIG_PATH) {}

void error_callback(int error, const char *description)
{
	fprintf(stderr, "Error: %s\n", description);
}

GLuint bindVAO()
{
	GLuint VertexArrayID;
	glGenVertexArrays(1, &VertexArrayID);
	glBindVertexArray(VertexArrayID);
	return VertexArrayID;
}

void disableAttributeArrays(int how_many)
{
	for (int i = 0; i < how_many; i++)
		glDisableVertexAttribArray(i);
}

void key_callback(GLFWwindow *window, int key, int scancode, int action, int mods)
{
	if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS)
		glfwSetWindowShouldClose(window, 1);
}

mat4 generateMVP(vec3 camPosition, vec3 camLookAt, vec3 upVector, float fov, float aspectRatio, float clippingRangeMin, float clippingRangeMax, const mat4 &modelTransform)
{
	mat4 ViewMatrix = lookAt(camPosition, camLookAt, upVector);
	mat4 ProjectionMatrix = perspective(fov, aspectRatio, clippingRangeMin, clippingRangeMax);
	return ProjectionMatrix * ViewMatrix * modelTransform;
}

Window::Window(int width, int height, const char *title)
{
    if (!glfwInit())
        throw SystemError("GLFW initialization failed");
	glfwWindowHint(GLFW_SAMPLES, 4);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
	this->width = width;
	this->height = height;
	this->aspectRatio = (float)width / (float)height;
	this->window = glfwCreateWindow(width, height, title, nullptr, nullptr);
	if (!this->window)
	{
		glfwTerminate();
		exit(2136);
	}
	glfwMakeContextCurrent(this->window);
	glfwGetFramebufferSize(window, &width, &height);
}

Window::Window(Resolution resolution, const char *title)
{
	glfwWindowHint(GLFW_SAMPLES, 4);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
	width = predefinedWidth(resolution);
	height = predefinedHeight(resolution);
	aspectRatio = (float)width / height;
	window = glfwCreateWindow(width, height, title, nullptr, nullptr);
	if (!window)
	{
		glfwTerminate();
		throw SystemError("GLFW window creation failed");
	}
	glfwMakeContextCurrent(window);
}

Window::~Window()
{
	destroy();
}



void Window::renderFramebufferToScreen()
{
	glfwSwapBuffers(this->window);
	glfwPollEvents();
}

void Window::showCursor()
{
	glfwSetInputMode(this->window, GLFW_CURSOR, GLFW_CURSOR_NORMAL);
}

void Window::disableCursor()
{
	glfwSetInputMode(this->window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
}

void Window::hideCursorWithinWindow()
{
	glfwSetInputMode(this->window, GLFW_CURSOR, GLFW_CURSOR_HIDDEN);
}

void Window::stickyKeys(bool sticky)
{
	if (sticky)
		glfwSetInputMode(this->window, GLFW_STICKY_KEYS, GL_TRUE);
	else
		glfwSetInputMode(this->window, GLFW_STICKY_KEYS, GL_FALSE);
}

void Window::stickyMouseButtons(bool sticky)
{
	if (sticky)
		glfwSetInputMode(this->window, GLFW_STICKY_MOUSE_BUTTONS, GL_TRUE);
	else
		glfwSetInputMode(this->window, GLFW_STICKY_MOUSE_BUTTONS, GL_FALSE);
}

void Window::setCallbacks(const GLFWkeyfun *keyCallback, const GLFWcharfun *charCallback, const GLFWmousebuttonfun *mouseButtonCallback, GLFWcursorposfun *cursorPosCallback, GLFWcursorenterfun *cursorEnterCallback, GLFWscrollfun *scrollCallback, GLFWdropfun *dropCallback)
{
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

bool Window::isOpen()
{
	return !glfwWindowShouldClose(this->window);
}

int Window::destroy() {
    glfwDestroyWindow(this->window);
    glfwTerminate();
    return 0;
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


	GLuint VertexShaderID =  vertexShader.getID();
	GLuint FragmentShaderID =  fragmentShader.getID();
	GLuint GeometryShaderID = geometryShader ? geometryShader->getID() : 0;

	programID = glCreateProgram();

	glAttachShader(programID, VertexShaderID);
	glAttachShader(programID, FragmentShaderID);
	if (geometryShader)
		glAttachShader(programID, GeometryShaderID);

	glLinkProgram(programID);
	glGetProgramiv(programID, GL_LINK_STATUS, &Result);
	glGetProgramiv(programID, GL_INFO_LOG_LENGTH, &InfoLogLength);

	if (InfoLogLength > 0)
	{
		vector<char> ProgramErrorMessage(InfoLogLength + 1);
		glGetProgramInfoLog(programID, InfoLogLength, NULL, &ProgramErrorMessage[0]);
		throw SystemError(&ProgramErrorMessage[0]);
	}

	glDetachShader(programID, VertexShaderID);
	if (geometryShader) glDetachShader(programID, GeometryShaderID);
	glDetachShader(programID, FragmentShaderID);

	glDeleteShader(VertexShaderID);
	if (geometryShader) glDeleteShader(GeometryShaderID);
	glDeleteShader(FragmentShaderID);

	shaderType = geometryShader ? GEOMETRY1 : CLASSIC;
}


ShaderProgram::ShaderProgram(const Shader &vertexShader, const Shader &fragmentShader, const Shader &geometryShader) :
	vertexShader(vertexShader),
	fragmentShader(fragmentShader),
	geometryShader(geometryShader){

	linkShaders();
}

ShaderProgram::ShaderProgram(const string &vertexPath, const string &fragPath) :
vertexShader(CodeFileDescriptor(vertexPath, false)),
fragmentShader(CodeFileDescriptor(fragPath, false)),
geometryShader(std::nullopt){
	linkShaders();
}

ShaderProgram::ShaderProgram(const Shader &vertexShader, const Shader &fragmentShader) :
	vertexShader(vertexShader),
	fragmentShader(fragmentShader),
	geometryShader(std::nullopt){

	linkShaders();
}


ShaderProgram::~ShaderProgram()
{
	glDeleteProgram(this->programID);
}

void ShaderProgram::use()
{
	glUseProgram(this->programID);
}

void ShaderProgram::initUniforms(const std::map<string, GLSLType> &uniforms)
{
	for (auto uni : uniforms)
		this->uniformTypes[uni.first] = uni.second;
	for (const auto &key: uniforms | std::views::keys)
	{
		GLuint uniformLocation = glGetUniformLocation(this->programID, key.c_str());
		this->uniformLocations[key] = uniformLocation;
	}
}

void ShaderProgram::setTextureSampler(const Texture *texture, int textureSlot) const
{
	GLuint samplerUniform = glGetUniformLocation(this->programID, texture->samplerName);
	glUniform1i(samplerUniform, textureSlot);
}

void ShaderProgram::setUniforms(const std::map<string, const GLfloat *> &uniformValues)
{
	for (auto const &uniform : uniformValues)
		setUniform(uniform.first, uniform.second);
	
}

void ShaderProgram::setUniform(const string &uniformName, const GLfloat *uniformValue)
{
	GLSLType uniformType = this->uniformTypes[uniformName];
	GLuint uniformLocation = this->uniformLocations[uniformName];

	switch (uniformType)
	{
	case FLOAT:
		glUniform1fv(uniformLocation, 1, uniformValue);
		break;
	case INT:
		glUniform1iv(uniformLocation, 1, (GLint *)uniformValue);
		break;
	case VEC2:
		glUniform2fv(uniformLocation, 1, uniformValue);
		break;
	case VEC3:
		glUniform3fv(uniformLocation, 1, uniformValue);
		break;
	case VEC4:
		glUniform4fv(uniformLocation, 1, uniformValue);
		break;
	case MAT2:
		glUniformMatrix2fv(uniformLocation, 1, GL_FALSE, uniformValue);
		break;
	case MAT3:
		glUniformMatrix3fv(uniformLocation, 1, GL_FALSE, uniformValue);
		break;
	case MAT4:
		glUniformMatrix4fv(uniformLocation, 1, GL_FALSE, uniformValue);
		break;

	default:
		throw std::invalid_argument("Uniform type not recognized");
	}
}

void ShaderProgram::setUniform(const string &uniformName, float uniformValue)
{
	if (this->uniformTypes[uniformName] != FLOAT)
		throw std::invalid_argument("Uniform type must be FLOAT");
	GLuint uniformLocation = this->uniformLocations[uniformName];
	glUniform1f(uniformLocation, uniformValue);
}

void ShaderProgram::setUniform(const string &uniformName, int uniformValue)
{
	if (this->uniformTypes[uniformName] != INT)
		throw std::invalid_argument("Uniform type must be INT");
	GLuint uniformLocation = this->uniformLocations[uniformName];
	glUniform1i(uniformLocation, uniformValue);
}

void ShaderProgram::setUniform(const string &uniformName, vec2 uniformValue)
{
	if (this->uniformTypes[uniformName] != VEC2)
		throw std::invalid_argument("Uniform type must be VEC2");
	GLuint uniformLocation = this->uniformLocations[uniformName];
	glUniform2f(uniformLocation, uniformValue.x, uniformValue.y);
}

void ShaderProgram::setUniform(const string &uniformName, vec3 uniformValue)
{
	if (this->uniformTypes[uniformName] != VEC3)
		throw std::invalid_argument("Uniform type must be VEC3");
	GLuint uniformLocation = this->uniformLocations[uniformName];
	glUniform3f(uniformLocation, uniformValue.x, uniformValue.y, uniformValue.z);
}

void ShaderProgram::setUniform(const string &uniformName, vec4 uniformValue)
{
	if (this->uniformTypes[uniformName] != VEC4)
		throw std::invalid_argument("Uniform type must be VEC4");
	GLuint uniformLocation = this->uniformLocations[uniformName];
	glUniform4f(uniformLocation, uniformValue.x, uniformValue.y, uniformValue.z, uniformValue.w);
}

void ShaderProgram::setUniform(const string &uniformName, mat2 uniformValue)
{
	if (this->uniformTypes[uniformName] != MAT2)
		throw std::invalid_argument("Uniform type must be MAT2");
	GLuint uniformLocation = this->uniformLocations[uniformName];
	glUniformMatrix2fv(uniformLocation, 1, GL_FALSE, &uniformValue[0][0]);
}

void ShaderProgram::setUniform(const string &uniformName, mat3 uniformValue)
{
	if (this->uniformTypes[uniformName] != MAT3)
		throw std::invalid_argument("Uniform type must be MAT3");
	GLuint uniformLocation = this->uniformLocations[uniformName];
	glUniformMatrix3fv(uniformLocation, 1, GL_FALSE, &uniformValue[0][0]);
}

void ShaderProgram::setUniform(const string& uniformName, mat4 uniformValue)
{
	if (this->uniformTypes[uniformName] != MAT4)
		throw std::invalid_argument("Uniform type must be MAT4");
	GLuint uniformLocation = this->uniformLocations[uniformName];
	glUniformMatrix4fv(uniformLocation, 1, GL_FALSE, &uniformValue[0][0]);
}

void ShaderProgram::setUniform(const string& uniformName, float x, float y)
{
	if (this->uniformTypes[uniformName] != VEC2)
		throw std::invalid_argument("Uniform type must be VEC2");
	GLuint uniformLocation = this->uniformLocations[uniformName];
	glUniform2f(uniformLocation, x, y);
}

void ShaderProgram::setUniform(const string &uniformName, float x, float y, float z)
{
	if (this->uniformTypes[uniformName] != VEC3)
		throw std::invalid_argument("Uniform type must be VEC3");
	GLuint uniformLocation = this->uniformLocations[uniformName];
	glUniform3f(uniformLocation, x, y, z);
}

void ShaderProgram::setUniform(const string &uniformName, float x, float y, float z, float w)
{
	if (this->uniformTypes[uniformName] != VEC4)
		throw std::invalid_argument("Uniform type must be VEC4");
	GLuint uniformLocation = this->uniformLocations[uniformName];
	glUniform4f(uniformLocation, x, y, z, w);
}


Camera::Camera()
{
	this->lookAtFunc =  std::make_shared<SmoothParametricCurve>(SmoothParametricCurve::constCurve(vec3(0)));
	this->up = [](float t) { return vec3(0, 0, 1); };
	this->fov_x = 45.0f;
	this->aspectRatio = 16.0f / 9.0f;
	this->clippingRangeMin = 0.1f;
	this->clippingRangeMax = 100.0f;
	this->trajectory = std::make_shared<SmoothParametricCurve>(SmoothParametricCurve::constCurve(vec3(2.0f, 3.0f, 1.0f)));
	this->moving = false;
	this->projectionMatrix = perspective(fov_x, aspectRatio, clippingRangeMin, clippingRangeMax);
}

Camera::Camera(vec3 position, vec3 lookAtPos, vec3 upVector, float fov_x, float aspectRatio, float clippingRangeMin, float clippingRangeMax)
{
	this->lookAtFunc = std::make_shared<SmoothParametricCurve>(SmoothParametricCurve::constCurve(lookAtPos));
	this->up = [upVector](float t) { return upVector; };
	this->fov_x = fov_x;
	this->aspectRatio = aspectRatio;
	this->clippingRangeMin = clippingRangeMin;
	this->clippingRangeMax = clippingRangeMax;
	this->trajectory = std::make_shared<SmoothParametricCurve>(SmoothParametricCurve::constCurve(position));
	this->moving = false;
	this->projectionMatrix = perspective(fov_x, aspectRatio, clippingRangeMin, clippingRangeMax);
}

Camera::Camera(float radius, float speed, float height, vec3 lookAtPos, vec3 upVector, float fov_x, float aspectRatio, float clippingRangeMin, float clippingRangeMax) :
	Camera(std::make_shared<SmoothParametricCurve>(
	[speed, radius, height](float t) {
		return vec3(radius * cos(speed * t), radius * sin(speed * t), height); }
	), lookAtPos, upVector, fov_x, aspectRatio, clippingRangeMin, clippingRangeMax) {}

Camera::Camera(const shared_ptr<SmoothParametricCurve> &trajectory, vec3 lookAtPos, vec3 upVector, float fov_x, float aspectRatio, float clippingRangeMin, float clippingRangeMax)
{
    this->lookAtFunc = std::make_shared<SmoothParametricCurve>(SmoothParametricCurve::constCurve(lookAtPos));
    this->up = [upVector](float t) { return upVector; };
	this->fov_x = fov_x;
	this->aspectRatio = aspectRatio;
	this->clippingRangeMin = clippingRangeMin;
	this->clippingRangeMax = clippingRangeMax;
	this->trajectory = trajectory;
	this->moving = true;
	this->projectionMatrix = perspective(fov_x, aspectRatio, clippingRangeMin, clippingRangeMax);
}

Camera::Camera(const std::shared_ptr<SmoothParametricCurve> &trajectory, const std::shared_ptr<SmoothParametricCurve> &lookAtPos,
        vec3 upVector, float fov_x, float aspectRatio, float clippingRangeMin, float clippingRangeMax) {

    this->lookAtFunc = lookAtPos;
    this->up = [upVector](float t) { return upVector; };
    this->fov_x = fov_x;
    this->aspectRatio = aspectRatio;
    this->clippingRangeMin = clippingRangeMin;
    this->clippingRangeMax = clippingRangeMax;
    this->trajectory = trajectory;
    this->moving = true;
    this->projectionMatrix = perspective(fov_x, aspectRatio, clippingRangeMin, clippingRangeMax);
}

Camera::Camera(const std::shared_ptr<SmoothParametricCurve> &trajectory, const std::shared_ptr<SmoothParametricCurve> &lookAtPos,
        const std::function<vec3(float)> &upVector, float fov_x, float aspectRatio, float clippingRangeMin, float clippingRangeMax) {

    this->lookAtFunc = lookAtPos;
    this->up = upVector;
    this->fov_x = fov_x;
    this->aspectRatio = aspectRatio;
    this->clippingRangeMin = clippingRangeMin;
    this->clippingRangeMax = clippingRangeMax;
    this->trajectory = trajectory;
    this->moving = true;
    this->projectionMatrix = perspective(fov_x, aspectRatio, clippingRangeMin, clippingRangeMax);
}

vec3 Camera::position(float t) const { return trajectory->operator()(t); }

vec3 Camera::lookAtPoint(float t) const { return lookAtFunc->operator()(t); }

vec3 Camera::upVector(float t) const { return up(t); }


mat4 Camera::viewMatrix(float t)
{
    return lookAt(position(t), lookAtPoint(t), upVector(t));
}

mat4 Camera::vp(float t)
{
	return projectionMatrix * viewMatrix(t);
}

mat4 Camera::mvp(float t, const mat4 &modelTransform)
{
	// return vp(t);
    return vp(t) * modelTransform;
}


Attribute::Attribute(const string &name, GLSLType type, int inputNumber)
{
	this->name = name;
	this->type = type;
	this->bufferAddress = 0;
	this->size = sizeOfGLSLType(type);
	this->enabled = false;
	this->bufferInitialized = false;
	this->inputNumber = inputNumber;
}

Attribute::~Attribute()
{
	if (enabled) Attribute::disable();
	if (bufferInitialized) Attribute::freeBuffer();
}

void Attribute::initBuffer()
{
	this->bufferInitialized = true;
	GLuint buffer;
	glGenBuffers(1, &buffer);
	this->bufferAddress = buffer;
    // glBufferData(GL_ARRAY_BUFFER, bufferLength * this->size, firstElementAdress, GL_STATIC_DRAW);
}

void Attribute::enable()
{
	glEnableVertexAttribArray(this->inputNumber);
	glBindBuffer(GL_ARRAY_BUFFER, this->bufferAddress);
	this->enabled = true;
	glVertexAttribPointer(this->inputNumber, lengthOfGLSLType(this->type), GL_FLOAT, GL_FALSE, 0, (void *)0);
}

void Attribute::disable()
{
	glDisableVertexAttribArray(this->inputNumber);
	this->enabled = false;
}

void Attribute::load(const void *firstElementAdress, int bufferLength)
{
	if (!bufferInitialized) {
	    initBuffer();
	    glBindBuffer(GL_ARRAY_BUFFER, this->bufferAddress);
	    glBufferData(GL_ARRAY_BUFFER, bufferLength * this->size, firstElementAdress, GL_STATIC_DRAW);
	    return;
	}
	glBindBuffer(GL_ARRAY_BUFFER, this->bufferAddress);
     glBufferData(GL_ARRAY_BUFFER, bufferLength * this->size, firstElementAdress, GL_STATIC_DRAW);
    // void *ptr = glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);
    // memcpy(ptr, firstElementAdress, bufferLength * this->size);
    // glUnmapBuffer(GL_ARRAY_BUFFER);

}

void Attribute::freeBuffer() {
    glDeleteBuffers(1, &this->bufferAddress);
    this->bufferAddress = -1;
    this->bufferInitialized = false;
    this->enabled = false;
}

RenderingStep::RenderingStep(const shared_ptr<ShaderProgram> &shader)
{
	this->shader = shader;
	this->attributes = vector<shared_ptr<Attribute>>();
	this->uniforms = std::map<string, GLSLType>();
	this->uniformSetters = std::map<string, shared_ptr<std::function<void(float, shared_ptr<ShaderProgram>)>>>();
	this->customStep = [](float t) {};
}

RenderingStep::RenderingStep(const RenderingStep& other)
{
	this->shader = other.shader;
	this->attributes = other.attributes;
	this->uniforms = other.uniforms;
	this->uniformSetters = other.uniformSetters;
	this->customStep = other.customStep;
}

RenderingStep::~RenderingStep()
{
	for (auto attribute : attributes)
		attribute.reset();
	shader.reset();
}

void RenderingStep::setWeakSuperMesh(const std::shared_ptr<IndexedMesh> &super) {
    this-> weak_super = super;
}

int RenderingStep::findAttributeByName(const string &name) {
	for (int i = 0; i < attributes.size(); i++)
		if (attributes[i]->name == name)
			return i;
	throw std::invalid_argument("Attribute not found");
}

void RenderingStep::initMaterialAttributes() {
    auto ambient = std::make_shared<Attribute>("ambientColor", VEC4, 4);
    auto diffuse = std::make_shared<Attribute>("diffuseColor", VEC4, 5);
    auto specular = std::make_shared<Attribute>("SpecularColor", VEC4, 6);
    auto intencities = std::make_shared<Attribute>("intencities", VEC4, 7);
    this->attributes.push_back(ambient);
    this->attributes.push_back(diffuse);
    this->attributes.push_back(specular);
    this->attributes.push_back(intencities);

    for (const auto &attribute: attributes) {
        attribute->initBuffer();
    }
}
void RenderingStep::initElementBuffer() {
    glGenBuffers(1, &elementBufferLoc);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, elementBufferLoc);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, weak_super->bufferIndexSize(), weak_super->bufferIndexLocation(), GL_STREAM_DRAW);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, elementBufferLoc);

}

void RenderingStep::loadElementBuffer() {
    if (!weakSuperLoaded())
        throw std::invalid_argument("Element buffer can only be loaded for weak super mesh");
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, elementBufferLoc);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, weak_super->bufferIndexSize(), weak_super->bufferIndexLocation(), GL_STATIC_DRAW);
}

void RenderingStep::initStdAttributes()
{
	auto positionAttribute = std::make_shared<Attribute>("position", VEC3, 0);
	auto normalAttribute = std::make_shared<Attribute>("normal", VEC3, 1);
	auto colorAttribute = std::make_shared<Attribute>("color", VEC4, 3);
	auto uvAttribute = std::make_shared<Attribute>("uv", VEC2, 2);
	this->attributes.push_back(positionAttribute);
	this->attributes.push_back(normalAttribute);
	this->attributes.push_back(uvAttribute);
	this->attributes.push_back(colorAttribute);
	if (weak_super->hasExtra0())
		this->attributes.push_back(std::make_shared<Attribute>("extra0", VEC4, 4));
	if (weak_super->hasExtra1())
		this->attributes.push_back(std::make_shared<Attribute>("extra1", VEC4, 5));
	if (weak_super->hasExtra2())
		this->attributes.push_back(std::make_shared<Attribute>("extra2", VEC4, 6));
	if (weak_super->hasExtra3())
		this->attributes.push_back(std::make_shared<Attribute>("extra3", VEC4, 7));
	if (weak_super->hasExtra4())
		this->attributes.push_back(std::make_shared<Attribute>("extra4", VEC4, 8));

	for (const auto& attribute : attributes)
		attribute->initBuffer();


}

void RenderingStep::resetAttributeBuffers()
{
	for (const auto& attribute : attributes)
		attribute->freeBuffer();
}

void RenderingStep::initUnusualAttributes(const std::vector<std::shared_ptr<Attribute>>& attributes) {
    this->attributes.insert(this->attributes.end(), attributes.begin(), attributes.end());
    for (const auto &attribute: attributes)
		attribute->initBuffer();
}


void RenderingStep::initExtraAttribute(int i) {
	attributes.push_back(std::make_shared<Attribute>("extra" + std::to_string(i), VEC4, 4 + i));
	attributes.back()->initBuffer();
}




void RenderingStep::loadMeshAttributes() {
    if (weakSuperLoaded())
    {
        for (auto i = 0; i < 4; i++)
            attributes[i]->load(weak_super->getBufferLocation(static_cast<CommonBufferType>(i)), weak_super->getBufferLength(static_cast<CommonBufferType>(i)));
    	if (weak_super->hasExtra0())
    		attributes[4]->load(weak_super->getBufferLocation(EXTRA0), weak_super->getBufferLength(EXTRA0));
    	if (weak_super->hasExtra1())
			attributes[5]->load(weak_super->getBufferLocation(EXTRA1), weak_super->getBufferLength(EXTRA1));
    	if (weak_super->hasExtra2())
    		attributes[6]->load(weak_super->getBufferLocation(EXTRA2), weak_super->getBufferLength(EXTRA2));
    	if (weak_super->hasExtra3())
			attributes[7]->load(weak_super->getBufferLocation(EXTRA3), weak_super->getBufferLength(EXTRA3));
    	if (weak_super->hasExtra4())
    		attributes[8]->load(weak_super->getBufferLocation(EXTRA4), weak_super->getBufferLength(EXTRA4));
    	return;
    }



	// for (int i = 0; i < 4; i++)
	// 	attributes[i]->load(model->mesh->bufferLocations[i], model->mesh->bufferSizes[i]);
}

void RenderingStep::initWeakMeshAttributes() {
	initStdAttributes();
	// for (int i = 0; i < 5; i++)
	// 	if (weak_super->hasExtra(i))
	// 		initExtraAttribute(i);
}

void RenderingStep::enableAttributes()
{
	for (const auto& attribute : attributes)
		attribute->enable();
}

void RenderingStep::disableAttributes()
{
	for (const auto& attribute : attributes)
		attribute->disable();
}

void RenderingStep::addUniform(string uniformName, GLSLType uniformType, shared_ptr<std::function<void(float, shared_ptr<ShaderProgram>)>> setter)
{
	this->uniforms[uniformName] = uniformType;
	this->uniformSetters[uniformName] = std::move(setter);
	this->shader->initUniforms({{uniformName, uniformType}});
}

void RenderingStep::addConstFloats(const std::map<string, float>& uniforms)
{
	for (auto uniform : uniforms)
		addUniform(uniform.first, FLOAT, std::make_shared<std::function<void(float, shared_ptr<ShaderProgram>)>>
				   ([uniform](float t, const shared_ptr<ShaderProgram>& shader)
				   {
					   shader->setUniform(uniform.first, uniform.second);
				   }));
}

void RenderingStep::addConstVec4(const string& uniformName, vec4 value)
{
	auto uniformSetter = [value, uniformName](float t, const shared_ptr<ShaderProgram> &shader) {
		shader->setUniform(uniformName, value);
	};
	addUniform(uniformName, VEC4,  std::make_shared<std::function<void(float, shared_ptr<ShaderProgram>)>>(uniformSetter));
}

void RenderingStep::addConstColor(const string &name, vec4 value) { addConstVec4(name, value); }

void RenderingStep::setUniforms(float t)
{
	for (const auto &key: uniforms | std::views::keys)
	{
		auto name = key;
		(*uniformSetters[name])(t, shader);
	}
}

void RenderingStep::addCameraUniforms(const std::shared_ptr<Camera>& camera)
{
	std::function<void(float, shared_ptr<ShaderProgram>)> MVPsetter;

	if (weakSuperLoaded())
		MVPsetter = [camera, this](float t, const shared_ptr<ShaderProgram> &shader) {
			mat4 mvp = camera->mvp(t, mat4(mat3(1)));
			shader->setUniform("mvp", mvp);
		};

	auto positionSetter = [camera, this](float t, const shared_ptr<ShaderProgram> &shader) {
		vec3 camPos = camera->position(t);
		shader->setUniform("camPosition", camPos);
	};
	addUniform("mvp", MAT4,  std::make_shared<std::function<void(float, shared_ptr<ShaderProgram>)>>(MVPsetter));
	addUniform("camPosition", VEC3,  std::make_shared<std::function<void(float, shared_ptr<ShaderProgram>)>>(positionSetter));
}

void RenderingStep::addLightUniform(const Light& pointLight, int lightIndex)
{
	string lightName = "light" + std::to_string(lightIndex);
	auto lightSetter = [pointLight, lightName](float t, const shared_ptr<ShaderProgram> &shader) {
		mat4 lightMat = pointLight.compressToMatrix();
		shader->setUniform(lightName, lightMat);
	};
	addUniform(lightName, MAT4, std::make_shared<std::function<void(float, shared_ptr<ShaderProgram>)>>(lightSetter));
}

void RenderingStep::addLightsUniforms(const std::vector<Light> &lights)
{
	for (int i = 1; i <= lights.size(); i++)
		addLightUniform(lights[i-1], i);
}




void RenderingStep::init(const shared_ptr<Camera> &cam, const std::vector<Light> &lights) {
    if (weakSuperLoaded()) {
        shader->use();
        initElementBuffer();
        initWeakMeshAttributes();
		addCameraUniforms(cam);
    	addLightsUniforms(lights);

        loadMeshAttributes();
        loadElementBuffer();
    }
}

void RenderingStep::addUniforms(const std::map<string, GLSLType> &uniforms, std::map<string, shared_ptr<std::function<void(float, shared_ptr<ShaderProgram>)>>> setters)
{
	for (const auto& uniform : uniforms)
		addUniform(uniform.first, uniform.second, setters[uniform.first]);
}

void RenderingStep::addCustomAction(const std::function<void(float)> &action)
{
	this->customStep = action;
}

bool RenderingStep::weakSuperLoaded() const {
	return weak_super != nullptr;
}


void RenderingStep::weakMeshRenderStep(float t) {
    // bindTextures();
    loadMeshAttributes();
    customStep(t);
    setUniforms(t);
    loadElementBuffer();

    enableAttributes();
    glDrawElements(GL_TRIANGLES, weak_super->bufferIndexSize(), GL_UNSIGNED_INT, 0);
    disableAttributes();
}


void RenderingStep::renderStep(float t)
{
    shader->use();
    if (weakSuperLoaded())
        weakMeshRenderStep(t);
    else
        throw ValueError("Mesh not loaded correctly");
}

RenderSettings::RenderSettings(vec4 bgColor, bool alphaBlending, bool depthTest, bool timeUniform, float speed, int maxFPS, bool takeScreenshots, Resolution resolution, float screenshotFrequency, const string &windowTitle):
bgColor(bgColor),
alphaBlending(alphaBlending),
depthTest(depthTest),
timeUniform(timeUniform),
speed(speed),
maxFPS(maxFPS),
takeScreenshots(takeScreenshots),
resolution(resolution),
screenshotFrequency(screenshotFrequency),
windowTitle(windowTitle) {

	if (takeScreenshots)
		throw NotImplementedError("Screenshot functionality not implemented yet."); // TODO screenshots
	this->screenshotFrequency = -1;

	screenshotDirectory = ConfigFile().getScreenshotsDir();
}

Renderer::Renderer(float animSpeed, vec4 bgColor, const string &screenshotDirectory, float screenshotFrequency)
: settings(bgColor,
		   true,
		   true,
		   true,
		   animSpeed,
		   120,
		   false,
		   UNKNOWN,
		   screenshotFrequency) {

	this->window = nullptr;
	this->vao = 0;
	this->camera = nullptr;
	this->time = 0;
	this->lights = std::vector<Light>();
	Logger::init();

	if (!glfwInit()) {
		throw SystemError("Failed to initialize GLFW");
	}
	this->animSpeed = [animSpeed](float t) { return animSpeed; };
	this->perFrameFunction = [](float t, float delta) {};

}

Renderer::Renderer(const RenderSettings &settings)
: Renderer(settings.speed, settings.bgColor, settings.screenshotDirectory.to_str(), settings.screenshotFrequency)
{
	this->settings = settings;
	this->animSpeed = [settings](float t) { return settings.speed; };

}


Renderer::~Renderer()
{
	if (window != nullptr)
		window->destroy();
	if (camera != nullptr)
		camera.reset();
}

void Renderer::initMainWindow(int width, int height, const char *title)
{
	this->window = std::make_unique<Window>(width, height, title);
	glewExperimental = true;
	if (glewInit() != GLEW_OK) {
		LOG_ERROR("Failed to initialize GLEW");
		throw SystemError("Failed to initialize GLEW");
	}
	this->vao = bindVAO();
	LOG("Window initialized with size: " + std::to_string(width) + "x" + std::to_string(height));
}

void Renderer::initMainWindow(Resolution resolution, const char *title) {
	initMainWindow(predefinedWidth(resolution), predefinedHeight(resolution), title);
}

void Renderer::initMainWindow() {
	const char *title = settings.windowTitle.c_str();
	initMainWindow(predefinedWidth(settings.resolution), predefinedHeight(settings.resolution), title);
}

void Renderer::resetTimer() {
	last_time_capture = glfwGetTime();
	this->time = 0;
}

void Renderer::addRenderingStep(std::shared_ptr<RenderingStep> renderingStep)
{
	this->renderingSteps.push_back(std::move(renderingStep));
}

void Renderer::addMeshStep(const ShaderProgram &shader, const std::shared_ptr<IndexedMesh> &model, const MaterialPhong &material) {
	auto renderingStep = make_shared<RenderingStep>(make_shared<ShaderProgram>(shader));
	renderingStep->setWeakSuperMesh(model);
	addRenderingStep(renderingStep);
}



void Renderer::setCamera(const std::shared_ptr<Camera> &camera)
{
	this->camera = camera;
}

void Renderer::setLights(const std::vector<Light> &lights)
{
	this->lights = lights;
}

void Renderer::setLightWithMesh(const Light &light, const MaterialPhong &material, const ShaderProgram &shader, float radius) {
	this->lights.push_back(light);
	addMeshStep(shader, std::make_shared<IndexedMesh>(icosphere(radius, 2, light.getPosition(), randomID())), material);
}

void Renderer::setLightsWithMesh(const std::vector<Light> &lights, const MaterialPhong &material, const ShaderProgram &shader, float radius) {
	for (const auto &light : lights)
		setLightWithMesh(light, material, shader, radius);
}

void Renderer::setLightWithMesh(const Light &light, float ambient, float diff, float spec, float shine, const ShaderProgram &shader, float radius) {
	setLightWithMesh(light, MaterialPhong(light.getColor(), light.getColor(), WHITE, ambient, diff, spec, shine), shader, radius);
}
void Renderer::setLightsWithMesh(const std::vector<Light> &lights, float ambient, float diff, float spec, float shine, const ShaderProgram &shader, float radius) {
	for (const auto &light : lights)
		setLightWithMesh(light, ambient, diff, spec, shine, shader, radius);
}

float Renderer::initFrame()
{
//	glViewport(0, 0, window->width, window->height);
	glClearColor(settings.bgColor.x, settings.bgColor.y, settings.bgColor.z, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

	float real_dt = 0;
	while (real_dt * settings.maxFPS < 1.0f) {
		auto t = glfwGetTime();
		real_dt += t - last_time_capture;
		last_time_capture = t;
	}
    dt = real_dt*animSpeed(time);
	time += dt;
	return time;
}



float Renderer::lastDeltaTime() const {
    return dt;
}

void Renderer::addPerFrameUniforms(const std::map<string, GLSLType> &uniforms, const std::map<string, shared_ptr<std::function<void(float, std::shared_ptr<ShaderProgram>)>>> &setters)
{
	for (const auto& renderingStep : renderingSteps)
		renderingStep->addUniforms(uniforms, setters);

}

void Renderer::addPerFrameUniform(const string &uniformName, GLSLType uniformType, const shared_ptr<std::function<void(float, std::shared_ptr<ShaderProgram>)>> &setter)
{
	for (const auto& renderingStep : renderingSteps)
		renderingStep->addUniform(uniformName, uniformType, setter);
}

void Renderer::addSurfaceFamilyDeformer(SurfaceParametricPencil &pencil, IndexedMesh &surface) {
	addCustomAction([&pencil, &surface](float t){
		surface.adjustToNewSurface(pencil(t));
	});
}

void Renderer::initRendering()
{

    window->initViewport();
    glBindVertexArray(vao);

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
	{
	    if (renderingStep->weakSuperLoaded())
	        renderingStep->init(camera, lights);

		else {
			throw ValueError("Rendering step does not have a weak super mesh set");
		}
	}

	LOG("OpenGL renderer initialized: " + string((const char *)glGetString(GL_RENDERER)));

}

void Renderer::addConstUniform(const string &uniformName, GLSLType uniformType, shared_ptr<std::function<void(std::shared_ptr<ShaderProgram>)>> setter)
{
	std::function<void(float, std::shared_ptr<ShaderProgram>)> setterWrapper = [setter](float t, const std::shared_ptr<ShaderProgram> &shader) { (*setter)(shader); };
	for (const auto& renderingStep : renderingSteps)
		renderingStep->addUniform(uniformName, uniformType, make_shared<std::function<void(float, std::shared_ptr<ShaderProgram>)>>(setterWrapper));
}

void Renderer::addTimeUniform()
{
	addPerFrameUniform("time", FLOAT, std::make_shared<std::function<void(float, std::shared_ptr<ShaderProgram>)>>(
	[](float t, const std::shared_ptr<ShaderProgram> &shader) {
		shader->setUniform("time", t);
	}));
}

void Renderer::addConstFloats(const std::map<string, float> &uniforms)
{
	for (const auto& renderingStep : renderingSteps)
		renderingStep->addConstFloats(uniforms);
}

void Renderer::addCustomAction(std::function<void(float)> action)
{
    perFrameFunction = [a=perFrameFunction, n=action](float t, float delta) {
		a(t, delta);
		n(t);
	};
}

void Renderer::addCustomAction(std::function<void(float, float)> action)
{
    perFrameFunction = [a=perFrameFunction, n=action](float t, float delta) {
        a(t, delta);
        n(t, delta);
    };
}

void Renderer::nonlinearSpeed(const std::function<float(float)> &speed) {
	animSpeed = speed;
}

void Renderer::addConstUniforms(const std::map<string, GLSLType>& uniforms, std::map<string, shared_ptr<std::function<void(std::shared_ptr<ShaderProgram>)>>> setters)
{	
	for (const auto& uniform : uniforms)
		addConstUniform(uniform.first, uniform.second,  setters[uniform.first]);
}

void Renderer::renderAllSteps()
{
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
