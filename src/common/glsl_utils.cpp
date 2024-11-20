
#include "glsl_utils.hpp"
#include <stdio.h>
#include <string>
#include <map>
#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <sstream>

#include <stdlib.h>
#include <string.h>

#include <GL/glew.h>

using namespace glm;
using std::vector, std::string, std::map, std::shared_ptr, std::unique_ptr, std::make_shared, std::make_unique;

GLuint LoadShaders(const char *vertex_file_path, const char *fragment_file_path)
{

	// Create the shaders
	GLuint VertexShaderID = glCreateShader(GL_VERTEX_SHADER);
	GLuint FragmentShaderID = glCreateShader(GL_FRAGMENT_SHADER);

	// Read the Vertex Shader code from the file
	string VertexShaderCode;
	std::ifstream VertexShaderStream(vertex_file_path, std::ios::in);
	if (VertexShaderStream.is_open())
	{
		std::stringstream sstr;
		sstr << VertexShaderStream.rdbuf();
		VertexShaderCode = sstr.str();
		VertexShaderStream.close();
	}
	else
	{
		printf("Impossible to open %s. Are you in the right directory ? Don't forget to read the FAQ !\n", vertex_file_path);
		getchar();
		return 0;
	}

	// Read the Fragment Shader code from the file
	string FragmentShaderCode;
	std::ifstream FragmentShaderStream(fragment_file_path, std::ios::in);
	if (FragmentShaderStream.is_open())
	{
		std::stringstream sstr;
		sstr << FragmentShaderStream.rdbuf();
		FragmentShaderCode = sstr.str();
		FragmentShaderStream.close();
	}

	GLint Result = GL_FALSE;
	int InfoLogLength;

	// Compile Vertex Shader
	printf("Compiling shader : %s\n", vertex_file_path);
	char const *VertexSourcePointer = VertexShaderCode.c_str();
	glShaderSource(VertexShaderID, 1, &VertexSourcePointer, NULL);
	glCompileShader(VertexShaderID);

	// Check Vertex Shader
	glGetShaderiv(VertexShaderID, GL_COMPILE_STATUS, &Result);
	glGetShaderiv(VertexShaderID, GL_INFO_LOG_LENGTH, &InfoLogLength);
	if (InfoLogLength > 0)
	{
		vector<char> VertexShaderErrorMessage(InfoLogLength + 1);
		glGetShaderInfoLog(VertexShaderID, InfoLogLength, NULL, &VertexShaderErrorMessage[0]);
		printf("%s\n", &VertexShaderErrorMessage[0]);
	}

	// Compile Fragment Shader
	printf("Compiling shader : %s\n", fragment_file_path);
	char const *FragmentSourcePointer = FragmentShaderCode.c_str();
	glShaderSource(FragmentShaderID, 1, &FragmentSourcePointer, NULL);
	glCompileShader(FragmentShaderID);

	// Check Fragment Shader
	glGetShaderiv(FragmentShaderID, GL_COMPILE_STATUS, &Result);
	glGetShaderiv(FragmentShaderID, GL_INFO_LOG_LENGTH, &InfoLogLength);
	if (InfoLogLength > 0)
	{
		vector<char> FragmentShaderErrorMessage(InfoLogLength + 1);
		glGetShaderInfoLog(FragmentShaderID, InfoLogLength, NULL, &FragmentShaderErrorMessage[0]);
		printf("%s\n", &FragmentShaderErrorMessage[0]);
	}

	// Link the program
	printf("Linking program\n");
	GLuint ProgramID = glCreateProgram();
	glAttachShader(ProgramID, VertexShaderID);
	glAttachShader(ProgramID, FragmentShaderID);
	glLinkProgram(ProgramID);

	// Check the program
	glGetProgramiv(ProgramID, GL_LINK_STATUS, &Result);
	glGetProgramiv(ProgramID, GL_INFO_LOG_LENGTH, &InfoLogLength);
	if (InfoLogLength > 0)
	{
		vector<char> ProgramErrorMessage(InfoLogLength + 1);
		glGetProgramInfoLog(ProgramID, InfoLogLength, NULL, &ProgramErrorMessage[0]);
		printf("%s\n", &ProgramErrorMessage[0]);
	}

	glDetachShader(ProgramID, VertexShaderID);
	glDetachShader(ProgramID, FragmentShaderID);

	glDeleteShader(VertexShaderID);
	glDeleteShader(FragmentShaderID);

	return ProgramID;
}

GLuint LoadShaders(const char *vertex_file_path, const char *geometry_file_path,const char *fragment_file_path)
{

	// Create the shaders
	GLuint VertexShaderID = glCreateShader(GL_VERTEX_SHADER);
	GLuint FragmentShaderID = glCreateShader(GL_FRAGMENT_SHADER);
	GLuint GeometryShaderID = glCreateShader(GL_GEOMETRY_SHADER);

	// Read the Vertex Shader code from the file
	string VertexShaderCode;
	std::ifstream VertexShaderStream(vertex_file_path, std::ios::in);
	if (VertexShaderStream.is_open())
	{
		std::stringstream sstr;
		sstr << VertexShaderStream.rdbuf();
		VertexShaderCode = sstr.str();
		VertexShaderStream.close();
	}
	else
	{
		printf("Impossible to open %s. Are you in the right directory ? Don't forget to read the FAQ !\n", vertex_file_path);
		getchar();
		return 0;
	}

	// Read the Fragment Shader code from the file
	string FragmentShaderCode;
	std::ifstream FragmentShaderStream(fragment_file_path, std::ios::in);
	if (FragmentShaderStream.is_open())
	{
		std::stringstream sstr;
		sstr << FragmentShaderStream.rdbuf();
		FragmentShaderCode = sstr.str();
		FragmentShaderStream.close();
	}

	// Read the Geometry Shader code from the file
	string GeometryShaderCode;
	std::ifstream GeometryShaderStream(geometry_file_path, std::ios::in);
	if (GeometryShaderStream.is_open())
	{
		std::stringstream sstr;
		sstr << GeometryShaderStream.rdbuf();
		GeometryShaderCode = sstr.str();
		GeometryShaderStream.close();
	}

	GLint Result = GL_FALSE;
	int InfoLogLength;

	// Compile Vertex Shader
	printf("Compiling vertex shader : %s\n", vertex_file_path);
	char const *VertexSourcePointer = VertexShaderCode.c_str();
	glShaderSource(VertexShaderID, 1, &VertexSourcePointer, NULL);
	glCompileShader(VertexShaderID);

	// Check Vertex Shader
	glGetShaderiv(VertexShaderID, GL_COMPILE_STATUS, &Result);
	glGetShaderiv(VertexShaderID, GL_INFO_LOG_LENGTH, &InfoLogLength);
	if (InfoLogLength > 0)
	{
		vector<char> VertexShaderErrorMessage(InfoLogLength + 1);
		glGetShaderInfoLog(VertexShaderID, InfoLogLength, NULL, &VertexShaderErrorMessage[0]);
		printf("%s\n", &VertexShaderErrorMessage[0]);
	}

	// Compile Geometry Shader
	printf("Compiling geometry shader : %s\n", geometry_file_path);
	char const *GeometrySourcePointer = GeometryShaderCode.c_str();
	glShaderSource(GeometryShaderID, 1, &GeometrySourcePointer, NULL);
	glCompileShader(GeometryShaderID);

	// Check Geometry Shader
	glGetShaderiv(GeometryShaderID, GL_COMPILE_STATUS, &Result);
	glGetShaderiv(GeometryShaderID, GL_INFO_LOG_LENGTH, &InfoLogLength);
	if (InfoLogLength > 0)
	{
		vector<char> GeometryShaderErrorMessage(InfoLogLength + 1);
		glGetShaderInfoLog(GeometryShaderID, InfoLogLength, NULL, &GeometryShaderErrorMessage[0]);
		printf("%s\n", &GeometryShaderErrorMessage[0]);
	}

	// Compile Fragment Shader
	printf("Compiling fragment shader : %s\n", fragment_file_path);
	char const *FragmentSourcePointer = FragmentShaderCode.c_str();
	glShaderSource(FragmentShaderID, 1, &FragmentSourcePointer, NULL);
	glCompileShader(FragmentShaderID);

	// Check Fragment Shader
	glGetShaderiv(FragmentShaderID, GL_COMPILE_STATUS, &Result);
	glGetShaderiv(FragmentShaderID, GL_INFO_LOG_LENGTH, &InfoLogLength);
	if (InfoLogLength > 0)
	{
		vector<char> FragmentShaderErrorMessage(InfoLogLength + 1);
		glGetShaderInfoLog(FragmentShaderID, InfoLogLength, NULL, &FragmentShaderErrorMessage[0]);
		printf("%s\n", &FragmentShaderErrorMessage[0]);
	}

	// Link the program
	printf("Linking program\n");
	GLuint ProgramID = glCreateProgram();
	glAttachShader(ProgramID, VertexShaderID);
	glAttachShader(ProgramID, FragmentShaderID);
	glAttachShader(ProgramID, GeometryShaderID);

	glLinkProgram(ProgramID);

	// Check the program
	glGetProgramiv(ProgramID, GL_LINK_STATUS, &Result);
	glGetProgramiv(ProgramID, GL_INFO_LOG_LENGTH, &InfoLogLength);
	if (InfoLogLength > 0)
	{
		vector<char> ProgramErrorMessage(InfoLogLength + 1);
		glGetProgramInfoLog(ProgramID, InfoLogLength, NULL, &ProgramErrorMessage[0]);
		printf("%s\n", &ProgramErrorMessage[0]);
	}

	glDetachShader(ProgramID, VertexShaderID);
	glDetachShader(ProgramID, GeometryShaderID);
	glDetachShader(ProgramID, FragmentShaderID);

	glDeleteShader(VertexShaderID);
	glDeleteShader(GeometryShaderID);
	glDeleteShader(FragmentShaderID);

	return ProgramID;
}

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
		default:
		throw std::invalid_argument("Resolution not recognized");
	}
	return -1;
}

int predefinedHeight(Resolution res) {
	switch (res) {
		case FHD:
			return 1080;
		case UHD:
			return 2160;
		default:
		throw std::invalid_argument("Resolution not recognized");
	}
	return -1;
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
	std::cout << "Error: unknown GLSLType" << std::endl;
	return -1;
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

mat4 generateMVP(vec3 camPosition, vec3 camLookAt, vec3 upVector, float fov, float aspectRatio, float clippingRangeMin, float clippingRangeMax, mat4 modelTransform)
{
	mat4 ViewMatrix = lookAt(camPosition, camLookAt, upVector);
	mat4 ProjectionMatrix = perspective(fov, aspectRatio, clippingRangeMin, clippingRangeMax);
	return ProjectionMatrix * ViewMatrix * modelTransform;
}

Window::Window(int width, int height, const char *title)
{
    if (!glfwInit())
        exit(2137);
	glfwWindowHint(GLFW_SAMPLES, 4);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
	this->width = width;
	this->height = height;
	this->aspectRatio = (float)width / (float)height;
	this->window = glfwCreateWindow(width, height, title, NULL, NULL);
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
	switch (resolution)
	{
	case FHD:
		this->width = 1920;
		this->height = 1080;
		break;
	case UHD:
		this->width = 3840;
		this->height = 2160;
		break;
	}
	this->aspectRatio = (float)width / (float)height;
	this->window = glfwCreateWindow(width, height, title, NULL, NULL);
	if (!this->window)
	{
		glfwTerminate();
		exit(2136);
	}
	glfwMakeContextCurrent(this->window);
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

void Window::setCallbacks(GLFWkeyfun *keyCallback, GLFWcharfun *charCallback, GLFWmousebuttonfun *mouseButtonCallback, GLFWcursorposfun *cursorPosCallback, GLFWcursorenterfun *cursorEnterCallback, GLFWscrollfun *scrollCallback, GLFWdropfun *dropCallback)
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
    exit(EXIT_SUCCESS);
    return 0;
}
void Window::initViewport() {
    glViewport(0, 0, this->width, this->height);
}

Shader::Shader(const char *vertex_file_path, const char *fragment_file_path)
{
	this->vertex_file_path = vertex_file_path;
	this->fragment_file_path = fragment_file_path;
	this->programID = LoadShaders(vertex_file_path, fragment_file_path);
	shaderType = CLASSIC;
}

Shader::Shader(const char *vertex_file_path, const char *fragment_file_path, const char *geometry_file_path)
{
	this->vertex_file_path = vertex_file_path;
	this->fragment_file_path = fragment_file_path;
	this->geometry_file_path = geometry_file_path;
	this->programID = LoadShaders(vertex_file_path, fragment_file_path, geometry_file_path);
	shaderType = GEOMETRY1;
}

Shader::Shader(string standard_file_path)
{
	string vertex_file_path_str = standard_file_path + ".vert";
	string fragment_file_path_str = standard_file_path + ".frag";
	this->vertex_file_path = vertex_file_path_str.c_str();
	this->fragment_file_path = fragment_file_path_str.c_str();
	this->programID = LoadShaders(vertex_file_path, fragment_file_path);
	this->uniformLocations = map<string, GLuint>();
	this->uniformTypes = map<string, GLSLType>();
	shaderType = CLASSIC;
}

Shader::~Shader()
{
	glDeleteProgram(this->programID);
}

void Shader::use()
{
	glUseProgram(this->programID);
}

void Shader::initUniforms(std::map<std::string, GLSLType> uniforms)
{
	for (auto uni : uniforms)
		this->uniformTypes[uni.first] = uni.second;
	for (auto const &uniform : uniforms)
	{
		GLuint uniformLocation = glGetUniformLocation(this->programID, uniform.first.c_str());
		this->uniformLocations[uniform.first] = uniformLocation;
	}
}

void Shader::setTextureSampler(Texture *texture, int slot)
{
	GLuint samplerUniform = glGetUniformLocation(this->programID, texture->samplerName);
	glUniform1i(samplerUniform, slot);
}

void Shader::setUniforms(map<string, const GLfloat *> uniformValues)
{
	for (auto const &uniform : uniformValues)
	{
		GLSLType uniformType = this->uniformTypes[uniform.first];
		setUniform(uniform.first, uniform.second);
	}
}

void Shader::setUniform(string uniformName, const GLfloat *uniformValue)
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

void Shader::setUniform(string uniformName, float uniformValue)
{
	if (this->uniformTypes[uniformName] != FLOAT)
		throw std::invalid_argument("Uniform type must be FLOAT");
	GLuint uniformLocation = this->uniformLocations[uniformName];
	glUniform1f(uniformLocation, uniformValue);
}

void Shader::setUniform(string uniformName, int uniformValue)
{
	if (this->uniformTypes[uniformName] != INT)
		throw std::invalid_argument("Uniform type must be INT");
	GLuint uniformLocation = this->uniformLocations[uniformName];
	glUniform1i(uniformLocation, uniformValue);
}

void Shader::setUniform(string uniformName, vec2 uniformValue)
{
	if (this->uniformTypes[uniformName] != VEC2)
		throw std::invalid_argument("Uniform type must be VEC2");
	GLuint uniformLocation = this->uniformLocations[uniformName];
	glUniform2f(uniformLocation, uniformValue.x, uniformValue.y);
}

void Shader::setUniform(string uniformName, vec3 uniformValue)
{
	if (this->uniformTypes[uniformName] != VEC3)
		throw std::invalid_argument("Uniform type must be VEC3");
	GLuint uniformLocation = this->uniformLocations[uniformName];
	glUniform3f(uniformLocation, uniformValue.x, uniformValue.y, uniformValue.z);
}

void Shader::setUniform(string uniformName, vec4 uniformValue)
{
	if (this->uniformTypes[uniformName] != VEC4)
		throw std::invalid_argument("Uniform type must be VEC4");
	GLuint uniformLocation = this->uniformLocations[uniformName];
	glUniform4f(uniformLocation, uniformValue.x, uniformValue.y, uniformValue.z, uniformValue.w);
}

void Shader::setUniform(string uniformName, mat2 uniformValue)
{
	if (this->uniformTypes[uniformName] != MAT2)
		throw std::invalid_argument("Uniform type must be MAT2");
	GLuint uniformLocation = this->uniformLocations[uniformName];
	glUniformMatrix2fv(uniformLocation, 1, GL_FALSE, &uniformValue[0][0]);
}

void Shader::setUniform(string uniformName, mat3 uniformValue)
{
	if (this->uniformTypes[uniformName] != MAT3)
		throw std::invalid_argument("Uniform type must be MAT3");
	GLuint uniformLocation = this->uniformLocations[uniformName];
	glUniformMatrix3fv(uniformLocation, 1, GL_FALSE, &uniformValue[0][0]);
}

void Shader::setUniform(string uniformName, mat4 uniformValue)
{
	if (this->uniformTypes[uniformName] != MAT4)
		throw std::invalid_argument("Uniform type must be MAT4");
	GLuint uniformLocation = this->uniformLocations[uniformName];
	glUniformMatrix4fv(uniformLocation, 1, GL_FALSE, &uniformValue[0][0]);
}

void Shader::setUniform(string uniformName, float x, float y)
{
	if (this->uniformTypes[uniformName] != VEC2)
		throw std::invalid_argument("Uniform type must be VEC2");
	GLuint uniformLocation = this->uniformLocations[uniformName];
	glUniform2f(uniformLocation, x, y);
}

void Shader::setUniform(string uniformName, float x, float y, float z)
{
	if (this->uniformTypes[uniformName] != VEC3)
		throw std::invalid_argument("Uniform type must be VEC3");
	GLuint uniformLocation = this->uniformLocations[uniformName];
	glUniform3f(uniformLocation, x, y, z);
}

void Shader::setUniform(string uniformName, float x, float y, float z, float w)
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
        glm::vec3 upVector, float fov_x, float aspectRatio, float clippingRangeMin, float clippingRangeMax) {

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
        const std::function<glm::vec3(float)> &upVector, float fov_x, float aspectRatio, float clippingRangeMin, float clippingRangeMax) {

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
    return projectionMatrix * viewMatrix(t) * modelTransform;	
}


Attribute::Attribute(std::string name, GLSLType type, int inputNumber)
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
	if (enabled)
		disable();
	if (bufferInitialized)
		freeBuffer();
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

RenderingStep::RenderingStep(const shared_ptr<Shader> &shader)
{
	this->shader = shader;
	this->attributes = vector<shared_ptr<Attribute>>();
	this->model = std::make_shared<Model3D>();
	this->uniforms = map<string, GLSLType>();
	this->uniformSetters = map<string, shared_ptr<std::function<void(float, shared_ptr<Shader>)>>>();
	this->customStep = [](float t) {};
}

RenderingStep::RenderingStep(const RenderingStep& other)
{
	this->shader = other.shader;
	this->attributes = other.attributes;
	this->model = other.model;
	this->uniforms = other.uniforms;
	this->uniformSetters = other.uniformSetters;
	this->customStep = other.customStep;
}

RenderingStep::~RenderingStep()
{
	for (auto attribute : attributes)
	{
		attribute.reset();
	}
	shader.reset();
	model.reset();
}

void RenderingStep::setModel(const std::shared_ptr<Model3D> &model)
{
	this->model = model;
}

void RenderingStep::setSuperMesh(const std::shared_ptr<SuperMesh> &super) {
	this->super = super;
	this-> model = nullptr;
    this-> weak_super = nullptr;
}

void RenderingStep::setWeakSuperMesh(const std::shared_ptr<WeakSuperMesh> &super) {
    this->super = nullptr;
    this-> model = nullptr;
    this-> weak_super = super;
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

	for (const auto& attribute : attributes)
	{
		attribute->initBuffer();
	}

}

void RenderingStep::resetAttributeBuffers()
{
	for (const auto& attribute : attributes)
	{
		attribute->freeBuffer();
	}
}

void RenderingStep::initUnusualAttributes(std::vector<std::shared_ptr<Attribute>> attributes) {
    this->attributes = attributes;
    for (const auto &attribute: attributes) {
        attribute->initBuffer();
    }
}



void RenderingStep::loadStandardAttributes() {
    if (weakSuperLoaded())
    {
        for (auto i = 0; i < 4; i++)
            attributes[i]->load(weak_super->getBufferLocation(CommonBufferType(i)), weak_super->getBufferLength(CommonBufferType(i)));
        return;
    }

	if (superLoaded())
	{
		for (int i = 0; i < 8; i++)
			attributes[i]->load(super->bufferLocations[i], super->bufferSizes[i]);
		return;
	}

	for (int i = 0; i < 4; i++)
		attributes[i]->load(model->mesh->bufferLocations[i], model->mesh->bufferSizes[i]);
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

void RenderingStep::addUniform(string uniformName, GLSLType uniformType, shared_ptr<std::function<void(float, shared_ptr<Shader>)>> setter)
{
	this->uniforms[uniformName] = uniformType;
	this->uniformSetters[uniformName] = std::move(setter);
	this->shader->initUniforms({{uniformName, uniformType}});
}

void RenderingStep::addConstFloats(const std::map<std::string, float>& uniforms)
{
	for (auto uniform : uniforms)
	{
		addUniform(uniform.first, FLOAT, std::make_shared<std::function<void(float, shared_ptr<Shader>)>>
			([uniform](float t, const shared_ptr<Shader>& shader)
				{
					shader->setUniform(uniform.first, uniform.second);	
				}));
	}
}

void RenderingStep::addConstVec4(const std::string& uniformName, vec4 value)
{
	auto uniformSetter = [value, uniformName](float t, const shared_ptr<Shader> &shader) {
		shader->setUniform(uniformName, value);
	};
	addUniform(uniformName, VEC4,  std::make_shared<std::function<void(float, shared_ptr<Shader>)>>(uniformSetter));
}

void RenderingStep::setUniforms(float t)
{
	for (const auto uniform_descriptor : uniforms)
	{
		auto name = uniform_descriptor.first;
		(*uniformSetters[name])(t, shader);
	}
}

void RenderingStep::addCameraUniforms(const std::shared_ptr<Camera>& camera)
{
	std::function<void(float, shared_ptr<Shader>)> MVPsetter;

	if (superLoaded() || weakSuperLoaded()) {
		MVPsetter = [camera, this](float t, const shared_ptr<Shader> &shader) {
			mat4 mvp = camera->mvp(t, mat4(mat3(1)));
			shader->setUniform("mvp", mvp);
		};
	}
	else {
		MVPsetter = [camera, this](float t, const shared_ptr<Shader> &shader) {
			mat4 mvp = camera->mvp(t, model->transform);
			shader->setUniform("mvp", mvp);
		};
	}
	auto positionSetter = [camera, this](float t, const shared_ptr<Shader> &shader) {
		vec3 camPos = camera->position(t);
		shader->setUniform("camPosition", camPos);
	};
	addUniform("mvp", MAT4,  std::make_shared<std::function<void(float, shared_ptr<Shader>)>>(MVPsetter));
	addUniform("camPosition", VEC3,  std::make_shared<std::function<void(float, shared_ptr<Shader>)>>(positionSetter));
}

void RenderingStep::addLightUniform(const std::shared_ptr<PointLight>& pointLight, int lightIndex)
{
	string lightName = "light" + std::to_string(lightIndex);
	auto lightSetter = [pointLight, this, lightName](float t, const shared_ptr<Shader> &shader) {
		mat4 lightMat = pointLight->compressToMatrix();
		shader->setUniform(lightName, lightMat);
	};
	addUniform(lightName, MAT4, std::make_shared<std::function<void(float, shared_ptr<Shader>)>>(lightSetter));
}

void RenderingStep::addLightsUniforms(const std::vector<std::shared_ptr<PointLight>> &lights)
{
	for (int i = 1; i <= lights.size(); i++)
		addLightUniform(lights[i-1], i);
}

void RenderingStep::addMaterialUniform() {
    if (weakSuperLoaded())
        throw std::invalid_argument("Material uniform can only be added for super mesh");
    mat4 lightMat = weak_super->getMaterial().compressToMatrix();
    auto materialSetter = [lightMat](float t, const shared_ptr<Shader> &s) { s->setUniform("material", lightMat); };
    addUniform("material", MAT4, std::make_shared<std::function<void(float, shared_ptr<Shader>)>>(materialSetter));
}
void RenderingStep::addTexturedMaterialUniforms() {
    if (!weakSuperLoaded())
        throw std::invalid_argument("Material uniform can only be added for super mesh");
    vec4 intenc = weak_super->getMaterial().compressIntencities();
    auto materialSetter = [intenc](float t, const shared_ptr<Shader> &s) { s->setUniform("intencities", intenc); };
    addUniform("intencities", VEC4, std::make_shared<std::function<void(float, shared_ptr<Shader>)>>(materialSetter));

    auto textureSetter = [m=weak_super->getMaterial().texture_ambient](float t, const shared_ptr<Shader> &s) {
        s->setTextureSampler(m.get(), m->textureSlot);
    };
    addUniform("texture_ambient", SAMPLER2D, std::make_shared<std::function<void(float, shared_ptr<Shader>)>>(textureSetter));

    auto textureSetter2 = [i=weak_super->getMaterial().texture_diffuse->textureSlot, tex=weak_super->getMaterial().texture_diffuse.get()](float t, const shared_ptr<Shader> &s) {
        setUniformTextureSampler(s->programID, tex, i);
    };
    addUniform("texture_diffuse", SAMPLER2D, std::make_shared<std::function<void(float, shared_ptr<Shader>)>>(textureSetter2));

    auto textureSetter3 = [this](float t, const shared_ptr<Shader> &s) {
        setUniformTextureSampler(s->programID, weak_super->getMaterial().texture_specular.get(),
                                 weak_super->getMaterial().texture_specular->textureSlot);
    };
    addUniform("texture_specular", SAMPLER2D, std::make_shared<std::function<void(float, shared_ptr<Shader>)>>(textureSetter3));
}

void RenderingStep::init(const shared_ptr<Camera> &cam, const std::vector<shared_ptr<PointLight>> &lights) {
    if (weakSuperLoaded()) {
        shader->use();
        initElementBuffer();
        initStdAttributes();
        addCameraUniforms(cam);
        addLightsUniforms(lights);
        initTextures();
        addTexturedMaterialUniforms();
        loadStandardAttributes();
        loadElementBuffer();
    }
}
void RenderingStep::initTextures() { weak_super->initGlobalTextures(); }

void RenderingStep::bindTextures() {
    weak_super->getMaterial().texture_ambient->bind();
    weak_super->getMaterial().texture_diffuse->bind();
    weak_super->getMaterial().texture_specular->bind();
}

void RenderingStep::addUniforms(const map<std::string, GLSLType> &uniforms, map<string, shared_ptr<std::function<void(float, shared_ptr<Shader>)>>> setters)
{
	for (const auto& uniform : uniforms)
		addUniform(uniform.first, uniform.second, setters[uniform.first]);
}

void RenderingStep::addCustomAction(const std::function<void(float)> &action)
{
	this->customStep = action;
}


void RenderingStep::weakMeshRenderStep(float t) {
    bindTextures();
    loadStandardAttributes();
    customStep(t);
    setUniforms(t);
    loadElementBuffer();

    enableAttributes();
    glDrawElements(GL_TRIANGLES, weak_super->bufferIndexSize(), GL_UNSIGNED_INT, 0);
    disableAttributes();
}

void RenderingStep::superMeshRenderStep(float t) {

    loadStandardAttributes();
    enableAttributes();
    customStep(t);
    setUniforms(t);
    glDrawArrays(GL_TRIANGLES, 0, super->bufferSizes[0]*3);
    disableAttributes();
}
void RenderingStep::modelRenderStep(float t) {
    loadStandardAttributes();
    enableAttributes();
    customStep(t);
    setUniforms(t);
    glDrawArrays(GL_TRIANGLES, 0, model->mesh->bufferSizes[0]*3);
    disableAttributes();
}

void RenderingStep::renderStep(float t)
{
    shader->use();
    if (weakSuperLoaded())
        weakMeshRenderStep(t);
    else if (superLoaded())
        superMeshRenderStep(t);
    else
        modelRenderStep(t);
}

Renderer::Renderer(float animSpeed, vec4 bgColor)
{
	this->window = nullptr;
	this->vao = 0;
	this->camera = nullptr;
	this->time = 0;
	this->lights = vector<shared_ptr<PointLight>>();
	if (!glfwInit())
		exit(2137);
	this->bgColor = bgColor;
	this->animSpeed = [animSpeed](float t) { return animSpeed; };
	this->perFrameFunction = make_unique<std::function<void(float, float)>>([](float t, float delta){});
}

Renderer::Renderer(int width, int height, const char *title,
                   const std::shared_ptr<Camera> &camera,
                   const std::vector<std::shared_ptr<PointLight>> &lights,
						const std::vector<std::shared_ptr<RenderingStep>>& renderingSteps,
						float animSpeed, vec4 bgColor) : Renderer(animSpeed, bgColor)
{
	initMainWindow(width, height, title);
	setCamera(camera);
	setLights(lights);
	for (const auto& renderingStep : renderingSteps)
		addRenderingStep(renderingStep);
}

Renderer::Renderer(Resolution resolution, const char *title,
                   const std::shared_ptr<Camera> &camera,
                   const std::vector<std::shared_ptr<PointLight>> &lights,
						const std::vector<std::shared_ptr<RenderingStep>>& renderingSteps,
						float animSpeed, vec4 bgColor) :
			Renderer(predefinedWidth(resolution), predefinedHeight(resolution), title, camera, lights,
				renderingSteps, animSpeed, bgColor) {}

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
	if (glewInit() != GLEW_OK)
		exit(2138);
	this->vao = bindVAO();
}

void Renderer::initMainWindow(Resolution resolution, const char *title) 
{
	this->window = std::make_unique<Window>(resolution, title);
	glewExperimental = true;
	if (glewInit() != GLEW_OK)
		exit(2138);
	this->vao = bindVAO();
}

void Renderer::addRenderingStep(std::shared_ptr<RenderingStep> renderingStep)
{

	this->renderingSteps.push_back(std::move(renderingStep));
}

void Renderer::setCamera(const std::shared_ptr<Camera> &camera)
{
	this->camera = camera;
}

void Renderer::setLights(const std::vector<std::shared_ptr<PointLight>> &lights)
{
	this->lights = lights;
}

float Renderer::initFrame()
{
	// glViewport(0, 0, window->width, window->height);
	glClearColor(this->bgColor.x, this->bgColor.y, this->bgColor.z, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);

    this->dt = (glfwGetTime() - this->frameOlderTimeThanThePublicOne)*animSpeed(this->time);
	this->time += dt;
	this->frameOlderTimeThanThePublicOne = glfwGetTime();
	return this->time;
}

float Renderer::lastDeltaTime() const {
	if (this->frameOlderTimeThanThePublicOne == NULL)
		return 0.f;
    return this->dt;
}

void Renderer::addPerFrameUniforms(std::map<std::string, GLSLType> uniforms, std::map<std::string, shared_ptr<std::function<void(float, std::shared_ptr<Shader>)>>> setters)
{
	for (const auto& renderingStep : renderingSteps)
		renderingStep->addUniforms(uniforms, setters);
	
}

void Renderer::addPerFrameUniform(std::string uniformName, GLSLType uniformType, shared_ptr<std::function<void(float, std::shared_ptr<Shader>)>> setter)
{
	for (const auto& renderingStep : renderingSteps)
		renderingStep->addUniform(uniformName, uniformType, setter);
}

void Renderer::initRendering()
{
    window->initViewport();
    glBindVertexArray(vao);
    // glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);
    glShadeModel(GL_SMOOTH);
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
 //    glEnable(GL_LIGHTING);
	// glDepthFunc(GL_3D);
	for (const auto& renderingStep : renderingSteps)
	{
	    if (renderingStep->weakSuperLoaded()) {
	        renderingStep->init(camera, lights);
	    }
		else if (renderingStep->superLoaded()) {
			renderingStep->initStdAttributes();
			renderingStep->initMaterialAttributes();
			renderingStep->addCameraUniforms(camera);
			renderingStep->addLightsUniforms(lights);
		}
		else {
			renderingStep->initStdAttributes();
			renderingStep->addCameraUniforms(camera);
			renderingStep->addLightsUniforms(lights);
			renderingStep->addMaterialUniform();
		}
	}
}

void Renderer::addConstUniform(std::string uniformName, GLSLType uniformType, shared_ptr<std::function<void(std::shared_ptr<Shader>)>> setter)
{
	std::function<void(float, std::shared_ptr<Shader>)> setterWrapper = [setter](float t, std::shared_ptr<Shader> shader) { (*setter)(shader); };
	for (const auto& renderingStep : renderingSteps)
		renderingStep->addUniform(uniformName, uniformType, make_shared<std::function<void(float, std::shared_ptr<Shader>)>>(setterWrapper));
}

void Renderer::addTimeUniform()
{
	addPerFrameUniform("time", FLOAT, std::make_shared<std::function<void(float, std::shared_ptr<Shader>)>>(
	[](float t, std::shared_ptr<Shader> shader) {
		shader->setUniform("time", t);
	}));
}

void Renderer::addConstFloats(std::map<std::string, float> uniforms)
{
	for (const auto& renderingStep : renderingSteps)
		renderingStep->addConstFloats(uniforms);
}

void Renderer::addCustomAction(std::function<void(float)> action)
{
	perFrameFunction = make_unique<std::function<void(float, float)>>([a=*this->perFrameFunction, n=action](float t, float delta) {
		a(t, delta);
		n(t);
	});
}

void Renderer::addCustomAction(std::function<void(float, float)> action)
{
    perFrameFunction = make_unique<std::function<void(float, float)>>([a=*this->perFrameFunction, n=action](float t, float delta) {
        a(t, delta);
        n(t, delta);
    });
}

void Renderer::addConstUniforms(std::map<std::string, GLSLType> uniforms, std::map<std::string, shared_ptr<std::function<void(std::shared_ptr<Shader>)>>> setters)
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
    while (window->isOpen()) {
        initFrame();
        (*perFrameFunction)(time, dt);
        renderAllSteps();
        window->renderFramebufferToScreen();
    }
    return window->destroy();
}
