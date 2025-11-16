#pragma once
#include <GL/glew.h>
#include <glfw/glfw3.h>
#include "shaderTypes.hpp"
#include "window.hpp"

using gl_id = GLuint;
using gl_uniform_loc = GLint;

class GLCommand {
	static gl_id currentProgram;
	static gl_id currentVAO;

public:
	static string formatGLenum(GLenum value);

	static void bindProgram(gl_id program);
	static void unbindProgram();
	static gl_id createProgram();
	static void attachShaderToProgram(gl_id programID, gl_id shaderID);
	static void linkProgram(gl_id programID);
	static void deleteProgram(gl_id programID);

	static gl_uniform_loc getUniformLocation(const string& name);
	static gl_uniform_loc getUniformLocation(const string& name, gl_id program);
	static void setUniform(gl_uniform_loc location, raw_data_ptr data, GLSLPrimitive type, array_len arraySize=1);
	static void setUniform(gl_uniform_loc location, int value);
	static void setUniform(gl_uniform_loc location, float value);
	static void setUniform(gl_uniform_loc location, uint value);
	static void setUniform(gl_uniform_loc location, bool value);

	static gl_id createShader(GLenum shaderType);
	static void compileShaderSource(gl_id shaderID, const string& source);
	static void deleteShader(gl_id shaderID);

	static void bindVAO(gl_id vao);
	static void unbindVAO();
	static void deleteVAO(gl_id vao);
	static void createVAO(gl_id* id);

	static void createTexture(gl_id* id);
	static void deleteTexture(gl_id id);
	static void bindTexture2D(gl_id id, unsigned int slot);
	static void unbindTexture2D();
	static void setTexture2DFilters(GLenum minFilter, GLenum magFilter, GLenum wrapS, GLenum wrapT);
	static void calculateMipmapsTexture2D();
	static void loadTexture2(GLenum internalFormat, array_len width, array_len height, GLenum format, data_ptr<uchar> data);
	static void setSampler2D(const string& samplerName, uint slot);

	static void createBuffer(gl_id* id);
	static void deleteBuffer(gl_id* id);
	static void loadBufferData(gl_id bufferID, raw_data_ptr firstElementAdress, byte_size bufferSize, GLenum usage);
	static void updateBufferData(gl_id bufferID, raw_data_ptr firstElementAdress, byte_size bufferSize);

	static void bindSSBO(gl_id bufferID, int bindingPoint);
	static void unbindSSBO();
	static void loadSSBOData(gl_id bufferID, raw_data_ptr firstElementAdress, byte_size bufferSize);
	static void updateSSBOData(gl_id bufferID, raw_data_ptr firstElementAdress, byte_size bufferSize);

	static void bindUBO(gl_id bufferID, int bindingPoint);
	static void unbindUBO();
	static void loadUBOData(gl_id bufferID, raw_data_ptr firstElementAdress, byte_size bufferSize);
	static void updateUBOData(gl_id bufferID, raw_data_ptr firstElementAdress, byte_size bufferSize);

	static void drawIndexedTriangles(array_len numberOfIndices);
	static void runComputeShader(int numGroupsX, int numGroupsY, int numGroupsZ);

	static void pointAtElementBuffer(gl_id buffer);
	static void pointAtAttributeBuffer(int inputNumber, gl_id bufferAddress, GLSLPrimitive type, byte_size stride,  raw_data_ptr offset);

	static void clearScreen(const vec4& color);
	static void enableBlending();
	static void enableDepth();
	static void init();
	static void initViewport(int width, int height);
};

class GLFWCommand {
public:
	static void init();
	static void terminate();

	static GLFWwindow* createWindow(int width, int height, const string& title);
	static void destroyWindow(GLFWwindow* window);
	static void pollEvents();
	static void swapFrameBuffers(GLFWwindow* window);

	static void setWindowData(GLFWwindow* window, data_ptr_mut<WindowSettings> data);
	static WindowSettings& getWindowData(GLFWwindow* window);
	static vec2 getCursorPosition(GLFWwindow* window);
};