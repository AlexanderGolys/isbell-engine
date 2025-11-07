#pragma once
#include <GL/glew.h>
#include "exceptions.hpp"
#include "logging.hpp"
#include "shaderTypes.hpp"


class GLCommand {
	static GLuint currentProgram;
	static GLuint currentVAO;
public:

	static void bindProgram(GLuint program) {
		if (currentProgram != program) {
			glUseProgram(program);
			currentProgram = program;
		}
	}

	static void unbindProgram() {
		glUseProgram(0);
		currentProgram = 0;
	}

	static GLint getUniformLocation(const string& name) {
		THROW_IF(currentProgram == 0, RuntimeError, "No shader program is currently bound.");
		GLint location = glGetUniformLocation(currentProgram, name.c_str());
		THROW_IF(location == -1, RuntimeError, "Uniform '" + name + "' not found in current shader program.");
		return location;
	}


	static GLint getUniformLocation(const string& name, GLuint program) {
		GLint location = glGetUniformLocation(program, name.c_str());
		THROW_IF(location == -1, RuntimeError, "Uniform '" + name + "' not found in given shader program.");
		return location;
	}

	static string formatGLenum(GLenum value) {
		switch (value) {
			case GL_VERTEX_SHADER: return "VERTEX SHADER";
			case GL_FRAGMENT_SHADER: return "FRAGMENT SHADER";
			case GL_GEOMETRY_SHADER: return "GEOMETRY SHADER";
			case GL_COMPUTE_SHADER: return "COMPUTE SHADER";
			default: return "UNKNOWN";
		}
	}

	static GLuint createShader(GLenum shaderType) {
		LOG(1, "Created shader of type " + formatGLenum(shaderType));
		return glCreateShader(shaderType);
	}

	static void compileShaderSource(GLuint shaderID, const string& source) {
		char const* sourcePtr = source.c_str();
		glShaderSource(shaderID, 1, &sourcePtr, nullptr);
		glCompileShader(shaderID);

		int InfoLogLength;
		glGetShaderiv(shaderID, GL_INFO_LOG_LENGTH, &InfoLogLength);
		if (InfoLogLength > 0) {
			vector<char> error(InfoLogLength + 1);
			glGetShaderInfoLog(shaderID, InfoLogLength, nullptr, &error[0]);
			THROW(SystemError, &error[0]);
		}
		LOG(0, "Shader compiled successfully.");
	}

	static void deleteShader(GLuint shaderID) {
		glDeleteShader(shaderID);
		LOG(1, "Shader deleted.");
	}

	static GLuint createProgram() {
		GLuint programID = glCreateProgram();
		LOG(1, "New shader program created.");
		return programID;
	}

	static void deleteProgram(GLuint programID) {
		if (currentProgram == programID)
			unbindProgram();

		glDeleteProgram(programID);
		LOG(1, "Shader program deleted.");
	}

	static void bindVAO(GLuint vao) {
		if (currentVAO != vao) {
			glBindVertexArray(vao);
			currentVAO = vao;
		}
	}

	static void unbindVAO() {
		glBindVertexArray(0);
		currentVAO = 0;
	}

	static void deleteVAO(GLuint vao) {
		if (currentVAO == vao)
			unbindVAO();

		glDeleteVertexArrays(1, &vao);
		LOG(1, "VAO deleted.");
	}

	static void createVAO(GLuint *id) {
		glGenVertexArrays(1, id);
		LOG(1, "VAO created.");
	}

	static void createBuffer(GLuint *id) {
		glCreateBuffers(1, id);
		LOG(1, "Buffer created.");
	}
	static void deleteBuffer(GLuint id) {
		glDeleteBuffers(1, &id);
		LOG(1, "Buffer deleted.");
	}

	static void drawIndexedTriangles(array_len numberOfIndices) {
		THROW_IF(currentVAO == 0, RuntimeError, "No VAO is currently bound.");
		THROW_IF(numberOfIndices == 0, ValueError, "Trying to draw 0 indices.");
		glDrawElements(GL_TRIANGLES, static_cast<GLsizei>(numberOfIndices), GL_UNSIGNED_INT, nullptr);
	}

	static void pointAtElementBuffer(GLuint buffer) {
		THROW_IF(currentVAO == 0, RuntimeError, "No VAO is currently bound.");
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, buffer);
	}

	static void pointAtAttributeBuffer(int inputNumber, GLuint bufferAddress, GLSLType type) {
		THROW_IF(currentVAO == 0, RuntimeError, "No VAO is currently bound.");
		THROW_IF(not supportedAttributeType(type), ValueError, "Unsupported GLSLType for attribute buffer (supported: float/double or vectors of these).");
		glBindBuffer(GL_ARRAY_BUFFER, bufferAddress);
		if (primitiveGLSLType(type) == GL_DOUBLE)
			glVertexAttribLPointer(inputNumber, lengthOfGLSLType(type), GL_DOUBLE, 0, (raw_data_ptr)0);
		else
			glVertexAttribPointer(inputNumber, lengthOfGLSLType(type), GL_FLOAT, GL_FALSE, 0, (raw_data_ptr)0);
		glEnableVertexAttribArray(inputNumber);
	}

	static void clearScreen(const vec4& color) {
		glClearColor(color.x, color.y, color.z, color.w);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	}


};
