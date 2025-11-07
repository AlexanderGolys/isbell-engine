#pragma once
#include <GL/glew.h>
#include "filesUtils.hpp"
#include "macroParsing.hpp"
#include "renderingUtils.hpp"
#include "shaderTypes.hpp"
#include "glCommand.hpp"

enum ShaderType {
	CLASSIC,
	GEOMETRY1,
};


class Shader {
protected:
	GLenum shaderType;
	GLuint shaderID;

	void compile(const string& code);
public:
	explicit Shader(const Path& file);
	Shader(const string& code, GLenum shaderType);
	~Shader();

	GLuint getID() const;
	GLenum getType() const;
	static GLenum getTypeFromExtension(const string& extension);
};

class VertexShader : public Shader {
public:
	explicit VertexShader(const Path& file);
	explicit VertexShader(const string& code);
};

class FragmentShader : public Shader {
public:
	explicit FragmentShader(const Path& file);
	explicit FragmentShader(const string& code);
};

class GeometryShader : public Shader {
public:
	explicit GeometryShader(const Path& file);
	explicit GeometryShader(const string& code);
};

class ComputeShader : public Shader {
public:
	explicit ComputeShader(const Path& file);
	explicit ComputeShader(const string& code);
};


class ShaderProgram {
	unordered_map<string, GLuint> cachedUniformLocations;
	GLuint programID = 0;

public:
	ShaderProgram(const Shader& vertexShader, const Shader& fragmentShader);
	ShaderProgram(const Shader& vertexShader, const Shader& fragmentShader, const Shader& geometryShader);
	~ShaderProgram();

	void bind() const;
	void unbind() const;
	GLuint getID() const;

	GLuint getUniformLocation(const string& uniformName);
};

