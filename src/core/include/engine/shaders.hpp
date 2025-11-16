#pragma once
#include <GL/glew.h>
#include "filesUtils.hpp"
#include "glCommand.hpp"


class Shader {
protected:
	CONST_PROPERTY(GLenum, shaderType);
	CONST_PROPERTY(gl_id, shaderID);

	void compile(const string& code);
public:
	explicit Shader(const Path& file);
	Shader(const string& code, GLenum shaderType);
	~Shader();

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
	CONST_PROPERTY(gl_id, programID);

public:
	ShaderProgram();
	ShaderProgram(const Shader& vertexShader, const Shader& fragmentShader);
	ShaderProgram(const Shader& vertexShader, const Shader& fragmentShader, const Shader& geometryShader);
	~ShaderProgram();

	void bind() const;
	void unbind() const;

	static sptr<ShaderProgram> standardShaderProgram(Path vertexShaderPath, Path fragmentShaderPath);
	static sptr<ShaderProgram> geometryShaderProgram(Path vertexShaderPath, Path fragmentShaderPath, Path geometryShaderPath);
};

class ComputeShaderProgram : public ShaderProgram {
public:
	explicit ComputeShaderProgram(const Shader& computeShader);
	void run(int numGroupsX, int numGroupsY=1, int numGroupsZ=1) const;

	static sptr<ComputeShaderProgram> standardComputeShaderProgram(Path computeShaderPath);
};

