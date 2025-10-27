#pragma once
#include "filesUtils.hpp"
#include "GL/glew.h"

class Shader {
	GLenum shaderType;
	GLuint shaderID;

public:
	explicit Shader(const Path& file);
	~Shader();

	GLuint getID() const;
	GLenum getType() const;
};


class ShaderProgram {
public:
	virtual ~ShaderProgram();
	virtual GLuint getProgramID() const = 0;
	virtual void bind() const = 0;
	virtual void unbind() const = 0;
};

class ClassicShaderProgram : public ShaderProgram {
	shared_ptr<Shader> vertexShader;
	shared_ptr<Shader> fragmentShader;
	GLuint programID;

public:
	ClassicShaderProgram(const shared_ptr<Shader>& vertexShader, const shared_ptr<Shader>& fragmentShader);
	~ClassicShaderProgram() override;
	GLuint getProgramID() const override;
	void bind() const override;
	void unbind() const override;
};

class GeometryShaderProgram : public ShaderProgram {
	shared_ptr<Shader> vertexShader;
	shared_ptr<Shader> fragmentShader;
	shared_ptr<Shader> geometryShader;
	GLuint programID;

public:
	GeometryShaderProgram(const shared_ptr<Shader>& vertexShader, const shared_ptr<Shader>& fragmentShader, const shared_ptr<Shader>& geometryShader);
	~GeometryShaderProgram() override;
	GLuint getProgramID() const override;
	void bind() const override;
	void unbind() const override;
};

class ComputeShaderProgram : public ShaderProgram {
	shared_ptr<Shader> computeShader;
	GLuint programID;

public:
	explicit ComputeShaderProgram(const shared_ptr<Shader>& computeShader);
	~ComputeShaderProgram() override;
	GLuint getProgramID() const override;
	void bind() const override;
	void unbind() const override;
	void run(int numGroupsX, int numGroupsY = 1, int numGroupsZ = 1) const;
};
