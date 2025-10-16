#pragma once
#include "logging.hpp"
#include "shaders.hpp"

#include "GL/glew.h"


GLenum gl_shader_type(ShaderType type);




class ShaderGL : public Shader {
	ShaderType type;
	GLuint shaderID = 0;
public:
	ShaderGL(CodeFileDescriptor &file);

	ShaderType getType() const override;
};


class ShaderProgramGL : public ShaderProgram {

public:
	ShaderProgramGL(const shared_ptr<Shader> &vertexShader, const shared_ptr<Shader> &fragmentShader, const shared_ptr<Shader> &geometryShader = nullptr) {

	}

	void bind() override;
	void unbind() override;
	void setUniform(const string &uniformName, std::any uniformValue, UniformType type) override;
	void setUniformArray(const string &uniformName, std::any uniformValue, UniformType type, int length) override;
};


class ComputeShaderProgramGL : public ComputeShaderProgram {

public:
	ComputeShaderProgramGL(const shared_ptr<Shader> &computeShader);

	void bind() override;
	void unbind() override;
	void setUniform(const string &uniformName, std::any uniformValue, UniformType type) override;
	void setUniformArray(const string &uniformName, std::any uniformValue, UniformType type, int length) override;
	void run(int numGroupsX, int numGroupsY=1, int numGroupsZ=1) override;
};
