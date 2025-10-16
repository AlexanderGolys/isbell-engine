#pragma once

#include "filesUtils.hpp"
#include "shaderTypes.hpp"

ShaderType identify_shader_extension(const string &extension);




class Shader {
public:
	virtual ~Shader() = default;

	virtual ShaderType getType() const = 0;

	static shared_ptr<Shader> compile_shader(CodeFileDescriptor &file);
};


class ShaderProgram {
public:
	virtual ~ShaderProgram() = default;
	virtual void bind() = 0;
	virtual void unbind() = 0;

	virtual void setUniform(const string &uniformName, std::any uniformValue, UniformType type) = 0;
	virtual void setUniformArray(const string &uniformName, std::any uniformValue, UniformType type, int length) = 0;

	static shared_ptr<ShaderProgram> create_program(const shared_ptr<Shader> &vertexShader, const shared_ptr<Shader> &fragmentShader, const shared_ptr<Shader> &geometryShader = nullptr);
};

class ComputeShaderProgram : public ShaderProgram {
public:
	virtual void run(int numGroupsX, int numGroupsY=1, int numGroupsZ=1) = 0;

	static shared_ptr<ComputeShaderProgram> create_compute_shader(const shared_ptr<Shader> &computeShader);
};
