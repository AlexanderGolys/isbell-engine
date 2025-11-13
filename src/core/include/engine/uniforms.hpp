#pragma once
#include "glCommand.hpp"
#include "renderLayers.hpp"
#include "shaderTypes.hpp"
#include "GL/glew.h"

class UniformComponent : public LayerComponent {
	PROPERTY(string, name);

public:
	explicit UniformComponent(const string& name);
	virtual void setUniformGL(gl_uniform_loc location) const = 0;
	void init() override {}
	void setDuringRender() const override;
};

template<typename T>
class ConstantUniform : public UniformComponent {
	PROPERTY(T, value);
public:
	ConstantUniform(const string& name, const T& value);
	void update(float t, float dt) override {}
	void setUniformGL(gl_uniform_loc location) const override;
};

template<typename T>
class ConstantUniformArray : public UniformComponent {
	PROPERTY(vector<T>, value);
public:
	ConstantUniformArray(const string& name, const vector<T>& value);
	void update(float t, float dt) override {}
	void setUniformGL(gl_uniform_loc location) const override;
	void set_value(array_len index, const T& newValue) { value[index] = newValue; }
};

class TimeUniform : public UniformComponent {
	float value;
public:
	TimeUniform() : UniformComponent("u_time"), value(0.0f) {}
	void update(float t, float dt) override { value = t;  }
	void setUniformGL(gl_uniform_loc location) const override;
};




// ------------------- Implementation ------------------ //


template <typename T>
ConstantUniform<T>::ConstantUniform(const string& name, const T& value)
: UniformComponent(name), value(value) {
	THROW_IF(not supportedUniformType(glslTypeVariable<T>), ValueError, "Type used in ConstantUniform is not a supported GLSL uniform type");
}

template <typename T>
void ConstantUniform<T>::setUniformGL(GLint location) const {
    if constexpr (std::is_same_v<T, float> or std::is_same_v<T, int> or std::is_same_v<T, unsigned int> or std::is_same_v<T, bool>)
    	GLCommand::setUniform(location, value);
	else if constexpr (std::is_same_v<T, vec2> or std::is_same_v<T, vec3> or std::is_same_v<T, vec4> or
	               std::is_same_v<T, ivec2> or std::is_same_v<T, ivec3> or std::is_same_v<T, ivec4>)
		GLCommand::setUniform(location, reinterpret_cast<raw_data_ptr>(&value[0]), glslTypeVariable<T>);
	else if constexpr (std::is_same_v<T, mat2> or std::is_same_v<T, mat3> or std::is_same_v<T, mat4> or
	               std::is_same_v<T, mat2x3> or std::is_same_v<T, mat2x4> or std::is_same_v<T, mat3x2> or
	               std::is_same_v<T, mat3x4> or std::is_same_v<T, mat4x2> or std::is_same_v<T, mat4x3>)
		GLCommand::setUniform(location, reinterpret_cast<raw_data_ptr>(&value[0][0]), glslTypeVariable<T>);
	else
		THROW(UnknownVariantError, "Unsupported GLSL type in ConstantUniform");
}

template <typename T>
ConstantUniformArray<T>::ConstantUniformArray(const string& name, const vector<T>& value): UniformComponent(name), value(value) {
	THROW_IF(not supportedUniformType(glslTypeVariable<T>), ValueError, "Type used in ConstantUniformArray is not a supported GLSL uniform type");
}

template <typename T>
void ConstantUniformArray<T>::setUniformGL(GLint location) const {
	if constexpr (std::is_same_v<T, float> or std::is_same_v<T, int> or std::is_same_v<T, unsigned int> or std::is_same_v<T, bool>)
		GLCommand::setUniform(location, reinterpret_cast<raw_data_ptr>(&value[0]), glslTypeVariable<T>, value.size());
	else if constexpr (std::is_same_v<T, vec2> or std::is_same_v<T, vec3> or std::is_same_v<T, vec4> or
				   std::is_same_v<T, ivec2> or std::is_same_v<T, ivec3> or std::is_same_v<T, ivec4>)
		GLCommand::setUniform(location, reinterpret_cast<raw_data_ptr>(&value[0][0]), glslTypeVariable<T>, value.size());
	else if constexpr (std::is_same_v<T, mat2> or std::is_same_v<T, mat3> or std::is_same_v<T, mat4> or
				   std::is_same_v<T, mat2x3> or std::is_same_v<T, mat2x4> or std::is_same_v<T, mat3x2> or
				   std::is_same_v<T, mat3x4> or std::is_same_v<T, mat4x2> or std::is_same_v<T, mat4x3>)
		GLCommand::setUniform(location, reinterpret_cast<raw_data_ptr>(&value[0][0][0]), glslTypeVariable<T>, value.size());
	else
		THROW(UnknownVariantError, "Unsupported GLSL type in ConstantUniform");
}


