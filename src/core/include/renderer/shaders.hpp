#pragma once
#include "exceptions.hpp"


class Shader {
public:
	virtual ~Shader() = default;
	virtual void bind();
	virtual void unbind();

	virtual void setUniform(const string &uniformName, float value);
};
