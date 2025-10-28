#pragma once

#include "exceptions.hpp"
#include "shaderDataTypes.hpp"
#include <GL/glew.h>

class Uniform {
	string name;
public:
	explicit Uniform(const string& name);
	virtual ~Uniform() = default;
	virtual void set(GLuint programID) const = 0;
	virtual void update(float time, float dt) = 0;

	string getName() const;
};

template<typename T>
class CachedUniform : public Uniform {
	T value;
	shared_ptr<HOM(float, T)> updater;
public:
	CachedUniform(const string& name, const T& value, shared_ptr<HOM(float, T)> updater=nullptr);

	void update(float time, float dt) override;
	T getValue() const;
	void setValue(const T& newValue);
};

class UniformFloat : public CachedUniform<float> {
public:
	using CachedUniform::CachedUniform;

	void set(GLuint programID) const override;
};

class UniformR2 : public CachedUniform<vec2> {
public:
	using CachedUniform::CachedUniform;
	void set(GLuint programID) const override;
};

class UniformR3 : public CachedUniform<vec3> {
public:
	using CachedUniform::CachedUniform;
	void set(GLuint programID) const override;
};

class UniformR4 : public CachedUniform<vec4> {
public:
	using CachedUniform::CachedUniform;
	void set(GLuint programID) const override;
};

class UniformMat2 : public CachedUniform<mat2> {
public:
	using CachedUniform::CachedUniform;
	void set(GLuint programID) const override;
};

class UniformMat3 : public CachedUniform<mat3> {
	public:
	using CachedUniform::CachedUniform;
	void set(GLuint programID) const override;
};

class UniformMat4 : public CachedUniform<mat4> {
public:
	using CachedUniform::CachedUniform;
	void set(GLuint programID) const override;
};

class UniformInt : public CachedUniform<int> {
public:
	using CachedUniform::CachedUniform;
	void set(GLuint programID) const override;
};




template <typename T>
CachedUniform<T>::CachedUniform(const string& name, const T& value, shared_ptr<std::function<T(float)>> updater): Uniform(name), value(value), updater(updater) {}

template <typename T>
void CachedUniform<T>::update(float time, float dt) {
	if (updater != nullptr)
		value = (*updater)(time);
}

template <typename T>
T CachedUniform<T>::getValue() const {
	return value;
}

template <typename T>
void CachedUniform<T>::setValue(const T& newValue) {
	value = newValue;
}

