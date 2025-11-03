#pragma once

#include "exceptions.hpp"
#include "shaders.hpp"
#include <GL/glew.h>
#include "shaderDataTypes.hpp"

class Uniform {
	string name;

	GLint getLocation(const ShaderProgram& prog) const;
	void setInt(const ShaderProgram& prog) const;
	void setUInt(const ShaderProgram& prog) const;
	void setFloat(const ShaderProgram& prog) const;
	void setMat(const ShaderProgram& prog, ShaderDataType matType) const;
	vs_dim getComponentLength() const;

public:
	explicit Uniform(const string& name) : name(name) {}
	virtual ~Uniform() = default;

	virtual void update(float time, float dt) = 0;
	virtual raw_data_ptr getDataPtr() const = 0;
	virtual ShaderDataType getType() const = 0;
	virtual array_len getArrayLength() const;

	void set(const ShaderProgram& prog) const;
	string getName() const;
};

template<typename T>
class SingleUniform : public Uniform {
	T data;
	ShaderDataType type;
	shared_ptr<BIHOM(float, float, T)> updater;
public:
	SingleUniform(const string& name, ShaderDataType type, const shared_ptr<BIHOM(float, float, T)>& updater);
	SingleUniform(const string& name, ShaderDataType type, const shared_ptr<HOM(float, T)>& updater);

	void update(float time, float dt) override;
	raw_data_ptr getDataPtr() const override;
	ShaderDataType getType() const override;
	void setData(const T& newData) { data = newData; }
};














template <typename T>
SingleUniform<T>::SingleUniform(const string& name, ShaderDataType type, const shared_ptr<std::function<T(float, float)>>& updater): Uniform(name), data((*updater)(0.0f, 0.0f)), type(type), updater(updater) {}

template <typename T>
SingleUniform<T>::SingleUniform(const string& name, ShaderDataType type, const shared_ptr<std::function<T(float)>>& updater): Uniform(name), data((*updater)(0.0f)), type(type), updater(nullptr) {
	updater = make_shared<BIHOM(float, float, T)>([updater](float time, float){ return (*updater)(time); });
}

template <typename T>
void SingleUniform<T>::update(float time, float dt) {
	data = (*updater)(time, dt);
}

template <typename T>
raw_data_ptr SingleUniform<T>::getDataPtr() const {
	ShaderDataTypeShape shape = dataTypeShape(type);
	switch (shape) {
	case SCALAR:
		return static_cast<raw_data_ptr>(&data);
	case VECTOR:
		return static_cast<raw_data_ptr>(&data[0]);
	case MATRIX:
		return static_cast<raw_data_ptr>(&data[0][0]);
	default:
		THROW(UnknownVariantError, "Data shape not recognized for uniform data retrieval");
	}
}

template <typename T>
ShaderDataType SingleUniform<T>::getType() const {
	return type;
}
