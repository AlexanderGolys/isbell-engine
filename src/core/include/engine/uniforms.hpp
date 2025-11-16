#pragma once
#include "glCommand.hpp"
#include "renderLayers.hpp"
#include "shaderTypes.hpp"
#include "GL/glew.h"


template<typename T>
concept uniform_struct = requires(T a) {
	{ a.data() } -> std::same_as<raw_data_ptr>;
	{ a.byteSize() } -> std::same_as<byte_size>;
};

template<typename T>
concept uniform_struct_dirty = uniform_struct<T> && requires(T a) {
	{ a.isDirty() } -> std::same_as<bool>;
	{ a.markClean() } -> std::same_as<void>;
	{ a.markDirty() } -> std::same_as<void>;
};

template<uniform_struct S>
struct arrayStruct {
	DIRTY_FLAG_STRUCT
	vector<S> dataArray;

	explicit arrayStruct(const vector<S>& dataArray) : dataArray(dataArray) {}
	raw_data_ptr data() const { return dataArray[0].data(); }
	byte_size byteSize() const { return dataArray.size() * dataArray[0].byteSize(); }
};


template<uniform_struct_dirty S>
class UniformBlockComponent : public LayerComponent {
	UniformBuffer uniformBuffer;
	S value;
public:
	UniformBlockComponent(const S& value, uint bindingPoint) : uniformBuffer(bindingPoint, value.byteSize()), value(value) {}
	void init() override { uniformBuffer.load(value.data()); }
	void setDuringRender() final;
	void setValue(const S& newValue) {
		value = newValue;
	}
};

class UniformComponent : public LayerComponent {
	CONST_PROPERTY(string, name);
	virtual void setUniformGL(gl_uniform_loc location) const = 0;

public:
	explicit UniformComponent(const string& name);
	void init() override {}
	void setDuringRender() final;
};

template<supported_uniform_type T>
class PrimitiveUniform : public UniformComponent {
	PROPERTY(T, value);
	void setUniformGL(gl_uniform_loc location) const final;
public:
	PrimitiveUniform(const string& name, const T& value) : UniformComponent(name), value(value) {}
};

template<supported_uniform_type T>
class ConstPrimitiveUniform : public PrimitiveUniform<T> {
public:
	using PrimitiveUniform<T>::PrimitiveUniform;
	void update(TimeStep) final {}
};

template<supported_uniform_type T>
class DynamicPrimitiveUniform : public PrimitiveUniform<T> {
	HOM(TimeStep, T) valueFunction;
public:
	DynamicPrimitiveUniform(const string& name, const HOM(TimeStep, T)& valueFunction)
	: PrimitiveUniform<T>(name, valueFunction(TimeStep())), valueFunction(valueFunction) {}
	void update(TimeStep t) final { this->set_value(valueFunction(t)); }
};

template<supported_uniform_type T>
class ArrayUniform : public UniformComponent {
	VECT_PROPERTY(T, value);
	void setUniformGL(gl_uniform_loc location) const final;

public:
	explicit ArrayUniform(const string& name) : UniformComponent(name) {}
};

template<supported_uniform_type T>
class ConstArrayUniform : public ArrayUniform<T> {
public:
	ConstArrayUniform(const string& name, const vector<T>& value) : ArrayUniform<T>(name) { this->set_value(value); }
	void update(TimeStep) final {}
};

template<supported_uniform_type T>
class DynamicArrayUniform : public ArrayUniform<T> {
	CONST_PROPERTY(BIHOM(TimeStep, array_index, T), valueFunction);
	CONST_PROPERTY(array_len, arraySize);
public:
	DynamicArrayUniform(const string& name, array_len size, const BIHOM(TimeStep, array_index, T)& valueFunction);
	void update(TimeStep t) final;
};




// ------------------- Implementation ------------------ //


template <uniform_struct_dirty S>
void UniformBlockComponent<S>::setDuringRender() {
	uniformBuffer.bind();
	if (value.isDirty()) {
		uniformBuffer.update(value.data());
		value.markClean();
	}
}

template <supported_uniform_type T>
void PrimitiveUniform<T>::setUniformGL(gl_uniform_loc location) const {
    if constexpr (typeShape<T> == GLSLDataTypeShape::SCALAR)
    	return GLCommand::setUniform(location, value);
	if constexpr (typeShape<T> == GLSLDataTypeShape::VECTOR)
		return GLCommand::setUniform(location, &value[0], glslTypeVariable<T>);
	if constexpr (typeShape<T> == GLSLDataTypeShape::MATRIX)
		return GLCommand::setUniform(location, &value[0][0], glslTypeVariable<T>);
	THROW(UnknownVariantError, "Unsupported GLSL type in ConstantUniform");
}

template <supported_uniform_type T>
void ArrayUniform<T>::setUniformGL(GLint location) const {
	if constexpr (typeShape<T> == GLSLDataTypeShape::SCALAR)
		return GLCommand::setUniform(location, &value[0], glslTypeVariable<T>, value.size());
	if constexpr (typeShape<T> == GLSLDataTypeShape::VECTOR)
		return GLCommand::setUniform(location, &value[0][0], glslTypeVariable<T>, value.size());
	if constexpr (typeShape<T> == GLSLDataTypeShape::MATRIX)
		return GLCommand::setUniform(location, &value[0][0][0], glslTypeVariable<T>, value.size());
	THROW(UnknownVariantError, "Unsupported GLSL type in ConstantUniform");
}

template <supported_uniform_type T>
DynamicArrayUniform<T>::DynamicArrayUniform(const string& name, array_len size, const std::function<T(TimeStep, array_index)>& valueFunction)
: ArrayUniform<T>(name), valueFunction(valueFunction), arraySize(size) {
	this->reserve_value(size);
	for (array_index i = 0; i < size; i++)
		this->append_value(valueFunction(TimeStep(), i));
}

template <supported_uniform_type T>
void DynamicArrayUniform<T>::update(TimeStep t) {
	for (array_index i = 0; i < arraySize; i++)
		this->set_value(i, valueFunction(t, i));
}


