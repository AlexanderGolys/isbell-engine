#pragma once
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include "exceptions.hpp"



enum class GLSLPrimitive {
	BYTE, BOOL,
	FLOAT, DOUBLE,
	SHORT, INT, UINT,
	IVEC2, IVEC3, IVEC4,
	VEC2, VEC3, VEC4,
	DVEC2, DVEC3, DVEC4,
	MAT2, MAT3, MAT4,
	MAT2x3, MAT2x4, MAT3x2, MAT3x4, MAT4x2, MAT4x3,
	SAMPLER1D, SAMPLER2D, SAMPLER3D,
	UNKNOWN
};

enum class GLSLDataTypeShape {
	SCALAR,
	VECTOR,
	MATRIX,
};

struct VertexBufferLayout {
	vector<GLSLPrimitive> types;
	byte_size stride;
	vector<byte_size> offsets;

	explicit VertexBufferLayout(const vector<GLSLPrimitive>& types);
	static VertexBufferLayout singleAttribute(GLSLPrimitive type);

	bool operator==(const VertexBufferLayout& other) const;
};



template<typename T>
constexpr GLSLPrimitive glslTypeVariable = GLSLPrimitive::UNKNOWN;

template<> constexpr GLSLPrimitive glslTypeVariable<bool> = GLSLPrimitive::BOOL;
template<> constexpr GLSLPrimitive glslTypeVariable<char> = GLSLPrimitive::BYTE;
template<> constexpr GLSLPrimitive glslTypeVariable<short> = GLSLPrimitive::SHORT;
template<> constexpr GLSLPrimitive glslTypeVariable<float> = GLSLPrimitive::FLOAT;
template<> constexpr GLSLPrimitive glslTypeVariable<double> = GLSLPrimitive::DOUBLE;
template<> constexpr GLSLPrimitive glslTypeVariable<int> = GLSLPrimitive::INT;
template<> constexpr GLSLPrimitive glslTypeVariable<uint> = GLSLPrimitive::UINT;
template<> constexpr GLSLPrimitive glslTypeVariable<array_len> = GLSLPrimitive::UINT;
template<> constexpr GLSLPrimitive glslTypeVariable<vec2> = GLSLPrimitive::VEC2;
template<> constexpr GLSLPrimitive glslTypeVariable<vec3> = GLSLPrimitive::VEC3;
template<> constexpr GLSLPrimitive glslTypeVariable<vec4> = GLSLPrimitive::VEC4;
template<> constexpr GLSLPrimitive glslTypeVariable<dvec2> = GLSLPrimitive::DVEC2;
template<> constexpr GLSLPrimitive glslTypeVariable<dvec3> = GLSLPrimitive::DVEC3;
template<> constexpr GLSLPrimitive glslTypeVariable<dvec4> = GLSLPrimitive::DVEC4;
template<> constexpr GLSLPrimitive glslTypeVariable<ivec2> = GLSLPrimitive::IVEC2;
template<> constexpr GLSLPrimitive glslTypeVariable<ivec3> = GLSLPrimitive::IVEC3;
template<> constexpr GLSLPrimitive glslTypeVariable<ivec4> = GLSLPrimitive::IVEC4;
template<> constexpr GLSLPrimitive glslTypeVariable<mat2> = GLSLPrimitive::MAT2;
template<> constexpr GLSLPrimitive glslTypeVariable<mat3> = GLSLPrimitive::MAT3;
template<> constexpr GLSLPrimitive glslTypeVariable<mat4> = GLSLPrimitive::MAT4;
template<> constexpr GLSLPrimitive glslTypeVariable<mat2x3> = GLSLPrimitive::MAT2x3;
template<> constexpr GLSLPrimitive glslTypeVariable<mat2x4> = GLSLPrimitive::MAT2x4;
template<> constexpr GLSLPrimitive glslTypeVariable<mat3x2> = GLSLPrimitive::MAT3x2;
template<> constexpr GLSLPrimitive glslTypeVariable<mat3x4> = GLSLPrimitive::MAT3x4;
template<> constexpr GLSLPrimitive glslTypeVariable<mat4x2> = GLSLPrimitive::MAT4x2;
template<> constexpr GLSLPrimitive glslTypeVariable<mat4x3> = GLSLPrimitive::MAT4x3;

template<typename T>
constexpr GLSLDataTypeShape typeShape = GLSLDataTypeShape::SCALAR;

template<> constexpr GLSLDataTypeShape typeShape<vec2> = GLSLDataTypeShape::VECTOR;
template<> constexpr GLSLDataTypeShape typeShape<vec3> = GLSLDataTypeShape::VECTOR;
template<> constexpr GLSLDataTypeShape typeShape<vec4> = GLSLDataTypeShape::VECTOR;
template<> constexpr GLSLDataTypeShape typeShape<ivec2> = GLSLDataTypeShape::VECTOR;
template<> constexpr GLSLDataTypeShape typeShape<ivec3> = GLSLDataTypeShape::VECTOR;
template<> constexpr GLSLDataTypeShape typeShape<ivec4> = GLSLDataTypeShape::VECTOR;
template<> constexpr GLSLDataTypeShape typeShape<mat2> = GLSLDataTypeShape::MATRIX;
template<> constexpr GLSLDataTypeShape typeShape<mat3> = GLSLDataTypeShape::MATRIX;
template<> constexpr GLSLDataTypeShape typeShape<mat4> = GLSLDataTypeShape::MATRIX;
template<> constexpr GLSLDataTypeShape typeShape<mat2x3> = GLSLDataTypeShape::MATRIX;
template<> constexpr GLSLDataTypeShape typeShape<mat2x4> = GLSLDataTypeShape::MATRIX;
template<> constexpr GLSLDataTypeShape typeShape<mat3x2> = GLSLDataTypeShape::MATRIX;
template<> constexpr GLSLDataTypeShape typeShape<mat3x4> = GLSLDataTypeShape::MATRIX;
template<> constexpr GLSLDataTypeShape typeShape<mat4x2> = GLSLDataTypeShape::MATRIX;
template<> constexpr GLSLDataTypeShape typeShape<mat4x3> = GLSLDataTypeShape::MATRIX;


inline int lengthOfGLSLType(GLSLPrimitive type) {
	switch (type) {
	case GLSLPrimitive::IVEC2:
	case GLSLPrimitive::VEC2:
	case GLSLPrimitive::DVEC2:
		return 2;
	case GLSLPrimitive::VEC3:
	case GLSLPrimitive::IVEC3:
	case GLSLPrimitive::DVEC3:
		return 3;
	case GLSLPrimitive::VEC4:
	case GLSLPrimitive::IVEC4:
	case GLSLPrimitive::MAT2:
	case GLSLPrimitive::DVEC4:
		return 4;
	case GLSLPrimitive::MAT2x3:
	case GLSLPrimitive::MAT3x2:
		return 6;
	case GLSLPrimitive::MAT2x4:
	case GLSLPrimitive::MAT4x2:
		return 8;
	case GLSLPrimitive::MAT3:
		return 9;
	case GLSLPrimitive::MAT3x4:
	case GLSLPrimitive::MAT4x3:
		return 12;
	case GLSLPrimitive::MAT4:
		return 16;
	case GLSLPrimitive::UNKNOWN:
		THROW(ValueError, "Unknown GLSLType has no length");
	}
	return 1;
}

inline GLenum primitiveTypeEnum(GLSLPrimitive type) {
	switch (type) {
	case GLSLPrimitive::INT:
	case GLSLPrimitive::IVEC2:
	case GLSLPrimitive::IVEC3:
	case GLSLPrimitive::IVEC4:
		return GL_INT;
	case GLSLPrimitive::UINT:
	case GLSLPrimitive::SAMPLER1D:
	case GLSLPrimitive::SAMPLER2D:
	case GLSLPrimitive::SAMPLER3D:
		return GL_UNSIGNED_INT;
	case GLSLPrimitive::BYTE:
		return GL_BYTE;
	case GLSLPrimitive::BOOL:
		return GL_BOOL;
	case GLSLPrimitive::SHORT:
		return GL_SHORT;
	case GLSLPrimitive::DOUBLE:
	case GLSLPrimitive::DVEC2:
	case GLSLPrimitive::DVEC3:
	case GLSLPrimitive::DVEC4:
		return GL_DOUBLE;
	case GLSLPrimitive::UNKNOWN:
		THROW(ValueError, "Unknown GLSLType has no primitive type");
	}
	return GL_FLOAT;
}

inline byte_size sizeOfGLSLType(GLSLPrimitive type) {
	byte_size SINGLE_PREC = 4;
	byte_size DOUBLE_PREC = 8;
	byte_size ONE_BYTE = 1;

	switch (type) {
	case GLSLPrimitive::FLOAT:
	case GLSLPrimitive::INT:
	case GLSLPrimitive::UINT:
	case GLSLPrimitive::SAMPLER1D:
	case GLSLPrimitive::SAMPLER2D:
	case GLSLPrimitive::SAMPLER3D:
		return SINGLE_PREC;
	case GLSLPrimitive::BYTE:
	case GLSLPrimitive::BOOL:
		return ONE_BYTE;
	case GLSLPrimitive::SHORT:
		return ONE_BYTE * 2;
	case GLSLPrimitive::DOUBLE:
		return DOUBLE_PREC;
	case GLSLPrimitive::DVEC2:
		return DOUBLE_PREC * 2;
	case GLSLPrimitive::DVEC3:
		return DOUBLE_PREC * 3;
	case GLSLPrimitive::DVEC4:
		return DOUBLE_PREC * 4;
	case GLSLPrimitive::VEC2:
	case GLSLPrimitive::IVEC2:
		return SINGLE_PREC * 2;
	case GLSLPrimitive::VEC3:
	case GLSLPrimitive::IVEC3:
		return SINGLE_PREC * 3;
	case GLSLPrimitive::VEC4:
	case GLSLPrimitive::IVEC4:
		return SINGLE_PREC * 4;
	case GLSLPrimitive::MAT2:
		return SINGLE_PREC * 2 * 2;
	case GLSLPrimitive::MAT3:
		return SINGLE_PREC * 3 * 3;
	case GLSLPrimitive::MAT4:
		return SINGLE_PREC * 4 * 4;
	case GLSLPrimitive::MAT2x3:
	case GLSLPrimitive::MAT3x2:
		return SINGLE_PREC * 2 * 3;
	case GLSLPrimitive::MAT2x4:
	case GLSLPrimitive::MAT4x2:
		return SINGLE_PREC * 2 * 4;
	case GLSLPrimitive::MAT3x4:
	case GLSLPrimitive::MAT4x3:
		return SINGLE_PREC * 3 * 4;
	}
	THROW(ValueError, "Unknown GLSLType");
}

inline bool supportedAttributeType(GLSLPrimitive type) {
	return	type == GLSLPrimitive::FLOAT or
			type == GLSLPrimitive::VEC2 or
			type == GLSLPrimitive::VEC3 or
			type == GLSLPrimitive::VEC4 or
			type == GLSLPrimitive::DOUBLE or
			type == GLSLPrimitive::DVEC2 or
			type == GLSLPrimitive::DVEC3 or
			type == GLSLPrimitive::DVEC4;
}

constexpr bool supportedUniformType(GLSLPrimitive type) {
	return
		type == GLSLPrimitive::BOOL or
		type == GLSLPrimitive::FLOAT or
		type == GLSLPrimitive::INT or
		type == GLSLPrimitive::UINT or
		type == GLSLPrimitive::IVEC2 or
		type == GLSLPrimitive::IVEC3 or
		type == GLSLPrimitive::IVEC4 or
		type == GLSLPrimitive::VEC2 or
		type == GLSLPrimitive::VEC3 or
		type == GLSLPrimitive::VEC4 or
		type == GLSLPrimitive::MAT2 or
		type == GLSLPrimitive::MAT3 or
		type == GLSLPrimitive::MAT4 or
		type == GLSLPrimitive::MAT2x3 or
		type == GLSLPrimitive::MAT2x4 or
		type == GLSLPrimitive::MAT3x2 or
		type == GLSLPrimitive::MAT3x4 or
		type == GLSLPrimitive::MAT4x2 or
		type == GLSLPrimitive::MAT4x3;
}

template<typename T>
concept supported_uniform_type = requires {
	supportedUniformType(glslTypeVariable<T>);
};

