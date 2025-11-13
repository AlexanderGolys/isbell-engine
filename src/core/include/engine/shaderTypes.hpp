#pragma once
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include "exceptions.hpp"



enum GLSLType {
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

struct VertexBufferLayout {
	vector<GLSLType> types;
	byte_size stride;
	vector<byte_size> offsets;

	explicit VertexBufferLayout(const vector<GLSLType>& types);
	static VertexBufferLayout singleAttribute(GLSLType type);

	bool operator==(const VertexBufferLayout& other) const;
};

template<typename T>
constexpr GLSLType glslTypeVariable = UNKNOWN;

template<> constexpr GLSLType glslTypeVariable<bool> = BOOL;
template<> constexpr GLSLType glslTypeVariable<char> = BYTE;
template<> constexpr GLSLType glslTypeVariable<short> = SHORT;
template<> constexpr GLSLType glslTypeVariable<float> = FLOAT;
template<> constexpr GLSLType glslTypeVariable<double> = DOUBLE;
template<> constexpr GLSLType glslTypeVariable<int> = INT;
template<> constexpr GLSLType glslTypeVariable<unsigned int> = UINT;
template<> constexpr GLSLType glslTypeVariable<array_len> = UINT;
template<> constexpr GLSLType glslTypeVariable<vec2> = VEC2;
template<> constexpr GLSLType glslTypeVariable<vec3> = VEC3;
template<> constexpr GLSLType glslTypeVariable<vec4> = VEC4;
template<> constexpr GLSLType glslTypeVariable<dvec2> = DVEC2;
template<> constexpr GLSLType glslTypeVariable<dvec3> = DVEC3;
template<> constexpr GLSLType glslTypeVariable<dvec4> = DVEC4;
template<> constexpr GLSLType glslTypeVariable<ivec2> = IVEC2;
template<> constexpr GLSLType glslTypeVariable<ivec3> = IVEC3;
template<> constexpr GLSLType glslTypeVariable<ivec4> = IVEC4;
template<> constexpr GLSLType glslTypeVariable<mat2> = MAT2;
template<> constexpr GLSLType glslTypeVariable<mat3> = MAT3;
template<> constexpr GLSLType glslTypeVariable<mat4> = MAT4;
template<> constexpr GLSLType glslTypeVariable<mat2x3> = MAT2x3;
template<> constexpr GLSLType glslTypeVariable<mat2x4> = MAT2x4;
template<> constexpr GLSLType glslTypeVariable<mat3x2> = MAT3x2;
template<> constexpr GLSLType glslTypeVariable<mat3x4> = MAT3x4;
template<> constexpr GLSLType glslTypeVariable<mat4x2> = MAT4x2;
template<> constexpr GLSLType glslTypeVariable<mat4x3> = MAT4x3;

inline int lengthOfGLSLType(GLSLType type) {
	switch (type) {
	case IVEC2:
	case VEC2:
	case DVEC2:
		return 2;
	case VEC3:
	case IVEC3:
	case DVEC3:
		return 3;
	case VEC4:
	case IVEC4:
	case MAT2:
	case DVEC4:
		return 4;
	case MAT2x3:
	case MAT3x2:
		return 6;
	case MAT2x4:
	case MAT4x2:
		return 8;
	case MAT3:
		return 9;
	case MAT3x4:
	case MAT4x3:
		return 12;
	case MAT4:
		return 16;
	case UNKNOWN:
		THROW(ValueError, "Unknown GLSLType has no length");
	}
	return 1;
}

inline GLenum primitiveTypeEnum(GLSLType type) {
	switch (type) {
	case INT:
	case IVEC2:
	case IVEC3:
	case IVEC4:
		return GL_INT;
	case UINT:
	case SAMPLER1D:
	case SAMPLER2D:
	case SAMPLER3D:
		return GL_UNSIGNED_INT;
	case BYTE:
		return GL_BYTE;
	case BOOL:
		return GL_BOOL;
	case SHORT:
		return GL_SHORT;
	case DOUBLE:
	case DVEC2:
	case DVEC3:
	case DVEC4:
		return GL_DOUBLE;
	case UNKNOWN:
		THROW(ValueError, "Unknown GLSLType has no primitive type");
	}
	return GL_FLOAT;
}

inline byte_size sizeOfGLSLType(GLSLType type) {
	byte_size SINGLE_PREC = 4;
	byte_size DOUBLE_PREC = 8;
	byte_size ONE_BYTE = 1;

	switch (type) {
	case FLOAT:
	case INT:
	case UINT:
	case SAMPLER1D:
	case SAMPLER2D:
	case SAMPLER3D:
		return SINGLE_PREC;
	case BYTE:
	case BOOL:
		return ONE_BYTE;
	case SHORT:
		return ONE_BYTE * 2;
	case DOUBLE:
		return DOUBLE_PREC;
	case DVEC2:
		return DOUBLE_PREC * 2;
	case DVEC3:
		return DOUBLE_PREC * 3;
	case DVEC4:
		return DOUBLE_PREC * 4;
	case VEC2:
	case IVEC2:
		return SINGLE_PREC * 2;
	case VEC3:
	case IVEC3:
		return SINGLE_PREC * 3;
	case VEC4:
	case IVEC4:
		return SINGLE_PREC * 4;
	case MAT2:
		return SINGLE_PREC * 2 * 2;
	case MAT3:
		return SINGLE_PREC * 3 * 3;
	case MAT4:
		return SINGLE_PREC * 4 * 4;
	case MAT2x3:
	case MAT3x2:
		return SINGLE_PREC * 2 * 3;
	case MAT2x4:
	case MAT4x2:
		return SINGLE_PREC * 2 * 4;
	case MAT3x4:
	case MAT4x3:
		return SINGLE_PREC * 3 * 4;
	}
	THROW(ValueError, "Unknown GLSLType");
}

inline bool supportedAttributeType(GLSLType type) {
	return type == FLOAT or type == VEC2 or type == VEC3 or type == VEC4
		or type == DOUBLE or type == DVEC2 or type == DVEC3 or type == DVEC4;
}

inline bool supportedUniformType(GLSLType type) {
	return
		type == BOOL or
		type == FLOAT or
		type == INT or
		type == UINT or
		type == IVEC2 or
		type == IVEC3 or
		type == IVEC4 or
		type == VEC2 or
		type == VEC3 or
		type == VEC4 or
		type == MAT2 or
		type == MAT3 or
		type == MAT4 or
		type == MAT2x3 or
		type == MAT2x4 or
		type == MAT3x2 or
		type == MAT3x4 or
		type == MAT4x2 or
		type == MAT4x3;
}