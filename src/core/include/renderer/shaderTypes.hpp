#pragma once
#include "func.hpp"

enum AttributeType {
	FLOAT,
	VEC2,
	VEC3,
	VEC4,
	INT,
	UINT,
	IVEC2,
	IVEC3,
	IVEC4,
	BYTE,
	UBYTE,
	SHORT,
	USHORT
};


size_t sizeOfAttributeType(AttributeType type);
int lengthOfAttributeType(AttributeType type);

enum UniformType {
	U_FLOAT,
	U_INT,
	U_UINT,
	U_VEC2,
	U_VEC3,
	U_VEC4,
	U_MAT2,
	U_MAT3,
	U_MAT4,
	U_MAT2x3,
	U_MAT3x2,
	U_MAT2x4,
	U_MAT4x2,
	U_MAT3x4,
	U_MAT4x3,
};

enum ShaderType {
	VERTEX_SHADER,
	GEOMETRY_SHADER,
	FRAGMENT_SHADER,
	COMPUTE_SHADER
};
