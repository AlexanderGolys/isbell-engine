#pragma once
#include <GL/glew.h>
#include "exceptions.hpp"



enum GLSLType {
	BYTE,
	FLOAT,
	INT,
	SHORT,
	UINT,
	DOUBLE,
	IVEC2,
	IVEC3,
	IVEC4,
	VEC2,
	VEC3,
	VEC4,
	DVEC2,
	DVEC3,
	DVEC4,
	MAT2,
	MAT3,
	MAT4,
	MAT2x3,
	MAT2x4,
	MAT3x2,
	MAT3x4,
	MAT4x2,
	MAT4x3,
	SAMPLER1D,
	SAMPLER2D,
	SAMPLER3D
};


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
	}
	return 1;
}

inline GLenum primitiveGLSLType(GLSLType type) {
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
	case SHORT:
		return GL_SHORT;
	case DOUBLE:
	case DVEC2:
	case DVEC3:
	case DVEC4:
		return GL_DOUBLE;
	}
	return GL_FLOAT;
}

inline bool supportedAttributeType(GLSLType type) {
	return type == FLOAT or type == VEC2 or type == VEC3 or type == VEC4
		or type == DOUBLE or type == DVEC2 or type == DVEC3 or type == DVEC4;
}