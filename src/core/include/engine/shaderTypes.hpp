#pragma once
#include <GL/glew.h>
#include "exceptions.hpp"

enum GLSLType {
	FLOAT,
	INT,
	VEC2,
	VEC3,
	VEC4,
	MAT2,
	MAT3,
	MAT4,
	SAMPLER1D,
	SAMPLER2D,
	SAMPLER3D
};


inline size_t sizeOfGLSLType(GLSLType type) {
	switch (type) {
	case FLOAT:
		return sizeof(float);
	case INT:
		return sizeof(int);
	case VEC2:
		return sizeof(vec2);
	case VEC3:
		return sizeof(vec3);
	case VEC4:
		return sizeof(vec4);
	case MAT2:
		return sizeof(mat2);
	case MAT3:
		return sizeof(mat3);
	case MAT4:
		return sizeof(mat4);
	case SAMPLER1D:
	case SAMPLER2D:
	case SAMPLER3D:
		return sizeof(GLuint);
	}
	throw UnknownVariantError("GLSL type not recognized", __FILE__, __LINE__);
}

inline int lengthOfGLSLType(GLSLType type) {
	switch (type) {
	case FLOAT:
	case INT:
	case SAMPLER1D:
	case SAMPLER2D:
	case SAMPLER3D:
		return 1;
	case VEC2:
		return 2;
	case VEC3:
		return 3;
	case VEC4:
	case MAT2:
		return 4;
	case MAT3:
		return 9;
	case MAT4:
		return 16;
	}
	std::cout << "Error: unknown GLSLType" << std::endl;
	return -1;
}

inline GLenum primitiveGLSLType(GLSLType type) {
	switch (type) {
	case FLOAT:
	case VEC2:
	case VEC3:
	case VEC4:
	case MAT2:
	case MAT3:
	case MAT4:
		return GL_FLOAT;
	case INT:
	case SAMPLER1D:
	case SAMPLER2D:
	case SAMPLER3D:
		return GL_INT;
	}
	std::cout << "Error: unknown GLSLType" << std::endl;
	return -1;
}
