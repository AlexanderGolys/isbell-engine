#pragma once
#include "exceptions.hpp"


enum ShaderDataType {
	BOOL,
	BYTE,
	INT,
	UINT,
	FLOAT,
	VEC2,
	VEC3,
	VEC4,
	IVEC2,
	IVEC3,
	IVEC4,
	MAT2,
	MAT3,
	MAT4,
	MAT2x3,
	MAT2x4,
	MAT3x2,
	MAT3x4,
	MAT4x2,
	MAT4x3,
};




inline size_t shaderDataTypeByteSize(ShaderDataType type) {
	switch (type) {
	case BOOL:
	case BYTE:
		return 1;

	case INT:
	case UINT:
	case FLOAT:
		return 4;

	case IVEC2:
	case VEC2:
		return 2 * 4;

	case IVEC3:
	case VEC3:
		return 3 * 4;

	case IVEC4:
	case VEC4:
		return 4 * 4;

	case MAT2:
		return 2 * 2 * 4;
	case MAT3:
		return 3 * 3 * 4;
	case MAT4:
		return 4 * 4 * 4;
	case MAT2x3:
	case MAT3x2:
		return 2 * 3 * 4;
	case MAT2x4:
	case MAT4x2:
		return 2 * 4 * 4;
	case MAT3x4:
	case MAT4x3:
		return 3 * 4 * 4;

	default:
		THROW(UnknownVariantError, "ShaderDataType not recognized");
	}
}

inline unsigned int shaderDataTypeComponentLength(ShaderDataType type) {
	switch (type) {
	case BOOL:
	case BYTE:
	case INT:
	case UINT:
	case FLOAT:
		return 1;

	case IVEC2:
	case VEC2:
		return 2;

	case IVEC3:
	case VEC3:
		return 3;

	case IVEC4:
	case VEC4:
		return 4;

	case MAT2:
		return 2 * 2;
	case MAT3:
		return 3 * 3;
	case MAT4:
		return 4 * 4;
	case MAT2x3:
	case MAT3x2:
		return 2 * 3;
	case MAT2x4:
	case MAT4x2:
		return 2 * 4;
	case MAT3x4:
	case MAT4x3:
		return 3 * 4;

	default:
		THROW(UnknownVariantError, "ShaderDataType not recognized");
	}
}
