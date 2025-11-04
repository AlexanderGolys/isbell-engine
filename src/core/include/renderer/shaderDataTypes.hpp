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

enum ShaderDataTypeShape {
	SCALAR = 1,
	VECTOR = 2,
	MATRIX = 3,
};

inline string to_string(ShaderDataType type) {
	switch (type) {
	case BOOL:    return "BOOL";
	case BYTE:    return "BYTE";
	case INT:     return "INT";
	case UINT:    return "UINT";
	case FLOAT:   return "FLOAT";
	case VEC2:    return "VEC2";
	case VEC3:    return "VEC3";
	case VEC4:    return "VEC4";
	case IVEC2:   return "IVEC2";
	case IVEC3:   return "IVEC3";
	case IVEC4:   return "IVEC4";
	case MAT2:    return "MAT2";
	case MAT3:    return "MAT3";
	case MAT4:    return "MAT4";
	case MAT2x3:  return "MAT2x3";
	case MAT2x4:  return "MAT2x4";
	case MAT3x2:  return "MAT3x2";
	case MAT3x4:  return "MAT3x4";
	case MAT4x2:  return "MAT4x2";
	case MAT4x3:  return "MAT4x3";
	default:
		return "UNKNOWN_DATA_TYPE";
	}
}


inline byte_size shaderDataTypeByteSize(ShaderDataType type) {
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

inline vs_dim shaderDataTypeComponentLength(ShaderDataType type) {
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

inline ShaderDataTypeShape dataTypeShape(ShaderDataType type) {
	switch (type) {
	case BOOL:
	case BYTE:
	case INT:
	case UINT:
	case FLOAT:
		return SCALAR;

	case IVEC2:
	case VEC2:
	case IVEC3:
	case VEC3:
	case IVEC4:
	case VEC4:
		return VECTOR;

	case MAT2:
	case MAT3:
	case MAT4:
	case MAT2x3:
	case MAT2x4:
	case MAT3x2:
	case MAT3x4:
	case MAT4x2:
	case MAT4x3:
		return MATRIX;

	default:
		THROW(UnknownVariantError, "ShaderDataType not recognized");
	}
}