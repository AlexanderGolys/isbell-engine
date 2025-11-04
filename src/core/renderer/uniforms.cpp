#include "../include/renderer/uniforms.hpp"


GLint Uniform::getLocation(const ShaderProgram& prog) const {
	return prog.getUniformLocation(name);
}

void Uniform::setInt(const ShaderProgram& prog) const {
	GLint location = getLocation(prog);
	THROW_IF(location == -1, ValueError, "Uniform '" + name + "' not found in shader program");
	vs_dim compLen = getComponentLength();
	data_ptr<GLint> dataPtr = static_cast<data_ptr<GLint>>(getDataPtr());
	array_len length = getArrayLength();
	switch (compLen) {
	case 1:
		glUniform1iv(location, length, dataPtr);
		break;
	case 2:
		glUniform2iv(location, length, dataPtr);
		break;
	case 3:
		glUniform3iv(location, length, dataPtr);
		break;
	case 4:
		glUniform4iv(location, length, dataPtr);
		break;
	default:
		THROW(UnknownVariantError, "Component length " + to_string(compLen) + " not supported for integer uniforms");
	}
}

void Uniform::setUInt(const ShaderProgram& prog) const {
	GLint location = getLocation(prog);
	THROW_IF(location == -1, ValueError, "Uniform '" + name + "' not found in shader program");
	vs_dim compLen = getComponentLength();
	data_ptr<GLuint> dataPtr = static_cast<data_ptr<GLuint>>(getDataPtr());
	array_len length = getArrayLength();
	switch (compLen) {
	case 1:
		glUniform1uiv(location, length, dataPtr);
		break;
	case 2:
		glUniform2uiv(location, length, dataPtr);
		break;
	case 3:
		glUniform3uiv(location, length, dataPtr);
		break;
	case 4:
		glUniform4uiv(location, length, dataPtr);
		break;
	default:
		THROW(UnknownVariantError, "Component length " + to_string(compLen) + " not supported for integer uniforms");
	}
}

void Uniform::setFloat(const ShaderProgram& prog) const {
	GLint location = getLocation(prog);
	THROW_IF(location == -1, ValueError, "Uniform '" + name + "' not found in shader program");
	vs_dim compLen = getComponentLength();
	data_ptr<GLfloat> dataPtr = static_cast<data_ptr<GLfloat>>(getDataPtr());
	array_len length = getArrayLength();
	switch (compLen) {
	case 1:
		glUniform1fv(location, length, dataPtr);
		break;
	case 2:
		glUniform2fv(location, length, dataPtr);
		break;
	case 3:
		glUniform3fv(location, length, dataPtr);
		break;
	case 4:
		glUniform4fv(location, length, dataPtr);
		break;
	default:
		THROW(UnknownVariantError, "Component length " + to_string(compLen) + " not supported for integer uniforms");
	}
}

void Uniform::setMat(const ShaderProgram& prog, ShaderDataType matType) const {
	GLint location = getLocation(prog);
	THROW_IF(location == -1, ValueError, "Uniform '" + name + "' not found in shader program");
	data_ptr<GLfloat> dataPtr = static_cast<data_ptr<GLfloat>>(getDataPtr());
	array_len length = getArrayLength();
	switch (matType) {
	case MAT2:
		glUniformMatrix2fv(location, length, GL_FALSE, dataPtr);
		break;
	case MAT3:
		glUniformMatrix3fv(location, length, GL_FALSE, dataPtr);
		break;
	case MAT4:
		glUniformMatrix4fv(location, length, GL_FALSE, dataPtr);
		break;
	case MAT2x3:
		glUniformMatrix2x3fv(location, length, GL_FALSE, dataPtr);
		break;
	case MAT3x2:
		glUniformMatrix3x2fv(location, length, GL_FALSE, dataPtr);
		break;
	case MAT2x4:
		glUniformMatrix2x4fv(location, length, GL_FALSE, dataPtr);
		break;
	case MAT4x2:
		glUniformMatrix4x2fv(location, length, GL_FALSE, dataPtr);
		break;
	case MAT3x4:
		glUniformMatrix3x4fv(location, length, GL_FALSE, dataPtr);
		break;
	case MAT4x3:
		glUniformMatrix4x3fv(location, length, GL_FALSE, dataPtr);
		break;
	default:
		THROW(UnknownVariantError, "ShaderDataType " + to_string(matType) + " is not a matrix type");
	}
}

array_len Uniform::getArrayLength() const {
	return 1;
}

void Uniform::set(const ShaderProgram& prog) const {
	ShaderDataType type = getType();
	switch (type) {
	case INT:
	case IVEC2:
	case IVEC3:
	case IVEC4:
		setInt(prog);
		break;
	case FLOAT:
	case VEC2:
	case VEC3:
	case VEC4:
		setFloat(prog);
		break;
	case UINT:
	case BYTE:
	case BOOL:
		setUInt(prog);
		break;
	case MAT2:
	case MAT3:
	case MAT4:
	case MAT2x3:
	case MAT2x4:
	case MAT3x2:
	case MAT3x4:
	case MAT4x2:
	case MAT4x3:
		setMat(prog, type);
		break;
	default:
		THROW(UnknownVariantError, "ShaderDataType " + to_string(type) + " not supported for uniform setting");
	}
}


vs_dim Uniform::getComponentLength() const {
	return shaderDataTypeComponentLength(getType());
}

string Uniform::getName() const {
	return name;
}

