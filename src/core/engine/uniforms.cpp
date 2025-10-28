#include "uniforms.hpp"

Uniform::Uniform(const string& name): name(name) {}

string Uniform::getName() const { return name;}

void UniformFloat::set(GLuint programID) const {
	GLint location = glGetUniformLocation(programID, getName().c_str());
	glUniform1f(location, getValue());
}

void UniformR2::set(GLuint programID) const {
	GLint location = glGetUniformLocation(programID, getName().c_str());
	vec2 value = getValue();
	glUniform2f(location, value.x, value.y);
}

void UniformR3::set(GLuint programID) const {
	GLint location = glGetUniformLocation(programID, getName().c_str());
	vec3 value = getValue();
	glUniform3f(location, value.x, value.y, value.z);
}

void UniformR4::set(GLuint programID) const {
	GLint location = glGetUniformLocation(programID, getName().c_str());
	vec4 value = getValue();
	glUniform4f(location, value.x, value.y, value.z, value.w);
}

void UniformMat2::set(GLuint programID) const {
	GLint location = glGetUniformLocation(programID, getName().c_str());
	mat2 value = getValue();
	glUniformMatrix2fv(location, 1, GL_FALSE, &value[0][0]);
}

void UniformMat3::set(GLuint programID) const {
	GLint location = glGetUniformLocation(programID, getName().c_str());
	mat3 value = getValue();
	glUniformMatrix3fv(location, 1, GL_FALSE, &value[0][0]);
}

void UniformMat4::set(GLuint programID) const {
	GLint location = glGetUniformLocation(programID, getName().c_str());
	mat4 value = getValue();
	glUniformMatrix4fv(location, 1, GL_FALSE, &value[0][0]);
}

void UniformInt::set(GLuint programID) const {
	GLint location = glGetUniformLocation(programID, getName().c_str());
	glUniform1i(location, getValue());
}
