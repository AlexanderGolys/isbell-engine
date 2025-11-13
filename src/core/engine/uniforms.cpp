#include "uniforms.hpp"


UniformComponent::UniformComponent(const string& name)
: name(name) {
	if (not name.starts_with("u_"))
		LOG_WARN("Uniform name '" + name + "' does not start with 'u_' prefix");
}


void TimeUniform::setUniformGL(GLint location) const {
	GLCommand::setUniform(location, value);
}


void UniformComponent::setDuringRender() const {
	setUniformGL(GLCommand::getUniformLocation(get_name()));
}
