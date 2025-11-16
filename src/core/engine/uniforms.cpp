#include "uniforms.hpp"


UniformComponent::UniformComponent(const string& name): name(name) {
	if (not name.starts_with("u_"))
		THROW(ValueError, "Uniform name must start with 'u_'. Given name: " + name);
}

void UniformComponent::setDuringRender() {
	setUniformGL(GLCommand::getUniformLocation(get_name()));
}
