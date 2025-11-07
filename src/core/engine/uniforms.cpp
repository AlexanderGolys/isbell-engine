#include "uniforms.hpp"


string UniformComponent::getName() const {
	return name;
}

void UniformComponent::setDuringRender() const {
	set(GLCommand::getUniformLocation(getName()));
}
