#include "buffers.hpp"



AttributeBuffer::AttributeBuffer(const string& name, GLSLType type, int inputNumber)
: name(name), type(type), inputNumber(inputNumber), bufferAddress(0) {
	glCreateBuffers(1, &bufferAddress);
}

AttributeBuffer::~AttributeBuffer() {
	glDeleteBuffers(1, &bufferAddress);
}


void AttributeBuffer::enable() const {
	glEnableVertexAttribArray(this->inputNumber);
	glBindBuffer(GL_ARRAY_BUFFER, this->bufferAddress);
	glVertexAttribPointer(this->inputNumber, lengthOfGLSLType(this->type), GL_FLOAT, GL_FALSE, 0, (void*)0);
}

void AttributeBuffer::disable() const {
	glDisableVertexAttribArray(this->inputNumber);
}

void AttributeBuffer::load(raw_data_ptr firstElementAdress, byte_size bufferSize) const {
	glBindBuffer(GL_ARRAY_BUFFER, bufferAddress);
	glBufferData(GL_ARRAY_BUFFER, bufferSize, firstElementAdress, GL_DYNAMIC_DRAW);
}

void AttributeBuffer::update(raw_data_ptr firstElementAdress, byte_size bufferSize) const {
	glBindBuffer(GL_ARRAY_BUFFER, bufferAddress);
	glBufferSubData(GL_ARRAY_BUFFER, 0, bufferSize, firstElementAdress);
}

ElementBuffer::ElementBuffer() {
	glCreateBuffers(1, &bufferAddress);
}

ElementBuffer::~ElementBuffer() {
	glDeleteBuffers(1, &bufferAddress);
}

void ElementBuffer::bind() const {
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, bufferAddress);
}

void ElementBuffer::unbind() const {
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
}


void ElementBuffer::load(raw_data_ptr firstElementAdress, byte_size bufferSize) const {
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, bufferAddress);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, bufferSize, firstElementAdress, GL_STATIC_DRAW);
}

GLuint ElementBuffer::getAddress() const {
	return bufferAddress;
}
