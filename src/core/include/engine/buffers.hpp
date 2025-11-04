#pragma once
#include "renderingUtils.hpp"

class AttributeBuffer {
public:
	string name;
	GLuint bufferAddress;
	GLSLType type;
	int inputNumber;

	AttributeBuffer(const string& name, GLSLType type, int inputNumber);
	~AttributeBuffer();

	void enable() const;
	void disable() const;
	void load(raw_data_ptr firstElementAdress, byte_size bufferSize) const;
	void update(raw_data_ptr firstElementAdress, byte_size bufferSize) const;
};

class ElementBuffer {
	GLuint bufferAddress;

public:
	ElementBuffer();
	~ElementBuffer();

	void bind() const;
	void unbind() const;
	void load(raw_data_ptr firstElementAdress, byte_size bufferSize) const;
	GLuint getAddress() const;
};

class VertexArray {
	GLuint vaoAddress = 0;
public:
	VertexArray() {
		glGenVertexArrays(1, &vaoAddress);
		glBindVertexArray(vaoAddress);
	}
	~VertexArray() {
		glDeleteVertexArrays(1, &vaoAddress);
	}
	void bind() const {
		glBindVertexArray(vaoAddress);
	}
	void unbind() const {
		glBindVertexArray(0);
	}
	void addElementBuffer(const ElementBuffer& elementBuffer) const {
		bind();
		elementBuffer.bind();
	}
	void addAttributeBuffer(const AttributeBuffer& attributeBuffer) const {
		bind();
		attributeBuffer.enable();
	}
};