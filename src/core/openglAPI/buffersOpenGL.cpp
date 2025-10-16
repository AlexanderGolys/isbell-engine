#include "buffersOpenGL.hpp"



size_t AttributeBufferGL::size() const {
	return length * attribute.size;
}

GLenum baseType(AttributeType type) {
	switch (type) {
		case FLOAT:
		case VEC2:
		case VEC3:
		case VEC4:
			return GL_FLOAT;
		case INT:
		case IVEC2:
		case IVEC3:
		case IVEC4:
			return GL_INT;
		case BYTE:
			return GL_BYTE;
		case SHORT:
			return GL_SHORT;
		case UBYTE:
			return GL_UNSIGNED_BYTE;
		case USHORT:
			return GL_UNSIGNED_SHORT;
		case UINT:
			return GL_UNSIGNED_INT;
	}
	THROW(UnknownVariantError, "Attribute type not recognized");
}

AttributeBufferGL::AttributeBufferGL(const ShaderAttribute &attribute, int length)
: attribute(attribute), length(length) {
	glCreateBuffers(1, &vboID);
	glBindBuffer(GL_ARRAY_BUFFER, vboID);
	glBufferData(GL_ARRAY_BUFFER, length * attribute.size, nullptr, GL_DYNAMIC_DRAW);
}

AttributeBufferGL::~AttributeBufferGL() {
	glDeleteBuffers(1, &vboID);
}

void AttributeBufferGL::bind() {
	glBindBuffer(GL_ARRAY_BUFFER, vboID);
}

void AttributeBufferGL::unbind() {
	glBindBuffer(GL_ARRAY_BUFFER, 0);
}

void AttributeBufferGL::load(void *firstElementAdress) {
	bind();
	glBufferSubData(GL_ARRAY_BUFFER, 0, size(), firstElementAdress);
}

void AttributeBufferGL::load(void *firstElementAdress, int length) {
	this->length = length;
	bind();
	glBufferData(GL_ARRAY_BUFFER, size(), firstElementAdress, GL_DYNAMIC_DRAW);
}

const ShaderAttribute & AttributeBufferGL::getAttribute() const {
	return attribute;
}

int AttributeBufferGL::getLength() const {
	return length;
}



size_t IndexBufferGL::size() const {
	return length * sizeof(ivec3);
}

VertexArrayObjectGL::VertexArrayObjectGL() {
	glCreateVertexArrays(1, &vaoID);
}
VertexArrayObjectGL::~VertexArrayObjectGL() {
	glDeleteVertexArrays(1, &vaoID);
}

void VertexArrayObjectGL::bind() {
	glBindVertexArray(vaoID);
}

void VertexArrayObjectGL::unbind() {
	glBindVertexArray(0);
}

void VertexArrayObjectGL::addAttributeBuffer(const shared_ptr<AttributeBuffer> &buffer) {
	ShaderAttribute attribute = buffer->getAttribute();
	bind();
	buffer->bind();
	glEnableVertexAttribArray(attributeBuffers.size());

	if (baseType(attribute.type) == GL_FLOAT)
		glVertexAttribPointer(attributeBuffers.size(),
							  attribute.length,
							  baseType(attribute.type),
							  GL_FALSE,
							  0,
							  (const void*)0);
	else:
	glVertexAttribIPointer(attributeBuffers.size(),
						   attribute.length,
						   baseType(attribute.type),
						   0,
						   (const void*)0);

	attributeBuffers.push_back(buffer);

}

void VertexArrayObjectGL::setIndexBuffer(const shared_ptr<IndexBuffer> &buffer) {
	bind();
	buffer->bind();
	this->indexBuffer = buffer;

}

const vector<shared_ptr<AttributeBuffer>> & VertexArrayObjectGL::getAttributeBuffers() const {
	return attributeBuffers;
}

const shared_ptr<IndexBuffer> & VertexArrayObjectGL::getIndexBuffer() const {
	return indexBuffer;
}

IndexBufferGL::IndexBufferGL(int length): length(length) {
	glCreateBuffers(1, &eboID);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, eboID);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, length * sizeof(ivec3), nullptr, GL_DYNAMIC_DRAW);
}

void IndexBufferGL::bind() {
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, eboID);
}

void IndexBufferGL::unbind() {
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
}

void IndexBufferGL::load(void *firstElementAdress) {
	bind();
	glBufferSubData(GL_ELEMENT_ARRAY_BUFFER, 0, size(), firstElementAdress);
}

void IndexBufferGL::load(void *firstElementAdress, int bufferLength) {
	length = bufferLength;
	bind();
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, size(), firstElementAdress, GL_DYNAMIC_DRAW);
}

int IndexBufferGL::getLength() const {
	return length;
}




ShaderStorageBufferObjectGL::ShaderStorageBufferObjectGL(size_t bufferSize, const void *data) {
	this->bufferSize = bufferSize;
	glGenBuffers(1, &ssbo);
	glBindBuffer(GL_SHADER_STORAGE_BUFFER, ssbo);
	glBufferData(GL_SHADER_STORAGE_BUFFER, bufferSize, data, GL_DYNAMIC_COPY);
	glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 3, ssbo);
}

ShaderStorageBufferObjectGL::~ShaderStorageBufferObjectGL() {
	glDeleteBuffers(1, &ssbo);
}

void ShaderStorageBufferObjectGL::load(void *data, size_t size) {
	bufferSize = size;
	load(data);
}

void ShaderStorageBufferObjectGL::bind() {
	glBindBuffer(GL_SHADER_STORAGE_BUFFER, ssbo);
}

void ShaderStorageBufferObjectGL::unbind() {
	glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);
}

void ShaderStorageBufferObjectGL::load(void *data) {
	glBindBuffer(GL_SHADER_STORAGE_BUFFER, ssbo);
	glGetBufferSubData(GL_SHADER_STORAGE_BUFFER, 0, (GLsizeiptr)bufferSize, data);
}

size_t ShaderStorageBufferObjectGL::getSize() const {
	return bufferSize;
}
