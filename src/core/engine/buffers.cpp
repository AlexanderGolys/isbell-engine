#include "buffers.hpp"
#include "exceptions.hpp"


VertexAttribute::VertexAttribute(const string& name, ShaderDataType type)
: name(name), type(type) {}

unsigned int VertexAttribute::elementLength() const {
	if (type == FLOAT)
		return 1;
	if (type == VEC2)
		return 2;
	if (type == VEC3)
		return 3;
	if (type == VEC4)
		return 4;
	THROW(UnknownVariantError, "VertexAttribute has invalid type " + to_string(type));
}

size_t VertexAttribute::elementSize() const {
	return elementLength() * 4;
}

VertexBufferLayout::VertexBufferLayout(initializer_list<pair<string, ShaderDataType>> attrs)
: strideBytes(0) {
	for (const auto& [name, type] : attrs) {
		attributes.emplace_back(name, type);
		attributeOffsets.push_back(strideBytes);
		strideBytes += attributes.back().elementSize();
	}
}

const void* VertexBufferLayout::offsetAt(int i) const {
	if (i < 0)
		i += attributeOffsets.size();
	THROW_IF(i < 0 || i >= attributeOffsets.size(), IndexOutOfBounds, i, attributeOffsets.size(), "VertexBufferLayout offsetAt");
	return reinterpret_cast<const void*>(attributeOffsets[i]);
}

unsigned int VertexBufferLayout::length() const {
	return attributes.size();
}

VertexBuffer::VertexBuffer(const VertexBufferLayout& layout)
: layout(layout) {
	glCreateBuffers(1, &id);
}

VertexBuffer::~VertexBuffer() {
	glDeleteBuffers(1, &id);
}

GLuint VertexBuffer::vboID() const {
	return id;
}

const VertexBufferLayout& VertexBuffer::getLayout() const {
	return layout;
}

void VertexBuffer::bind() const {
	glBindBuffer(GL_ARRAY_BUFFER, id);
}

void VertexBuffer::unbind() const {
	glBindBuffer(GL_ARRAY_BUFFER, 0);
}

void VertexBuffer::uploadData(const void* data, unsigned int size) const {
	glNamedBufferData(id, size, data, GL_STATIC_DRAW);
}

void VertexBuffer::updateData(const void* data, unsigned int size) const {
	glNamedBufferSubData(id, 0, size, data);
}

unsigned int VertexBuffer::length() const {
	return layout.length();
}


IndexBuffer::IndexBuffer(unsigned int trianglesCount) : count(trianglesCount * 3) {
	glCreateBuffers(1, &id);
	glNamedBufferData(id, trianglesCount * sizeof(uint32_t), nullptr, GL_DYNAMIC_DRAW);
}

IndexBuffer::~IndexBuffer() {
	glDeleteBuffers(1, &id);
}

GLuint IndexBuffer::eboID() const {
	return id;
}

void IndexBuffer::bind() const {
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, id);
}

void IndexBuffer::unbind() const {
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
}

void IndexBuffer::uploadData(const vector<ivec3>& indices) {
	count = indices.size() * 3;
	bind();
	glNamedBufferData(id, count*sizeof(uint32_t), indices.data(), GL_STATIC_DRAW);
}

unsigned int IndexBuffer::getCount() const {
	return count;
}

VertexArray::VertexArray() {
	glGenVertexArrays(1, &id);
}

VertexArray::VertexArray(const shared_ptr<IndexBuffer>& indexBuffer, const vector<shared_ptr<VertexBuffer>>& vertexBuffers): VertexArray() {
	setIndexBuffer(indexBuffer);
	for (const auto& vb : vertexBuffers) {
		addVertexBuffer(vb);
	}
}

VertexArray::~VertexArray() {
	glDeleteVertexArrays(1, &id);
}

GLuint VertexArray::vaoID() const {
	return id;
}

void VertexArray::bind() const {
	glBindVertexArray(id);
}

void VertexArray::unbind() const {
	glBindVertexArray(0);
}

void VertexArray::addVertexBuffer(const shared_ptr<VertexBuffer>& vertexBuffer) {
	bind();
	vertexBuffer->bind();
	const VertexBufferLayout& layout = vertexBuffer->getLayout();
	unsigned int stride = layout.strideBytes;
	for (int i = 0; i < layout.length(); i++) {
		VertexAttribute attr = layout.attributes[i];
		glEnableVertexAttribArray(attributeCount);
		glVertexAttribPointer(
			attributeCount,
			attr.elementLength(),
			GL_FLOAT,
			GL_FALSE,
			stride,
			layout.offsetAt(i));
		attributeCount++;
	}
}

void VertexArray::setIndexBuffer(const shared_ptr<IndexBuffer>& indexBuffer) {
	bind();
	indexBuffer->bind();
	ebo = indexBuffer;
}

unsigned int VertexArray::getIndexCount() const {
	return ebo->getCount();
}

