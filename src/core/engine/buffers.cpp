#include "buffers.hpp"

#include "formatters.hpp"
#include "glCommand.hpp"
#include "logging.hpp"


AttributeBuffer::AttributeBuffer(const string& name, GLSLType type, int inputNumber)
: name(name), bufferID(0), type(type), inputNumber(inputNumber), bufferedSize(0) {
	glCreateBuffers(1, &bufferID);
}

AttributeBuffer::~AttributeBuffer() {
	glDeleteBuffers(1, &bufferID);
}


void AttributeBuffer::pointAtAttributeBuffer() const {
	GLCommand::pointAtAttributeBuffer(inputNumber, bufferID, type);
}

void AttributeBuffer::load(raw_data_ptr firstElementAdress, byte_size bufferSize) {
	if (bufferedSize != 0)
		LOG_WARN("Loading the attribute buffer again instead of updating it");
	glNamedBufferData(bufferID, bufferSize, firstElementAdress, GL_DYNAMIC_DRAW);
	bufferedSize = bufferSize;
	LOG(1, "Loaded attribute buffer " + name + " of size " + formatByteSize(bufferSize));
}

void AttributeBuffer::update(raw_data_ptr firstElementAdress, byte_size bufferSize) {
	glNamedBufferSubData(bufferID, 0, bufferSize, firstElementAdress);
	if (bufferSize != bufferedSize)
		LOG(1, "Resized attribute buffer " + name + " from " + formatByteSize(bufferedSize) + " to " + formatByteSize(bufferSize));
	bufferedSize = bufferSize;
	LOG(0, "Updated attribute buffer " + name + " with new data of size " + formatByteSize(bufferSize));

}

ElementBuffer::ElementBuffer() {
	GLCommand::createBuffer(&bufferID);
}

ElementBuffer::~ElementBuffer() {
	GLCommand::deleteBuffer(bufferID);
}


void ElementBuffer::load(raw_data_ptr firstElementAdress, byte_size bufferSize) {
	if (bufferedSize != 0)
		LOG_WARN("Loading the attribute buffer again instead of updating it");
	glNamedBufferData(bufferID, bufferSize, firstElementAdress, GL_STATIC_DRAW);
	bufferedSize = bufferSize;
	LOG(1, "Loaded element buffer with indices of size " + formatByteSize(bufferSize));
}

GLuint ElementBuffer::getID() const {
	return bufferID;
}

array_len ElementBuffer::getNumberOfIndices() const {
	return bufferedSize / sizeof(GLuint);
}

VertexArray::VertexArray() {
	GLCommand::createVAO(&vaoID);
	GLCommand::bindVAO(vaoID);
}

VertexArray::~VertexArray() {
	GLCommand::deleteVAO(vaoID);
}

void VertexArray::bind() const {
	GLCommand::bindVAO(vaoID);
}

void VertexArray::unbind() const {
	GLCommand::unbindVAO();
}

void VertexArray::addElementBuffer(sptr<ElementBuffer> elementBuffer) {
	THROW_IF(elementBuffer != nullptr, "Element buffer is already set for this VAO");
	GLCommand::bindVAO(vaoID);
	GLCommand::pointAtElementBuffer(elementBuffer->getID());
	this->elementBuffer = elementBuffer;
}

void VertexArray::addAttributeBuffer(sptr<AttributeBuffer> attributeBuffer) {
	GLCommand::bindVAO(vaoID);
	attributeBuffer->pointAtAttributeBuffer();
	attributeBuffers.push_back(attributeBuffer);
}

sptr<ElementBuffer> VertexArray::getElementBuffer() const {
	return elementBuffer;
}

vector<sptr<AttributeBuffer>>::iterator VertexArray::begin() {
	return attributeBuffers.begin();
}

vector<sptr<AttributeBuffer>>::iterator VertexArray::end() {
	return attributeBuffers.end();
}

void VertexArray::draw() const {
	THROW_IF(not elementBuffer, "No element buffer set for this VAO");
	GLCommand::drawIndexedTriangles(elementBuffer->getNumberOfIndices());
}
