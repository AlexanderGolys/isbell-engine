#include "buffers.hpp"
#include "formatters.hpp"
#include "glCommand.hpp"
#include "logging.hpp"
#include <GL/glew.h>


VertexBufferLayout::VertexBufferLayout(const vector<GLSLPrimitive>& types): types(types) {
	stride = 0;
	offsets.reserve(types.size());
	for (const auto& t : types) {
		offsets.emplace_back(stride);
		stride += sizeOfGLSLType(t);
	}
}

VertexBufferLayout VertexBufferLayout::singleAttribute(GLSLPrimitive type) {
	return VertexBufferLayout({type});
}

bool VertexBufferLayout::operator==(const VertexBufferLayout& other) const {
	return types == other.types;
}

VertexBuffer::VertexBuffer(const VertexBufferLayout& layout, int firstInputNumber)
: layout(layout), bufferedSize(0), firstInputNumber(firstInputNumber) {
	GLCommand::createBuffer(bufferID);
}

VertexBuffer::~VertexBuffer() {
	GLCommand::deleteBuffer(bufferID);
}

void VertexBuffer::pointAtAttributes() const {
	for (int i = 0; i < layout.types.size(); i++) {
		int inputNumber = firstInputNumber + i;
		byte_size offset = layout.offsets[i];
		GLCommand::pointAtAttributeBuffer(inputNumber, bufferID, layout.types[i], layout.stride, reinterpret_cast<raw_data_ptr>(offset));
	}
}

void VertexBuffer::load(raw_data_ptr firstElementAdress, byte_size bufferSize) {
	if (bufferedSize != 0)
		LOG_WARN("Loading the vertex buffer again instead of updating it");
	GLCommand::loadBufferData(bufferID, firstElementAdress, bufferSize, GL_DYNAMIC_DRAW);
	bufferedSize = bufferSize;
}

void VertexBuffer::update(raw_data_ptr firstElementAdress, byte_size bufferSize) {
	GLCommand::updateBufferData(bufferID, firstElementAdress, bufferSize);
	if (bufferSize != bufferedSize)
		LOG(1, "Resized vertex buffer from " + formatByteSize(bufferedSize) + " to " + formatByteSize(bufferSize));
	bufferedSize = bufferSize;
}

void VertexBuffer::update(raw_data_ptr firstElementAdress) const {
	GLCommand::updateBufferData(bufferID, firstElementAdress, bufferedSize);
}

int VertexBuffer::getFirstInputNumber() const {
	return firstInputNumber;
}

int VertexBuffer::getNumberOfAttributes() const {
	return layout.types.size();
}

AttributeBuffer::AttributeBuffer(const string& name, GLSLPrimitive type, int inputNumber)
: VertexBuffer(VertexBufferLayout::singleAttribute(type), inputNumber), name(name), type(type) {
}

ElementBuffer::ElementBuffer() {
	GLCommand::createBuffer(bufferID);
}

ElementBuffer::~ElementBuffer() {
	GLCommand::deleteBuffer(bufferID);
}


void ElementBuffer::load(raw_data_ptr firstElementAdress, byte_size bufferSize) {
	if (bufferedSize != 0)
		LOG_WARN("Loading the element buffer again instead of updating it");
	GLCommand::loadBufferData(bufferID, firstElementAdress, bufferSize, GL_STATIC_DRAW);
	bufferedSize = bufferSize;
}

GLuint ElementBuffer::getID() const {
	return bufferID;
}

array_len ElementBuffer::getNumberOfIndices() const {
	return bufferedSize / sizeof(GLuint);
}

VertexArray::VertexArray() {
	GLCommand::createVAO(vaoID);
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
	THROW_IF(this->elementBuffer != nullptr, ValueError, "Element buffer is already set for this VAO");
	GLCommand::bindVAO(vaoID);
	GLCommand::pointAtElementBuffer(elementBuffer->getID());
	this->elementBuffer = elementBuffer;
}

void VertexArray::addElementBuffer() {
	addElementBuffer(make_shared<ElementBuffer>());
}

void VertexArray::addVertexBuffer(sptr<VertexBuffer> vertexBuffer) {
	GLCommand::bindVAO(vaoID);
	vertexBuffer->pointAtAttributes();
	vertexBuffers.push_back(vertexBuffer);
}

bool VertexArray::empty() const {
	return vertexBuffers.empty() and elementBuffer == nullptr;
}

vector<sptr<VertexBuffer>>::iterator VertexArray::begin() {
	return vertexBuffers.begin();
}

vector<sptr<VertexBuffer>>::iterator VertexArray::end() {
	return vertexBuffers.end();
}

void VertexArray::loadIndexedMesh(sptr<IndexedMesh3D> mesh) {
	THROW_IF(not empty(), ValueError, "Vertex array is already populated");
	sptr<ElementBuffer> elementBuffer = make_shared<ElementBuffer>();
	elementBuffer->load(mesh->bufferIndexLocation(), mesh->faceIndicesDataSize());
	addElementBuffer(elementBuffer);

	sptr<AttributeBuffer> positionBuffer = make_shared<AttributeBuffer>("position", GLSLPrimitive::VEC3, 0);
	positionBuffer->load(mesh->getBufferLocation("position"), mesh->getAttributeDataSize("position"));
	addVertexBuffer(positionBuffer);

	sptr<AttributeBuffer> normalBuffer = make_shared<AttributeBuffer>("normal", GLSLPrimitive::VEC3, 1);
	normalBuffer->load(mesh->getBufferLocation("normal"), mesh->getAttributeDataSize("normal"));
	addVertexBuffer(normalBuffer);

	sptr<AttributeBuffer> uvBuffer = make_shared<AttributeBuffer>("uv", GLSLPrimitive::VEC2, 2);
	uvBuffer->load(mesh->getBufferLocation("uv"), mesh->getAttributeDataSize("uv"));
	addVertexBuffer(uvBuffer);

	sptr<AttributeBuffer> colorBuffer = make_shared<AttributeBuffer>("color", GLSLPrimitive::VEC4, 3);
	colorBuffer->load(mesh->getBufferLocation("color"), mesh->getAttributeDataSize("color"));
	addVertexBuffer(colorBuffer);

	int attribIndex = 3;
	for (string attrName : mesh->getActiveExtraBuffers()) {
		sptr<AttributeBuffer> attr = make_shared<AttributeBuffer>(attrName, GLSLPrimitive::VEC4, ++attribIndex);
		colorBuffer->load(mesh->getBufferLocation(attrName), mesh->getAttributeDataSize(attrName));
		addVertexBuffer(attr);
	}
	LOG("Loaded mesh of total size " + formatByteSize(mesh->totalByteSize()) + ".");
}

void VertexArray::updateIndexedMesh(sptr<IndexedMesh3D> mesh) const {
	for (auto& vertexBuffer : vertexBuffers) {
		auto attrBuffer = std::dynamic_pointer_cast<AttributeBuffer>(vertexBuffer);
		attrBuffer->update(mesh->getBufferLocation(attrBuffer->get_name()));
	}
	LOG("Updated buffered mesh attributes of size " + formatByteSize(mesh->totalAttributeDataSize()) + ".");
}

void VertexArray::loadGeometricData(sptr<GeometricData> data) {
	THROW_IF(not empty(), ValueError, "VertexArray is not empty");

	sptr<ElementBuffer> eb = make_shared<ElementBuffer>();
	eb->load(data->indexBufferData(), data->indexBufferSize());
	addElementBuffer(eb);

	sptr<VertexBuffer> vb = make_shared<VertexBuffer>(data->layout());
	vb->load(data->vertexBufferData(), data->vertexBufferSize());
	addVertexBuffer(vb);

	data->markClean();
	LOG("Loaded geometric data with total size " + formatByteSize(data->totalSize()) + ".");
}

void VertexArray::updateGeometricData(sptr<GeometricData> data) const {
	THROW_IF(not elementBuffer, ValueError, "Element buffer is not set");
	THROW_IF(vertexBuffers.size() != 1, ValueError, "Vertex array must have exactly one vertex buffer to update geometric data");
	THROW_IF(vertexBuffers[0]->get_layout() != data->layout(), ValueError, "Vertex buffer layout does not match geometric data layout");

	if (not data->isDirty())
		return;

	vertexBuffers[0]->update(data->vertexBufferData(), data->vertexBufferSize());
	data->markClean();
	LOG("Updated geometric data vertex buffer of size " + formatByteSize(data->vertexBufferSize()) + ".");

}

void VertexArray::draw() const {
	THROW_IF(not elementBuffer, ValueError, "No element buffer set for this VAO");
	GLCommand::drawIndexedTriangles(elementBuffer->getNumberOfIndices());
}

ShaderStorageBuffer::ShaderStorageBuffer() {
	GLCommand::createBuffer(bufferID);
}

ShaderStorageBuffer::~ShaderStorageBuffer() {
	GLCommand::deleteBuffer(bufferID);
}

void ShaderStorageBuffer::bind(int bindingPoint) const {
	GLCommand::bindSSBO(bufferID, bindingPoint);
}

void ShaderStorageBuffer::unbind() const {
	GLCommand::unbindSSBO();
}

void ShaderStorageBuffer::load(byte_size bufferSize, raw_data_ptr data) const {
	GLCommand::loadSSBOData(bufferID, data, bufferSize);
	LOG("Loaded SSBO of size " + formatByteSize(bufferSize) + ".");

}

void ShaderStorageBuffer::update(byte_size bufferSize, raw_data_ptr data) const {
	GLCommand::updateSSBOData(bufferID, data, bufferSize);
	LOG("Updated SSBO of size " + formatByteSize(bufferSize) + ".");
}

UniformBufferObject::UniformBufferObject(uint bindingPoint, size_t bufferSize, string_cr name)
: bindingPoint(bindingPoint), bufferSize(bufferSize) , name(name) {
	GLCommand::createBuffer(bufferID);
}

UniformBufferObject::~UniformBufferObject() {
	GLCommand::deleteBuffer(bufferID);
}

void UniformBufferObject::initPerShader(gl_id currentProgram) const {
	GLCommand::initUBO(bufferID, bindingPoint, name, currentProgram);
}

void UniformBufferObject::unbind() const {
	GLCommand::unbindUBO();
}

void UniformBufferObject::load(raw_data_ptr firstElementAdress) const {
	GLCommand::loadUBOData(bufferID, firstElementAdress, bufferSize);
}

void UniformBufferObject::update(raw_data_ptr firstElementAdress) const {
	GLCommand::updateUBOData(bufferID, firstElementAdress, bufferSize);
}
