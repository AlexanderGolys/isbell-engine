#include "buffers.hpp"
#include "exceptions.hpp"

AttributeDataArray::AttributeDataArray(const string& name, unsigned char elementLength)
: name(name), elementLen(elementLength) {
	THROW_IF(elementLength < 1, ValueError, "Attribute type must be real vector space of dimension at least 1");
	THROW_IF(elementLength > 4, ValueError, "Attribute type must be real vector space of dimension at most 4, got" + to_string(elementLength));
}

unsigned char AttributeDataArray::elementLength() const {
	return elementLen;
}

size_t AttributeDataArray::elementSize() const {
	return elementLen * sizeof(float);
}

unsigned int AttributeDataArray::length() const {
	return data.size() / elementLen;
}

size_t AttributeDataArray::size() const {
	return data.size() * sizeof(float);
}

const string& AttributeDataArray::getName() const {
	return name;
}

const float* AttributeDataArray::getDataPointer() const {
	return data.data();
}

void AttributeDataArray::setData(const vector<float>& new_data) {
	THROW_IF(elementLen != 1, ValueError, "to set data with vector<float>, element type must be float");
	this->data = new_data;
	markDirty();
}

void AttributeDataArray::setData(const vector<vec2>& new_data) {
	THROW_IF(elementLen != 2, ValueError, "to set data with vector<vec2>, element type must be vec2");
	if (new_data.size() == 0) {
		this->data = {};
		markDirty();
		return;
	}
	const float* ptr = reinterpret_cast<const float*>(new_data.data());
	this->data = vector(ptr, ptr + 2 * new_data.size());
	markDirty();
}

void AttributeDataArray::setData(const vector<vec3>& new_data) {
	THROW_IF(elementLen != 3, ValueError, "to set data with vector<vec3>, element type must be vec3");
	if (new_data.size() == 0) {
		this->data = {};
		markDirty();
		return;
	}
	const float* ptr = reinterpret_cast<const float*>(new_data.data());
	this->data = vector(ptr, ptr + 3 * new_data.size());
	markDirty();
}

void AttributeDataArray::setData(const vector<vec4>& new_data) {
	THROW_IF(elementLen != 4, ValueError, "to set data with vector<vec4>, element type must be vec4");
	if (new_data.size() == 0) {
		this->data = {};
		markDirty();
		return;
	}
	const float* ptr = reinterpret_cast<const float*>(new_data.data());
	this->data = vector(ptr, ptr + 4 * new_data.size());
	markDirty();
}

void AttributeDataArray::append(float value) {
	THROW_IF(elementLen != 1, ValueError, "to append float, element type must have len 1, not " + to_string(elementLen));
	data.push_back(value);
	markDirty();
}

void AttributeDataArray::append(vec2 value) {
	THROW_IF(elementLen != 2, ValueError, "to append vec2, element type must have len 2, not " + to_string(elementLen));
	data.push_back(value.x);
	data.push_back(value.y);
	markDirty();
}

void AttributeDataArray::append(vec3 value) {
	THROW_IF(elementLen != 3, ValueError, "to append vec3, element type must have len 3, not " + to_string(elementLen));
	data.push_back(value.x);
	data.push_back(value.y);
	data.push_back(value.z);
	markDirty();
}

void AttributeDataArray::append(vec4 value) {
	THROW_IF(elementLen != 4, ValueError, "to append vec4, element type must have len 4, not " + to_string(elementLen));
	data.push_back(value.x);
	data.push_back(value.y);
	data.push_back(value.z);
	data.push_back(value.w);
	markDirty();
}

void AttributeDataArray::setElement(int index, float value) {
	if (index < 0)
		index += length();
	THROW_IF(index < 0 || index >= length(), IndexOutOfBounds, index, length(), "AttributeDataArray " + name);
	THROW_IF(elementLen != 1, ValueError, "to set element with float, element type must have len 1, not " + to_string(elementLen));
	data[index] = value;
	markDirty();
}

void AttributeDataArray::setElement(int index, vec2 value) {
	if (index < 0)
		index += length();
	THROW_IF(index < 0 || index >= length(), IndexOutOfBounds, index, length(), "AttributeDataArray " + name);
	THROW_IF(elementLen != 2, ValueError, "to set element with vec2, element type must have len 2, not " + to_string(elementLen));
	data[index * 2] = value.x;
	data[index * 2 + 1] = value.y;
	markDirty();
}

void AttributeDataArray::setElement(int index, vec3 value) {
	if (index < 0)
		index += length();
	THROW_IF(index < 0 || index >= length(), IndexOutOfBounds, index, length(), "AttributeDataArray " + name);
	THROW_IF(elementLen != 3, ValueError, "to set element with vec3, element type must have len 3, not " + to_string(elementLen));
	data[index * 3] = value.x;
	data[index * 3 + 1] = value.y;
	data[index * 3 + 2] = value.z;
	markDirty();
}

void AttributeDataArray::setElement(int index, vec4 value) {
	if (index < 0)
		index += length();
	THROW_IF(index < 0 || index >= length(), IndexOutOfBounds, index, length(), "AttributeDataArray " + name);
	THROW_IF(elementLen != 4, ValueError, "to set element with vec4, element type must have len 4, not " + to_string(elementLen));
	data[index * 4] = value.x;
	data[index * 4 + 1] = value.y;
	data[index * 4 + 2] = value.z;
	data[index * 4 + 3] = value.w;
	markDirty();
}

float AttributeDataArray::getElement1f(int index) const {
	if (index < 0)
		index += length();
	THROW_IF(index < 0 || index >= length(), IndexOutOfBounds, index, length(), "AttributeDataArray " + name);
	THROW_IF(elementLen != 1, ValueError, "to get element as float, element type must be float");
	return data[index];
}

vec2 AttributeDataArray::getElement2f(int index) const {
	if (index < 0)
		index += length();
	THROW_IF(index < 0 || index >= length(), IndexOutOfBounds, index, length(), "AttributeDataArray " + name);
	THROW_IF(elementLen != 2, ValueError, "to get element as vec2, element type must be vec2");
	return vec2(data[index * 2], data[index * 2 + 1]);
}

vec3 AttributeDataArray::getElement3f(int index) const {
	if (index < 0)
		index += length();
	THROW_IF(index < 0 || index >= length(), IndexOutOfBounds, index, length(), "AttributeDataArray " + name);
	THROW_IF(elementLen != 3, ValueError, "to get element as vec3, element type must be vec3");
	return vec3(data[index * 3], data[index * 3 + 1], data[index * 3 + 2]);
}

vec4 AttributeDataArray::getElement4f(int index) const {
	if (index < 0)
		index += length();
	THROW_IF(index < 0 || index >= length(), IndexOutOfBounds, index, length(), "AttributeDataArray " + name);
	THROW_IF(elementLen != 4, ValueError, "to get element as vec4, element type must be vec4");
	return vec4(data[index * 4], data[index * 4 + 1], data[index * 4 + 2], data[index * 4 + 3]);
}

ShaderDataType AttributeDataArray::getBaseType() const {
	if (elementLen == 1)
		return FLOAT;
	if (elementLen == 2)
		return VEC2;
	if (elementLen == 3)
		return VEC3;
	if (elementLen == 4)
		return VEC4;
	THROW(UnknownVariantError, "AttributeDataArray has invalid element length " + to_string(elementLen));
}

void AttributeDataArray::markDirty() {
	dirty = true;
}

void AttributeDataArray::markClean() {
	dirty = false;
}

bool AttributeDataArray::isDirty() const {
	return dirty;
}

void AttributeDataArray::reserveSpace(unsigned int newLength) {
	data.reserve(newLength * elementLen);
}

VertexData::VertexData(std::initializer_list<pair<string, unsigned char>> attributeInfos) {
	for (const auto& attrInfo : attributeInfos)
		attributes.emplace_back(attrInfo.first, attrInfo.second);
}

void VertexData::reserve(unsigned int numberOfVertices, unsigned int numberOfTriangles) {
	for (auto& attr : attributes)
		attr.reserveSpace(numberOfVertices);
	indices.reserve(numberOfTriangles);
}

void VertexData::addAttribute(const string& name, unsigned char elementLength) {
	THROW_IF(vertexCount() > 0, ValueError, "Cannot add attribute after vertices have been added");
	attributes.emplace_back(name, elementLength);
}

void VertexData::setVertexAttribute(int vertexIndex, const string& attributeName, float value) {
	if (vertexIndex < 0)
		vertexIndex += vertexCount();
	THROW_IF(vertexIndex < 0 || vertexIndex >= vertexCount(), IndexOutOfBounds, vertexIndex, vertexCount(), "VertexData setVertexAttribute float");
	for (auto& attr : attributes) {
		if (attr.getName() == attributeName) {
			THROW_IF(attr.elementLength() != 1, ValueError, "Attribute " + attributeName + " is not of length 1");
			attr.setElement(vertexIndex, value);
			return;
		}
	}
	THROW(UnknownVariantError, "Attribute " + attributeName + " not found");
}

void VertexData::setVertexAttribute(int vertexIndex, const string& attributeName, const vec2& value) {
	if (vertexIndex < 0)
		vertexIndex += vertexCount();
	THROW_IF(vertexIndex < 0 || vertexIndex >= vertexCount(), IndexOutOfBounds, vertexIndex, vertexCount(), "VertexData setVertexAttribute vec2");
	for (auto& attr : attributes) {
		if (attr.getName() == attributeName) {
			THROW_IF(attr.elementLength() != 2, ValueError, "Attribute " + attributeName + " is not of length 2");
			attr.setElement(vertexIndex, value);
			return;
		}
	}
	THROW(UnknownVariantError, "Attribute " + attributeName + " not found");
}

void VertexData::setVertexAttribute(int vertexIndex, const string& attributeName, const vec3& value) {
	if (vertexIndex < 0)
		vertexIndex += vertexCount();
	THROW_IF(vertexIndex < 0 || vertexIndex >= vertexCount(), IndexOutOfBounds, vertexIndex, vertexCount(), "VertexData setVertexAttribute vec3");
	for (auto& attr : attributes) {
		if (attr.getName() == attributeName) {
			THROW_IF(attr.elementLength() != 3, ValueError, "Attribute " + attributeName + " is not of length 3");
			attr.setElement(vertexIndex, value);
			return;
		}
	}
	THROW(UnknownVariantError, "Attribute " + attributeName + " not found");
}

void VertexData::setVertexAttribute(int vertexIndex, const string& attributeName, const vec4& value) {
	if (vertexIndex < 0)
		vertexIndex += vertexCount();
	THROW_IF(vertexIndex < 0 || vertexIndex >= vertexCount(), IndexOutOfBounds, vertexIndex, vertexCount(), "VertexData setVertexAttribute vec4");
	for (auto& attr : attributes) {
		if (attr.getName() == attributeName) {
			THROW_IF(attr.elementLength() != 4, ValueError, "Attribute " + attributeName + " is not of length 4");
			attr.setElement(vertexIndex, value);
			return;
		}
	}
	THROW(UnknownVariantError, "Attribute " + attributeName + " not found");
}

void VertexData::addTriangle(int v1, int v2, int v3) {
	indices.emplace_back(v1, v2, v3);
}

void VertexData::addTriangle(const ivec3& triangle) {
	addTriangle(triangle.x, triangle.y, triangle.z);
}

unsigned int VertexData::vertexCount() const {
	if (attributes.empty())
		return 0;
	return attributes[0].length();
}

unsigned int VertexData::triangleCount() const {
	return indices.size();
}


AttributeBuffer::AttributeBuffer(const string& name, ShaderDataType type, int inputNumber)
: name(name), type(type), inputNumber(inputNumber) {}

AttributeBuffer::~AttributeBuffer() {
	if (initialized)
		free();
}

const string& AttributeBuffer::getName() const {
	return name;
}

GLuint AttributeBuffer::getVBO() const {
	return vboID;
}

ShaderDataType AttributeBuffer::getType() const {
	return type;
}

int AttributeBuffer::getInputNumber() const {
	return inputNumber;
}

bool AttributeBuffer::isInitialized() const {
	return initialized;
}

void AttributeBuffer::init() {
	if (initialized)
		return;

	glCreateBuffers(1, &vboID);
	initialized = true;
}

void AttributeBuffer::uploadData(const AttributeDataArray& dataArray) {
	if (!initialized)
		init();

	THROW_IF(dataArray.getBaseType() != type, ValueError,
			 "AttributeBuffer " + name + " type mismatch: expected " + to_string(type) + " but got " + to_string( dataArray.getBaseType()));

	glNamedBufferData(vboID, dataArray.size(), dataArray.getDataPointer(), GL_DYNAMIC_DRAW);
}

void AttributeBuffer::uploadDataIfDirty(AttributeDataArray& dataArray) {
	if (!dataArray.isDirty())
		return;

	uploadData(dataArray);
	dataArray.markClean();
}

void AttributeBuffer::free() {
	if (!initialized)
		return;

	glDeleteBuffers(1, &vboID);
	vboID = 0;
	initialized = false;
}


VAO::VAO() {}

VAO::~VAO() {
	if (initialized)
		free();
}

VAO::VAO(VAO&& other) noexcept
: vaoID(other.vaoID), eboID(other.eboID), attributes(std::move(other.attributes)), initialized(other.initialized), indexCount(other.indexCount) {
	other.vaoID = 0;
	other.eboID = 0;
	other.initialized = false;
	other.indexCount = 0;
}

VAO& VAO::operator=(VAO&& other) noexcept {
	if (this != &other) {
		if (initialized)
			free();

		vaoID = other.vaoID;
		eboID = other.eboID;
		attributes = std::move(other.attributes);
		initialized = other.initialized;
		indexCount = other.indexCount;

		other.vaoID = 0;
		other.eboID = 0;
		other.initialized = false;
		other.indexCount = 0;
	}
	return *this;
}

void VAO::init() {
	if (initialized)
		return;

	glCreateVertexArrays(1, &vaoID);
	glCreateBuffers(1, &eboID);
	glVertexArrayElementBuffer(vaoID, eboID);

	initialized = true;
}

void VAO::addAttribute(const shared_ptr<AttributeBuffer>& buffer, int components, GLenum type) {
	if (!initialized)
		init();

	if (!buffer->isInitialized())
		buffer->init();

	int location = attributes.size();

	glEnableVertexArrayAttrib(vaoID, location);
	glVertexArrayAttribFormat(vaoID, location, components, type, GL_FALSE, 0);
	glVertexArrayAttribBinding(vaoID, location, location);
	glVertexArrayVertexBuffer(vaoID, location, buffer->getVBO(), 0, 0);

	attributes.push_back(buffer);
}

void VAO::uploadIndices(const vector<ivec3>& indices) {
	if (!initialized)
		init();

	indexCount = indices.size() * 3;
	glNamedBufferData(eboID, indices.size() * sizeof(ivec3), indices.data(), GL_STATIC_DRAW);
}

void VAO::uploadIndices(const void* data, unsigned int count) {
	if (!initialized)
		init();

	indexCount = count;
	glNamedBufferData(eboID, count * sizeof(unsigned int), data, GL_STATIC_DRAW);
}

void VAO::updateDirtyBuffers(VertexData& vertexData) {
	for (auto& buffer : attributes) {
		for (auto& attrData : vertexData.attributes) {
			if (attrData.getName() == buffer->getName()) {
				buffer->uploadDataIfDirty(attrData);
				break;
			}
		}
	}
}

void VAO::bind() const {
	THROW_IF(!initialized, ValueError, "Cannot bind uninitialized VAO");
	glBindVertexArray(vaoID);
}

void VAO::unbind() const {
	glBindVertexArray(0);
}

void VAO::free() {
	if (!initialized)
		return;

	glDeleteVertexArrays(1, &vaoID);
	glDeleteBuffers(1, &eboID);

	vaoID = 0;
	eboID = 0;
	initialized = false;
	indexCount = 0;
}

GLuint VAO::getVAO() const {
	return vaoID;
}

GLuint VAO::getEBO() const {
	return eboID;
}

unsigned int VAO::getIndexCount() const {
	return indexCount;
}

bool VAO::isInitialized() const {
	return initialized;
}
