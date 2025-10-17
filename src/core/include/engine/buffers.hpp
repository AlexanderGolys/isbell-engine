#pragma once
#include "renderingUtils.hpp"
#include "shaderDataTypes.hpp"

// Modern OpenGL buffer management using Direct State Access (DSA) API
//
// Usage example:
//
//   VertexData vertexData({{"position", 3}, {"normal", 3}, {"uv", 2}});
//   vertexData.reserve(100, 50);
//
//   auto vao = make_shared<VAO>();
//   auto posBuffer = make_shared<AttributeBuffer>("position", VEC3, 0);
//   auto normalBuffer = make_shared<AttributeBuffer>("normal", VEC3, 1);
//
//   vao->addAttribute(0, posBuffer, 3, GL_FLOAT, sizeof(float) * 3, 0);
//   vao->addAttribute(1, normalBuffer, 3, GL_FLOAT, sizeof(float) * 3, 0);
//   vao->uploadIndices(vertexData.indices);
//
//   In render loop:
//   vao->updateDirtyBuffers(vertexData);
//   vao->bind();
//   glDrawElements(GL_TRIANGLES, vao->getIndexCount(), GL_UNSIGNED_INT, 0);
//



class AttributeDataArray {
	string name;
	vector<float> data;
	unsigned char elementLen;
	bool dirty = true;

public:
	AttributeDataArray(const string &name, unsigned char elementLength);

	unsigned char elementLength() const;
	size_t elementSize() const;
	unsigned int length() const;
	size_t size() const;
	const string &getName() const;
	const float *getDataPointer() const;

	void setData(const vector<float> &new_data);
	void setData(const vector<vec2> &new_data);
	void setData(const vector<vec3> &new_data);
	void setData(const vector<vec4> &new_data);

	void append(float value);
	void append(vec2 value);
	void append(vec3 value);
	void append(vec4 value);

	void setElement(int index, float value);
	void setElement(int index, vec2 value);
	void setElement(int index, vec3 value);
	void setElement(int index, vec4 value);

	float getElement1f(int index) const;
	vec2 getElement2f(int index) const;
	vec3 getElement3f(int index) const;
	vec4 getElement4f(int index) const;

	ShaderDataType getBaseType() const;

	void markDirty();
	void markClean();
	bool isDirty() const;

	void reserveSpace(unsigned int newLength) {
		data.reserve(newLength * elementLen);
	}
};


class VertexData {
public:
	vector<AttributeDataArray> attributes;
	vector<ivec3> indices;

	VertexData(std::initializer_list<pair<string, unsigned char>> attributeInfos) {
		for (const auto &attrInfo: attributeInfos)
			attributes.emplace_back(attrInfo.first, attrInfo.second);
	}

	void reserve(unsigned int numberOfVertices, unsigned int numberOfTriangles) {
		for (auto &attr: attributes)
			attr.reserveSpace(numberOfVertices);
		indices.reserve(numberOfTriangles);
	}

	template<typename VertexStruct>
	vector<VertexStruct> toVertexArray() const;

	template<typename VertexStruct>
	void addVertex(const VertexStruct &vertex) {
		for (auto &attr: attributes) {
			if (attr.elementLength() == 1)
				attr.append(vertex.getAttribute1f(attr.getName()));
			else if (attr.elementLength() == 2)
				attr.append(vertex.getAttribute2f(attr.getName()));
			else if (attr.elementLength() == 3)
				attr.append(vertex.getAttribute3f(attr.getName()));
			else if (attr.elementLength() == 4)
				attr.append(vertex.getAttribute4f(attr.getName()));
			else
				THROW(UnknownVariantError, "Unsupported attribute length in addVertex");
		}
	}

	template<typename VertexStruct>
	void fromVertexArray(const vector<VertexStruct> &vertexArray) {
		for (const auto &vertex: vertexArray) {
			addVertex(vertex);
		}
	}

	template<typename VertexStruct>
	void setVertex(int index, const VertexStruct &vertex) {
		if (index < 0)
			index += vertexCount();
		THROW_IF(index < 0 || index >= vertexCount(), IndexOutOfBounds, index, vertexCount(), "VertexData setVertex");
		for (auto &attr: attributes) {
			if (attr.elementLength() == 1)
				attr.setElement(index, vertex.getAttribute1f(attr.getName()));
			else if (attr.elementLength() == 2)
				attr.setElement(index, vertex.getAttribute2f(attr.getName()));
			else if (attr.elementLength() == 3)
				attr.setElement(index, vertex.getAttribute3f(attr.getName()));
			else if (attr.elementLength() == 4)
				attr.setElement(index, vertex.getAttribute4f(attr.getName()));
			else
				THROW(UnknownVariantError, "Unsupported attribute length in setVertex");
		}
	}

	void setVertexAttribute(int vertexIndex, const string &attributeName, float value) {
		if (vertexIndex < 0)
			vertexIndex += vertexCount();
		THROW_IF(vertexIndex < 0 || vertexIndex >= vertexCount(), IndexOutOfBounds, vertexIndex, vertexCount(), "VertexData setVertexAttribute float");
		for (auto &attr: attributes) {
			if (attr.getName() == attributeName) {
				THROW_IF(attr.elementLength() != 1, ValueError, "Attribute " + attributeName + " is not of length 1");
				attr.setElement(vertexIndex, value);
				return;
			}
		}
		THROW(UnknownVariantError, "Attribute " + attributeName + " not found");
	}

	void setVertexAttribute(int vertexIndex, const string &attributeName, const vec2 &value) {
		if (vertexIndex < 0)
			vertexIndex += vertexCount();
		THROW_IF(vertexIndex < 0 || vertexIndex >= vertexCount(), IndexOutOfBounds, vertexIndex, vertexCount(), "VertexData setVertexAttribute vec2");
		for (auto &attr: attributes) {
			if (attr.getName() == attributeName) {
				THROW_IF(attr.elementLength() != 2, ValueError, "Attribute " + attributeName + " is not of length 2");
				attr.setElement(vertexIndex, value);
				return;
			}
		}
		THROW(UnknownVariantError, "Attribute " + attributeName + " not found");
	}

	void setVertexAttribute(int vertexIndex, const string &attributeName, const vec3 &value) {
		if (vertexIndex < 0)
			vertexIndex += vertexCount();
		THROW_IF(vertexIndex < 0 || vertexIndex >= vertexCount(), IndexOutOfBounds, vertexIndex, vertexCount(), "VertexData setVertexAttribute vec3");
		for (auto &attr: attributes) {
			if (attr.getName() == attributeName) {
				THROW_IF(attr.elementLength() != 3, ValueError, "Attribute " + attributeName + " is not of length 3");
				attr.setElement(vertexIndex, value);
				return;
			}
		}
		THROW(UnknownVariantError, "Attribute " + attributeName + " not found");
	}

	void setVertexAttribute(int vertexIndex, const string &attributeName, const vec4 &value) {
		if (vertexIndex < 0)
			vertexIndex += vertexCount();
		THROW_IF(vertexIndex < 0 || vertexIndex >= vertexCount(), IndexOutOfBounds, vertexIndex, vertexCount(), "VertexData setVertexAttribute vec4");
		for (auto &attr: attributes) {
			if (attr.getName() == attributeName) {
				THROW_IF(attr.elementLength() != 4, ValueError, "Attribute " + attributeName + " is not of length 4");
				attr.setElement(vertexIndex, value);
				return;
			}
		}
		THROW(UnknownVariantError, "Attribute " + attributeName + " not found");
	}



	template<typename VertexStruct>
	VertexStruct getVertex(int index) const {
		if (index < 0)
			index += vertexCount();
		THROW_IF(index < 0 || index >= vertexCount(), IndexOutOfBounds, index, vertexCount(), "VertexData getVertex");

		VertexStruct vertex;
		for (const auto &attr: attributes) {
			if (attr.elementLength() == 1)
				vertex.setAttribute1f(attr.getName(), attr.getElement1f(index));
			else if (attr.elementLength() == 2)
				vertex.setAttribute2f(attr.getName(), attr.getElement2f(index));
			else if (attr.elementLength() == 3)
				vertex.setAttribute3f(attr.getName(), attr.getElement3f(index));
			else if (attr.elementLength() == 4)
				vertex.setAttribute4f(attr.getName(), attr.getElement4f(index));
			else
				THROW(UnknownVariantError, "Unsupported attribute length in getVertex");
		}
		return vertex;
	}

	void addTriangle(int v1, int v2, int v3) {
		int vertCount = vertexCount();
		indices.emplace_back(v1, v2, v3);
	}

	void addTriangle(const ivec3 &triangle) {
		addTriangle(triangle.x, triangle.y, triangle.z);
	}

	unsigned int vertexCount() const {
		if (attributes.empty())
			return 0;
		return attributes[0].length();
	}

	unsigned int triangleCount() const {
		return indices.size();
	}
};



class AttributeBuffer {
	string name;
	GLuint vboID = 0;
	ShaderDataType type;
	int inputNumber;
	bool initialized = false;

public:
	AttributeBuffer(const string &name, ShaderDataType type, int inputNumber);
	~AttributeBuffer();

	const string &getName() const;
	GLuint getVBO() const;
	ShaderDataType getType() const;
	int getInputNumber() const;
	bool isInitialized() const;

	void init();
	void uploadData(const AttributeDataArray &dataArray);
	void uploadDataIfDirty(AttributeDataArray &dataArray);
	void free();
};


class VAO {
	GLuint vaoID = 0;
	GLuint eboID = 0;
	vector<shared_ptr<AttributeBuffer>> attributes;
	bool initialized = false;
	unsigned int indexCount = 0;

public:
	VAO();
	~VAO();

	VAO(const VAO &) = delete;
	VAO &operator=(const VAO &) = delete;
	VAO(VAO &&other) noexcept;
	VAO &operator=(VAO &&other) noexcept;

	void init();
	void addAttribute(int location, shared_ptr<AttributeBuffer> buffer, int components, GLenum type, int stride, int offset);
	void uploadIndices(const vector<ivec3> &indices);
	void uploadIndices(const void *data, unsigned int count);
	void updateDirtyBuffers(VertexData &vertexData);
	void bind() const;
	void unbind() const;
	void free();

	GLuint getVAO() const;
	GLuint getEBO() const;
	unsigned int getIndexCount() const;
	bool isInitialized() const;
};
