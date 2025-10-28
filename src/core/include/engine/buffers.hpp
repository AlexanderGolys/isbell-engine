#pragma once
#include "shaderDataTypes.hpp"
#include <GL/glew.h>


struct VertexAttribute {
	string name;
	ShaderDataType type;

	VertexAttribute(const string& name, ShaderDataType type);
	unsigned int elementLength() const;
	size_t elementSize() const;
};

struct VertexBufferLayout {
	vector<VertexAttribute> attributes;
	vector<unsigned int> attributeOffsets;
	unsigned int strideBytes;

	VertexBufferLayout(initializer_list<pair<string, ShaderDataType>> attrs);
	const void* offsetAt(int i) const;
	unsigned int length() const;
};


class VertexBuffer {
	GLuint id = 0;
	VertexBufferLayout layout;

public:
	explicit VertexBuffer(const VertexBufferLayout& layout);
	~VertexBuffer();

	GLuint vboID() const;
	const VertexBufferLayout& getLayout() const;
	void bind() const;
	void unbind() const;
	void uploadData(const void* data, unsigned int size) const;
	void updateData(const void* data, unsigned int size) const;
	unsigned int length() const;
};


class IndexBuffer {
	GLuint id = 0;
	unsigned int count;

public:
	explicit IndexBuffer(unsigned int trianglesCount);
	explicit IndexBuffer(const vector<ivec3>& indices);
	~IndexBuffer();

	GLuint eboID() const;
	void bind() const;
	void unbind() const;
	void uploadData(const vector<ivec3>& indices);
	unsigned int getCount() const;
};

class VertexArray {
	GLuint id = 0;
	shared_ptr<IndexBuffer> ebo;
	vector<shared_ptr<VertexBuffer>> vbos;
	unsigned int attributeCount = 0;

public:
	VertexArray();
	VertexArray(const shared_ptr<IndexBuffer>& indexBuffer, const vector<shared_ptr<VertexBuffer>>& vertexBuffers);
	~VertexArray();

	GLuint vaoID() const;
	void bind() const;
	void unbind() const;
	void addVertexBuffer(const shared_ptr<VertexBuffer>& vertexBuffer);
	void setIndexBuffer(const shared_ptr<IndexBuffer>& indexBuffer);
	unsigned int getIndexCount() const;
};
