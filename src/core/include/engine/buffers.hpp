#pragma once
#include "shaderDataTypes.hpp"
#include <GL/glew.h>


struct VertexAttribute {
	string name;
	ShaderDataType type;

	VertexAttribute(const string& name, ShaderDataType type);
	vs_dim elementLength() const;
	byte_size elementSize() const;
};

struct VertexBufferLayout {
	vector<VertexAttribute> attributes;
	vector<byte_size> attributeOffsets;
	byte_size strideBytes;

	VertexBufferLayout(initializer_list<pair<string, ShaderDataType>> attrs);
	byte_size offsetAt(int i) const;
	unsigned int numberOfAttributes() const;
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
	void uploadData(raw_data_ptr data, byte_size size) const;
	void updateData(raw_data_ptr data, byte_size size) const;
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
