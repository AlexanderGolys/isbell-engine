#pragma once
#include "shaderTypes.hpp"


struct ShaderAttribute {
	string name;
	AttributeType type;
	size_t size;
	int length;

	ShaderAttribute(const string &name, AttributeType type);
};



class AttributeBuffer {
public:
	virtual ~AttributeBuffer() = default;

	virtual void bind() = 0;
	virtual void unbind() = 0;
	virtual void load(void* firstElementAdress) = 0;
	virtual void load(void* firstElementAdress, int length) = 0;


	virtual const ShaderAttribute& getAttribute() const = 0;
	virtual int getLength() const = 0;
	virtual size_t size() const = 0;

	static shared_ptr<AttributeBuffer> init(const ShaderAttribute &attribute, int length);

};


class IndexBuffer {
public:
	virtual ~IndexBuffer() = default;

	virtual void bind() = 0;
	virtual void unbind() = 0;
	virtual void load(void* firstElementAdress) = 0;
	virtual void load(void* firstElementAdress, int bufferLength) = 0;

	virtual int getLength() const = 0;
	virtual size_t size() const = 0;

	static shared_ptr<IndexBuffer> init(int length);
};

class VertexArray {
public:
	virtual ~VertexArray() = default;

	virtual void bind() = 0;
	virtual void unbind() = 0;

	virtual void addAttributeBuffer(const shared_ptr<AttributeBuffer> &buffer) = 0;
	virtual void setIndexBuffer(const shared_ptr<IndexBuffer> &buffer) = 0;

	virtual const vector<shared_ptr<AttributeBuffer>>& getAttributeBuffers() const = 0;
	virtual const shared_ptr<IndexBuffer>& getIndexBuffer() const = 0;

	static shared_ptr<VertexArray> init();
};

class StorageBuffer {
public:
	virtual ~StorageBuffer() = default;

	virtual void bind() = 0;
	virtual void unbind() = 0;

	virtual void load(void* data, size_t size) = 0;
	virtual void load(void* data) = 0;

	virtual size_t getSize() const = 0;

	static shared_ptr<StorageBuffer> init(size_t size, const void* data);
};
