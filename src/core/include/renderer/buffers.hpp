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

	virtual void bind();
	virtual void unbind();
	virtual void load(void* firstElementAdress);
	virtual void load(void* firstElementAdress, int length);


	virtual const ShaderAttribute& getAttribute() const;
	virtual int getLength() const;
	virtual size_t size() const;

	static shared_ptr<AttributeBuffer> init(ShaderAttribute, int length);
};


class IndexBuffer {
public:
	virtual ~IndexBuffer() = default;

	virtual void bind();
	virtual void unbind();
	virtual void load(void* firstElementAdress);
	virtual void load(void* firstElementAdress, int bufferLength);

	virtual int getLength() const;
	virtual size_t size() const;

	static shared_ptr<AttributeBuffer> init(int length);
};

class VertexArray {
public:
	virtual ~VertexArray() = default;

	virtual void bind();
	virtual void unbind();

	virtual void addAttributeBuffer(const shared_ptr<AttributeBuffer> &buffer);
	virtual void setIndexBuffer(const shared_ptr<IndexBuffer> &buffer);

	virtual const vector<shared_ptr<AttributeBuffer>>& getAttributeBuffers() const;
	virtual const shared_ptr<IndexBuffer>& getIndexBuffer() const;

	static shared_ptr<VertexArray> init();
};

class StorageBuffer {
public:
	virtual ~StorageBuffer() = default;

	virtual void bind();
	virtual void unbind();

	virtual void load(void* data, size_t size);
	virtual void load(void* data);

	virtual size_t getSize() const;

	static shared_ptr<StorageBuffer> init(size_t size, const void* data);
};
