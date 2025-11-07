#pragma once
#include "renderingUtils.hpp"

class AttributeBuffer {
	string name;
	GLuint bufferID;
	GLSLType type;
	int inputNumber;
	byte_size bufferedSize;

public:
	AttributeBuffer(const string& name, GLSLType type, int inputNumber);
	~AttributeBuffer();

	void pointAtAttributeBuffer() const;
	void load(raw_data_ptr firstElementAdress, byte_size bufferSize);
	void update(raw_data_ptr firstElementAdress, byte_size bufferSize);
};

class ElementBuffer {
	GLuint bufferID;
	byte_size bufferedSize = 0;

public:
	ElementBuffer();
	~ElementBuffer();

	void load(raw_data_ptr firstElementAdress, byte_size bufferSize);
	GLuint getID() const;
	array_len getNumberOfIndices() const;
};

class VertexArray {
	GLuint vaoID = 0;
	sptr<ElementBuffer> elementBuffer = nullptr;
	vector<sptr<AttributeBuffer>> attributeBuffers;
public:
	VertexArray();
	~VertexArray();

	void bind() const;
	void unbind() const;
	void addElementBuffer(sptr<ElementBuffer> elementBuffer);
	void addAttributeBuffer(sptr<AttributeBuffer> attributeBuffer);

	sptr<ElementBuffer> getElementBuffer() const;
	vector<sptr<AttributeBuffer>>::iterator begin();
	vector<sptr<AttributeBuffer>>::iterator end();

	void draw() const;
};