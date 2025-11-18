#pragma once
#include "glCommand.hpp"
#include "indexedRendering.hpp"
#include "indexedMesh.hpp"



class VertexBuffer {
	CONST_PROPERTY(VertexBufferLayout, layout);
	GLuint bufferID;
	byte_size bufferedSize;
	int firstInputNumber;

public:
	explicit VertexBuffer(const VertexBufferLayout& layout, int firstInputNumber=0);
	virtual ~VertexBuffer();
	VertexBuffer(const VertexBuffer&) = delete;
	VertexBuffer& operator=(const VertexBuffer&) = delete;

	void pointAtAttributes() const;
	void load(raw_data_ptr firstElementAdress, byte_size bufferSize);
	void update(raw_data_ptr firstElementAdress, byte_size bufferSize);
	void update(raw_data_ptr firstElementAdress) const;

	int getFirstInputNumber() const;
	int getNumberOfAttributes() const;
};

class AttributeBuffer : public VertexBuffer {
	CONST_PROPERTY(string, name);
	GLSLPrimitive type;

public:
	AttributeBuffer(const string& name, GLSLPrimitive type, int inputNumber);

};

class ElementBuffer {
	GLuint bufferID;
	byte_size bufferedSize = 0;

public:
	ElementBuffer();
	~ElementBuffer();
	ElementBuffer(const ElementBuffer&) = delete;
	ElementBuffer& operator=(const ElementBuffer&) = delete;

	void load(raw_data_ptr firstElementAdress, byte_size bufferSize);
	GLuint getID() const;
	array_len getNumberOfIndices() const;
};

class VertexArray {
	GLuint vaoID = 0;
	sptr<ElementBuffer> elementBuffer = nullptr;
	vector<sptr<VertexBuffer>> vertexBuffers;
public:
	VertexArray();
	~VertexArray();
	VertexArray(const VertexArray&) = delete;
	VertexArray& operator=(const VertexArray&) = delete;

	void bind() const;
	void unbind() const;
	void addElementBuffer(sptr<ElementBuffer> elementBuffer);
	void addElementBuffer();
	void addVertexBuffer(sptr<VertexBuffer> vertexBuffer);
	bool empty() const;

	vector<sptr<VertexBuffer>>::iterator begin();
	vector<sptr<VertexBuffer>>::iterator end();

	void loadIndexedMesh(sptr<IndexedMesh3D> mesh);
	void updateIndexedMesh(sptr<IndexedMesh3D> mesh) const;

	void loadGeometricData(sptr<GeometricData> data);
	void updateGeometricData(sptr<GeometricData> data) const;

	void draw() const;
};

class ShaderStorageBuffer {
	gl_id bufferID = 0;
public:
	ShaderStorageBuffer();
	~ShaderStorageBuffer();
	ShaderStorageBuffer(const ShaderStorageBuffer&) = delete;
	ShaderStorageBuffer& operator=(const ShaderStorageBuffer&) = delete;

	void bind(int bindingPoint) const;
	void unbind() const;
	void load(byte_size bufferSize, raw_data_ptr data) const;
	void update(byte_size bufferSize, raw_data_ptr data) const;
};

class UniformBufferObject {
	gl_id bufferID = 0;
	uint bindingPoint;
	byte_size bufferSize;
	string name;
public:
	UniformBufferObject(uint bindingPoint, byte_size bufferSize, string_cr name);
	~UniformBufferObject();
	UniformBufferObject(const UniformBufferObject&) = delete;
	UniformBufferObject& operator=(const UniformBufferObject&) = delete;

	void initPerShader(gl_id currentProgram) const;
	void unbind() const;
	void load(raw_data_ptr firstElementAdress) const;
	void update(raw_data_ptr firstElementAdress) const;
};