#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include "buffers.hpp"


GLenum baseType(AttributeType type);




class AttributeBufferGL : public AttributeBuffer {
	GLuint vboID = 0;
	ShaderAttribute attribute;
	int length;

public:
	AttributeBufferGL(const ShaderAttribute &attribute, int length);
	~AttributeBufferGL() override;

	void bind() override;
	void unbind() override;
	void load(void* firstElementAdress) override;
	void load(void* firstElementAdress, int length) override;
	const ShaderAttribute& getAttribute() const override;
	int getLength() const override;
	size_t size() const override;
};





class IndexBufferGL : public IndexBuffer {
	GLuint eboID = 0;
	int length;

public:
	IndexBufferGL(int length);

	void bind() override;
	void unbind() override;
	void load(void *firstElementAdress) override;
	void load(void *firstElementAdress, int bufferLength) override;
	int getLength() const override;
	size_t size() const override;
};



class VertexArrayObjectGL : public VertexArray {
	GLuint vaoID = 0;
	vector<shared_ptr<AttributeBuffer>> attributeBuffers;
	shared_ptr<IndexBuffer> indexBuffer;

public:
	VertexArrayObjectGL();
	~VertexArrayObjectGL() override;

	void bind() override;
	void unbind() override;
	void addAttributeBuffer(const shared_ptr<AttributeBuffer> &buffer) override;
	void setIndexBuffer(const shared_ptr<IndexBuffer> &buffer) override;
	const vector<shared_ptr<AttributeBuffer>>& getAttributeBuffers() const override;
	const shared_ptr<IndexBuffer>& getIndexBuffer() const override;
};

class ShaderStorageBufferObjectGL : public StorageBuffer {
	GLuint ssbo = 0;
	size_t bufferSize;

public:
	ShaderStorageBufferObjectGL(size_t bufferSize, const void* data);
	~ShaderStorageBufferObjectGL() override;

	void load(void *data, size_t size) override;
	void load(void *data) override;
	size_t getSize() const override;

	void bind() override;
	void unbind() override;
};
