#include "buffers.hpp"
#include "renderer.hpp"
#include "exceptions.hpp"
#include "buffersOpenGL.hpp"

size_t sizeOfAttributeType(AttributeType type) {
	switch (type) {
		case FLOAT:
		case UINT:
		case INT: return 4;

		case VEC2:
		case IVEC2: return 4 * 2;

		case VEC3:
		case IVEC3: return 4 * 3;

		case VEC4:
		case IVEC4: return 4 * 4;

		case UBYTE:
		case BYTE: return 1;

		case USHORT:
		case SHORT: return 2;
	}
	THROW(UnknownVariantError, "Attribute type not recognized");
}

int lengthOfAttributeType(AttributeType type) {
	switch (type) {
		case INT:
		case FLOAT:
		case BYTE:
		case SHORT:
		case UINT:
		case UBYTE:
		case USHORT:
			return 1;

		case VEC2:
		case IVEC2:
			return 2;

		case VEC3:
		case IVEC3:
			return 3;

		case VEC4:
		case IVEC4:
			return 4;
	}
	THROW(UnknownVariantError, "Attribute type not recognized");
}

ShaderAttribute::ShaderAttribute(const string &name, AttributeType type):
name(name), type(type), size(sizeOfAttributeType(type)), length(lengthOfAttributeType(type)) {}

shared_ptr<AttributeBuffer> AttributeBuffer::init(const ShaderAttribute &attribute, int length) {
	RenderAPI api = Renderer::getAPI();
	if (api == NO_RENDERING)
		THROW(RendererError, "Rendering API not set.");
	if (api == OPENGL)
		return make_shared<AttributeBufferGL>(attribute, length);
	if (api == VULKAN) {
		THROW(NotImplementedVariantError, "Vulkan API", "Rendering API");
	}
	THROW(UnknownVariantError, "Rendering API not recognized.");
}

shared_ptr<IndexBuffer> IndexBuffer::init(int length)
{
	RenderAPI api = Renderer::getAPI();
	if (api == NO_RENDERING)
		THROW(RendererError, "Rendering API not set.");
	if (api == OPENGL)
		return make_shared<IndexBufferGL>(length);
	if (api == VULKAN)
		THROW(NotImplementedVariantError, "Vulkan API", "Rendering API");
	THROW(UnknownVariantError, "Rendering API not recognized.");
}

shared_ptr<VertexArray> VertexArray::init()
{
	RenderAPI api = Renderer::getAPI();
	if (api == NO_RENDERING)
		THROW(RendererError, "Rendering API not set.");
	if (api == OPENGL)
		return make_shared<VertexArrayObjectGL>();
	if (api == VULKAN)
		THROW(NotImplementedVariantError, "Vulkan API", "Rendering API");
	THROW(UnknownVariantError, "Rendering API not recognized.");
}

shared_ptr<StorageBuffer> StorageBuffer::init(size_t size, const void* data)
{
	RenderAPI api = Renderer::getAPI();
	if (api == NO_RENDERING)
		THROW(RendererError, "Rendering API not set.");
	if (api == OPENGL)
		return make_shared<ShaderStorageBufferObjectGL>(size, data);
	if (api == VULKAN)
		THROW(NotImplementedVariantError, "Vulkan API", "Rendering API");
	THROW(UnknownVariantError, "Rendering API not recognized.");
}
