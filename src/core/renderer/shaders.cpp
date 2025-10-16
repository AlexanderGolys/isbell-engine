#include "shaders.hpp"
#include "renderer.hpp"
#include "exceptions.hpp"
#include "shadersOpenGL.hpp"


ShaderType identify_shader_extension(const string &extension) {
	string ext = remove(extension, ".");
	if (ext == "vert")
		return VERTEX_SHADER;
	if (ext == "frag")
		return FRAGMENT_SHADER;
	if (ext == "geom" or ext == "geo")
		return GEOMETRY_SHADER;
	if (ext == "comp")
		return COMPUTE_SHADER;
	THROW(UnknownVariantError, "Shader extension not recognized: " + extension);
}

shared_ptr<Shader> Shader::compile_shader(CodeFileDescriptor &file)
{
	RenderAPI api = Renderer::getAPI();
	if (api == NO_RENDERING)
		THROW(RendererError, "Rendering API not set.");
	if (api == OPENGL)
		return make_shared<ShaderGL>(file);
	if (api == VULKAN)
		THROW(NotImplementedVariantError, "Vulkan API", "Rendering API");
	THROW(UnknownVariantError, "Rendering API not recognized.");
}

shared_ptr<ShaderProgram> ShaderProgram::create_program(const shared_ptr<Shader> &vertexShader, const shared_ptr<Shader> &fragmentShader, const shared_ptr<Shader> &geometryShader)
{
	RenderAPI api = Renderer::getAPI();
	if (api == NO_RENDERING)
		THROW(RendererError, "Rendering API not set.");
	if (api == OPENGL)
		return make_shared<ShaderProgramGL>(vertexShader, fragmentShader, geometryShader);
	if (api == VULKAN)
		THROW(NotImplementedVariantError, "Vulkan API", "Rendering API");
	THROW(UnknownVariantError, "Rendering API not recognized.");
}

shared_ptr<ComputeShaderProgram> ComputeShaderProgram::create_compute_shader(const shared_ptr<Shader> &computeShader)
{
	RenderAPI api = Renderer::getAPI();
	if (api == NO_RENDERING)
		THROW(RendererError, "Rendering API not set.");
	if (api == OPENGL)
		return make_shared<ComputeShaderProgramGL>(computeShader);
	if (api == VULKAN)
		THROW(NotImplementedVariantError, "Vulkan API", "Rendering API");
	THROW(UnknownVariantError, "Rendering API not recognized.");
}
