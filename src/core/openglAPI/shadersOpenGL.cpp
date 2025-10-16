#include "shadersOpenGL.hpp"
#include "exceptions.hpp"


GLenum gl_shader_type(ShaderType type) {
	switch (type) {
		case VERTEX_SHADER: return GL_VERTEX_SHADER;
		case FRAGMENT_SHADER: return GL_FRAGMENT_SHADER;
		case GEOMETRY_SHADER: return GL_GEOMETRY_SHADER;
		case COMPUTE_SHADER: return GL_COMPUTE_SHADER;
	}
	THROW(UnknownVariantError, "ShaderType not recognized");
}

ShaderGL::ShaderGL(CodeFileDescriptor &file) {
	type = identify_shader_extension(file.extension());
	shaderID = glCreateShader(gl_shader_type(type));
	std::string code = file.getCode();
	const char* sourcePtr = code.c_str();

	int InfoLogLength;

	glShaderSource(shaderID, 1, &sourcePtr, nullptr);
	glCompileShader(shaderID);

	glGetShaderiv(shaderID, GL_INFO_LOG_LENGTH, &InfoLogLength);
	if (InfoLogLength > 0) {
		vector<char> error(InfoLogLength + 1);
		glGetShaderInfoLog(shaderID, InfoLogLength, nullptr, &error[0]);
		THROW(CompilationError, &error[0]);
	}
	LOG("Shader compiled: " + file.getPath().string());

}

ShaderType ShaderGL::getType() const {
	return type;
}
