#include "shaders.hpp"
#include "glCommand.hpp"


GLenum Shader::getTypeFromExtension(const string& extension) {
	string nodot = extension;
	if (nodot[0] == '.')
		nodot = nodot.substr(1);

	if (nodot == "vert")
		return GL_VERTEX_SHADER;
	if (nodot == "frag")
		return GL_FRAGMENT_SHADER;
	if (nodot == "geom" || nodot == "geo")
		return GL_GEOMETRY_SHADER;
	if (nodot == "comp")
		return GL_COMPUTE_SHADER;
	THROW(SystemError, "Unknown shader getExtension: ." + nodot);
}

VertexShader::VertexShader(const Path& file)
: Shader(file) {
	THROW_IF(get_shaderType() != GL_VERTEX_SHADER, ValueError, "Provided shader is not a vertex shader.");
}

VertexShader::VertexShader(const string& code)
: Shader(code, GL_VERTEX_SHADER) {}

FragmentShader::FragmentShader(const Path& file)
: Shader(file) {
	THROW_IF(get_shaderType() != GL_FRAGMENT_SHADER, ValueError, "Provided shader is not a fragment shader.");
}

FragmentShader::FragmentShader(const string& code)
: Shader(code, GL_FRAGMENT_SHADER) {}

GeometryShader::GeometryShader(const Path& file)
: Shader(file) {
	THROW_IF(get_shaderType() != GL_GEOMETRY_SHADER, ValueError, "Provided shader is not a geometry shader.");
}

GeometryShader::GeometryShader(const string& code)
: Shader(code, GL_FRAGMENT_SHADER) {}

ComputeShader::ComputeShader(const Path& file)
: Shader(file) {
	THROW_IF(get_shaderType() != GL_COMPUTE_SHADER, ValueError, "Provided shader is not a compute shader.");
}

ComputeShader::ComputeShader(const string& code)
: Shader(code, GL_COMPUTE_SHADER) {}


void Shader::compile(const string& code) {
	shaderID = GLCommand::createShader(shaderType);
	GLCommand::compileShaderSource(shaderID, code);
}

Shader::Shader(const string& code, GLenum shaderType)
: shaderType(shaderType){
	compile(code);
}

Shader::~Shader() {
	GLCommand::deleteShader(shaderID);
}


Shader::Shader(const Path& file) {
	CodeFileDescriptor fd(file);
	shaderType = getTypeFromExtension(fd.extension());
	string code = fd.getCode();
	compile(code);
}

ShaderProgram::ShaderProgram(const Shader& vertexShader, const Shader& fragmentShader, const Shader& geometryShader) : ShaderProgram() {
	GLCommand::attachShaderToProgram(programID, vertexShader.get_shaderID());
	GLCommand::attachShaderToProgram(programID, fragmentShader.get_shaderID());
	GLCommand::attachShaderToProgram(programID, geometryShader.get_shaderID());

	GLCommand::linkProgram(programID);
}

ShaderProgram::ShaderProgram() {
	programID = GLCommand::createProgram();
}

ShaderProgram::ShaderProgram(const Shader& vertexShader, const Shader& fragmentShader) : ShaderProgram() {
	GLCommand::attachShaderToProgram(programID, vertexShader.get_shaderID());
	GLCommand::attachShaderToProgram(programID, fragmentShader.get_shaderID());
	GLCommand::linkProgram(programID);
}


ShaderProgram::~ShaderProgram() {
	GLCommand::deleteProgram(programID);
}

void ShaderProgram::bind() const {
	GLCommand::bindProgram(programID);
}

void ShaderProgram::unbind() const {
	GLCommand::unbindProgram();
}

sptr<ShaderProgram> ShaderProgram::standardShaderProgram(Path vertexShaderPath, Path fragmentShaderPath) {
	return make_shared<ShaderProgram>(Shader(vertexShaderPath), Shader(fragmentShaderPath));
}

sptr<ShaderProgram> ShaderProgram::geometryShaderProgram(Path vertexShaderPath, Path fragmentShaderPath, Path geometryShaderPath) {
	return make_shared<ShaderProgram>(Shader(vertexShaderPath), Shader(fragmentShaderPath), Shader(geometryShaderPath));
}

ComputeShaderProgram::ComputeShaderProgram(const Shader& computeShader) : ShaderProgram() {
	GLCommand::attachShaderToProgram(get_programID(), computeShader.get_shaderID());
	GLCommand::linkProgram(get_programID());
}

void ComputeShaderProgram::run(int numGroupsX, int numGroupsY, int numGroupsZ) const {
	bind();
	GLCommand::runComputeShader(numGroupsX, numGroupsY, numGroupsZ);
}

sptr<ComputeShaderProgram> ComputeShaderProgram::standardComputeShaderProgram(Path computeShaderPath) {
	return make_shared<ComputeShaderProgram>(Shader(computeShaderPath));
}
