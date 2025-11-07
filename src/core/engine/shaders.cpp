#include "shaders.hpp"

#include <ranges>


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
	THROW_IF(getType() != GL_VERTEX_SHADER, ValueError, "Provided shader is not a vertex shader.");
}

VertexShader::VertexShader(const string& code)
: Shader(code, GL_VERTEX_SHADER) {}

FragmentShader::FragmentShader(const Path& file)
: Shader(file) {
	THROW_IF(getType() != GL_FRAGMENT_SHADER, ValueError, "Provided shader is not a fragment shader.");
}

FragmentShader::FragmentShader(const string& code)
: Shader(code, GL_FRAGMENT_SHADER) {}

GeometryShader::GeometryShader(const Path& file)
: Shader(file) {
	THROW_IF(getType() != GL_GEOMETRY_SHADER, ValueError, "Provided shader is not a geometry shader.");
}

GeometryShader::GeometryShader(const string& code)
: Shader(code, GL_FRAGMENT_SHADER) {}

ComputeShader::ComputeShader(const Path& file)
: Shader(file) {
	THROW_IF(getType() != GL_COMPUTE_SHADER, ValueError, "Provided shader is not a compute shader.");
}

ComputeShader::ComputeShader(const string& code)
: Shader(code, GL_COMPUTE_SHADER) {}


GLuint Shader::getID() const {
	return shaderID;
}

GLenum Shader::getType() const {
	return shaderType;
}

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

// TODO commands
ShaderProgram::ShaderProgram(const Shader& vertexShader, const Shader& fragmentShader, const Shader& geometryShader) {
	GLint Result = GL_FALSE;
	int InfoLogLength;

	GLuint VertexShaderID = vertexShader.getID();
	GLuint FragmentShaderID = fragmentShader.getID();
	GLuint GeometryShaderID = geometryShader.getID();

	programID = glCreateProgram();

	glAttachShader(programID, VertexShaderID);
	glAttachShader(programID, FragmentShaderID);
	glAttachShader(programID, GeometryShaderID);

	glLinkProgram(programID);
	glGetProgramiv(programID, GL_LINK_STATUS, &Result);
	glGetProgramiv(programID, GL_INFO_LOG_LENGTH, &InfoLogLength);

	if (InfoLogLength > 0) {
		vector<char> ProgramErrorMessage(InfoLogLength + 1);
		glGetProgramInfoLog(programID, InfoLogLength, nullptr, &ProgramErrorMessage[0]);
		THROW(SystemError, &ProgramErrorMessage[0]);
	}

	glDetachShader(programID, VertexShaderID);
	glDetachShader(programID, GeometryShaderID);
	glDetachShader(programID, FragmentShaderID);
}

ShaderProgram::ShaderProgram(const Shader& vertexShader, const Shader& fragmentShader) {
	GLint Result = GL_FALSE;
	int InfoLogLength;

	GLuint VertexShaderID = vertexShader.getID();
	GLuint FragmentShaderID = fragmentShader.getID();

	programID = GLCommand::createProgram();

	glAttachShader(programID, VertexShaderID);
	glAttachShader(programID, FragmentShaderID);

	glLinkProgram(programID);
	glGetProgramiv(programID, GL_LINK_STATUS, &Result);
	glGetProgramiv(programID, GL_INFO_LOG_LENGTH, &InfoLogLength);

	if (InfoLogLength > 0) {
		vector<char> ProgramErrorMessage(InfoLogLength + 1);
		glGetProgramInfoLog(programID, InfoLogLength, nullptr, &ProgramErrorMessage[0]);
		THROW(SystemError, &ProgramErrorMessage[0]);
	}

	glDetachShader(programID, VertexShaderID);
	glDetachShader(programID, FragmentShaderID);
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

GLuint ShaderProgram::getID() const {
	return programID;
}

GLuint ShaderProgram::getUniformLocation(const string& uniformName) {
	if (cachedUniformLocations.contains(uniformName))
		return cachedUniformLocations.at(uniformName);
	GLuint uniformLocation = GLCommand::getUniformLocation(uniformName, programID);
	cachedUniformLocations[uniformName] = uniformLocation;
	return uniformLocation;
}
