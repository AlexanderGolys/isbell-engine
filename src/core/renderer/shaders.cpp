#include "../include/renderer/shaders.hpp"

#include "logging.hpp"


Shader::Shader(const Path& file) {
	CodeFileDescriptor fileDesc(file, false);
	string extension = fileDesc.extension();
	string code = fileDesc.getCode();

	if (extension[0] == '.')
		extension = extension.substr(1);

	switch (extension) {
	case "vert":
		shaderType = GL_VERTEX_SHADER;
		break;
	case "frag":
		shaderType = GL_FRAGMENT_SHADER;
		break;
	case "geom":
	case "geo":
		shaderType = GL_GEOMETRY_SHADER;
		break;
	case "comp":
		shaderType = GL_COMPUTE_SHADER;
		break;
	default: THROW(UnknownVariantError, "Unknown shader extension: ." + extension);
	}

	shaderID = glCreateShader(shaderType);
	GLint Result = GL_FALSE;
	int InfoLogLength;
	char const* sourcePtr = code.c_str();

	glShaderSource(shaderID, 1, &sourcePtr, nullptr);
	glCompileShader(shaderID);
	glGetShaderiv(shaderID, GL_COMPILE_STATUS, &Result);
	glGetShaderiv(shaderID, GL_INFO_LOG_LENGTH, &InfoLogLength);
	if (InfoLogLength > 0) {
		vector<char> error(InfoLogLength + 1);
		glGetShaderInfoLog(shaderID, InfoLogLength, nullptr, &error[0]);
		THROW(SystemError, &error[0]);
	}
}

Shader::~Shader() {
	glDeleteShader(shaderID);
}

GLuint Shader::getID() const {
	return shaderID;
}

GLenum Shader::getType() const {
	return shaderType;
}

GLint ShaderProgram::getUniformLocation(const string& name) {
	if (uniformLocations.contains(name))
		return uniformLocations[name];
	GLint location = glGetUniformLocation(getProgramID(), name.c_str());
	uniformLocations[name] = location;
	return location;
}

ClassicShaderProgram::ClassicShaderProgram(const shared_ptr<Shader>& vertexShader, const shared_ptr<Shader>& fragmentShader) {
	this->vertexShader = vertexShader;
	THROW_IF(vertexShader == nullptr, SystemError, "Vertex shader cannot be null");
	THROW_IF(vertexShader->getType() != GL_VERTEX_SHADER, SystemError, "Provided vertex shader is not of type GL_VERTEX_SHADER");
	this->fragmentShader = fragmentShader;
	THROW_IF(fragmentShader == nullptr, SystemError, "Fragment shader cannot be null");
	THROW_IF(fragmentShader->getType() != GL_FRAGMENT_SHADER, SystemError, "Provided fragment shader is not of type GL_FRAGMENT_SHADER");

	GLint Result = GL_FALSE;
	int InfoLogLength;

	GLuint VertexShaderID = vertexShader->getID();
	GLuint FragmentShaderID = fragmentShader->getID();

	programID = glCreateProgram();

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
	glDetachShader(this->programID, vertexShader->getID());
	glDetachShader(this->programID, fragmentShader->getID());
}

ClassicShaderProgram::~ClassicShaderProgram() {
	glDeleteProgram(this->programID);
}

GLuint ClassicShaderProgram::getProgramID() const {
	return programID;
}

void ClassicShaderProgram::bind() const {
	glUseProgram(programID);
}

void ClassicShaderProgram::unbind() const {
	glUseProgram(0);
}

GeometryShaderProgram::GeometryShaderProgram(const shared_ptr<Shader>& vertexShader, const shared_ptr<Shader>& fragmentShader, const shared_ptr<Shader>& geometryShader) {
	this->vertexShader = vertexShader;
	THROW_IF(vertexShader == nullptr, SystemError, "Vertex shader cannot be null");
	THROW_IF(vertexShader->getType() != GL_VERTEX_SHADER, SystemError, "Provided vertex shader is not of type GL_VERTEX_SHADER");
	this->fragmentShader = fragmentShader;
	THROW_IF(fragmentShader == nullptr, SystemError, "Fragment shader cannot be null");
	THROW_IF(fragmentShader->getType() != GL_FRAGMENT_SHADER, SystemError, "Provided fragment shader is not of type GL_FRAGMENT_SHADER");
	this->geometryShader = geometryShader;
	THROW_IF(geometryShader == nullptr, SystemError, "Geometry shader cannot be null");
	THROW_IF(geometryShader->getType() != GL_GEOMETRY_SHADER, SystemError, "Provided geometry shader is not of type GL_GEOMETRY_SHADER");

	GLint Result = GL_FALSE;
	int InfoLogLength;

	GLuint VertexShaderID = vertexShader->getID();
	GLuint FragmentShaderID = fragmentShader->getID();
	GLuint GeometryShaderID = geometryShader->getID();

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
	glDetachShader(this->programID, vertexShader->getID());
	glDetachShader(this->programID, fragmentShader->getID());
	glDetachShader(this->programID, geometryShader->getID());
}

GeometryShaderProgram::~GeometryShaderProgram() {
	glDeleteProgram(this->programID);
}

GLuint GeometryShaderProgram::getProgramID() const {
	return programID;
}

void GeometryShaderProgram::bind() const {
	glUseProgram(programID);
}

void GeometryShaderProgram::unbind() const {
	glUseProgram(0);
}

ComputeShaderProgram::ComputeShaderProgram(const shared_ptr<Shader>& computeShader) {
	this->computeShader = computeShader;
	THROW_IF(computeShader == nullptr, SystemError, "Compute shader cannot be null");
	THROW_IF(computeShader->getType() != GL_COMPUTE_SHADER, SystemError, "Provided compute shader is not of type GL_COMPUTE_SHADER");

	GLint Result = GL_FALSE;
	int InfoLogLength;

	GLuint ComputeShaderID = computeShader->getID();

	programID = glCreateProgram();

	glAttachShader(programID, ComputeShaderID);

	glLinkProgram(programID);
	glGetProgramiv(programID, GL_LINK_STATUS, &Result);
	glGetProgramiv(programID, GL_INFO_LOG_LENGTH, &InfoLogLength);

	if (InfoLogLength > 0) {
		vector<char> ProgramErrorMessage(InfoLogLength + 1);
		glGetProgramInfoLog(programID, InfoLogLength, nullptr, &ProgramErrorMessage[0]);
		THROW(SystemError, &ProgramErrorMessage[0]);
	}
}

ComputeShaderProgram::~ComputeShaderProgram() {
	glDetachShader(this->programID, computeShader->getID());
	glDeleteProgram(this->programID);
}

GLuint ComputeShaderProgram::getProgramID() const {
	return programID;
}

void ComputeShaderProgram::bind() const {
	glUseProgram(programID);
}

void ComputeShaderProgram::unbind() const {
	glUseProgram(0);
}

void ComputeShaderProgram::run(int numGroupsX, int numGroupsY, int numGroupsZ) const {
	bind();
	glDispatchCompute(numGroupsX, numGroupsY, numGroupsZ);
	glMemoryBarrier(GL_SHADER_IMAGE_ACCESS_BARRIER_BIT);
	glFinish();
}
