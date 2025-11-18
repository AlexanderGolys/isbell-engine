#include "glCommand.hpp"

#include "formatters.hpp"
#include "logging.hpp"

GLuint GLCommand::currentVAO = 0;
bool GLFWCommand::initialized = false;

void GLCommand::bindProgram(gl_id program) {
	glUseProgram(program);
}

void GLCommand::unbindProgram() {
	bindProgram(0);
}

string GLCommand::formatGLenum(GLenum value) {
	switch (value) {
		case GL_VERTEX_SHADER: return "VERTEX SHADER";
		case GL_FRAGMENT_SHADER: return "FRAGMENT SHADER";
		case GL_GEOMETRY_SHADER: return "GEOMETRY SHADER";
		case GL_COMPUTE_SHADER: return "COMPUTE SHADER";
	}
	return "UNKNOWN";
}

void GLCommand::createShader(gl_id& shaderID, GLenum shaderType) {
	shaderID = glCreateShader(shaderType);
}

void GLCommand::compileShaderSource(gl_id shaderID, const string& source) {
	char const* sourcePtr = source.c_str();
	glShaderSource(shaderID, 1, &sourcePtr, nullptr);
	glCompileShader(shaderID);

	int InfoLogLength;
	glGetShaderiv(shaderID, GL_INFO_LOG_LENGTH, &InfoLogLength);
	if (InfoLogLength > 0) {
		vector<char> error(InfoLogLength + 1);
		glGetShaderInfoLog(shaderID, InfoLogLength, nullptr, &error[0]);
		THROW(SystemError, &error[0]);
	}
}

void GLCommand::deleteShader(gl_id& shaderID) {
	glDeleteShader(shaderID);
}

void GLCommand::createProgram(gl_id& shaderID) {
	shaderID = glCreateProgram();
}

void GLCommand::attachShaderToProgram(gl_id programID, gl_id shaderID) {
	glAttachShader(programID, shaderID);
}

void GLCommand::linkProgram(gl_id programID) {
	glLinkProgram(programID);
	GLint program_linked;
	glGetProgramiv(programID, GL_LINK_STATUS, &program_linked);
	if (program_linked != GL_TRUE)
	{
		GLsizei log_length = 0;
		GLchar message[1024];
		glGetProgramInfoLog(programID, 1024, &log_length, message);
		THROW(SystemError, "Error linking program: " + string(message));
	}
	LOG(2, "Shader program linked successfully.");
}

void GLCommand::deleteProgram(gl_id& programID) {
	glDeleteProgram(programID);
	programID = 0;
}

void GLCommand::bindVAO(gl_id vao) {
	glBindVertexArray(vao);
	currentVAO = vao;
}

void GLCommand::unbindVAO() {
	glBindVertexArray(0);
	currentVAO = 0;
}

void GLCommand::deleteVAO(gl_id vao) {
	if (currentVAO == vao)
		unbindVAO();

	glDeleteVertexArrays(1, &vao);
}

void GLCommand::createVAO(gl_id& voaID) {
	glCreateVertexArrays(1, &voaID);
}

void GLCommand::createTexture(gl_id& id) {
	glGenTextures(1, &id);
}

void GLCommand::deleteTexture(gl_id& id) {
	glDeleteTextures(1, &id);
	id = 0;
}

void GLCommand::bindTexture2D(gl_id id, uint slot) {
	THROW_IF(slot >= 32, ValueError, "Texture slot out of range (0-31).");
	glActiveTexture(GL_TEXTURE0 + slot);
	glBindTexture(GL_TEXTURE_2D, id);
}

void GLCommand::unbindTexture2D() {
	glBindTexture(GL_TEXTURE_2D, 0);
}

void GLCommand::setTexture2DFilters(GLenum minFilter, GLenum magFilter, GLenum wrapS, GLenum wrapT) {
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, minFilter);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, magFilter);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, wrapS);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, wrapT);
}

void GLCommand::calculateMipmapsTexture2D() {
	glGenerateMipmap(GL_TEXTURE_2D);
}

void GLCommand::loadTexture2(GLenum internalFormat, array_len width, array_len height, GLenum format, raw_data_ptr<uchar> data) {
	glTexImage2D(GL_TEXTURE_2D, 0, internalFormat, width, height, 0, format, GL_UNSIGNED_BYTE, data);
	LOG(1, "Texture2D loaded with size " + to_string(width) + "x" + to_string(height) + ".");
}

void GLCommand::setSampler2D(const string& samplerName, uint slot, gl_id programID) {
	gl_uniform_loc location = glGetUniformLocation(programID, samplerName.c_str());
	glUniform1i(location, slot);
}


void GLCommand::createBuffer(gl_id& id) {
	glCreateBuffers(1, &id);
}

void GLCommand::deleteBuffer(gl_id& id) {
	glDeleteBuffers(1, &id);
	id = 0;
}

void GLCommand::loadBufferData(gl_id bufferID, raw_data_ptr firstElementAdress, byte_size bufferSize, GLenum usage) {
	glNamedBufferData(bufferID, bufferSize, firstElementAdress, usage);
}

void GLCommand::updateBufferData(gl_id bufferID, raw_data_ptr firstElementAdress, byte_size bufferSize) {
	glNamedBufferSubData(bufferID, 0, bufferSize, firstElementAdress);
}

void GLCommand::bindSSBO(gl_id bufferID, int bindingPoint) {
	glBindBuffer(GL_SHADER_STORAGE_BUFFER, bufferID);
	glBindBufferBase(GL_SHADER_STORAGE_BUFFER, bindingPoint, bufferID);
}

void GLCommand::unbindSSBO() {
	glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);
}

void GLCommand::loadSSBOData(gl_id bufferID, raw_data_ptr firstElementAdress, byte_size bufferSize) {
	glBindBuffer(GL_SHADER_STORAGE_BUFFER, bufferID);
	glBufferData(GL_SHADER_STORAGE_BUFFER, bufferSize, firstElementAdress, GL_DYNAMIC_DRAW);
}

void GLCommand::updateSSBOData(gl_id bufferID, raw_data_ptr firstElementAdress, byte_size bufferSize) {
	glBindBuffer(GL_SHADER_STORAGE_BUFFER, bufferID);
	glBufferSubData(GL_SHADER_STORAGE_BUFFER, 0, bufferSize, firstElementAdress);
}

void GLCommand::initUBO(gl_id bufferID, int bindingPoint, const string& blockName, gl_id currentProgram) {
	gl_id blockIndex = glGetUniformBlockIndex(currentProgram, blockName.c_str());
	glUniformBlockBinding(currentProgram, blockIndex, bindingPoint);

	glBindBuffer(GL_UNIFORM_BUFFER, bufferID);
	glBindBufferBase(GL_UNIFORM_BUFFER, bindingPoint, bufferID);
}

void GLCommand::unbindUBO() {
	glBindBuffer(GL_UNIFORM_BUFFER, 0);
}




void GLCommand::loadUBOData(gl_id bufferID, raw_data_ptr firstElementAdress, byte_size bufferSize) {
	glBindBuffer(GL_UNIFORM_BUFFER, bufferID);
	glBufferData(GL_UNIFORM_BUFFER, bufferSize, firstElementAdress, GL_DYNAMIC_DRAW);
}

void GLCommand::updateUBOData(gl_id bufferID, raw_data_ptr firstElementAdress, byte_size bufferSize) {
	glBindBuffer(GL_UNIFORM_BUFFER, bufferID);
	glBufferSubData(GL_UNIFORM_BUFFER, 0, bufferSize, firstElementAdress);
}


void GLCommand::drawIndexedTriangles(array_len numberOfIndices) {
	THROW_IF(currentVAO == 0, RuntimeError, "No VAO is currently bound.");
	THROW_IF(numberOfIndices == 0, ValueError, "Trying to draw 0 indices.");
	glDrawElements(GL_TRIANGLES, static_cast<GLsizei>(numberOfIndices), GL_UNSIGNED_INT, nullptr);
}

void GLCommand::runComputeShader(int numGroupsX, int numGroupsY, int numGroupsZ) {
	glDispatchCompute(numGroupsX, numGroupsY, numGroupsZ);
	glMemoryBarrier(GL_SHADER_IMAGE_ACCESS_BARRIER_BIT);
	glFinish();
}

void GLCommand::pointAtElementBuffer(gl_id buffer) {
	THROW_IF(currentVAO == 0, RuntimeError, "No VAO is currently bound.");
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, buffer);
}

void GLCommand::pointAtAttributeBuffer(int inputNumber, gl_id bufferAddress, GLSLPrimitive type, byte_size stride,  raw_data_ptr offset) {
	THROW_IF(currentVAO == 0, RuntimeError, "No VAO is currently bound.");
	THROW_IF(not supportedAttributeType(type), ValueError, "Unsupported GLSLType for attribute buffer (supported: float/double or vectors of these).");
	glBindBuffer(GL_ARRAY_BUFFER, bufferAddress);
	if (primitiveTypeEnum(type) == GL_DOUBLE)
		glVertexAttribLPointer(inputNumber, lengthOfGLSLType(type), GL_DOUBLE, stride, offset);
	else
		glVertexAttribPointer(inputNumber, lengthOfGLSLType(type), GL_FLOAT, GL_FALSE, stride, offset);
	glEnableVertexAttribArray(inputNumber);
}

void GLCommand::clearScreen(const vec4& color) {
	glClearColor(color.x, color.y, color.z, color.w);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
}

void GLCommand::enableBlending() {
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
}

void GLCommand::enableDepth() {
	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LESS);
}

void GLCommand::init() {
	glewExperimental = true;
	if (glewInit() != GLEW_OK)
		THROW(SystemError, "GLEW initialization failed.");
	LOG("GLEW initialized.");
}

void GLCommand::initViewport(int width, int height) {
	glViewport(0, 0, width,  height);
}

void GLFWCommand::init() {
	if (not glfwInit())
		THROW(SystemError, "GLFW initialization failed.");

	glfwWindowHint(GLFW_SAMPLES, 4);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
	glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);

	LOG("GLFW initialized.");
	initialized = true;
}

GLFWwindow* GLFWCommand::createWindow(int width, int height, const string& title) {
	GLFWwindow* w = glfwCreateWindow(width, height, title.c_str(), nullptr, nullptr);
	THROW_IF(not w, SystemError, "GLFW window creation failed.");
	glfwMakeContextCurrent(w);
	return w;
}

void GLFWCommand::destroyWindow(GLFWwindow* window) {
	glfwDestroyWindow(window);
}

void GLFWCommand::pollEvents() {
	glfwPollEvents();
}

void GLFWCommand::swapFrameBuffers(GLFWwindow* window) {
	glfwSwapBuffers(window);
}

void GLFWCommand::setWindowData(GLFWwindow* window, raw_data_mut data) {
	glfwSetWindowUserPointer(window, data);
}

WindowSettings& GLFWCommand::getWindowData(GLFWwindow* window) {
	return *(static_cast<WindowSettings*>(glfwGetWindowUserPointer(window)));
}

vec2 GLFWCommand::getCursorPosition(GLFWwindow* window) {
	double xPos, yPos;
	glfwGetCursorPos(window, &xPos, &yPos);
	return vec2(xPos, yPos);
}

bool GLFWCommand::isWindowOpen(GLFWwindow* window) {
	if (not window)
		return false;
	auto winData = getWindowData(window);
	return winData.open and not glfwWindowShouldClose(window);
}

void GLFWCommand::terminate() {
	glfwTerminate();
	initialized = false;
}


