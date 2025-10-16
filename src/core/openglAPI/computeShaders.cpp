#include "computeShaders.hpp"

using namespace openglAPI;


ComputeShaderProgram::ComputeShaderProgram(const shared_ptr<Shader> &computeShader): ShaderProgram(), computeShader(computeShader) {
	shaderType = COMPUTE;
	computeShader->compile();
	linkComputeShader();
}

void ComputeShaderProgram::linkComputeShader() {
	programID = glCreateProgram();
	glAttachShader(programID, computeShader->getID());
	glLinkProgram(programID);
}

void ComputeShaderProgram::run(int numGroupsX, int numGroupsY, int numGroupsZ) {
	use();
	glDispatchCompute(numGroupsX, numGroupsY, numGroupsZ);
	glMemoryBarrier(GL_SHADER_IMAGE_ACCESS_BARRIER_BIT);
	glFinish();
}
