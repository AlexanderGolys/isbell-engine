#pragma once
#include "openglRenderer.hpp"

namespace openglAPI {

	class ComputeShaderProgram : public ShaderProgram {
		shared_ptr<Shader> computeShader;
	public:
		ComputeShaderProgram(const shared_ptr<Shader> &computeShader);
		void linkComputeShader();

		void run(ivec3 numGroups) {
			glDispatchCompute(numGroups.x, numGroups.y, numGroups.z);
			glMemoryBarrier(GL_SHADER_IMAGE_ACCESS_BARRIER_BIT);
		}
	};

}
