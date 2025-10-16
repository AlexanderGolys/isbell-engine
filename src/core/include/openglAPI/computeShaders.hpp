#pragma once
#include "openglRenderer.hpp"

namespace openglAPI {

	class ComputeShaderProgram : public ShaderProgram {
		shared_ptr<Shader> computeShader;
	public:
		ComputeShaderProgram(const shared_ptr<Shader> &computeShader);
		void linkComputeShader();

		void run(int numGroupsX, int numGroupsY=1, int numGroupsZ=1);
	};

}
