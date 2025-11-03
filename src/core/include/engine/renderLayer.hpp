#pragma once

#include "buffers.hpp"
#include "shaders.hpp"
#include "layers.hpp"

class RenderLayer : public Layer {
protected:
	shared_ptr<ShaderProgram> shaderProgram;
	shared_ptr<VertexArray> vertexArray;
	vector<shared_ptr<Uniform>> uniforms;

public:
	RenderLayer(const shared_ptr<ShaderProgram>& shaderProgram, const shared_ptr<VertexArray>& vertexArray, const vector<shared_ptr<Uniform>>& uniforms);

	virtual void customUpdate(float time, float dt);
	void addUniform(const shared_ptr<Uniform>& uniform);
	void onRender() override;
	void onUpdate(float time, float dt) override;
};
