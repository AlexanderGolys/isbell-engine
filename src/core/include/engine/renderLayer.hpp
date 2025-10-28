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
	RenderLayer(const shared_ptr<ShaderProgram>& shaderProgram, const shared_ptr<VertexArray>& vertexArray, const vector<shared_ptr<Uniform>>& uniforms)
	: shaderProgram(shaderProgram), vertexArray(vertexArray), uniforms(uniforms) {}

	virtual void customUpdate(float time, float dt) {}

	void addUniform(const shared_ptr<Uniform>& uniform) {
		uniforms.push_back(uniform);
	}

	void onRender() override {
		shaderProgram->bind();
		vertexArray->bind();
		for (const auto& uniform : uniforms)
			shaderProgram->setUniform(*uniform);
		glDrawElements(GL_TRIANGLES, vertexArray->getIndexCount(), GL_UNSIGNED_INT, nullptr);

	}

	void onUpdate(float time, float dt) override {
		for (const auto& uniform : uniforms)
			uniform->update(time, dt);
		customUpdate(time, dt);
	}
};
