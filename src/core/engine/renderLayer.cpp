#include "renderLayer.hpp"

RenderLayer::RenderLayer(const shared_ptr<ShaderProgram>& shaderProgram, const shared_ptr<VertexArray>& vertexArray, const vector<shared_ptr<Uniform>>& uniforms): shaderProgram(shaderProgram), vertexArray(vertexArray), uniforms(uniforms) {}

void RenderLayer::customUpdate(float time, float dt) {}

void RenderLayer::addUniform(const shared_ptr<Uniform>& uniform) {
	uniforms.push_back(uniform);
}

void RenderLayer::onRender() {
	shaderProgram->bind();
	vertexArray->bind();
	for (const auto& uniform : uniforms)
		uniform->set(*shaderProgram);
	glDrawElements(GL_TRIANGLES, vertexArray->getIndexCount(), GL_UNSIGNED_INT, nullptr);

}

void RenderLayer::onUpdate(float time, float dt) {
	for (const auto& uniform : uniforms)
		uniform->update(time, dt);
	customUpdate(time, dt);
}
