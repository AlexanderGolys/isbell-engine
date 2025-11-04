#include "../include/renderer/renderLayer.hpp"

RenderLayer::RenderLayer(const shared_ptr<ShaderProgram>& shaderProgram, const shared_ptr<VertexArray>& vertexArray, const vector<shared_ptr<Uniform>>& uniforms): shaderProgram(shaderProgram), vertexArray(vertexArray), uniforms(uniforms) {}


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

DynamicMeshLayer::DynamicMeshLayer(const shared_ptr<ShaderProgram>& shaderProgram, const shared_ptr<DynamicMeshInterface>& mesh, const VertexBufferLayout& layout,
	const vector<shared_ptr<Uniform>>& uniforms): RenderLayer(shaderProgram, nullptr, uniforms), mesh(mesh), layout(layout), vertexBuffer(nullptr) {
	vertexBuffer = make_shared<VertexBuffer>(layout);
	auto indexBuffer = make_shared<IndexBuffer>(mesh->indexDataSize());
	vertexBuffer->uploadData(mesh->getVertexData(), mesh->vertexDataSize());
	indexBuffer->uploadData(mesh->getIndexData(), mesh->indexDataSize());
	vertexArray = make_shared<VertexArray>(indexBuffer, vector{vertexBuffer});
}

void DynamicMeshLayer::customUpdate(float time, float dt) {
	if (mesh->isDirty()) {
		vertexBuffer->updateData(mesh->getVertexData(), mesh->vertexDataSize());
		mesh->markClean();
	}
}
