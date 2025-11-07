#include "renderLayers.hpp"


IndexedDrawLayer::IndexedDrawLayer(sptr<ShaderProgram> shader, sptr<VertexArray> vao, const vector<sptr<LayerComponent>>& components)
: shader(shader), vao(vao), components(components) {}

void IndexedDrawLayer::addComponent(sptr<LayerComponent> comp) {
	components.push_back(comp);
}

void IndexedDrawLayer::init() {
	for (auto& comp : components)
		comp->init();
}

void IndexedDrawLayer::renderStep() {
	vao->bind();
	shader->bind();
	for (auto& uniform : components)
		uniform->setDuringRender();
	vao->draw();
	vao->unbind();
	shader->unbind();
}

void IndexedDrawLayer::updateStep(float t, float dt) {
	vao->bind();
	shader->bind();
	for (auto& comp : components)
		comp->update(t, dt);
	shader->unbind();
	vao->unbind();
}
