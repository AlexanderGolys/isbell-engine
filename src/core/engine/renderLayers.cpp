#include "renderLayers.hpp"


DrawLayer::DrawLayer(sptr<ShaderProgram> shader, sptr<VertexArray> vao, const vector<sptr<LayerComponent>>& components)
: vao(vao), shader(shader), components(components) {}

void DrawLayer::addComponent(sptr<LayerComponent> comp) {
	components.push_back(comp);
}

void DrawLayer::addEventListener(sptr<EventListener> listener) {
	eventListeners.push_back(listener);
}

void DrawLayer::init() {
	for (auto& comp : components)
		comp->init();
	for (auto& listener : eventListeners)
		listener->init();
}

void DrawLayer::renderStep() {
	vao->bind();
	shader->bind();
	for (auto& uniform : components)
		uniform->setDuringRender();
	customRenderStep();
	vao->draw();
	shader->unbind();
	vao->unbind();
}

void DrawLayer::updateStep(float t, float dt) {
	for (auto& comp : components)
		comp->update(t, dt);
}

void DrawLayer::onEvent(const Event& event, float t, float dt) const {
	for (auto& listener : eventListeners)
		if (listener->listensToEventType(event.getEventType()))
			listener->onEvent(event, t, dt);
}
