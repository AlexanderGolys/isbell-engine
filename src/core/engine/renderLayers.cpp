#include "renderLayers.hpp"


void Layer::onEvent(sptr<Event> event, TimeStep timeStep) const {
	for (auto& listener : eventListeners)
		if (listener->listensToEventType(event->getEventType()))
			listener->onEvent(event, timeStep);
}

void Layer::addEventListener(sptr<EventListener> listener) {
	eventListeners.push_back(listener);
}

CombinedLayerComponent::CombinedLayerComponent(const vector<sptr<LayerComponent>>& components): components(components) {}

void CombinedLayerComponent::addComponent(sptr<LayerComponent> comp) { components.push_back(comp); }

void CombinedLayerComponent::init() {
	for (auto& comp : components)
		comp->init();
}

void CombinedLayerComponent::update(TimeStep timeStep) {
	for (auto& comp : components)
		comp->update(timeStep);
}

void CombinedLayerComponent::setDuringRender() {
	for (auto& comp : components)
		comp->setDuringRender();
}

DrawLayer::DrawLayer(sptr<ShaderProgram> shader, sptr<VertexArray> vao, const vector<sptr<LayerComponent>>& components)
: vao(vao), shader(shader), components(components) {}

void DrawLayer::addComponent(sptr<LayerComponent> comp) {
	components.push_back(comp);
}

void DrawLayer::init() {
	for (auto& comp : components)
		comp->init();

}

void DrawLayer::renderStep() {
	if (not prerenderStep())
		return;
	vao->bind();
	shader->bind();
	for (auto& uniform : components)
		uniform->setDuringRender();
	customRenderStep();
	vao->draw();
	shader->unbind();
	vao->unbind();
}

void DrawLayer::updateStep(TimeStep timeStep) {
	for (auto& comp : components)
		comp->update(timeStep);
}

GenericMeshLayer::GenericMeshLayer(sptr<ShaderProgram> shader, sptr<GeometricData> mesh)
: DrawLayer(shader, make_shared<VertexArray>(), {}), mesh(mesh) {}

void GenericMeshLayer::init() {
	DrawLayer::init();
	get_vao()->loadGeometricData(mesh);
}

void GenericMeshLayer::customRenderStep() {
	get_vao()->updateGeometricData(mesh);
}
