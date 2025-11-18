#include "renderLayers.hpp"


void DrawLayer::addComponent(const sptr<LayerComponent>& comp) {
	return components.push_back(comp);
}

void DrawLayer::init() {
	shader->bind();
	for (auto& comp : components)
		comp->initPerShader(shader->get_programID());
	shader->unbind();

}

void DrawLayer::renderStep() {
	vao.bind();
	shader->bind();
	for (auto& uniform : components)
		uniform->setDuringRender();
	customRenderStep();
	vao.draw();
	shader->unbind();
	vao.unbind();
}

void DrawLayer::updateStep(TimeStep timeStep) {
	for (auto& comp : components)
		comp->update(timeStep);
}

void DrawLayer::finalize() {
	for (auto& comp : components)
		comp->finalize();
}

gl_id DrawLayer::get_shaderProgramID() const {
	return shader->get_programID();
}


void GenericMeshLayer::init() {
	DrawLayer::init();
	vao.loadGeometricData(mesh);
}

void GenericMeshLayer::customRenderStep() {
	vao.updateGeometricData(mesh);
}
