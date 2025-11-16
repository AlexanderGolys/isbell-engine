#include "sceneRendering.hpp"
#include "glm/gtc/matrix_transform.hpp"


CameraSettings::CameraSettings(float fov_x, float aspectRatio, float clippingRangeMin, float clippingRangeMax)
: fov_x(fov_x), aspectRatio(aspectRatio), clippingRangeMin(clippingRangeMin), clippingRangeMax(clippingRangeMax) {
	projMatrix = glm::perspective(fov_x, aspectRatio, clippingRangeMin, clippingRangeMax);
}

CameraTransform::CameraTransform(vec3 position, vec3 lookAtPos, vec3 upVector):	position(position), lookAtPos(lookAtPos), upVector(upVector) {}

vec3 CameraTransform::getDirection() const {
	return normalize(lookAtPos - position);
}

vec3 CameraTransform::getRightVector() const {
	return normalize(cross(getDirection(), upVector));
}

Camera::Camera(const CameraSettings& settings, vec3 position, vec3 lookAtPos, vec3 upVector)
: position(position), lookAtPos(lookAtPos), upVector(upVector), settings(settings){
}

mat4 Camera::v() const {
	return glm::lookAt(position, lookAtPos, upVector);
}

mat4 Camera::p() const {
	return settings.projMatrix;
}

void Camera::set_transform(vec3 newPosition, vec3 newLookAt, vec3 newUp) {
	set_position(newPosition);
	set_lookAtPos(newLookAt);
	set_upVector(newUp);
}

void Camera::set_transform(const CameraTransform& transform) {
	set_transform(transform.position, transform.lookAtPos, transform.upVector);
}

CameraTransform Camera::get_transform() const {
	return CameraTransform(position, lookAtPos, upVector);
}

mat4 Camera::vp() const {
	return p() * v();
}


mat4 Camera::mvp(const mat4& modelTransform) const {
	return vp() * modelTransform;
}

CameraUniform::CameraUniform(sptr<Camera> camera, sptr<Object3D> attachedObject){
	addComponent(make_shared<DynamicPrimitiveUniform<mat4>>("u_mvp", [camera, attachedObject](TimeStep){
		if (attachedObject)
			return camera->mvp(attachedObject->get_transform().toMat4());
		return camera->vp();
	}));
	addComponent(make_shared<DynamicPrimitiveUniform<vec3>>("u_cameraPos", [camera](TimeStep){
		return camera->get_position();
	}));
}

MaterialModelPhong::MaterialModelPhong(sptr<TextureData> ambientTextureData, sptr<TextureData> diffuseTextureData, sptr<TextureData> specularTextureData, float ambientIntensity, float diffuseIntensity, float specularIntensity, float shininess)
:	ambientTextureData(ambientTextureData),
	diffuseTextureData(diffuseTextureData),
	specularTextureData(specularTextureData),
	ambientIntensity(ambientIntensity),
	diffuseIntensity(diffuseIntensity),
	specularIntensity(specularIntensity),
	shininess(shininess) {}

MaterialModelPhong::MaterialModelPhong(Color ambientColor, Color diffuseColor, Color specularColor, float ambientIntensity, float diffuseIntensity, float specularIntensity, float shininess)
:	ambientTextureData(make_shared<ConstColorTextureData>(ambientColor)),
	diffuseTextureData(make_shared<ConstColorTextureData>(diffuseColor)),
	specularTextureData(make_shared<ConstColorTextureData>(specularColor)),
	ambientIntensity(ambientIntensity),
	diffuseIntensity(diffuseIntensity),
	specularIntensity(specularIntensity),
	shininess(shininess) {}

sptr<TextureData> MaterialModelPhong::getAmbientTextureData() const {
	return ambientTextureData;
}

sptr<TextureData> MaterialModelPhong::getDiffuseTextureData() const {
	return diffuseTextureData;
}

sptr<TextureData> MaterialModelPhong::getSpecularTextureData() const {
	return specularTextureData;
}

vec4 MaterialModelPhong::getIntensities() const {
	return vec4(ambientIntensity, diffuseIntensity, specularIntensity, shininess);
}

MaterialModelPhongComponent::MaterialModelPhongComponent(sptr<MaterialModelPhong> materialModel):
materialModel(materialModel),
ambientTexture(materialModel->getAmbientTextureData(), 0, "u_texture_ambient"),
diffuseTexture(materialModel->getDiffuseTextureData(), 1, "u_texture_diffuse"),
specularTexture(materialModel->getSpecularTextureData(), 2, "u_texture_specular"),
intensitiesUni("u_intencities", materialModel->getIntensities()) {}

void MaterialModelPhongComponent::init() {
	ambientTexture.init();
	diffuseTexture.init();
	specularTexture.init();
}

void MaterialModelPhongComponent::update(TimeStep timeStep) {}

void MaterialModelPhongComponent::setDuringRender() {
	ambientTexture.setDuringRender();
	diffuseTexture.setDuringRender();
	specularTexture.setDuringRender();
	intensitiesUni.setDuringRender();
}

SimplePointLight::SimplePointLight(vec3 position, Color color, vec3 intensities):
position(position),
color(color),
intensities(intensities) {}

raw_data_ptr SimplePointLight::data() const {
	return &position[0];
}

byte_size SimplePointLight::byteSize() const {
	return 3*sizeof(vec4);
}


PointLightsComponent::PointLightsComponent(sptr<arrayStruct<SimplePointLight>> lights): UniformBlockComponent(*lights, 0), lights(lights) {}

void PointLightsComponent::update(TimeStep timeStep) {
	if (lights->isDirty()) {
		setValue(*lights);
		lights->markClean();
	}
}

OldMesh3DLayer::OldMesh3DLayer(sptr<ShaderProgram> shader, sptr<IndexedMesh3D> mesh, sptr<Camera> camera, sptr<MaterialModelPhong> materialModel, sptr<arrayStruct<SimplePointLight>> lights)
	: DrawLayer(shader, make_shared<VertexArray>(), {}), mesh(mesh)
{
	addComponent(make_shared<MaterialModelPhongComponent>(materialModel));
	addComponent(make_shared<PointLightsComponent>(lights));
	addComponent(make_shared<DynamicPrimitiveUniform<float>>("u_time", [](TimeStep t){ return t.t; }));
	addComponent(make_shared<CameraUniform>(camera));
}

void OldMesh3DLayer::init() {
	DrawLayer::init();
	get_vao()->loadIndexedMesh(mesh);
}

Mesh3DLayer::Mesh3DLayer(sptr<ShaderProgram> shader, sptr<Mesh3D> mesh, sptr<Camera> camera, sptr<MaterialModelPhong> materialModel)
: GenericMeshLayer(shader, mesh) {
	addComponent(make_shared<MaterialModelPhongComponent>(materialModel));
	addComponent(make_shared<CameraUniform>(camera, mesh));
	addComponent(make_shared<DynamicPrimitiveUniform<float>>("u_time", [](TimeStep t){ return t.t; }));
}

Scene3DLayer::Scene3DLayer(sptr<Camera> camera, sptr<arrayStruct<SimplePointLight>> lights)
: pointLightsComponent(make_shared<PointLightsComponent>(lights)), camera(camera) {}

void Scene3DLayer::addMeshLayer(sptr<ShaderProgram> shader, sptr<Mesh3D> mesh, sptr<MaterialModelPhong> materialModel) {
	meshLayers.emplace_back(shader, mesh, camera, materialModel);
	meshLayers.back().addComponent(pointLightsComponent);
}

void Scene3DLayer::init() {
	for (auto& layer : meshLayers)
		layer.init();
}

void Scene3DLayer::renderStep() {
	for (auto& layer : meshLayers)
		layer.renderStep();
}

void Scene3DLayer::updateStep(TimeStep timeStep) {
	for (auto& layer : meshLayers)
		layer.updateStep(timeStep);
	for (auto& comp : updateComponents)
		comp->update(timeStep);
}

void Scene3DLayer::onEvent(sptr<Event> event, TimeStep timeStep) const {
	for (auto& layer : meshLayers)
		layer.onEvent(event, timeStep);
	Layer::onEvent(event, timeStep);
}

void Scene3DLayer::addComponentToEachMesh(sptr<LayerComponent> comp) {
	for (auto& layer : meshLayers)
		layer.addComponent(comp);
}

void Scene3DLayer::addEventListenerToEachMesh(sptr<EventListener> listener) {
	for (auto& layer : meshLayers)
		layer.addEventListener(listener);
}

void Scene3DLayer::addUpdateComponent(sptr<UpdateComponent> comp) {
	updateComponents.push_back(comp);
}


