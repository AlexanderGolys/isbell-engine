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

CameraUniformData::CameraUniformData(vec3 position, const mat4& vp)
: position(position), vp(vp) {}

raw_data_ptr CameraUniformData::data() const {
	return &position[0];
}

byte_size CameraUniformData::blockSize() const { return sizeof(vec3) + sizeof(mat4); }

Camera::Camera(const CameraSettings& settings, vec3 position, vec3 lookAtPos, vec3 upVector)
: transform(position, lookAtPos, upVector), settings(settings), uniformData(position, vp()) {}

vec3 Camera::get_position() const {
	return transform.position;
}

vec3 Camera::get_lookAtPos() const {
	return transform.lookAtPos;
}

vec3 Camera::get_upVector() const {
	return transform.upVector;
}

void Camera::set_position(const vec3& position) {
	transform.position = position;
	setUniformData();
	markDirty();
}

void Camera::set_lookAtPos(const vec3& lookAtPos) {
	transform.lookAtPos = lookAtPos;
	setUniformData();
	markDirty();
}

void Camera::set_upVector(const vec3& upVector) {
	transform.upVector = upVector;
	setUniformData();
	markDirty();
}

raw_data_ptr Camera::data() const {
	return uniformData.data();
}

byte_size Camera::blockSize() const {
	return sizeof(mat4) + sizeof(vec4);
}

void Camera::set_transform(const CameraTransform& transform) {
	this->transform = transform;
	setUniformData();
	markDirty();
}

const CameraTransform& Camera::get_transform() const {
	return transform;
}

mat4 Camera::v() const {
	return glm::lookAt(get_position(), get_lookAtPos(), get_upVector());
}

mat4 Camera::p() const {
	return settings.projMatrix;
}

void Camera::setUniformData() {
	uniformData.position = transform.position;
	uniformData.vp = vp();
}


mat4 Camera::vp() const {
	return p() * v();
}


MaterialModelPhong::MaterialModelPhong(sptr<TextureData> ambientTextureData, sptr<TextureData> diffuseTextureData, sptr<TextureData> specularTextureData, float ambientIntensity, float diffuseIntensity, float specularIntensity, float shininess)
:	ambientTextureData(ambientTextureData),
	diffuseTextureData(diffuseTextureData),
	specularTextureData(specularTextureData),
	intensities(vec4(ambientIntensity, diffuseIntensity, specularIntensity, shininess))
{}

MaterialModelPhong::MaterialModelPhong(Color ambientColor, Color diffuseColor, Color specularColor, float ambientIntensity, float diffuseIntensity, float specularIntensity, float shininess)
:	ambientTextureData(make_shared<ConstColorTextureData>(ambientColor)),
	diffuseTextureData(make_shared<ConstColorTextureData>(diffuseColor)),
	specularTextureData(make_shared<ConstColorTextureData>(specularColor)),
	intensities(vec4(ambientIntensity, diffuseIntensity, specularIntensity, shininess)){}


sptr<TextureData> MaterialModelPhong::getAmbientTextureData() const {
	return ambientTextureData;
}

sptr<TextureData> MaterialModelPhong::getDiffuseTextureData() const {
	return diffuseTextureData;
}

sptr<TextureData> MaterialModelPhong::getSpecularTextureData() const {
	return specularTextureData;
}



MaterialModelPhongComponent::MaterialModelPhongComponent(sptr<MaterialModelPhong> materialModel):
ambientTexture(materialModel->getAmbientTextureData(), 0, "u_texture_ambient"),
diffuseTexture(materialModel->getDiffuseTextureData(), 1, "u_texture_diffuse"),
specularTexture(materialModel->getSpecularTextureData(), 2, "u_texture_specular"){}

void MaterialModelPhongComponent::init() {
	ambientTexture.init();
	diffuseTexture.init();
	specularTexture.init();
}


void MaterialModelPhongComponent::setDuringRender() {
	ambientTexture.setDuringRender();
	diffuseTexture.setDuringRender();
	specularTexture.setDuringRender();
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



OldMesh3DLayer::OldMesh3DLayer(sptr<ShaderProgram> shader, sptr<IndexedMesh3D> mesh, sptr<Camera> camera,
							   sptr<MaterialModelPhong> materialModel, sptr<arrayStruct<SimplePointLight>> lights)
	: DrawLayer(shader, make_shared<VertexArray>(), {}), mesh(mesh)
{
	addComponent(make_shared<UniformBlock<MaterialModelPhong>>(materialModel, 4, "Material"));
	addComponent(make_shared<UniformArrayBlock<SimplePointLight>>(lights, 0, "Lights"));
	addComponent(make_shared<UniformVecInnerBlock<vec2>>(vec2(0), 3, "Time", [](TimeStep t){ return vec2(t.t, t.dt); }));
	addComponent(make_shared<UniformBlock<Camera>>(camera, 1, "Camera"));
	addComponent(make_shared<MaterialModelPhongComponent>(materialModel));

}

void OldMesh3DLayer::init() {
	DrawLayer::init();
	get_vao()->loadIndexedMesh(mesh);
}











Mesh3DLayer::Mesh3DLayer(sptr<ShaderProgram> shader, sptr<Mesh3D> mesh, sptr<MaterialModelPhong> materialModel)
: GenericMeshLayer(shader, mesh) {
	addComponent(make_shared<MaterialModelPhongComponent>(materialModel));
	addComponent(make_shared<UniformBlock<Mesh3D>>(mesh, 2, "Model"));
	addComponent(make_shared<UniformVecInnerBlock<vec2>>(vec2(0), 3, "Time", [](TimeStep t){ return vec2(t.t, t.dt); }));
	addComponent(make_shared<UniformBlock<MaterialModelPhong>>(materialModel, 4, "Material"));
}




Scene3DLayer::Scene3DLayer(sptr<Camera> camera, sptr<arrayStruct<SimplePointLight>> lights){
	addLocalComponent(make_shared<UniformArrayBlock<SimplePointLight>>(lights, 0, "Lights"));
	addLocalComponent(make_shared<UniformBlock<Camera>>(camera, 1, "Camera"));
}

void Scene3DLayer::addMeshLayer(sptr<ShaderProgram> shader, sptr<Mesh3D> mesh, sptr<MaterialModelPhong> materialModel) {
	meshLayers.emplace_back(shader, mesh, materialModel);
	for (auto& comp : localComponents)
		meshLayers.back().addComponent(comp);
}

void Scene3DLayer::init() {
	for (auto& layer : meshLayers)
		layer.init();
	for (auto& comp : globalComponents)
		comp->init();
}

void Scene3DLayer::renderStep() {
	for (auto& layer : meshLayers)
		layer.renderStep();
}

void Scene3DLayer::updateStep(TimeStep timeStep) {
	for (auto& layer : meshLayers)
		layer.updateStep(timeStep);
	for (auto& comp : globalComponents)
		comp->update(timeStep);
}

void Scene3DLayer::onEvent(sptr<Event> event, TimeStep timeStep) const {
	for (auto& layer : meshLayers)
		layer.onEvent(event, timeStep);
	Layer::onEvent(event, timeStep);
}

void Scene3DLayer::addLocalComponent(sptr<LayerComponent> comp) {
	localComponents.push_back(comp);
	for (auto& layer : meshLayers)
		layer.addComponent(comp);
}

void Scene3DLayer::addEventListenerToEachMesh(sptr<EventListener> listener) {
	for (auto& layer : meshLayers)
		layer.addEventListener(listener);
}

void Scene3DLayer::addGlobalComponent(sptr<UpdateComponent> comp) {
	globalComponents.push_back(comp);
}

void Scene3DLayer::finalize() {
	for (auto& layer : meshLayers)
		layer.finalize();
	for (auto& comp : globalComponents)
		comp->finalize();
}


