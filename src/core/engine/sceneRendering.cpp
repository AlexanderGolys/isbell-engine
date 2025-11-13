#include "sceneRendering.hpp"
#include "glm/gtc/matrix_transform.hpp"


CameraSettings::CameraSettings(float fov_x, float aspectRatio, float clippingRangeMin, float clippingRangeMax)
: fov_x(fov_x), aspectRatio(aspectRatio), clippingRangeMin(clippingRangeMin), clippingRangeMax(clippingRangeMax) {
	projMatrix = glm::perspective(fov_x, aspectRatio, clippingRangeMin, clippingRangeMax);
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

mat4 Camera::vp() const {
	return p() * v();
}


mat4 Camera::mvp(const mat4& modelTransform) const {
	return vp() * modelTransform;
}

CameraUniform::CameraUniform(sptr<Camera> camera)
:	camera(camera),
	vp_uni("u_mvp", camera->vp()),
	campos_uni("u_camPosition", camera->get_position()) {}

void CameraUniform::setDuringRender() const {
	vp_uni.setDuringRender();
	campos_uni.setDuringRender();
}

void CameraUniform::update(float t, float dt) {
	vp_uni.set_value(camera->vp());
	campos_uni.set_value(camera->get_position());
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

void MaterialModelPhongComponent::update(float t, float dt) {
	intensitiesUni.set_value(materialModel->getIntensities());
}

void MaterialModelPhongComponent::setDuringRender() const {
	ambientTexture.setDuringRender();
	diffuseTexture.setDuringRender();
	specularTexture.setDuringRender();
	intensitiesUni.setDuringRender();
}

PointLightsComponent::PointLightsComponent(const vector<std::shared_ptr<SimplePointLight>>& lights)
:	lights(lights),
	lightsPositionsUni("u_lightPos", vector(lights.size(), vec3(0))),
	lightsColorsUni("u_lightColor", vector(lights.size(), vec4(0))),
	lightsIntensitiesUni("u_lightIntencities", vector(lights.size(), vec3(0))) {}

void PointLightsComponent::update(float t, float dt) {
	for (array_index i = 0; i < lights.size(); i++) {
		if (not lights[i]->isDirty())
			continue;
		lightsPositionsUni.set_value(i, lights[i]->get_position());
		lightsColorsUni.set_value(i, lights[i]->get_color());
		lightsIntensitiesUni.set_value(i, lights[i]->get_intensities());
		lights[i]->markClean();
	}
}

void PointLightsComponent::setDuringRender() const {
	lightsPositionsUni.setDuringRender();
	lightsColorsUni.setDuringRender();
	lightsIntensitiesUni.setDuringRender();
}

Mesh3DLayer::Mesh3DLayer(sptr<ShaderProgram> shader, sptr<IndexedMesh3D> mesh, sptr<Camera> camera, sptr<MaterialModelPhong> materialModel, const vector<sptr<SimplePointLight>>& lights)
	: DrawLayer(shader, make_shared<VertexArray>(), {}), mesh(mesh)
{
	addComponent(make_shared<MaterialModelPhongComponent>(materialModel));
	addComponent(make_shared<PointLightsComponent>(lights));
	addComponent(make_shared<TimeUniform>());
	addComponent(make_shared<CameraUniform>(camera));
}

void Mesh3DLayer::init() {
	DrawLayer::init();
	get_vao()->loadIndexedMesh(mesh);
}

void GenericMeshLayer::customRenderStep() {
	get_vao()->updateGeometricData(mesh);
}

GenericMeshLayer::GenericMeshLayer(sptr<ShaderProgram> shader, sptr<GeometricData> mesh)
: DrawLayer(shader, make_shared<VertexArray>(), {}), mesh(mesh) {}

void GenericMeshLayer::init() {
	DrawLayer::init();
	get_vao()->loadGeometricData(mesh);
}
