#pragma once
#include "colors.hpp"
#include "indexedRendering.hpp"
#include "textures.hpp"
#include "uniforms.hpp"
#include "indexedMesh.hpp"
#include "mesh3D.hpp"

struct CameraSettings {
	float fov_x;
	float aspectRatio;
	float clippingRangeMin;
	float clippingRangeMax;
	mat4 projMatrix;

	explicit CameraSettings(float fov_x = PI/4, float aspectRatio = 16/9.f, float clippingRangeMin = .01f, float clippingRangeMax = 100.f);
};

struct CameraTransform {
	vec3 position;
	vec3 lookAtPos;
	vec3 upVector;

	CameraTransform(vec3 position, vec3 lookAtPos, vec3 upVector = vec3(0,0,1));
	vec3 getDirection() const;
	vec3 getRightVector() const;
};

class Camera {
	DIRTY_FLAG;
	PROPERTY_DIRTY(vec3, position);
	PROPERTY_DIRTY(vec3, lookAtPos);
	PROPERTY_DIRTY(vec3, upVector);
	CameraSettings settings;

public:
	Camera(const CameraSettings& settings, vec3 position, vec3 lookAtPos, vec3 upVector = vec3(0, 0, 1));
	mat4 mvp(const mat4& modelTransform) const;
	mat4 vp() const;
	mat4 v() const;
	mat4 p() const;
	void set_transform(vec3 newPosition, vec3 newLookAt, vec3 newUp);
	void set_transform(const CameraTransform& transform);
	CameraTransform get_transform() const;
};


class CameraUniform : public CombinedLayerComponent {
public:
	explicit CameraUniform(sptr<Camera> camera, sptr<Object3D> attachedObject = nullptr);
};

class MaterialModelPhong {
	sptr<TextureData> ambientTextureData;
	sptr<TextureData> diffuseTextureData;
	sptr<TextureData> specularTextureData;
	float ambientIntensity;
	float diffuseIntensity;
	float specularIntensity;
	float shininess;

public:
	MaterialModelPhong(sptr<TextureData> ambientTextureData, sptr<TextureData> diffuseTextureData, sptr<TextureData> specularTextureData, float ambientIntensity, float diffuseIntensity, float specularIntensity, float shininess);
	MaterialModelPhong(Color ambientColor, Color diffuseColor, Color specularColor, float ambientIntensity, float diffuseIntensity, float specularIntensity, float shininess);

	sptr<TextureData> getAmbientTextureData() const;
	sptr<TextureData> getDiffuseTextureData() const;
	sptr<TextureData> getSpecularTextureData() const;
	vec4 getIntensities() const;
};

class MaterialModelPhongComponent : public LayerComponent {
	sptr<MaterialModelPhong> materialModel;
	Texture2D ambientTexture;
	Texture2D diffuseTexture;
	Texture2D specularTexture;
	ConstPrimitiveUniform<vec4> intensitiesUni;

public:
	explicit MaterialModelPhongComponent(sptr<MaterialModelPhong> materialModel);
	void init() override;
	void update(TimeStep timeStep) override;
	void setDuringRender() override;
};

struct SimplePointLight {
	alignas(16) vec3 position;
	vec4 color;
	alignas(16) vec3 intensities;

	SimplePointLight(vec3 position, Color color, vec3 intensities);
	raw_data_ptr data() const;
	byte_size byteSize() const;
};

class PointLightsComponent : public UniformBlockComponent<arrayStruct<SimplePointLight>> {
	sptr<arrayStruct<SimplePointLight>> lights;
public:
	explicit PointLightsComponent(sptr<arrayStruct<SimplePointLight>> lights);
	void update(TimeStep timeStep) override;
};

class OldMesh3DLayer : public DrawLayer {
	sptr<IndexedMesh3D> mesh;
public:
	OldMesh3DLayer(sptr<ShaderProgram> shader, sptr<IndexedMesh3D> mesh, sptr<Camera> camera, sptr<MaterialModelPhong> materialModel, sptr<arrayStruct<SimplePointLight>> lights);
	void init() override;
};



class Mesh3DLayer : public GenericMeshLayer {
public:
	Mesh3DLayer(sptr<ShaderProgram> shader, sptr<Mesh3D> mesh, sptr<Camera> camera, sptr<MaterialModelPhong> materialModel);
};

class Scene3DLayer : public Layer {
	vector<Mesh3DLayer> meshLayers;
	sptr<PointLightsComponent> pointLightsComponent;
	sptr<Camera> camera;
	vector<sptr<UpdateComponent>> updateComponents;
public:
	Scene3DLayer(sptr<Camera> camera, sptr<arrayStruct<SimplePointLight>> lights);

	void init() override;
	void renderStep() override;
	void updateStep(TimeStep timeStep) override;
	void onEvent(sptr<Event> event, TimeStep timeStep) const override;

	void addMeshLayer(sptr<ShaderProgram> shader, sptr<Mesh3D> mesh, sptr<MaterialModelPhong> materialModel);
	void addComponentToEachMesh(sptr<LayerComponent> comp);
	void addEventListenerToEachMesh(sptr<EventListener> listener);
	void addUpdateComponent(sptr<UpdateComponent> comp);

};
