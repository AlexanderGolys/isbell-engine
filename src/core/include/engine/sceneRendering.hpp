#pragma once
#include "colors.hpp"
#include "indexedRendering.hpp"
#include "textures.hpp"
#include "uniforms.hpp"
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

struct CameraUniformData : IDataBlock {
	alignas(16) vec3 position;
	mat4 vp;

	CameraUniformData(vec3 position, const mat4& vp);
	raw_data_ptr data() const override;
	byte_size blockSize() const override;
};

class Camera : public IDirtyDataBlock {
	CameraTransform transform;
	CameraSettings settings;
	CameraUniformData uniformData;

	void setUniformData();
	mat4 vp() const;
	mat4 v() const;
	mat4 p() const;
public:
	Camera(const CameraSettings& settings, vec3 position, vec3 lookAtPos, vec3 upVector = vec3(0, 0, 1));

	vec3 get_position() const;
	vec3 get_lookAtPos() const;
	vec3 get_upVector() const;
	void set_position(const vec3& position);
	void set_lookAtPos(const vec3& lookAtPos);
	void set_upVector(const vec3& upVector);
	void set_transform(const CameraTransform& transform);
	const CameraTransform& get_transform() const;
	raw_data_ptr data() const override;
	byte_size blockSize() const override;
};


class MaterialModelPhong {
	DIRTY_FLAG;
	PROPERTY_DIRTY(vec4, intensities)
	sptr<TextureData> ambientTextureData;
	sptr<TextureData> diffuseTextureData;
	sptr<TextureData> specularTextureData;

public:
	MaterialModelPhong(sptr<TextureData> ambientTextureData, sptr<TextureData> diffuseTextureData, sptr<TextureData> specularTextureData, float ambientIntensity, float diffuseIntensity, float specularIntensity, float shininess);
	MaterialModelPhong(Color ambientColor, Color diffuseColor, Color specularColor, float ambientIntensity, float diffuseIntensity, float specularIntensity, float shininess);

	sptr<TextureData> getAmbientTextureData() const;
	sptr<TextureData> getDiffuseTextureData() const;
	sptr<TextureData> getSpecularTextureData() const;

	raw_data_ptr data() const { return &intensities[0];}
	byte_size byteSize() const { return sizeof(vec4); }
};

class MaterialModelPhongComponent : public LayerComponent {
	Texture2D ambientTexture;
	Texture2D diffuseTexture;
	Texture2D specularTexture;

public:
	explicit MaterialModelPhongComponent(sptr<MaterialModelPhong> materialModel);
	void init() override;
	void update(TimeStep timeStep) override {}
	void setDuringRender() override;
};

struct SimplePointLight : IDataBlock {
	alignas(16) vec3 position;
	vec4 color;
	alignas(16) vec3 intensities;

	SimplePointLight(vec3 position, Color color, vec3 intensities);
	SimplePointLight(const SimplePointLight& other) = delete;
	SimplePointLight& operator=(const SimplePointLight& other) = delete;

	raw_data_ptr data() const;
	byte_size byteSize() const;
};


class OldMesh3DLayer : public DrawLayer {
	sptr<IndexedMesh3D> mesh;
public:
	OldMesh3DLayer(sptr<ShaderProgram> shader, sptr<IndexedMesh3D> mesh, sptr<Camera> camera, sptr<MaterialModelPhong> materialModel, sptr<arrayStruct<SimplePointLight>> lights);
	void init() override;
};



class Mesh3DLayer : public GenericMeshLayer {
public:
	Mesh3DLayer(sptr<ShaderProgram> shader, sptr<Mesh3D> mesh, sptr<MaterialModelPhong> materialModel);
};


class Scene3DLayer : public Layer {
	vector<Mesh3DLayer> meshLayers;
	vector<sptr<UpdateComponent>> globalComponents;
	vector<sptr<LayerComponent>> localComponents;
public:
	Scene3DLayer(sptr<Camera> camera, sptr<arrayStruct<SimplePointLight>> lights);

	void init() override;
	void renderStep() override;
	void updateStep(TimeStep timeStep) override;
	void onEvent(sptr<Event> event, TimeStep timeStep) const override;
	void finalize() override;


	void addMeshLayer(sptr<ShaderProgram> shader, sptr<Mesh3D> mesh, sptr<MaterialModelPhong> materialModel);
	void addLocalComponent(sptr<LayerComponent> comp);
	void addEventListenerToEachMesh(sptr<EventListener> listener);
	void addGlobalComponent(sptr<UpdateComponent> comp);
};
