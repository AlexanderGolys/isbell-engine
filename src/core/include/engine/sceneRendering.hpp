#pragma once
#include "colors.hpp"
#include "indexedRendering.hpp"
#include "textures.hpp"
#include "uniforms.hpp"
#include "indexedMesh.hpp"

struct CameraSettings {
	float fov_x;
	float aspectRatio;
	float clippingRangeMin;
	float clippingRangeMax;
	mat4 projMatrix;

	explicit CameraSettings(float fov_x = PI/4, float aspectRatio = 16/9.f, float clippingRangeMin = .01f, float clippingRangeMax = 100.f);
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
};


class CameraUniform : public LayerComponent {
	sptr<Camera> camera;
	ConstantUniform<mat4> vp_uni;
	ConstantUniform<vec3> campos_uni;
public:
	explicit CameraUniform(sptr<Camera> camera);
	void setDuringRender() const override;
	void update(float t, float dt) override;
	void init() override {}
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
	ConstantUniform<vec4> intensitiesUni;

public:
	explicit MaterialModelPhongComponent(sptr<MaterialModelPhong> materialModel);
	void init() override;
	void update(float t, float dt) override;
	void setDuringRender() const override;
};

class SimplePointLight {
	DIRTY_FLAG;
	PROPERTY_DIRTY(vec3, position);
	PROPERTY_DIRTY(Color, color);
	PROPERTY_DIRTY(vec3, intensities);
public:
	SimplePointLight(vec3 position, Color color, vec3 intensities) : position(position), color(color), intensities(intensities) {}
};

class PointLightsComponent : public LayerComponent {
	vector<sptr<SimplePointLight>> lights;
	ConstantUniformArray<vec3> lightsPositionsUni;
	ConstantUniformArray<vec4> lightsColorsUni;
	ConstantUniformArray<vec3> lightsIntensitiesUni;

public:
	explicit PointLightsComponent(const vector<sptr<SimplePointLight>>& lights);
	void init() override {}
	void update(float t, float dt) override;
	void setDuringRender() const override;
};

class Mesh3DLayer : public DrawLayer {
	sptr<IndexedMesh3D> mesh;
public:
	Mesh3DLayer(sptr<ShaderProgram> shader, sptr<IndexedMesh3D> mesh, sptr<Camera> camera, sptr<MaterialModelPhong> materialModel, const vector<sptr<SimplePointLight>>& lights);
	void init() override;
};

class GenericMeshLayer : public DrawLayer {
	sptr<GeometricData> mesh;
public:
	GenericMeshLayer(sptr<ShaderProgram> shader, sptr<GeometricData> mesh);
	void init() override;
	void customRenderStep() override;
};

