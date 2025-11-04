#pragma once
#include "dynamicMesh.hpp"
#include "exceptions.hpp"
#include "layers.hpp"
#include "mat.hpp"
#include "renderLayer.hpp"

struct Camera {
	vec3 position;
	vec3 lookAt;
	vec3 upVector;
	float fov;
	float aspectRatio;
	float nearClip;
	float farClip;

	Camera(vec3 position, vec3 lookAt, vec3 upVector = vec3(0, 0, 1), float fov = 45.0f, float aspectRatio = 16.0f / 9.0f, float nearClip = 0.1f, float farClip = 100.0f);

	mat4 projectionMatrix() const;
	mat4 viewMatrix() const;
	mat4 vp() const;
	mat4 mvp(const mat4& modelTransform) const;
};

class Light {
public:
	virtual ~Light() = default;
	virtual mat4 computeMatrix() const = 0;
};

class PointLight : public Light {
	vec3 position;
	vec3 color;
	float intensityConst, intensityLinear, intensityQuad;

public:
	PointLight(vec3 position, vec3 color, float intensityConst, float intensityLinear, float intensityQuad);
	PointLight(vec3 position, vec3 color, vec3 attenuation);
};

class Mesh3DObject {
	SE3 transform;
	shared_ptr<DynamicMesh3D> mesh;
	bool dirty=false;

public:
	explicit Mesh3DObject(const shared_ptr<DynamicMesh3D>& mesh, const SE3& transform = SE3::identity()) : mesh(mesh), transform(transform) {}
	void setTransform(const SE3& newTransform) { transform = newTransform; }
	const SE3& getTransform() const { return transform; }
	const shared_ptr<DynamicMesh3D>& getMesh() const { return mesh; }

	void markDirty() { dirty = true; }
	void markClean() { dirty = false; }
	bool isDirty() const { return dirty; }

};

class Mesh3DLayer : public RenderLayer {
	shared_ptr<Mesh3DObject> meshObject;
	shared_ptr<VertexBuffer> vertexBuffer;
	shared_ptr<IndexBuffer> indexBuffer;
public:
	Mesh3DLayer(const shared_ptr<ShaderProgram>& shaderProgram, const shared_ptr<Mesh3DObject>& meshObject)
	: RenderLayer(shaderProgram, nullptr, {}), meshObject(meshObject), vertexBuffer(nullptr), indexBuffer(nullptr)
	{
		indexBuffer = make_shared<IndexBuffer>(meshObject->getMesh()->getIndexArray());
		vertexBuffer = make_shared<VertexBuffer>(VertexBufferLayout{
		  {"position", VEC3},
		  {"normal", VEC3},
		  {"uv", VEC2},
		  {"color", VEC4}
	  });
		vertexBuffer->uploadData(meshObject->getMesh()->getVertexData(), meshObject->getMesh()->size());
		vertexArray = make_shared<VertexArray>(indexBuffer, vector{vertexBuffer});
	};
};

class Scene3D : public Layer {
	vector<shared_ptr<Mesh3DObject>> meshObjects;
	vector<Light> lights;
	Camera camera;

	int activeCameraIndex = 0;
};
