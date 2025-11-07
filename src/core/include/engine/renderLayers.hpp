#pragma once
#include "buffers.hpp"
#include "logging.hpp"
#include "shaders.hpp"
#include "uniforms.hpp"


class Layer {
public:
	virtual ~Layer() = default;
	virtual void init() = 0;
	virtual void renderStep() = 0;
	virtual void updateStep(float t, float dt) = 0;
};

class LayerComponent {
public:
	virtual ~LayerComponent() = default;
	virtual void init() = 0;
	virtual void update(float t, float dt) = 0;
	virtual void setDuringRender() const = 0;
};


class IndexedDrawLayer : public Layer {
	sptr<ShaderProgram> shader;
	sptr<VertexArray> vao;
	vector<sptr<LayerComponent>> components;

public:
	IndexedDrawLayer(sptr<ShaderProgram> shader, sptr<VertexArray> vao, const vector<sptr<LayerComponent>>& components = {});
	void addComponent(sptr<LayerComponent> comp);

	void init() override;
	void renderStep() override;
	void updateStep(float t, float dt) override;
};

class MeshLayer : public IndexedDrawLayer {
	sptr<IndexedMesh> mesh;
	sptr<MaterialPhong> material;
};