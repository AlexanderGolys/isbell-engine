#pragma once
#include "buffers.hpp"
#include "clock.hpp"
#include "event.hpp"
#include "listeners.hpp"
#include "logging.hpp"
#include "shaders.hpp"


class Layer {
	vector<sptr<EventListener>> eventListeners;

public:
	virtual ~Layer() = default;
	virtual void init() = 0;
	virtual void renderStep() = 0;
	virtual void updateStep(TimeStep timeStep) {}
	virtual void onEvent(sptr<Event> event, TimeStep timeStep) const;

	void addEventListener(sptr<EventListener> listener);
};

class LayerComponent {
public:
	virtual ~LayerComponent() = default;
	virtual void init() = 0;
	virtual void update(TimeStep timeStep) = 0;
	virtual void setDuringRender() = 0;
};

class CombinedLayerComponent : public LayerComponent {
	vector<sptr<LayerComponent>> components;
public:
	explicit CombinedLayerComponent(const vector<sptr<LayerComponent>>& components={});

	void addComponent(sptr<LayerComponent> comp);

	void init() override;
	void update(TimeStep timeStep) override;
	void setDuringRender() override;
};

class UpdateComponent : public LayerComponent {
public:
	void setDuringRender() final {}
};


class DrawLayer : public Layer {
	CONST_PROPERTY(sptr<VertexArray>, vao);
	sptr<ShaderProgram> shader;
	vector<sptr<LayerComponent>> components;

public:
	DrawLayer(sptr<ShaderProgram> shader, sptr<VertexArray> vao, const vector<sptr<LayerComponent>>& components = {});
	void addComponent(sptr<LayerComponent> comp);

	void init() override;
	void renderStep() final;
	virtual void customRenderStep() {}
	virtual bool prerenderStep() { return true; }
	void updateStep(TimeStep timeStep) override;
};

class GenericMeshLayer : public DrawLayer {
	sptr<GeometricData> mesh;
public:
	GenericMeshLayer(sptr<ShaderProgram> shader, sptr<GeometricData> mesh);
	void init() override;
	void customRenderStep() override;
};
