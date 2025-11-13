#pragma once
#include "buffers.hpp"
#include "event.hpp"
#include "logging.hpp"
#include "shaders.hpp"


class Layer {
public:
	virtual ~Layer() = default;
	virtual void init() = 0;
	virtual void renderStep() = 0;
	virtual void updateStep(float t, float dt) = 0;
	virtual void onEvent(const Event& event, float t, float delta) const {}
};

class LayerComponent {
public:
	virtual ~LayerComponent() = default;
	virtual void init() = 0;
	virtual void update(float t, float dt) = 0;
	virtual void setDuringRender() const = 0;
};

class EventListener {
public:
	virtual ~EventListener() = default;
	virtual void init() {}
	virtual void onEvent(const Event& event, float t, float dt) = 0;
	virtual bool listensToEventType(EventType type) const = 0;
};

class DrawLayer : public Layer {
	CONST_PROPERTY(sptr<VertexArray>, vao);
	sptr<ShaderProgram> shader;
	vector<sptr<LayerComponent>> components;
	vector<sptr<EventListener>> eventListeners;

public:
	DrawLayer(sptr<ShaderProgram> shader, sptr<VertexArray> vao, const vector<sptr<LayerComponent>>& components = {});
	void addComponent(sptr<LayerComponent> comp);
	void addEventListener(sptr<EventListener> listener);

	void init() override;
	void renderStep() final;
	virtual void customRenderStep() {}
	void updateStep(float t, float dt) override;
	void onEvent(const Event& event, float t, float dt) const override;
};

