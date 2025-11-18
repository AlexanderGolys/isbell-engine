#pragma once
#include "buffers.hpp"
#include "clock.hpp"
#include "event.hpp"
#include "listeners.hpp"
#include "logging.hpp"
#include "shaders.hpp"




class LayerComponent {
public:
	LayerComponent(const LayerComponent&) = delete;
	virtual ~LayerComponent();

	virtual void init() = 0;
	virtual void initPerShader(gl_id currentProgram) = 0;
	virtual void update(TimeStep timeStep) = 0;
	virtual void setDuringRender() = 0;
	virtual void finalize() = 0;
};

template<typename LayerCompType>
concept layer_component = derived_from<LayerCompType, LayerComponent>;

class Layer {
public:
	Layer(const Layer&) = delete;
	Layer& operator=(const Layer&) = delete;
	virtual ~Layer();

	virtual void init() = 0;
	virtual void renderStep() = 0;
	virtual void updateStep(TimeStep timeStep) = 0;
	virtual void finalize() = 0;

	virtual void addComponent(const sptr<LayerComponent>& comp) = 0;

	template<layer_component LayerCompType, typename... Args>
	void emplaceComponent(Args&&... args);
};

template <layer_component LayerCompType, typename ... Args>
void Layer::emplaceComponent(Args&&... args) {
	addComponent(make_shared<LayerCompType>(std::forward<Args>(args)...));
}


class DrawLayer : public Layer {
	sptr<ShaderProgram> shader;
	vector<sptr<LayerComponent>> components;

protected:
	VertexArray vao;
	gl_id getProgramID() const;

public:
	explicit DrawLayer(sptr<ShaderProgram> shader);

	void init() override;
	virtual void customRenderStep() {}
	void updateStep(TimeStep timeStep) override;
	void finalize() override;
	void addComponent(const sptr<LayerComponent>& comp) override;
	void renderStep() final;
};

class GenericMeshLayer : public DrawLayer {
	sptr<GeometricData> mesh;
public:
	GenericMeshLayer(sptr<ShaderProgram> shader, sptr<GeometricData> mesh);
	void init() override;
	void customRenderStep() override;
};
