#pragma once
#include "renderLayers.hpp"
#include "sceneRendering.hpp"
#include "GL/glew.h"

class UniformComponent : public LayerComponent {
	string name;
public:
	explicit UniformComponent(const string& name) : name(name) {}
	string getName() const;

	virtual void set(GLint location) const = 0;
	void init() override {}
	void setDuringRender() const override;
};

template<typename T, GLSLType type>
class ConstantUniform : public UniformComponent {
	T value;
};


class CameraUniform : public LayerComponent {
	sptr<Camera> camera;
public:
	explicit CameraUniform(sptr<Camera> camera) : camera(camera) {}
	void setDuringRender() const override {

	}

};