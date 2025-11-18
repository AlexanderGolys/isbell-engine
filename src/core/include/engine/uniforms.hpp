#pragma once
#include "controllers.hpp"
#include "renderLayers.hpp"
#include "shaderTypes.hpp"





template<uniform_struct_dirty S>
class UniformBlock : public LayerComponent {
	UniformBufferObject uniformBuffer;
	sptr<S> value;
	uptr<Controller<S>> controller = nullptr;

public:
	UniformBlock(sptr<S> value, uint bindingPoint, const string& name);
	void initPerShader(gl_id currentProgram) final;
	void setDuringRender() final;
	void update(TimeStep t) final;
	void finalize() final;

	template<controller_function<S> Hom>
	void addController(const Hom& setterFunction);
};

template<vec_float_type VecType>
class UniformVecInnerBlock : public UniformBlock<primitiveVecStruct<VecType>> {
public:
	UniformVecInnerBlock(VecType value, uint bindingPoint, const string& name);
	UniformVecInnerBlock(VecType value, uint bindingPoint, const string& name, const HOM(TimeStep, VecType)& setterFunction);

	void addController(const HOM(TimeStep, VecType)& setterFunction);
};

template<uniform_struct S>
using UniformArrayBlock = UniformBlock<arrayStruct<S>>;

template<vec_float_type VecType>
using UniformVecBlock = UniformBlock<primitiveVecStruct<VecType>>;

using UniformFloatInnerBlock = UniformVecInnerBlock<float>;

// ------------------- Implementation ------------------ //


template <uniform_struct_dirty S>
UniformBlock<S>::UniformBlock(sptr<S> value, uint bindingPoint, const string& name)
: uniformBuffer(bindingPoint, value->byteSize(), name),value(value) {}

template <uniform_struct_dirty S>
void UniformBlock<S>::initPerShader(gl_id currentProgram) {
	uniformBuffer.init(currentProgram);
	uniformBuffer.load(value->data());
	if (controller)
		controller->init();
}

template <uniform_struct_dirty S>
void UniformBlock<S>::setDuringRender() {
	if (not value->isDirty())
		return;
	uniformBuffer.update(value->data());
}

template <uniform_struct_dirty S>
void UniformBlock<S>::update(TimeStep t) {
	if (controller)
		controller->update(t);
}

template <uniform_struct_dirty S>
void UniformBlock<S>::finalize() {
	if (controller)
		controller->finalize();
}

template <uniform_struct_dirty S>
template <controller_function<S> Hom>
void UniformBlock<S>::addController(const Hom& setterFunction) {
	if (controller)
		THROW(ValueError, "Controller already exists for this UniformBlock");
	controller = make_unique<Controller<S>>(value, setterFunction);
}

template <vec_float_type VecType>
UniformVecInnerBlock<VecType>::UniformVecInnerBlock(VecType value, uint bindingPoint, const string& name)
: UniformBlock<primitiveVecStruct<VecType>>(make_shared<primitiveVecStruct<VecType>>(value), bindingPoint, name) {}

template <vec_float_type VecType>
UniformVecInnerBlock<VecType>::UniformVecInnerBlock(VecType value, uint bindingPoint, const string& name, const HOM(TimeStep, VecType)& setterFunction)
: UniformVecInnerBlock(value, bindingPoint, name) {
	addController(setterFunction);
}

template <vec_float_type VecType>
void UniformVecInnerBlock<VecType>::addController(const HOM(TimeStep, VecType)& setterFunction) {
	BIHOM(sptr<primitiveVecStruct<VecType>>&, TimeStep, void) controlFunction =
		[setterFunction](sptr<primitiveVecStruct<VecType>>& obj, TimeStep t) {
			obj->set_value(setterFunction(t));
	};
	UniformBlock<primitiveVecStruct<VecType>>::addController(controlFunction);
}
