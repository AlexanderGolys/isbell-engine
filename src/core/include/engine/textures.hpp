#pragma once
#include "colors.hpp"
#include "renderLayers.hpp"
#include "GL/glew.h"

class TextureData {
public:
	virtual ~TextureData() = default;
	virtual array_len width() const = 0;
	virtual array_len height() const = 0;
	virtual GLenum format() const = 0;
	virtual GLenum targetFormat() const = 0;
	virtual raw_data_ptr<uchar> dataPtr() const = 0;
};

class ConstColorTextureData : public TextureData {
	array<uchar, 4> data;
public:
	explicit ConstColorTextureData(const Color& color);
	array_len width() const override;
	array_len height() const override;
	GLenum format() const override;
	GLenum targetFormat() const override;
	raw_data_ptr<uchar> dataPtr() const override;
};

class Texture2D : public LayerComponent {
	sptr<TextureData> imageData;
	gl_id id = 0;
	uint slot;
	string samplerName;
public:
	Texture2D(sptr<TextureData> imageData, uint slot, string_cr samplerName);
	void init() override;
	void update(TimeStep timeStep) override {}
	void setDuringRender() override;

	static sptr<Texture2D> constColorTexture(const Color& color, uint slot, string_cr samplerName);
};