#include "textures.hpp"

#include "glCommand.hpp"

ConstColorTextureData::ConstColorTextureData(const Color& color) {
	data[0] = static_cast<uchar>(round(color.r() * 255));
	data[1] = static_cast<uchar>(round(color.g() * 255));
	data[2] = static_cast<uchar>(round(color.b() * 255));
	data[3] = static_cast<uchar>(round(color.a() * 255));
}

array_len ConstColorTextureData::width() const {
	return 1;
}

array_len ConstColorTextureData::height() const {
	return 1;
}

GLenum ConstColorTextureData::format() const {
	return GL_RGBA;
}

GLenum ConstColorTextureData::targetFormat() const {
	return GL_RGBA8;
}

raw_data_ptr<uchar> ConstColorTextureData::dataPtr() const {
	return data.data();
}

Texture2D::Texture2D(sptr<TextureData> imageData, uint slot, string_cr samplerName)
: imageData(imageData), slot(slot), samplerName(samplerName) {}

void Texture2D::init() {
	GLCommand::createTexture(id);
	GLCommand::bindTexture2D(id, slot);
	GLCommand::setTexture2DFilters(GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR, GL_REPEAT, GL_REPEAT);
	GLCommand::loadTexture2(imageData->targetFormat(), imageData->width(), imageData->height(), imageData->format(), imageData->dataPtr());
	GLCommand::calculateMipmapsTexture2D();
}

void Texture2D::setDuringRender() {
	GLCommand::bindTexture2D(id, slot);
	GLCommand::setSampler2D(samplerName, slot);
}

sptr<Texture2D> Texture2D::constColorTexture(const Color& color, uint slot, string_cr samplerName) {
	return make_shared<Texture2D>(make_shared<ConstColorTextureData>(color), slot, samplerName);
}
