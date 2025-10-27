#pragma once

#include "buffers.hpp"


class MeshData : public VertexData {
public:
	MeshData(std::initializer_list<pair<string, unsigned char>> attributeInfos)
	: VertexData({{"position", 3}, {"normal", 3}, {"uv", 2}, {"color", 4}}) {
		for (const auto& info : attributeInfos)
			addAttribute(info.first, info.second);
	}
};
