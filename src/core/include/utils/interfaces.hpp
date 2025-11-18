#pragma once
#include "macros.hpp"

class DirtyFlag {
	bool dirty = true;
public:
	void markDirty() { dirty = true; }
	void markClean() { dirty = false; }
	bool isDirty() const { return dirty; }
};

class IDataBlock {
public:
	virtual ~IDataBlock() = default;
	virtual raw_data_ptr data() const = 0;
	virtual byte_size blockSize() const = 0;
};

class IDirtyDataBlock : public IDataBlock, public DirtyFlag {};