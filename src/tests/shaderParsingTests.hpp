#pragma once
#include "unittests.hpp"
#include "shaderGenerator.hpp"


inline UnitTestResult shaderParsingTests__all() {
	UnitTestResult res;
	ShaderFunctionModule mod(Path("../../../../src/core/shaders/glsl-macro-modules/lightTools.glsl"), "lightTools");

	res.runTest([&]() {
		return assertTrue_UT_(mod.containsFunction("parsePointLight"));
	});

	res.runTest([&]() {
		return assertTrue_UT_(mod.containsFunction("rayTraceColor"));
	});

	res.runTest([&]() {
		return assertFalse_UT_(mod.containsStruct("PointLight"));
	});

	res.runTest([&]() {
		return assertMore_UT_(mod.functionNames.size(), 5);
	});

	return res;
}
