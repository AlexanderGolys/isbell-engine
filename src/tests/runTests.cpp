#include "hyperbolicTests.hpp"
#include "discreteFuncTests.hpp"
#include "formattersTests.hpp"
#include "experimentalTests.hpp"
#include "filesystemTests.hpp"
#include "quatGLSLModuleTests.hpp"

#include "../utils/logging.hpp"

int main()
  {
	logging::Logger::init();
	UnitTestResult total_result;

	runTest("Hyperbolic Tests", hyperbolicTests__all, total_result);
	runTest("Discrete Function Tests", discreteFuncTests__all, total_result);
	runTest("Formatters Tests", formattersTests__all, total_result);
	runTest("Experimental Tests", experimentalTests__all, total_result);
	runTest("Filesystem Tests", filesystemTests__all, total_result);
	runTest("Quaternion GLSL Module Tests", quatGLSLModuleTests__all, total_result);
	LOG_PURE("--------------------------------");
	printTestResult("All Tests", total_result);
  }
