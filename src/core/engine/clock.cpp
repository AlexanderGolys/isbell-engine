#include "clock.hpp"

#include "GLFW/glfw3.h"


FPSStatsResults::FPSStatsResults(float averageFPS, float worstFPS, float window): averageFPS(averageFPS), worstFPS(worstFPS), window(window) {}

float FPSClock::getTime() const {
	return glfwGetTime() - timeZero;
}

FPSClock::FPSClock(float avgWindowSeconds) {
	if (avgWindowSeconds > 0) {
		measureFPS = true;
		avgWindow = avgWindowSeconds;
	} else {
		measureFPS = false;
		avgWindow = 0;
	}
	timeZero = glfwGetTime();
	time = 0;
}

void FPSClock::reset() {
	timeZero = glfwGetTime();
	time = 0;
}

void FPSClock::measureFPSStats(const FPSStatsResults& result) {
	LOG(0, "FPS Stats (last " + to_string(avgWindow) + "s): Average FPS: " + to_string(result.averageFPS) + ", Worst FPS: " + to_string(result.worstFPS));
	statsHistory.push_back(result);
}

TimeStep FPSClock::tick() {
	float dt = getTime() - time;
	time += dt;
	if (measureFPS) {
		avgAccumulator += dt;
		frameCount++;
		if (dt > worstDelta)
			worstDelta = dt;
		if (avgAccumulator >= avgWindow) {
			float averageFPS = frameCount / avgAccumulator;
			float worstFPS = 1000.0f / worstDelta;
			measureFPSStats(FPSStatsResults(averageFPS, worstFPS, avgWindow));
			avgAccumulator = 0;
			frameCount = 0;
			worstDelta = 0;
		}
	}
	return TimeStep{time, dt};
}