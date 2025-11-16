#include "clock.hpp"

#include "GLFW/glfw3.h"


FPSStatsResults::FPSStatsResults(uint averageFPS, uint worstFPS, float window): averageFPS(averageFPS), worstFPS(worstFPS), window(window) {}

void FPSStatsResults::log() const {
	LOG(format("FPS: Average FPS = {}, Worst FPS = {} (over {:.1f}s)", averageFPS, worstFPS, window));
}

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

void FPSClock::measureFPSStats(uint averageFPS, uint worstFPS, float window) {
	statsHistory.emplace_back(averageFPS, worstFPS, window);
	statsHistory.back().log();
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
			uint averageFPS = frameCount / avgAccumulator;
			worstDelta = std::max(worstDelta, .00001f);
			uint worstFPS = (uint) 1.0f / worstDelta;
			measureFPSStats(averageFPS, worstFPS, avgWindow);
			avgAccumulator = 0;
			frameCount = 0;
			worstDelta = 0;
		}
	}
	return TimeStep{time, dt};
}

WaitTimer::WaitTimer(float waitDuration): lastTime(-waitDuration), waitDuration(waitDuration) {}

bool WaitTimer::available(float currentTime) {
	if (currentTime - lastTime < waitDuration)
		return false;
	lastTime = currentTime;
	return true;
}
