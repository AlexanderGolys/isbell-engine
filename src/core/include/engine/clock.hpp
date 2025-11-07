#pragma once
#include "exceptions.hpp"
#include "logging.hpp"


struct FPSStatsResults {
	float averageFPS;
	float worstFPS;
	float window;

	FPSStatsResults(float averageFPS, float worstFPS, float window);
};

struct TimeStep {
	float time;
	float deltaTime;
};




class FPSClock {
	float time = 0;
	float timeZero = 0;
	float avgWindow;
	bool measureFPS;
	float avgAccumulator = 0;
	int frameCount = 0;
	float worstDelta = 0;
	vector<FPSStatsResults> statsHistory = {};
public:
	float getTime() const;
	explicit FPSClock(float avgWindowSeconds = -1.f);
	void reset();
	void measureFPSStats(const FPSStatsResults &result);
	TimeStep tick();
};
