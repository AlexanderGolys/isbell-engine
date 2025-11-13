#pragma once
#include "exceptions.hpp"
#include "logging.hpp"


struct FPSStatsResults {
	uint averageFPS;
	uint worstFPS;
	float window;

	FPSStatsResults(uint averageFPS, uint worstFPS, float window);
	void log() const;
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
	explicit FPSClock(float avgWindowSeconds = -1.f);

	float getTime() const;
	void reset();
	void measureFPSStats(uint averageFPS, uint worstFPS, float window);
	TimeStep tick();
};
