#pragma once

#include "spdlog/spdlog.h"
#include "spdlog/sinks/stdout_color_sinks.h"

#include <memory>



class Logger {
	static std::shared_ptr<spdlog::logger> engine_logger;
	static std::shared_ptr<spdlog::logger> external_logger;
	static std::shared_ptr<spdlog::logger> pure_logger;


public:
	static void init();
	static std::shared_ptr<spdlog::logger>& getEngineLogger();
	static std::shared_ptr<spdlog::logger>& getExternalLogger();
	static std::shared_ptr<spdlog::logger>& getPureLogger();
};

#define LOG(...) ::Logger::getEngineLogger()->info(__VA_ARGS__)
#define LOG_ERROR(...) Logger::getEngineLogger()->error(__VA_ARGS__)
#define LOG_WARN(...) Logger::getEngineLogger()->warn(__VA_ARGS__)

#define APP_LOG(...) Logger::getExternalLogger()->info(__VA_ARGS__)
#define APP_LOG_ERROR(...) Logger::getExternalLogger()->error(__VA_ARGS__)
#define APP_LOG_WARN(...) Logger::getExternalLogger()->warn(__VA_ARGS__)

#define PURE_LOG(...) Logger::getPureLogger()->info(__VA_ARGS__)
#define PURE_LOG_ERROR(...) Logger::getPureLogger()->error(__VA_ARGS__)
#define PURE_LOG_WARN(...) Logger::getPureLogger()->warn(__VA_ARGS__)
