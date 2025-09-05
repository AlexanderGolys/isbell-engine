#include "macros.hpp"
#include "logging.hpp"
#include "spdlog/spdlog.h"
#include "spdlog/sinks/stdout_color_sinks.h"
#include <cstdarg>

namespace logging {

	// Opaque logger handles
	void* Logger::engine_logger = nullptr;
	void* Logger::external_logger = nullptr;
	void* Logger::pure_logger = nullptr;

	void Logger::init() {
		Logger::engine_logger = (void*)spdlog::stdout_color_mt("Engine").get();
		Logger::external_logger = (void*)spdlog::stdout_color_mt("External").get();
		Logger::pure_logger = (void*)spdlog::stdout_color_mt("Pure").get();
	}

	void* Logger::getEngineLogger() { return engine_logger; }
	void* Logger::getExternalLogger() { return external_logger; }
	void* Logger::getPureLogger() { return pure_logger; }

	static spdlog::logger* get_logger(void* handle) {
		return reinterpret_cast<spdlog::logger*>(handle);
	}

	void log_info(const char* fmt, ...) {
		char buf[1024];
		va_list args;
		va_start(args, fmt);
		vsnprintf(buf, sizeof(buf), fmt, args);
		va_end(args);
		get_logger(Logger::engine_logger)->info(buf);
	}
	void log_error(const char* fmt, ...) {
		char buf[1024];
		va_list args;
		va_start(args, fmt);
		vsnprintf(buf, sizeof(buf), fmt, args);
		va_end(args);
		get_logger(Logger::engine_logger)->error(buf);
	}
	void log_warn(const char* fmt, ...) {
		char buf[1024];
		va_list args;
		va_start(args, fmt);
		vsnprintf(buf, sizeof(buf), fmt, args);
		va_end(args);
		get_logger(Logger::engine_logger)->warn(buf);
	}

	void log_info(const string& msg) {
		get_logger(Logger::engine_logger)->info(msg);
	}
	void log_error(const string& msg) {
		get_logger(Logger::engine_logger)->error(msg);
	}
	void log_warn(const string& msg) {
		get_logger(Logger::engine_logger)->warn(msg);
	}

} // namespace logging
