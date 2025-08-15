#include "logging.hpp"

std::shared_ptr<spdlog::logger> Logger::engine_logger;
std::shared_ptr<spdlog::logger> Logger::external_logger;
std::shared_ptr<spdlog::logger> Logger::pure_logger;

void Logger::init() {
	engine_logger = std::make_shared<spdlog::logger>("Engine", std::make_shared<spdlog::sinks::stdout_color_sink_mt>());
	spdlog::register_logger(engine_logger);
	engine_logger->set_pattern("%^[%T:%e] %v");
	engine_logger->set_level(spdlog::level::trace);
	engine_logger->flush_on(spdlog::level::trace);

	external_logger = std::make_shared<spdlog::logger>("External", std::make_shared<spdlog::sinks::stdout_color_sink_mt>());
	spdlog::register_logger(external_logger);
	external_logger->set_pattern("%^[%T:%e] [SCENE] %v");
	external_logger->set_level(spdlog::level::trace);
	external_logger->flush_on(spdlog::level::trace);

	pure_logger = std::make_shared<spdlog::logger>("Pure", std::make_shared<spdlog::sinks::stdout_color_sink_mt>());
	spdlog::register_logger(pure_logger);
	pure_logger->set_pattern("%v");
	pure_logger->set_level(spdlog::level::trace);
	pure_logger->flush_on(spdlog::level::trace);

	pure_logger->info("------------------------------------");
}

std::shared_ptr<spdlog::logger>& Logger::getEngineLogger() {
	return engine_logger;
}

std::shared_ptr<spdlog::logger>& Logger::getExternalLogger() {
	return external_logger;
}

std::shared_ptr<spdlog::logger> & Logger::getPureLogger() {
	return pure_logger;
}
