#pragma once
#include <format>
#include "func.hpp"
// ReSharper disable CppParameterMayBeConstPtrOrRef

using std::formatter, std::format, std::to_string;

namespace std
{
	template<>
	struct formatter<Complex>
	{
		constexpr format_parse_context::const_iterator parse(format_parse_context &ctx) {
			return ctx.begin();
		}

		template<typename FormatContext>
		std::format_context::iterator format(const Complex& c, FormatContext& ctx) const{
			if (nearlyEqual(c.imag(), 0.f))
				return format_to(ctx.out(), "{}", c.real());
			if (nearlyEqual(c.real(), 0.f))
				return format_to(ctx.out(), "{}i", c.imag());
			if (c.imag() < 0.f)
				return format_to(ctx.out(), "{} - {}i", c.real(), -c.imag());
			return format_to(ctx.out(), "{} + {}i", c.real(), c.imag());
		}
	};

	template<>
	struct formatter<vec3>
	{
		constexpr auto parse(format_parse_context& ctx) {
			return ctx.begin();
		}

		template<typename FormatContext>
		std::format_context::iterator format(const vec3& v, FormatContext& ctx) const{
			return format_to(ctx.out(), "({}, {}, {})",
				static_cast<double>(v.x), static_cast<double>(v.y), static_cast<double>(v.z));
		}
	};

	template<>
	struct formatter<vec2>
	{
		constexpr auto parse(format_parse_context& ctx) { return ctx.begin(); }

		template<typename FormatContext>
		std::format_context::iterator format(const vec2& v, FormatContext& ctx) const{
			return format_to(ctx.out(), "({}, {})",
				static_cast<double>(v.x), static_cast<double>(v.y));
		}
	};

	template<>
	struct formatter<vec4>
	{
		constexpr auto parse(format_parse_context& ctx) { return ctx.begin(); }

		template<typename FormatContext>
		std::format_context::iterator format(const vec4& v, FormatContext& ctx) const{
			return format_to(ctx.out(), "({}, {}, {}, {})",
				static_cast<double>(v.x), static_cast<double>(v.y), static_cast<double>(v.z), static_cast<double>(v.w));
		}
	};

	template<>
	struct formatter<vec5>
	{
		constexpr auto parse(format_parse_context& ctx) { return ctx.begin(); }

		template<typename FormatContext>
		std::format_context::iterator format(const vec5& v, FormatContext& ctx) const{
			return format_to(ctx.out(), "({}, {}, {}, {}, {})", v.x, v.y, v.z, v.w, v.v);
		}
	};

	template<>
	struct formatter<mat2>
	{
		constexpr auto parse(format_parse_context& ctx) { return ctx.begin(); }

		template<typename FormatContext>
		std::format_context::iterator format(const mat2& m, FormatContext& ctx) const{
			return format_to(ctx.out(), "|{} {}|\n|{} {}|",
				static_cast<double>(m[0][0]), static_cast<double>(m[0][1]),
				static_cast<double>(m[1][0]), static_cast<double>(m[1][1]));
		}
	};

	template<>
	struct formatter<mat3>
	{
		constexpr auto parse(format_parse_context& ctx) { return ctx.begin(); }

		template<typename FormatContext>
		std::format_context::iterator format(const mat3& m, FormatContext& ctx) const {
			return format_to(
				ctx.out(),
				"|{} {} {}|\n|{} {} {}|\n|{} {} {}|)",
				static_cast<double>(m[0][0]), static_cast<double>(m[0][1]), static_cast<double>(m[0][2]),
				static_cast<double>(m[1][0]), static_cast<double>(m[1][1]), static_cast<double>(m[1][2]),
				static_cast<double>(m[2][0]), static_cast<double>(m[2][1]), static_cast<double>(m[2][2])
			);
		}
	};

	template<>
	struct formatter<mat4>
	{
		constexpr auto parse(format_parse_context& ctx) { return ctx.begin(); }

		template<typename FormatContext>
		std::format_context::iterator format(const mat4& m, FormatContext& ctx) const
		{
			return format_to(
				ctx.out(),
				"|{} {} {} {}|\n|{} {} {} {}|\n|{} {} {} {}|\n|{} {} {} {}|",
				static_cast<double>(m[0][0]), static_cast<double>(m[0][1]), static_cast<double>(m[0][2]), static_cast<double>(m[0][3]),
				static_cast<double>(m[1][0]), static_cast<double>(m[1][1]), static_cast<double>(m[1][2]), static_cast<double>(m[1][3]),
				static_cast<double>(m[2][0]), static_cast<double>(m[2][1]), static_cast<double>(m[2][2]), static_cast<double>(m[2][3]),
				static_cast<double>(m[3][0]), static_cast<double>(m[3][1]), static_cast<double>(m[3][2]), static_cast<double>(m[3][3])
			);
		}
	};
}

inline string formatByteSize(byte_size s) {
	if (s < 1024)
		return to_string(s) + " B";
	if (s < 1024 * 1024)
		return format("{:.1f} KB", static_cast<double>(s) / 1024.0);
	if (s < 1024 * 1024 * 1024)
		return format("{:.1f} MB", static_cast<double>(s) / (1024.0 * 1024.0));
	return format("{:.1f} GB", static_cast<double>(s) / (1024.0 * 1024.0 * 1024.0));
}
