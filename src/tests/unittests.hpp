#pragma once
#include <cassert>

template <typename... Args>
void assertEqual_(const Args&... args) {
	assert(args == ...);
};


template <typename T, typename... Args>
void assertNearlyEqual_(const T& A0, const Args&... args) {
	for (const auto& A : {args...})
		assert(nearlyEqual(A0, A));
};

#define assertEqual(...) assertEqual_(__VA_ARGS__)
#define assertNearlyEqual(...) assertNearlyEqual_(__VA_ARGS__)
