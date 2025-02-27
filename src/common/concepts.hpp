#pragma once
#include "macros.hpp"


template<typename A>
constexpr A one = A(1);

template<typename A>
constexpr A zero = A(0);

template<typename T>
    concept AbelianSemigroup = requires(T a, T b) {
		{ a + b } -> std::convertible_to<T>;
	};

template<typename T>
    concept Semigroup = requires(T a, T b) {
		{ a*b } -> std::convertible_to<T>;
	};

template<typename T>
    concept TotalOrder = requires(T a, T b) {
		{ a < b } -> std::convertible_to<bool>;
		{ a == b } -> std::convertible_to<bool>;
	};

template<typename T>
    concept TotallyOrderedAbelianSemigroup = AbelianSemigroup<T> && TotalOrder<T>;

template<typename T>
    concept AbelianMonoid = AbelianSemigroup<T> && requires {
		{ zero<T> } -> std::convertible_to<T>; };

template<typename G>
    concept Monoid = Semigroup<G> && requires {
		{ one<G> } -> std::convertible_to<G>; };

template<typename G>
    concept GroupConcept = Monoid<G> && (requires(G g) {
		{ ~g } -> std::convertible_to<G>; } || requires(G g) {
		{ inverse(g) } -> std::convertible_to<G>; } || requires(G g) {
		{ one<G>/g } -> std::convertible_to<G>; });

template<typename G>
    concept AbelianGroupConcept = AbelianMonoid<G> && requires(G g) {
		{ -g } -> std::convertible_to<G>; };

template<typename T> concept Rng = AbelianGroupConcept<T> && Semigroup<T>;


template<typename T> concept RingConcept = Rng<T> && Monoid<T>;


template<typename T> concept DivisionRing = RingConcept<T> && GroupConcept<T>;

template<typename A, typename R> concept ModuleConcept = Rng<R> && AbelianGroupConcept<A> && requires(A a, R r) {
	{ a*r } -> std::convertible_to<A>;
};

template<typename A, typename R> concept Algebra = ModuleConcept<A, R> && Rng<A>;


template<typename V, typename K> concept VectorSpaceConcept = ModuleConcept<V, K> && DivisionRing<K>;

template<typename T> concept Normed = requires(T a, float c) {
	{norm(a)} -> std::convertible_to<float>;
	{a/c} -> std::convertible_to<T>;
};

template<typename V> concept EuclideanSpaceConcept = VectorSpaceConcept<V, float> && requires(V a, V b) {
	{ dot(a, b) } -> std::convertible_to<float>;
};
