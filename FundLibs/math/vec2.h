#pragma once
#include <cstdint>
#include <cmath>

template <typename T, uint32_t vecN>
struct vec {
public:
	T x[vecN];
	inline vec() { for (uint32_t i(vecN); i--;) x[i] = 0; }
	inline vec(T a) { for (uint32_t i(vecN); i--;) x[i] = a; }
	inline vec(T x, T y) { this->x[0] = x; this->x[1] = y; }
	template <typename O> inline vec<T, vecN> operator=(const vec<O, vecN>& o) { for (uint32_t i(vecN); i--;) x[i] = (T)(o.get(i)); return *this; }
	inline T& operator()(const uint32_t& i) { return x[i]; }
	inline T get(const uint32_t& i) const { return x[i]; }

	//c
	inline T abs2() const { T l2(0); for (uint32_t i(vecN); i--;) l2 += x[i] * x[i]; return l2; }
	inline T abs() const { return (T)sqrt(abs2()); }
	inline T abs(T abs2) const { return (T)sqrt(abs2); }
	template <typename O> T dot(vec<O, vecN> o) const { T dp(0); for (uint32_t i(vecN); i--;) dp += x[i] * o.get(i); return dp; }
	//logical
	inline uint32_t operator()() const { bool b(false); for (uint32_t i(vecN); i--;) b = b || x[i]; return b; }
	inline uint32_t operator!() const { return !operator()(); }
	template <typename O> inline uint32_t operator==(const vec<O, vecN>& o) const { bool b(true); for (uint32_t i(vecN); i--;) b = b && (x[i] == o.get(i)); return b; }
	template <typename O> inline uint32_t operator!=(const vec<O, vecN>& o) const { return !(operator==(o)); }
	//v-v
	inline vec<T, vecN> norm(T r) { T n = abs2(); if (n) { n = r / abs(n); operator*=(n); } return *this; }
	inline vec<T, vecN> get_norm(T r) { T n = abs2(); if (n) n = r / abs(n); return (*this) * n; }
	inline vec<T, vecN> operator-() const { vec<T, vecN> v = (*this); for (uint32_t i(vecN); i--;) v(i) = -(v(i)); return v; }
	template <typename O> inline vec<T, vecN> operator+=(const vec<O, vecN>& o) { for (uint32_t i(vecN); i--;) x[i] += o.get(i); return *this; }
	template <typename O> inline vec<T, vecN> operator-=(const vec<O, vecN>& o) { for (uint32_t i(vecN); i--;) x[i] -= o.get(i); return *this; }
	template <typename O> inline vec<T, vecN> operator+ (const vec<O, vecN>& o) const { vec<T, vecN> v = (*this); v += o; return v; }
	template <typename O> inline vec<T, vecN> operator- (const vec<O, vecN>& o) const { return (*this) + (-o); }
	template <typename O> inline vec<T, vecN> operator*=(const vec<O, vecN>& o) { for (uint32_t i(vecN); i--;) x[i] *= o.get(i); return *this; }
	template <typename O> inline vec<T, vecN> operator/=(const vec<O, vecN>& o) { for (uint32_t i(vecN); i--;) x[i] /= o.get(i); return *this; }
	template <typename O> inline vec<T, vecN> operator* (const vec<O, vecN>& o) const { vec<T, vecN> v = (*this); v *= o; return v; }
	template <typename O> inline vec<T, vecN> operator/ (const vec<O, vecN>& o) const { vec<T, vecN> v = (*this); v /= o; return v; }
	//v-c
	template <typename O> inline vec<T, vecN> operator*=(const O& o) { for (uint32_t i(vecN); i--;) x[i] *= o; return *this; }
	template <typename O> inline vec<T, vecN> operator/=(const O& o) { double w = 1.0 / o; operator*=(w); return *this; }
	template <typename O> inline vec<T, vecN> operator* (const O& o) const { vec<T, vecN> v = (*this); v *= o; return v; }
	template <typename O> inline vec<T, vecN> operator/ (const O& o) const { O w = 1.0 / o; return (*this) * w; }
};

typedef struct vec<float, 2> vec2f;
typedef struct vec<double, 2> vec2d;
typedef struct vec<int8_t, 2>  vec2i8;
typedef struct vec<int16_t, 2> vec2i16;
typedef struct vec<int32_t, 2> vec2i;
typedef struct vec<int64_t, 2> vec2i64;