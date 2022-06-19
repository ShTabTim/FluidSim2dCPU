#pragma once
#include <cstdint>
#include <cmath>

template <typename T>
struct vec2 {
public:
	T x = 0, y = 0;
	inline vec2() {}
	inline vec2(T x, T y) { this->x = x; this->y = y; }
	template <typename O> inline vec2<T> operator=(const vec2<O>& o) { x = T(o.x); y = T(o.y); return *this; }

	//c
	inline T abs2() const { return x * x + y * y; }
	inline T abs() const { return (T) sqrt(x * x + y * y); }
	template <typename O> T dot(vec2<O> o) const { return (T) (x * o.x + y * o.y); }
	//logical
	inline uint32_t operator()() const { return x || y; }
	inline uint32_t operator!() const { return !x && !y; }
	template <typename O> inline uint32_t operator==(const vec2<O>& o) const { return x == o.x && y == o.y; }
	template <typename O> inline uint32_t operator!=(const vec2<O>& o) const { return x != o.x || y != o.y; }
	//v-v
	inline vec2<T> norm(T r) { T n = abs(); if (n) { n = r / n; x *= n; y *= n; } return *this; }
	inline vec2<T> get_norm(T r) { T n = abs(); if (n) n = r / n; return vec2(x * n, y * n); }
	inline vec2<T> operator-() const { return vec2<T>(-x, -y); }
	template <typename O> inline vec2<T> operator+=(const vec2<O>& o) { x += o.x; y += o.y; return *this; }
	template <typename O> inline vec2<T> operator-=(const vec2<O>& o) { x -= o.x; y -= o.y; return *this; }
	template <typename O> inline vec2<T> operator+(const vec2<O>& o) const { return vec2<T>(x + o.x, y + o.y); }
	template <typename O> inline vec2<T> operator-(const vec2<O>& o) const { return vec2<T>(x - o.x, y - o.y); }
	template <typename O> inline vec2<T> operator*=(const vec2<O>& o) { x *= o.x; y *= o.y; return *this; }
	template <typename O> inline vec2<T> operator/=(const vec2<O>& o) { x /= o.x; y /= o.y; return *this; }
	template <typename O> inline vec2<T> operator*(const vec2<O>& o) const { return vec2<T>(x * o.x, y * o.y); }
	template <typename O> inline vec2<T> operator/(const vec2<O>& o) const { return vec2<T>(x / o.x, y / o.y); }
	//v-c
	template <typename O> inline vec2<T> operator*=(const O& o) { x *= o; y *= o; return *this; }
	template <typename O> inline vec2<T> operator/=(const O& o) { double w = 1.0 / o; x *= w; y *= w; return *this; }
	template <typename O> inline vec2<T> operator*(const O& o) const { return vec2<T>(x * o, y * o); }
	template <typename O> inline vec2<T> operator/(const O& o) const { double w = 1.0 / o; return vec2<T>(x * w, y * w); }
};

typedef struct vec2<float> vec2f;
typedef struct vec2<double> vec2d;
typedef struct vec2<int8_t>  vec2i8;
typedef struct vec2<int16_t> vec2i16;
typedef struct vec2<int32_t> vec2i;
typedef struct vec2<int64_t> vec2i64;