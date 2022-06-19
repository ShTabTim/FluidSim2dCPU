#pragma once
#include <cstdint>
#include <cmath>

template <typename T>
struct vec3 {
public:
	T x = 0, y = 0, z = 0;
	inline vec3() {}
	inline vec3(T x, T y, T z) { this->x = x; this->y = y; this->z = z; }
	template <typename O> inline vec3<T> operator=(const vec3<O>& o) { x = T(o.x); y = T(o.y); z = T(o.z); return *this; }

	//c
	inline T abs2() const { return x * x + y * y + z * z; }
	inline T abs() const { return (T)sqrt((double)(x * x + y * y + z * z)); }
	template <typename O> inline T dot(vec3<O>& o) const { return x * o.x + y * o.y + z * o.z; }
	//logical
	inline uint32_t operator()() const { return x || y || z; }
	inline uint32_t operator!() const { return !x && !y && !z; }
	template <typename O> inline uint32_t operator==(const vec3<O>& o) const { return x == o.x && y == o.y && z == o.z; }
	template <typename O> inline uint32_t operator!=(const vec3<O>& o) const { return x != o.x || y != o.y || z != o.z; }
	//v-v
	inline vec3<T> norm(T r) { T n = abs(); if (n) { n = r / n; x *= n; y *= n; z *= n; } return *this; }
	inline vec3<T> get_norm(T r) { T n = abs(); if (n) n = r / n; return vec3(x * n, y * n, z * n); }
	inline vec3<T> operator-() const { return vec3<T>(-x, -y, -z); }
	template <typename O> inline vec3<T> operator+=(const vec3<O>& o) { x += o.x; y += o.y; z += o.z; return *this; }
	template <typename O> inline vec3<T> operator-=(const vec3<O>& o) { x -= o.x; y -= o.y; z -= o.z; return *this; }
	template <typename O> inline vec3<T> operator+(const vec3<O>& o) const { return vec3<T>(x + o.x, y + o.y, z + o.z); }
	template <typename O> inline vec3<T> operator-(const vec3<O>& o) const { return vec3<T>(x - o.x, y - o.y, z - o.z); }
	template <typename O> inline vec3<T> operator*=(const vec3<O>& o) { x *= o.x; y *= o.y; z *= o.z; return *this; }
	template <typename O> inline vec3<T> operator/=(const vec3<O>& o) { x /= o.x; y /= o.y; z /= o.z; return *this; }
	template <typename O> inline vec3<T> operator*(const vec3<O>& o) const { return vec3<T>(x * o.x, y * o.y, z * o.z); }
	template <typename O> inline vec3<T> operator/(const vec3<O>& o) const { return vec3<T>(x / o.x, y / o.y, z / o.z); }
	//v-c
	template <typename O> inline vec3<T> operator*=(const O& o) { x *= o.x; y *= o.y; z *= o.z; return *this; }
	template <typename O> inline vec3<T> operator/=(const O& o) { double w = 1.0 / o; x *= w; y *= w; z *= w; return *this; }
	template <typename O> inline vec3<T> operator*(const O& o) const { return vec3<T>(x * o, y * o, z * o); }
	template <typename O> inline vec3<T> operator/(const O& o) const { double w = 1.0 / o; return vec3<T>(x * w, y * w, z * w); }
};

typedef struct vec3<float> vec3f;
typedef struct vec3<double> vec3d;
typedef struct vec3<int8_t>  vec3i8;
typedef struct vec3<int16_t> vec3i16;
typedef struct vec3<int32_t> vec3i;
typedef struct vec3<int64_t> vec3i64;