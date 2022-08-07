#include <Windows.h>
#include "FundLibs/sh_rend_cpu/sh_win.h"
#include "FundLibs/math/vec2.h"
#include "FundLibs/key/keys.h"
#include <iostream>
#include <cassert>

#define CeilSizeX 4
#define CeilSizeY 4

#define W 128
#define H 128

#define ScrW W*CeilSizeX
#define ScrH H*CeilSizeY

#define Force 50

#define DeltaX 1.0f
#define Sigma 512.0f
#define WaterViscosity 1.006f * DeltaX * DeltaX
#define AirViscosity   0.151f * DeltaX * DeltaX
#define N 2

#define sum2d(i, j, operation) for (uint32_t i(W); i--;) for (uint32_t j(H); j--;) operation;
#define sum1d(i, operation) for (uint32_t i(W*H); i--;) operation;

typedef struct dynamic_data {
	void* data;
	bool is_operated;
	bool is_used;
	dynamic_data() {
		set_all_zero();
	}
	void init(bool is_operated, uint32_t bit_size) {
		if (!data)
			data = malloc(bit_size);
		else
			data = realloc(data, bit_size);
		this->is_operated = is_operated;
		is_used = true;
	}
	void finit() {
		if (data)
			free(data);
		set_all_zero();
	}
	void inline set_all_zero() {
		data = 0;
		is_operated = 0;
		is_used = 0;
	}
} dynamic_data;

typedef struct dynamic_data_hendler {
	dynamic_data* data_cell_set;
	uint32_t dynamic_set_size;
	dynamic_data_hendler(uint32_t dynamic_set_size = 128) {
		data_cell_set = 0;
		init(dynamic_set_size);
		for (uint32_t i(dynamic_set_size); i--;)
			operator()(i).set_all_zero();
	}
	void init(uint32_t dynamic_set_size) {
		this->dynamic_set_size = dynamic_set_size;
		if (!data_cell_set)
			data_cell_set = (dynamic_data*)malloc(dynamic_set_size * sizeof(dynamic_data));
		else
			data_cell_set = (dynamic_data*)realloc(data_cell_set, dynamic_set_size * sizeof(dynamic_data));
	}
	dynamic_data& operator()(uint32_t i) {
		assert(i < dynamic_set_size);
		return data_cell_set[i];
	}
	dynamic_data get(uint32_t i) {
		assert(i < dynamic_set_size);
		return data_cell_set[i];
	}
	uint32_t get_free_data_cell() {
		for (uint32_t i(dynamic_set_size); i--;)
			if (!data_cell_set[i].is_used)
				return i;
		return dynamic_set_size;
	}
	void clear_not_operated() {
		for (uint32_t i(dynamic_set_size); i--;)
			if (data_cell_set[i].is_used && !data_cell_set[i].is_operated)
				data_cell_set[i].finit();
	}
	void finit() {
		for (uint32_t i(dynamic_set_size); i--;)
			if (data_cell_set[i].is_used)
				data_cell_set[i].finit();
	}
} dynamic_data_hendler;
dynamic_data_hendler DDH(256);

template<typename T>
class sheet {
	uint32_t data_id;
	uint32_t sheet_amount = 256;
	sheet() {
		data_id = DDH.get_free_data_cell();
		DDH.data_cell_set[data_id].init(false, sizeof(T) * sheet_amount);
	}
	sheet(uint32_t data_id, bool is_operated = true) {
		DDH(this->data_id).finit();
		this->data_id = data_id;
		DDH.data_cell_set[data_id].init(is_operated, sizeof(T) * sheet_amount);
	}
	T& operator() (const uint32_t& i) {
		return ((T*)(DDH(data_id).data))[i];
	}
	T get(const uint32_t& i) {
		return ((T*)(DDH.get(data_id).data))[i];
	}
};

template<typename T>
class field {
public:
	uint32_t data_id;
	T zero;
	field() {
		zero = T();
		data_id = DDH.get_free_data_cell();
		DDH.data_cell_set[data_id].init(false, sizeof(T) * W * H);
	}
	field(uint32_t data_id, bool is_operated = true) {
		zero = T();
		this->data_id = data_id;
		DDH.data_cell_set[data_id].init(is_operated, sizeof(T) * W * H);
	}
	void init(uint32_t data_id, bool is_operated = true) {
		DDH(this->data_id).finit();
		this->data_id = data_id;
		DDH(data_id).init(is_operated, sizeof(T) * W * H);
	}



	T& operator()(const uint32_t& i) {
		return ((T*)(DDH(data_id).data))[i];
	}
	T& operator()(const uint32_t& x, const uint32_t& y) {
		return operator()(H * x + y);
	}
	T inline get(const uint32_t& i) const {
		return ((T*)(DDH.get(data_id).data))[i];
	}
	T inline get(const uint32_t& x, const uint32_t& y) const {
		if (x >= 0 && x < W && y >= 0 && y < H) return get(H * x + y);
		return zero;
	}



	field<T> inline operator=(const field<T>& f) {
		zero = f.zero;
		sum1d(i, operator()(i) = f.get(i));
		return (*this);
	}



	template<typename M>
	field<T> inline operator*(const M& o) const {
		field<T> f;
		f.zero = zero * o;
		sum1d(i, f(i) = get(i) * o);
		return f;
	}
	template<typename M>
	field<T> inline operator*(const field<M>& o) const {
		field<T> f;
		f.zero = zero * o.zero;
		sum1d(i, f(i) = get(i) * o.get(i));
		return f;
	}
	template<typename M>
	field<T> inline operator/(const M& o) const {
		field<T> f;
		f.zero = zero / o;
		sum1d(i, f(i) = get(i) / o);
		return f;
	}
	template<typename M>
	field<T> inline operator/(const field<M>& o) const {
		field<T> f;
		f.zero = zero / o.zero;
		sum1d(i, f(i) = get(i) / o.get(i));
		return f;
	}



	field<T> inline operator+(const field<T>& o) const {
		field<T> f;
		f.zero = zero + o.zero;
		sum1d(i, f(i) = o.get(i) + get(i));
		return f;
	}
	field<T> inline operator+=(const field<T>& o) {
		zero += o.zero;
		sum1d(i, operator()(i) += o.get(i));
		return *this;
	}

	field<T> inline operator+(const T& o) const {
		field<T> f;
		f.zero = zero + o;
		sum1d(i, f(i) = o + get(i));
		return f;
	}
	field<T> inline operator+=(const T& o) {
		zero += o;
		sum1d(i, operator()(i) += o);
		return *this;
	}



	field<T> inline operator-() const {
		field<T> f;
		f.zero = -zero;
		sum1d(i, f(i) = -get(i));
		return f;
	}
	field<T> inline operator-(const field<T>& o) const {
		return (*this) + (-o);
	}
	field<T> operator-=(const field<T>& o) {
		zero -= o.zero;
		sum1d(i, operator()(i) -= o.get(i));
		return *this;
	}
	field<T> inline operator-(const T& o) const {
		return (*this) + (-o);
	}
	field<T> operator-=(const T& o) {
		zero -= o;
		sum1d(i, operator()(i) -= o);
		return *this;
	}
};

//template<typename T>
//class mat2 {
//	T aij[2 * 2];
//	mat() {
//		for (uint32_t i(2 * 2); i--;)aij[i] = 0;
//	}
//	mat(T aii) {
//		for (uint32_t i(2); i--;)
//			for (uint32_t j(2); j--;)
//				operator()(i, j) = aii * (i == j);
//	}
//	T& operator()(const uint32_t& i, const uint32_t& j) const {
//		return aij[i * 2 + j];
//	}
//
//	mat operator-(const mat& o) {}
//};


template<typename T>
class dnfield_dtn {
public:
	field<T>* dnf_dtn = 0;
	dnfield_dtn(uint32_t data_id, uint32_t K) {
		if (!dnf_dtn)
			dnf_dtn = new field<T>[K];
		for (uint32_t i(K); i--;)
			dnf_dtn[i].init(data_id + i);
	}
	field<T>& operator()(const uint32_t& n) {
		return dnf_dtn[n];
	}
};



template<typename T>
T inline delta_x(const field<T>& f, uint32_t x, uint32_t y) {
	return (f.get(x + 1, y) - f.get(x - 1, y)) * 0.5f / DeltaX;
}

template<typename T>
T inline delta_y(const field<T>& f, uint32_t x, uint32_t y) {
	return (f.get(x, y + 1) - f.get(x, y - 1)) * 0.5f / DeltaX;
}

template<typename T>
field<T> inline delta_x(const field<T>& _f) {
	field<T> f;
	f.zero = _f.zero;
	for (uint32_t i(W); i--;)
		for (uint32_t j(H); j--;)
			f(i, j) = delta_x(_f, i, j);
	return f;
}

template<typename T>
field<T> inline delta_y(const field<T>& _f) {
	field<T> f;
	f.zero = _f.zero;
	for (uint32_t i(W); i--;)
		for (uint32_t j(H); j--;)
			f(i, j) = delta_y(_f, i, j);
	return f;
}

template<typename T>
T inline delta2_x(const field<T>& f, uint32_t x, uint32_t y) {
	return (f.get(x + 1, y) + f.get(x - 1, y) - f.get(x, y) * 2) / (DeltaX * DeltaX);
}
template<typename T>
T inline delta2_y(const field<T>& f, uint32_t x, uint32_t y) {
	return (f.get(x, y + 1) + f.get(x, y - 1) - f.get(x, y) * 2) / (DeltaX * DeltaX);
}
template<typename T>
field<T> inline delta2_x(const field<T>& _f) {
	field<T> f;
	f.zero = _f.zero;
	for (uint32_t i(W); i--;)
		for (uint32_t j(H); j--;)
			f(i, j) = delta2_x(_f, i, j);
	return f;
}
template<typename T>
field<T> inline delta2_y(const field<T>& _f) {
	field<T> f;
	f.zero = _f.zero;
	for (uint32_t i(W); i--;)
		for (uint32_t j(H); j--;)
			f(i, j) = delta2_y(_f, i, j);
	return f;
}

//template<typename T>
//matn<T> inline nab_tenz_prod_vec(const field<vec2<T>>& f_v, uint32_t x, uint32_t y) {
//	return matn<T>();// delta_x(f_b, x, y)* (f_v.get(x, y)).x + delta_y(f_b, x, y) * (f_v.get(x, y)).y;
//}

template<typename T, typename M>
field<T> inline dot(const field<vec<T, N>>& f_v, const vec<M, N>& v) {
	field<T> f;
	f.zero = f_v.zero.dot(v);
	sum1d(i, f(i) = f_v.get(i).dot(v))
		return f;
}

template<typename T, typename M>
field<T> inline dot(const field<vec<T, N>>& f_v0, const field<vec<M, N>>& f_v1) {
	field<T> f;
	f.zero = f_v0.zero.dot(f_v1.zero);
	sum1d(i, f(i) = f_v0.get(i).dot(f_v1.get(i)))
		return f;
}

template<typename T>
T inline vnab_b(const field<vec2f>& f_v, const field<T>& f_b, uint32_t x, uint32_t y) {
	return delta_x(f_b, x, y) * (f_v.get(x, y)).get(0) + delta_y(f_b, x, y) * (f_v.get(x, y)).get(1);
}

template<typename T>
field<T> inline vnab_b(const field<vec2f>& f_v, const field<T>& f_b) {
	field<T> f;
	f.zero = f_b.zero;
	for (uint32_t i(W); i--;)
		for (uint32_t j(H); j--;)
			f(i, j) = vnab_b(f_v, f_b, i, j);
	return f;
}

template<typename T>
T inline lapl(const field<T>& f, uint32_t x, uint32_t y) {
	return delta2_x(f, x, y) + delta2_y(f, x, y);
}

template<typename T>
field<T> inline lapl(const field<T>& _f) {
	field<T> f;
	f.zero = _f.zero;
	for (uint32_t i(W); i--;)
		for (uint32_t j(H); j--;)
			f(i, j) = lapl(_f, i, j);
	return f;
}

template<typename T>
T inline div(const field<vec<T, N>>& v, uint32_t x, uint32_t y) {
	return delta_x(v, x, y).get(0) + delta_y(v, x, y).get(1);
}

template<typename T>
field<T> inline div(const field<vec<T, N>>& v) {
	field<T> f;
	f.zero = v.zero.get(0) + v.zero.get(1);
	for (uint32_t i(W); i--;)
		for (uint32_t j(H); j--;)
			f(i, j) = div(v, i, j);
	return f;
}

template<typename T>
vec<T, N> inline grad(const field<T>& s, uint32_t x, uint32_t y) {
	return vec<T, N>(delta_x(s, x, y), delta_y(s, x, y));
}

template<typename T>
field<vec<T, N>> inline grad(const field<T>& s) {
	field<vec<T, N>> f;
	f.zero = vec<T, N>(s.zero);
	for (uint32_t i(W); i--;)
		for (uint32_t j(H); j--;)
			f(i, j) = grad(s, i, j);
	return f;
}

field<vec2f> force(0);

field<float> pressure_W(1);
field<float> new_pressure(2);
field<float> pressure(3);

field<float> old_xi(4);
field<float> new_xi(5);
dnfield_dtn<float> xi_field(6, 4);

field<vec2f> old_velosity(10);
field<vec2f> new_velosity(11);
dnfield_dtn<vec2f> vel_field(12, 4);


field<float> inline sqrt(const field<float>& x) {
	field<float> f;
	sum1d(i, f(i) = sqrt(x.get(i)));
	return f;
}

template<typename T>
field<T> inline alpha(const T& WaterAlpha, const T& AirAlpha) {
	return xi_field(0) * (WaterAlpha - AirAlpha) + AirAlpha;
}

template<typename T>
field<T> inline dalpha_dt(const T& WaterAlpha, const T& AirAlpha) {
	return xi_field(1) * (WaterAlpha - AirAlpha);
}

template<typename T>
field<T> inline d2alpha_dt2(const T& WaterAlpha, const T& AirAlpha) {
	return xi_field(2) * (WaterAlpha - AirAlpha);
}

field<vec2f> fg(field<vec2f> v) {
	return v * div(v);
}

field<vec2f> dfg_dt(field<vec2f> v, field<vec2f> dv_dt) {
	return v * div(dv_dt) + dv_dt * div(v);
}

field<vec2f> d2fg_dt2(field<vec2f> v, field<vec2f> dv_dt, field<vec2f> d2v_dt2) {
	return v * div(d2v_dt2) + dv_dt * div(dv_dt) * 2 + d2v_dt2 * div(v);
}

float pressure_mute(uint32_t i, uint32_t j) {
	return (
		pressure.get(i + 1, j) + pressure.get(i - 1, j) +
		pressure.get(i, j + 1) + pressure.get(i, j - 1) +
		-pressure_W.get(i, j) / (DeltaX * DeltaX)) * 0.5f / N;
}

void comp_pressure(uint32_t k, float t) {
	pressure_W = div(vel_field(0));
	uint32_t f;
	for (uint32_t gh(k); gh--;) {
		for (uint32_t i(W); i--;)
			for (uint32_t j(H); j--;) {
				new_pressure(i, j) = pressure(i, j) * t + (1 - t) * pressure_mute(i, j);
			}
		f = pressure.data_id;
		pressure.data_id = new_pressure.data_id;
		new_pressure.data_id = f;
	}
}

class MainRenderer : public sh_dwaw_win_cpu {
	bool sh_init() {
		AppName = L"CPURenderer";
		DDH.clear_not_operated();

		xi_field(0).zero = 1;
		for (uint32_t i(W); i--;)
			for (uint32_t j(H); j--;) {
				pressure(i, j) = 0;
				vel_field(0)(i, j) = vec2f(0, 0);
				xi_field(0)(i, j) = 1;
			}
		return 1;
	}

	float mami(float x) {
		return ((x < 0) - (x >= 0)) * (exp(x * ((x < 0) - (x >= 0))) - 1);
	}

	float prot(float x) {
		return (x > 0 && x < 2) * (2 - x) * x;
	}

	float sqr(float d) { return d * d; }

	bool sh_loop(double dt) {
		key_loop(get_hwnd());
		static vec2f ptn = vec2f(uint32_t(get_x() / CeilSizeX), uint32_t(get_y() / CeilSizeY)), ptold = ptn;
		ptn = vec2f(uint32_t(get_x() / CeilSizeX), uint32_t(get_y() / CeilSizeY));
		if (get_key(VK_LBUTTON).held && (get_x() >= 0 && get_x() < ScrW && get_y() >= 0 && get_y() < ScrH))
			sum2d(i, j, force(i, j) += (ptn - ptold).get_norm(Force * exp(-0.02f * (sqr(ptn.get(0) - i) + sqr(ptn.get(1) - j)))));
		if (get_key(VK_RBUTTON).held && (get_x() >= 0 && get_x() < ScrW && get_y() >= 0 && get_y() < ScrH))
			sum2d(i, j, xi_field(0)(i, j) += ((sqr(ptn.get(0) - i) + sqr(ptn.get(1) - j)) <= 36) * (1 - 2 * get_key('A').held));
		if (get_key('Q').held && (get_x() >= 0 && get_x() < ScrW && get_y() >= 0 && get_y() < ScrH))
			sum2d(i, j, pressure(i, j) += 300 * pow(0.05f, max(0, -10 + (sqr(ptn.get(0) - i) + sqr(ptn.get(1) - j)))));
		if (get_key('W').held && (get_x() >= 0 && get_x() < ScrW && get_y() >= 0 && get_y() < ScrH))
			sum2d(i, j, pressure(i, j) -= 300 * pow(0.05f, max(0, -10 + (sqr(ptn.get(0) - i) + sqr(ptn.get(1) - j)))));
		ptold = ptn;

		//////////////////////////////////////////////////////////////

		for (uint32_t i(W); i--;)
			for (uint32_t j(H); j--;) {
				xi_field(0)(i, j) = max(0, min(1, xi_field(0)(i, j)));
			}

		old_xi = xi_field(0);
		old_velosity = vel_field(0);
		field<float> inv_abs_grad_xi;
		sum2d(i, j, inv_abs_grad_xi(i, j) = (grad(old_xi, i, j).abs() * 0xFFFFFFFF > 1) ? 1.0f / grad(old_xi, i, j).abs() : 0xFFFFFFFF);
		field<vec2f> surf_force;
		surf_force = grad(old_xi) * div(grad(old_xi) * inv_abs_grad_xi);
		//surf_force = fg(grad(old_xi));
		surf_force = -(surf_force)*Sigma;
		sum2d(i, j, surf_force(i, j).norm(min(256, surf_force.get(i, j).abs())));





		field<float> viscosity = alpha(WaterViscosity, AirViscosity);

		new_velosity = vel_field(0);
		vel_field(1) = -vnab_b(old_velosity, old_velosity) + force + lapl(old_velosity) * viscosity;
		vel_field(0) += vel_field(1) * dt;
		vel_field(0) += surf_force * dt;
		comp_pressure(40, 0.05f);
		vel_field(0) -= grad(pressure);
		vel_field(1) = (vel_field(0) - new_velosity) * (1.0f / dt);

		xi_field(1) = -dot(old_velosity, grad(old_xi));
		xi_field(0) += xi_field(1) * dt;

		for (uint32_t i(W); i--;)
			for (uint32_t j(H); j--;) {
				vel_field(0)(i, j).norm((vel_field(0)(i, j).abs2() > 0.0f) * min(20, vel_field(0)(i, j).abs()));
				xi_field(0)(i, j) = max(0, min(1, xi_field(0)(i, j)));
			}


		sum2d(i, j, inv_abs_grad_xi(i, j) = (inv_abs_grad_xi(i, j) < 0xFFFFFF) ? inv_abs_grad_xi(i, j) : 0xFFFFFF);

		field<float> dinv_abs_grad_xi_dt;
		dinv_abs_grad_xi_dt = -dot(grad(xi_field(0)), grad(xi_field(1))) * inv_abs_grad_xi * inv_abs_grad_xi * inv_abs_grad_xi;




		field<float> dviscosity_dt = dalpha_dt(WaterViscosity, AirViscosity);

		field<vec2f> dsurf_force_dt;
		//dsurf_force_dt = grad(xi_field(1)) * div(grad(old_xi) * inv_abs_grad_xi) + grad(old_xi) * div((grad(xi_field(1)) - grad(old_xi) * dot(grad(old_xi), grad(xi_field(1))) * inv_abs_grad_xi * inv_abs_grad_xi) * inv_abs_grad_xi);
		dsurf_force_dt =
			grad(xi_field(1)) * div(grad(old_xi) * inv_abs_grad_xi) +
			grad(old_xi) * div(
				grad(xi_field(1)) * inv_abs_grad_xi +
				grad(old_xi) * dinv_abs_grad_xi_dt
			);
		//dsurf_force_dt = grad(xi_field(1)) * dfg_dt(grad(old_xi), grad(xi_field(1)));
		dsurf_force_dt = -(dsurf_force_dt)*Sigma;

		sum2d(i, j, dsurf_force_dt(i, j).norm(min(256, dsurf_force_dt.get(i, j).abs())));

		new_velosity = vel_field(0);
		vel_field(2) = -(vnab_b(vel_field(1), old_velosity) + vnab_b(old_velosity, vel_field(1))) + lapl(vel_field(1)) * viscosity + lapl(old_velosity) * dviscosity_dt;
		vel_field(0) += vel_field(2) * (0.5f * dt * dt);
		vel_field(0) += dsurf_force_dt * (0.5f * dt * dt);
		comp_pressure(40, 0.05f);
		vel_field(0) -= grad(pressure);
		vel_field(2) = (vel_field(0) - new_velosity) * (2.0f / (dt * dt));

		xi_field(2) = -(dot(old_velosity, grad(xi_field(1))) + dot(vel_field(1), grad(old_xi)));
		xi_field(0) += xi_field(2) * (0.5f * dt * dt);

		for (uint32_t i(W); i--;)
			for (uint32_t j(H); j--;) {
				vel_field(0)(i, j).norm((vel_field(0)(i, j).abs2() > 0.0f) * min(20, vel_field(0)(i, j).abs()));
				xi_field(0)(i, j) = max(0, min(1, xi_field(0)(i, j)));
			}

		field<float> d2inv_abs_grad_xi_dt2;
		d2inv_abs_grad_xi_dt2 = -(dot(grad(xi_field(0)), grad(xi_field(2))) + dot(grad(xi_field(1)), grad(xi_field(1)))) * dinv_abs_grad_xi_dt * 3;





		field<float> d2viscosity_dt2 = d2alpha_dt2(WaterViscosity, AirViscosity);

		field<vec2f> d2surf_force_dt2;
		//d2surf_force_dt2 = grad(xi_field(1)) * d2fg_dt2(grad(old_xi), grad(xi_field(1)), grad(xi_field(2))) + grad(xi_field(2)) * dfg_dt(grad(old_xi), grad(xi_field(1)));
		d2surf_force_dt2 =
			grad(xi_field(2)) * div(grad(old_xi) * inv_abs_grad_xi) +
			grad(xi_field(1)) * div(
				grad(xi_field(1)) * inv_abs_grad_xi +
				grad(old_xi) * dinv_abs_grad_xi_dt
			) * 2 +
			grad(old_xi) * div(
				grad(xi_field(2)) * inv_abs_grad_xi +
				grad(xi_field(1)) * dinv_abs_grad_xi_dt * 2 +
				grad(old_xi) * d2inv_abs_grad_xi_dt2
			);
		d2surf_force_dt2 = -(d2surf_force_dt2)*Sigma;

		sum2d(i, j, d2surf_force_dt2(i, j).norm(min(256, d2surf_force_dt2.get(i, j).abs())));

		new_velosity = vel_field(0);
		vel_field(3) = -(vnab_b(vel_field(2), old_velosity) + vnab_b(vel_field(1), vel_field(1)) * 2 + vnab_b(old_velosity, vel_field(2))) + lapl(vel_field(2)) * viscosity + lapl(vel_field(1)) * dviscosity_dt * 2 + lapl(old_velosity) * d2viscosity_dt2;
		vel_field(0) += vel_field(3) * (dt * dt * dt / 6.0f);
		vel_field(0) += d2surf_force_dt2 * (dt * dt * dt / 6.0f);
		comp_pressure(40, 0.05f);
		vel_field(0) -= grad(pressure);
		vel_field(3) = (vel_field(0) - new_velosity) * (6.0f / (dt * dt * dt));

		xi_field(3) = -(dot(vel_field(2), grad(old_xi)) + dot(vel_field(1), grad(xi_field(1))) * 2 + dot(old_velosity, grad(xi_field(2))));
		xi_field(0) += xi_field(3) * (dt * dt * dt / 6.0f);

		for (uint32_t i(W); i--;)
			for (uint32_t j(H); j--;) {
				vel_field(0)(i, j).norm((vel_field(0)(i, j).abs2() > 0.0f) * min(20, vel_field(0)(i, j).abs()));
				xi_field(0)(i, j) = max(0, min(1, xi_field(0)(i, j)));
			}

		//////////////////////////////////////////////////////////////
		for (uint32_t i(W); i--;)
			for (uint32_t j(H); j--;) {
				if (force(i, j).abs() > 1)
					force(i, j) *= 0.2f;
				else force(i, j) = vec2f(0, 0);
			}

		static uint32_t dr = 0;
		for (uint8_t k(0); k < 10; k++)
			if (get_key('0' + k).held) dr = k;
		if (dr == 2 || dr == 3 || dr == 4) fill_rect(0, 0, get_dr_w(), get_dr_h(), 100, 0, 0);
		for (uint32_t i(W); i--;)
			for (uint32_t j(H); j--;) {
				switch (dr) {
				case 0:
				{
					float fhj = 2 * xi_field(0)(i, j) * ((xi_field(0)(i, j) > 0.6f) * (!get_key('Z').held) + get_key('Z').held);
					fill_rect(i * CeilSizeX, j * CeilSizeY, (i + 1) * CeilSizeX, (j + 1) * CeilSizeY, 255 * prot(fhj), 255 * prot(fhj - 1), 0);// 255 * prot(fhj - 2));
				}
				break;
				case 1:
				{
					float pp = 1.5f + 1.5f * mami(0.1f * pressure(i, j));
					fill_rect(i * CeilSizeX, j * CeilSizeY, (i + 1) * CeilSizeX, (j + 1) * CeilSizeY, 255 * prot(pp), 255 * prot(pp - 1), 255 * prot(pp - 2));
					break;
				}
				case 2:
					draw_line((i + 0.5f) * CeilSizeX, (j + 0.5f) * CeilSizeY, (i + 0.5f) * CeilSizeX + vel_field(0)(i, j)(0), (j + 0.5f) * CeilSizeY + vel_field(0)(i, j)(1), 0xFF, 0xFF, 0xFF);
					break;
				case 3:
					draw_line((i + 0.5f) * CeilSizeX, (j + 0.5f) * CeilSizeY, (i + 0.5f) * CeilSizeX + force(i, j)(0), (j + 0.5f) * CeilSizeY + force(i, j)(1), 0xFF, 0xFF, 0xFF);
					break;
				case 4:
					draw_line((i + 0.5f) * CeilSizeX, (j + 0.5f) * CeilSizeY, (i + 0.5f) * CeilSizeX + surf_force(i, j)(0), (j + 0.5f) * CeilSizeY + surf_force(i, j)(1), 0xFF, 0xFF, 0xFF);
					break;
				}
			}

		DDH.clear_not_operated();

		return 1;
	}
	bool sh_finit() {
		DDH.finit();
		return 1;
	}
};

int main() {
	MainRenderer simulation;
	if (simulation.init(0, ScrW, ScrH, ScrW, ScrH))
		simulation.run();
	return 0;
}