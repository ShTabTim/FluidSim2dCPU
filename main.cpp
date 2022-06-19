#include <Windows.h>
#include "FundLibs/sh_rend_cpu/sh_win.h"
#include "FundLibs/math/vec2.h"
#include "FundLibs/key/keys.h"
#include <iostream>

#define CeilSizeX 4
#define CeilSizeY 4

#define W 256
#define H 64

#define ScrW W*CeilSizeX
#define ScrH H*CeilSizeY

#define Force 300

#define Pressure 2
#define Density 2
#define Viscosity 2.1f*(W+H)/128
#define Viscosity_conc 0.3f*(W+H)/128
#define Accuracy 4
#define Gravity 100

template<typename T>
class field {
public:
	T* _field;
	T zero = T();
	field() {
		_field = new T[W * H];
	}
	T& operator()(uint32_t x, uint32_t y) {
		if (x >= 0 && x < W && y >= 0 && y < H) return _field[H * x + y];
		return zero;
	}
	void set_zero(T zero) { this->zero = zero; }
};

template<typename T>
T inline delta_x(field<T> f, uint32_t x, uint32_t y) {
	return (f(x + 1, y) - f(x - 1, y)) * 0.5f;
	return f(x + 1, y) - f(x, y);
}

template<typename T>
T inline delta_y(field<T> f, uint32_t x, uint32_t y) {
	return (f(x, y + 1) - f(x, y - 1)) * 0.5f;
	return f(x, y + 1) - f(x, y);
}

template<typename T>
T inline delta2_x(field<T> f, uint32_t x, uint32_t y) {
	return f(x + 1, y) + f(x - 1, y) - f(x, y) * 2;
}
template<typename T>
T inline delta2_y(field<T> f, uint32_t x, uint32_t y) {
	return f(x, y + 1) + f(x, y - 1) - f(x, y) * 2;
}

template<typename T>
T inline vnab_b(field<vec2f> f_v, field<T> f_b, uint32_t x, uint32_t y) {
	return delta_x(f_b, x, y) * f_v(x, y).x + delta_y(f_b, x, y) * f_v(x, y).y ;
}

template<typename T>
T inline lapl(field<T> f, uint32_t x, uint32_t y) {
	return delta2_x(f, x, y) + delta2_y(f, x, y);
}


field<uint32_t> box;
field<float> pressure;
field<float> col_conc;
field<vec2f> velosity;
field<vec2f> force;

field<float> new_pressure;
field<float> new_col_conc;
field<vec2f> new_velosity;

float div(field<vec2f> v, uint32_t x, uint32_t y) {
	return delta_x(v, x, y).x + delta_y(v, x, y).y;
}
vec2f grad(field<float> s, uint32_t x, uint32_t y) {
	return vec2f(delta_x(s, x, y), delta_y(s, x, y));
}


void inline jacobi_pressure(uint32_t x, uint32_t y, float alpha, float beta) {
#define GET_M_BOX_PRESSURE(x0, y0) ((!box(x0, y0))?(pressure(x0, y0)):( pressure(x, y)))
#define GET_M_BOX_VELOSITY(x0, y0) ((!box(x0, y0))?(velosity(x0, y0)):(-velosity(x, y)))
	new_pressure(x, y) = (
		GET_M_BOX_PRESSURE(x+1, y  ) + 
		GET_M_BOX_PRESSURE(x-1, y  ) + 
		GET_M_BOX_PRESSURE(x,   y+1) + 
		GET_M_BOX_PRESSURE(x,   y-1) + (
			GET_M_BOX_VELOSITY(x+1, y  ).x -
			GET_M_BOX_VELOSITY(x-1, y  ).x +
			GET_M_BOX_VELOSITY(x,   y+1).y -
			GET_M_BOX_VELOSITY(x,   y-1).y) * 0.5f * alpha) / beta;
#undef GET_M_BOX_PRESSURE
#undef GET_M_BOX_VELOSITY
}

void compute_pressure(uint32_t N) {
	for (uint32_t i(N); i--;) {
		for (uint32_t i(W); i--;)
			for (uint32_t j(H); j--;) {
				if (!box(i, j)) jacobi_pressure(i, j, -Pressure, 4);
				if(get_key('A').held)pressure(i, j) = new_pressure(i, j);
			}
		for(uint32_t i(W);i--;)
			for(uint32_t j(H);j--;)
				pressure(i, j) = new_pressure(i, j);
	}
}

class MainRenderer : public sh_dwaw_win_cpu {
    bool sh_init() {
        AppName = L"CPURenderer";

		for (uint32_t i(W); i--;)
			for (uint32_t j(H); j--;)
				pressure(i, j) = 1;
		box.set_zero(255);
		for(uint32_t i(W);i--;)
			for(uint32_t j(H);j--;)
				box(i, j) = 0;
		for (uint32_t i(W); i--;) {
			box(i, 0) = 255;
			box(i, H - 1) = 255;
		}
		for (uint32_t j(H); j--;) {
			box(0, j) = 255;
			box(W - 1, j) = 255;
		}
        return 1;
    }
float sqr(float d) { return d * d; }
    bool sh_loop(double dt) {
		key_loop(get_hwnd());
		static vec2f ptn = vec2f(uint32_t(get_x() / CeilSizeX), uint32_t(get_y() / CeilSizeY)), ptold = ptn;
		ptn = vec2f(uint32_t(get_x() / CeilSizeX), uint32_t(get_y() / CeilSizeY));
		static uint8_t block = 255;
		for (uint8_t i(10); i--;)
			if (get_key(VK_NUMPAD0 + i).held)
				block = i*28;
		if(get_key('V').held && (get_x()>=0 && get_x() < ScrW && get_y()>=0 && get_y()<ScrH)){
			box(ptn.x, ptn.y) = block;
		}
		if (get_key('N').held && (get_x()>=0 && get_x() < ScrW && get_y()>=0 && get_y()<ScrH))
			for (uint32_t i(1); i < W - 1; i++)
				for (uint32_t j(1); j < H - 1; j++) {
					if (!box(i, j))force(i, j) += (ptn - ptold).get_norm(Force) * exp(- 0.02f*(sqr(ptn.x - i) + sqr(ptn.y - j)));
				}
		if (get_key('B').held && (get_x()>=0 && get_x() < ScrW && get_y()>=0 && get_y()<ScrH))
			for (uint32_t i(1); i < W - 1; i++)
				for (uint32_t j(1); j < H - 1; j++) {
					if (!box(i, j))col_conc(i, j) += 100*pow(0.05f, max(0, -10 + (sqr(ptn.x - i) + sqr(ptn.y - j))));
				}
		if (get_key('Q').held && (get_x()>=0 && get_x() < ScrW && get_y()>=0 && get_y()<ScrH))
			for (uint32_t i(1); i < W - 1; i++)
				for (uint32_t j(1); j < H - 1; j++) {
					if (!box(i, j))pressure(i, j) += 1000*pow(0.05f, max(0, -10 + (sqr(ptn.x - i) + sqr(ptn.y - j))));
				}
		if (get_key('W').held && (get_x()>=0 && get_x() < ScrW && get_y()>=0 && get_y()<ScrH))
			for (uint32_t i(1); i < W - 1; i++)
				for (uint32_t j(1); j < H - 1; j++) {
					if (!box(i, j))pressure(i, j) -= 1000*pow(0.05f, max(0, -10 + (sqr(ptn.x - i) + sqr(ptn.y - j))));
				}
		ptold = ptn;

		//dt /= Accuracy;
		//for (uint32_t k(Accuracy); k--;) {
			for (uint32_t i(W); i--;)
				for (uint32_t j(H); j--;) {
					if (!box(i, j)) {
						new_velosity(i, j) = velosity(i, j);
						new_col_conc(i, j) = col_conc(i, j);
					} else {
						new_velosity(i, j) = vec2f(0, 0);
						new_col_conc(i, j) = 0;
					}
				}
			for (uint32_t i(W); i--;)
				for (uint32_t j(H); j--;) {
					if (!box(i, j))new_velosity(i, j) += (-vnab_b(velosity, velosity, i, j) + lapl(velosity, i, j) * Viscosity + force(i, j)) * dt;
				}
			for (uint32_t i(W); i--;)
				for (uint32_t j(H); j--;)
					velosity(i, j) = new_velosity(i, j);

			compute_pressure(30);

			for (uint32_t i(W); i--;)
				for (uint32_t j(H); j--;) {
					if (!box(i, j)) {
#define GET_PRESSURE(x0, y0, x1, y1) ((!box(x0, y0))?pressure(x0, y0):pressure(x1, y1))
						new_velosity(i, j) -= vec2f(
							GET_PRESSURE(i+1, j, i, j) - 
							GET_PRESSURE(i-1, j, i, j), 
							GET_PRESSURE(i, j+1, i, j) - 
							GET_PRESSURE(i, j-1, i, j)) * 0.5f / Density;
#undef GET_PRESSURE
						new_col_conc(i, j) += (-vnab_b(velosity, col_conc, i, j) + Viscosity_conc * lapl(col_conc, i, j)) * dt;
						if (force(i, j).abs() > 1) force(i, j) *= 0.9f;
						else force(i, j) = vec2f(0, 0);
					}
				}
			for (uint32_t i(W); i--;)
				for (uint32_t j(H); j--;)
					col_conc(i, j) = max(0, min(240, new_col_conc(i, j)));
			for (uint32_t i(W); i--;)
				for (uint32_t j(H); j--;)
					if (!box(i, j))velosity(i, j) = new_velosity(i, j).get_norm(min(20, new_velosity(i, j).abs()));
		//}
		//dt *= Accuracy;
		static uint32_t dr = 0;
		for (uint8_t k(0); k < 10; k++)
			if (get_key('0' + k).held) dr = k;
		if(dr == 2 || dr == 3) fill_rect(0, 0, get_dr_w(), get_dr_h(), 100, 0, 0);
		for (uint32_t i(W); i--;)
			for (uint32_t j(H); j--;) {
				switch (dr) {
				case 0:
					fill_rect(i*CeilSizeX, j*CeilSizeY, (i + 1)*CeilSizeX, (j + 1)*CeilSizeY, 0xFF & uint32_t(col_conc(i, j)), (0xFF00 & uint32_t(col_conc(i, j))) >> 8, (0xFF0000 & uint32_t(col_conc(i, j))) >> 16);
					break;
				case 1:
					fill_rect(i*CeilSizeX, j*CeilSizeY, (i + 1)*CeilSizeX, (j + 1)*CeilSizeY, min(255, uint32_t(10*pressure(i, j))), min(255, uint32_t(10*pressure(i, j)) >> 8), min(255, uint32_t(10*pressure(i, j)) >> 16));
					break;
				case 2:
					draw_line((i+0.5f)*CeilSizeX, (j+0.5f)*CeilSizeY, (i+0.5f)*CeilSizeX + velosity(i, j).x, (j + 0.5f)*CeilSizeY + velosity(i, j).y, 0xFF, 0xFF, 0xFF);
					break;
				case 3:
					draw_line((i+0.5f)*CeilSizeX, (j+0.5f)*CeilSizeY, (i+0.5f)*CeilSizeX + force(i, j).x, (j + 0.5f)*CeilSizeY + force(i, j).y, 0xFF, 0xFF, 0xFF);
					break;
				}
				if(box(i, j))fill_rect(i*CeilSizeX, j*CeilSizeY, (i + 1)*CeilSizeX, (j + 1)*CeilSizeY, 0xFF&box(i, j), 0xFF&(box(i, j)>>8), 0xFF&(box(i, j)>>16));
			}
		
        return 1;
    }
    bool sh_finit() {
        return 1;
    }
};

int main() {
	MainRenderer simulation;
	if (simulation.init(0, ScrW, ScrH, ScrW, ScrH))
		simulation.run();
	return 0;
}