#include "Simulation.h"
#include <algorithm>
#include <memory>
#include <string>
#include <Windows.h>

Simulation::Simulation(Graphics& gfx, Float viscosity, Float density,
	const Vec2I& dim, Float delta_time)
	:
	viscosity_(viscosity),
	density_(density),
	dim_(dim),
	nx(gfx.ScreenWidth),
	ny(gfx.ScreenHeight),
	nc(nx * ny),
	nbf(nc * sizeof(Float)),
	nbb(nc * sizeof(bool)),
	dx(dim_.x / (Float)nx),
	dy(dim_.y / (Float)ny),
	dt(delta_time),
	viscosity_qf(mm_set(viscosity, viscosity, viscosity, viscosity)),
	density_qf(mm_set(density, density, density, density)),
	dx_qf(mm_set(dx, dx, dx, dx)),
	dy_qf(mm_set(dy, dy, dy, dy)),
	dx2_qf(mm_set(dx * dx, dx * dx, dx * dx, dx * dx)),
	dy2_qf(mm_set(dy * dy, dy * dy, dy * dy, dy * dy)),
	dt_qf(mm_set(dt, dt, dt, dt)),
	gfx_(gfx)
{
	p = (Float*)_aligned_malloc(nbf, 64);
	pn = (Float*)_aligned_malloc(nbf, 64);
	u = (Float*)_aligned_malloc(nbf, 64);
	un = (Float*)_aligned_malloc(nbf, 64);
	v = (Float*)_aligned_malloc(nbf, 64);
	vn = (Float*)_aligned_malloc(nbf, 64);
	b_ = (Float*)_aligned_malloc(nbf, 64);
	is_solid = (bool*)_aligned_malloc(nbb, 64);

	memset(p, 0, nbf);
	memset(pn, 0, nbf);
	memset(u, 0, nbf);
	memset(un, 0, nbf);
	memset(v, 0, nbf);
	memset(vn, 0, nbf);
	memset(b_, 0, nbf);
	memset(is_solid, 0, nbb);

	// initialize the field
	InitField();
}

Simulation::~Simulation()
{
	_aligned_free(p);
	_aligned_free(pn);
	_aligned_free(u);
	_aligned_free(un);
	_aligned_free(v);
	_aligned_free(vn);
	_aligned_free(b_);
	_aligned_free(is_solid);
}

void Simulation::Step()
{
	//#pragma omp parallel
	{
		// reset boundary conditions
		ResetBoundaryConditions();

		// precalculations
		const Float dx2 = dx * dx;
		const Float dy2 = dy * dy;

		// update old_cells
		memcpy(un, u, nc * sizeof(Float));
		memcpy(vn, v, nc * sizeof(Float));
		memcpy(pn, p, nc * sizeof(Float));

		/////////////////////////////////////////////////////////////////////////////
		// precalculate const values in pressure calculation, store in a buffer
#pragma omp parallel for
		for (int y = 1; y < ny - 1; ++y)
		{
			for (int x = 1; x < nx - 1; x += 4)
			{
				const int idx = y * nx + x;

				if (is_solid[idx] == true)
				{
					continue;
				}

				const QF u_left = *(QF*)&un[idx - 1];
				const QF u_right = *(QF*)&un[idx + 1];
				const QF u_down = *(QF*)&un[idx - nx];
				const QF u_up = *(QF*)&un[idx + nx];

				const QF v_left = *(QF*)&vn[idx - 1];
				const QF v_right = *(QF*)&vn[idx + 1];
				const QF v_down = *(QF*)&vn[idx - nx];
				const QF v_up = *(QF*)&vn[idx + nx];

				*(QF*)&b_[idx] = (density_qf * dx2_qf * dy2_qf) / (mm_set1(kTwoF) * (dx2_qf + dy2_qf))
					* ((mm_set1(kOneF) / dt_qf) * ((u_right - u_left) / (mm_set1(kTwoF) * dx_qf) + (v_up - v_down) / (mm_set1(kTwoF) * dy_qf))
						- ((u_right - u_left) / (mm_set1(kTwoF) * dx_qf)) * ((u_right - u_left) / (mm_set1(kTwoF) * dx_qf))
						- mm_set1(kTwoF) * (((u_up - u_down) / (mm_set1(kTwoF) * dy_qf)) * ((v_right - v_left) / (mm_set1(kTwoF) * dx_qf)))
						- ((v_up - v_down) / (mm_set1(kTwoF) * dy_qf)) * ((v_up - v_down) / (mm_set1(kTwoF) * dy_qf)));
			}
		}

		// iteratively solve for pressure using the poisson equation
		Float l1norm = l1norm_target + kOneF;
		while (l1norm > l1norm_target)
		{
			// update old pressure
			memcpy(pn, p, nbf);

			/////////////////////////////////////////////////////////////////////////////
#pragma omp parallel for
			for (int y = 1; y < ny - 1; ++y)
			{
				for (int x = 1; x < nx - 1; x += 4)
				{
					const int idx = y * nx + x;

					//if (is_solid[idx] == true)
					//	continue;
					if (is_solid[idx] && is_solid[idx + 1] && is_solid[idx + 2] && is_solid[idx + 3])
						continue;

					const QF p_left = *(QF*)&pn[idx - 1];
					const QF p_right = *(QF*)&pn[idx + 1];
					const QF p_down = *(QF*)&pn[idx - nx];
					const QF p_up = *(QF*)&pn[idx + nx];

					*(QF*)&p[idx] = ((p_right + p_left) * dy2_qf + (p_up + p_down) * dx2_qf) / (mm_set1(kTwoF) * (dx2_qf + dy2_qf)) - *(QF*)&b_[idx];
				}
			}

			// reset boundary conditions
			ResetBoundaryConditions();

			// update the l1norm
			l1norm = GetL1Norm();
		}

		/////////////////////////////////////////////////////////////////////////////
		// update all inner cells
#pragma omp parallel for
		for (int y = 1; y < ny - 1; ++y)
		{
			for (int x = 1; x < nx - 1; ++x)
			{
				const int idx = y * nx + x;

				const Float u_left = un[idx - 1];
				const Float u_right = un[idx + 1];
				const Float u_down = un[idx - nx];
				const Float u_up = un[idx + nx];

				const Float v_left = vn[idx - 1];
				const Float v_right = vn[idx + 1];
				const Float v_down = vn[idx - nx];
				const Float v_up = vn[idx + nx];

				const Float p_left = is_solid[idx - 1] ? p[idx] : p[idx - 1];
				const Float p_right = is_solid[idx + 1] ? p[idx] : p[idx + 1];
				const Float p_down = is_solid[idx - nx] ? p[idx] : p[idx - nx];
				const Float p_up = is_solid[idx + nx] ? p[idx] : p[idx + nx];

				// x velocity
				u[idx] = un[idx]
					// convection
					- un[idx] * (dt / dx) * (un[idx] - u_left)
					- vn[idx] * (dt / dy) * (un[idx] - u_down)
					// diffusion
					+ (viscosity_ * dt / dx2) * (u_right - kTwoF * un[idx] + u_left)
					+ (viscosity_ * dt / dy2) * (u_up - kTwoF * un[idx] + u_down)
					// pressure
					- dt / (density_ * kTwoF * dx) * (p_right - p_left);

				// y velocity
				v[idx] = vn[idx]
					// convection
					- un[idx] * (dt / dx) * (vn[idx] - v_left)
					- vn[idx] * (dt / dy) * (vn[idx] - v_down)
					// diffusion
					+ (viscosity_ * dt / dx2) * (v_right - kTwoF * vn[idx] + v_left)
					+ (viscosity_ * dt / dy2) * (v_up - kTwoF * vn[idx] + v_down)
					// pressure
					- dt / (density_ * kTwoF * dy) * (p_up - p_down);
			}
		}
	}
}

void Simulation::Draw()
{
	// plot magnitude of u
	Float min_mag = Vec2(u[0], v[0]).Magnitude();
	Float max_mag = min_mag;

	for (int i = 0; i < nc; ++i)
	{
		min_mag = std::min(min_mag, Vec2(u[i], v[i]).Magnitude());
		max_mag = std::max(max_mag, Vec2(u[i], v[i]).Magnitude());
	}

	for (int y = 0; y < ny; y += 2)
	{
		for (int x = 0; x < nx; x += 2)
		{
			int idx = y * nx + x;
			Float inv_delta = 1 / (max_mag - min_mag);
			Float mag[4] = {
			 Vec2(u[idx], v[idx]).Magnitude(),
			 Vec2(u[idx + 1], v[idx + 1]).Magnitude(),
			 Vec2(u[idx + nx], v[idx + nx]).Magnitude(),
			 Vec2(u[idx + nx + 1], v[idx + nx + 1]).Magnitude()
			};
			Color res =
				(Cell::mc1 * ((max_mag - mag[0]) * inv_delta) + Cell::mc2 * ((mag[0] - min_mag) * inv_delta)) * 0.25f +
				(Cell::mc1 * ((max_mag - mag[1]) * inv_delta) + Cell::mc2 * ((mag[1] - min_mag) * inv_delta)) * 0.25f +
				(Cell::mc1 * ((max_mag - mag[2]) * inv_delta) + Cell::mc2 * ((mag[2] - min_mag) * inv_delta)) * 0.25f +
				(Cell::mc1 * ((max_mag - mag[3]) * inv_delta) + Cell::mc2 * ((mag[3] - min_mag) * inv_delta)) * 0.25f;
			gfx_.PutPixel(x / 2, ny - y / 2 - 1, res); // left top
		}
	}

	// plot pressure
	Float min_p = p[0];
	Float max_p = min_p;
	for (int i = 0; i < nc; ++i)
	{
		min_p = std::min(min_p, p[i]);
		max_p = std::max(max_p, p[i]);
	}

	for (int y = 0; y < ny; y += 2)
	{
		for (int x = 0; x < nx; x += 2)
		{
			int idx = y * nx + x;
			Float inv_delta = 1 / (max_p - min_p);
			Float pressure[4] = {
			 p[idx],
			 p[idx + 1],
			 p[idx + nx],
			 p[idx + nx + 1]
			};
			Color res =
				(Cell::pc1 * (max_p - pressure[0]) * inv_delta) + Cell::pc2 * ((pressure[0] - min_p) * inv_delta) * 0.25f +
				(Cell::pc1 * (max_p - pressure[1]) * inv_delta) + Cell::pc2 * ((pressure[1] - min_p) * inv_delta) * 0.25f +
				(Cell::pc1 * (max_p - pressure[2]) * inv_delta) + Cell::pc2 * ((pressure[2] - min_p) * inv_delta) * 0.25f +
				(Cell::pc1 * (max_p - pressure[3]) * inv_delta) + Cell::pc2 * ((pressure[3] - min_p) * inv_delta) * 0.25f;
			gfx_.PutPixel(nx / 2 + x / 2, ny - y / 2 - 1, res); // right top
		}
	}

	// plot b
	Float min_cor = b_[0];
	Float max_cor = min_cor;

	for (int i = 0; i < nc; ++i)
	{
		min_cor = std::min(min_cor, b_[i]);
		max_cor = std::max(max_cor, b_[i]);
	}

	for (int y = 0; y < ny; y += 2)
	{
		for (int x = 0; x < nx; x += 2)
		{
			int idx = y * nx + x;
			Float inv_delta = 1 / (max_cor - min_cor);
			Float cor[4] = {
			 b_[idx],
			 b_[idx + 1],
			 b_[idx + nx],
			 b_[idx + nx + 1]
			};
			Color res =
				(Cell::bc1 * ((max_cor - cor[0]) * inv_delta) + Cell::bc2 * ((cor[0] - min_cor) * inv_delta)) * 0.25f +
				(Cell::bc1 * ((max_cor - cor[1]) * inv_delta) + Cell::bc2 * ((cor[1] - min_cor) * inv_delta)) * 0.25f +
				(Cell::bc1 * ((max_cor - cor[2]) * inv_delta) + Cell::bc2 * ((cor[2] - min_cor) * inv_delta)) * 0.25f +
				(Cell::bc1 * ((max_cor - cor[3]) * inv_delta) + Cell::bc2 * ((cor[3] - min_cor) * inv_delta)) * 0.25f;
			gfx_.PutPixel(x / 2, ny - ny / 2 - y / 2 - 1, res); // right top
		}
	}

	OutputDebugStringA(("Min p = " + std::to_string(min_p) + ", max p = " + std::to_string(max_p) + "\n").c_str());
	OutputDebugStringA(("Min cor = " + std::to_string(min_cor) + ", max cor = " + std::to_string(max_cor) + "\n").c_str());
	OutputDebugStringA(("Min vel = " + std::to_string(min_mag) + ", max vel = " + std::to_string(max_mag) + "\n").c_str());
	OutputDebugStringA("\n");
}

void Simulation::InitField()
{
	for (int x = 0; x < nx; ++x)
	{
		p[x] = p[x + nx];
		u[x] = kZeroF;
		v[x] = kZeroF;
		is_solid[x] = true;

		p[x + (ny - 1) * nx] = p[x + (ny - 2) * nx];
		u[x + (ny - 1) * nx] = kZeroF;
		v[x + (ny - 1) * nx] = kZeroF;
		is_solid[x + (ny - 1) * nx] = true;
	}
	for (int y = 1; y < ny - 1; ++y)
	{
		p[y * nx] = const_pressure;
		u[y * nx] = u[y * nx + 1];
		v[y * nx] = v[y * nx + 1];
		is_solid[y * nx] = true;
		p[y * nx + nx - 1] = -const_pressure;
		u[y * nx + nx - 1] = u[y * nx + nx - 2];
		v[y * nx + nx - 1] = v[y * nx + nx - 2];
		is_solid[y * nx + nx - 1] = true;
	}

	int
		min_y = int(ny * 0.1f),
		max_y = int(ny * 0.8f),
		min_x = int(nx * 0.1f),
		max_x = int(nx * 0.125f);

	for (int y = min_y; y < max_y; ++y)
	{
		for (int x = min_x; x < max_x; ++x)
		{
			u[y * nx + x] = kZeroF;
			v[y * nx + x] = kZeroF;
			is_solid[y * nx + x] = true;
		}
	}

	ResetBoundaryConditions();
}

void Simulation::ResetBoundaryConditions()
{
	for (int x = 1; x < nx - 1; x += 4)
	{
		*(QF*)&p[x] = *(QF*)&p[x + nx];
		*(QF*)&p[x + (ny - 1) * nx] = *(QF*)&p[x + (ny - 2) * nx];
	}

	for (int y = 1; y < ny - 1; ++y)
	{
		p[y * nx] = const_pressure;
		u[y * nx] = u[y * nx + 1];
		v[y * nx] = v[y * nx + 1];

		u[y * nx + nx - 1] = u[y * nx + nx - 2];
		v[y * nx + nx - 1] = v[y * nx + nx - 2];
		p[y * nx + nx - 1] = -const_pressure;
	}

#pragma omp parallel for
	for (int y = 1; y < ny - 1; ++y)
	{
		for (int x = 1; x < nx - 1; ++x)
		{
			int idx = y * nx + x;

			if (is_solid[idx] == false)
				continue;
			//if (!is_solid[idx] && !is_solid[idx + 1] && !is_solid[idx + 2] && !is_solid[idx + 3])
			//	continue;

			int count = 0;
			Float avg_pressure = kZeroF;

			// left
			count += (int)is_solid[idx - 1];
			avg_pressure += is_solid[idx - 1] ? kZeroF : p[idx - 1];

			// right
			count += (int)is_solid[idx + 1];
			avg_pressure += is_solid[idx + 1] ? kZeroF : p[idx + 1];

			// down
			count += (int)is_solid[idx - nx];
			avg_pressure += is_solid[idx - nx] ? kZeroF : p[idx - nx];

			// up
			count += (int)is_solid[idx + nx];
			avg_pressure += is_solid[idx + nx] ? kZeroF : p[idx + nx];

			if (count == 0)
			{
				continue;
			}

			u[idx] = kZeroF;
			v[idx] = kZeroF;
			p[idx] = avg_pressure / (Float)count;
		}
	}
}

Float Simulation::GetL1Norm()
{
	//*
	Float sum_diff = kZeroF;
	Float sum_old = kZeroF;
	
#pragma omp parallel for default(shared) reduction(+:sum_old, sum_diff)
	for (int i = 0; i < nc; i += 1)
	{
		//if (is_solid[i] == true)
		//	continue;
		//if (is_solid[i] && is_solid[i + 1] && is_solid[i + 2] && is_solid[i + 3])
		//	continue;
	
		const Float old = Abs(pn[i]);
		const Float cur = Abs(p[i]);
		sum_diff += cur - old;
		sum_old += old;
	}
	
	return sum_diff / sum_old;
	/*/
	QF sum_diff = mm_set1(kZeroF);
	QF sum_old = mm_set1(kZeroF);
	
//#pragma omp parallel for default(shared) reduction(+:sum_old, sum_diff)
	for (int i = 0; i < nc; i += 4)
	{
		//if (is_solid[i] == true)
		//	continue;
		//if (is_solid[i] && is_solid[i + 1] && is_solid[i + 2] && is_solid[i + 3])
		//	continue;
	
		const QF old = mm_abs(*(QF*)&pn[i]);
		const QF cur = mm_abs(*(QF*)&p[i]);
		sum_diff = sum_diff + cur - old;
		sum_old = sum_old + old;
	}
	
	return mm_sum_of_elements(sum_diff / sum_old);
	//*/
}
