#include "Simulation.h"
#include <algorithm>
#include <memory>
#include <Windows.h>


///////////////////////////////////////////////////////////////////////////////
// Simulation
///////////////////////////////////////////////////////////////////////////////

Simulation::Simulation(Graphics& gfx, Float viscosity, Float density,
	const Vec2I& dim, Float delta_time)
	:
	viscosity_(viscosity),
	density_(density),
	force_u_(2.0f),
	force_v_(0.0f),
	dim_(dim),
	nx(gfx.ScreenWidth),
	ny(gfx.ScreenHeight),
	nc(nx * ny),
	nbf(nc * sizeof(Float)),
	nbb(nc * sizeof(bool)),
	dx(dim_.x / (Float)nx),
	dy(dim_.y / (Float)ny),
	dt(delta_time),
	viscosity_qf(mm_set(viscosity)),
	density_qf(mm_set(density)),
	force_u_qf(mm_set(force_u_)),
	force_v_qf(mm_set(force_v_)),
	dx_qf(mm_set(dx, dx, dx, dx)),
	dy_qf(mm_set(dy, dy, dy, dy)),
	dx2_qf(mm_set(dx * dx, dx * dx, dx * dx, dx * dx)),
	dy2_qf(mm_set(dy * dy, dy * dy, dy * dy, dy * dy)),
	dt_qf(mm_set(dt, dt, dt, dt)),
	ones(mm_set(Int(-1))),
	zeros(mm_set(Int(0))),
	kTwoQF(mm_set(kTwoF)),
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
	// reset boundary conditions
	ResetBoundaryConditions();

	// precalculations
	const Float dx2 = dx * dx;
	const Float dy2 = dy * dy;

	// update old_cells
	memcpy(un, u, nc * sizeof(Float));
	memcpy(vn, v, nc * sizeof(Float));
	memcpy(pn, p, nc * sizeof(Float));

	//////////////////////////////////////////////////////////////////////////////
	// Update velocities
	UpdateVelocities();

	//////////////////////////////////////////////////////////////////////////////
	// solve the Poisson Pressure equation using the iterative jacobi method
	SolveForPressure();

	//////////////////////////////////////////////////////////////////////////////
	// subtract the pressure gradient
	SubtractPressureGradient();
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
			if (is_solid[idx])
			{
				gfx_.PutPixel(x / 2, ny - y / 2 - 1, Colors::Green * 0.5f);
				continue;
			}
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
			if (is_solid[idx])
			{
				gfx_.PutPixel(nx / 2 + x / 2, ny - y / 2 - 1, Colors::Gray);
				continue;
			}
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


///////////////////////////////////////////////////////////////////////////////
// Boundary control
///////////////////////////////////////////////////////////////////////////////

void Simulation::InitField()
{
	for (int y = 1; y < ny - 1; ++y)
	{
		Float height_fraction = Float(y) / Float(ny);
		for (int x = 1; x < nx - 1; ++x)
		{
			int idx = y * nx + x;

			is_solid[idx] = true;


			// control variables
			const Float periods = kOneF;
			const Float width = 0.2f;
			const Float offset = 0.1f;
			const int num_sines = 15;

			for (int i = 0; i < num_sines; ++i)
			{
				Float fraction = Float(x + i - num_sines / 2) / nx;
				Float height = Sin(fraction * kTwoPi * periods);
				height *= height;
				height *= kOneF - kTwoF * width;
				height += width;

				if (int((height + offset) * ny) >= y && int((height - offset) * ny) <= y)
				{
					is_solid[idx] = false;
					u[idx] = kZeroF;
					v[idx] = kZeroF;
				}
			}
		}
	}

	ResetBoundaryConditions();
}

void Simulation::ResetBoundaryConditions()
{
	// reset the edges
	ResetEdges();

#pragma omp parallel for
	for (int y = 1; y < ny - 1; ++y)
	{
		for (int x = 1; x < nx - 1; ++x)
		{
			int idx = y * nx + x;

			if (*(int*)&is_solid[idx] == 0)
				continue;

			QI mask = mm_set(
				(Int)is_solid[idx + 3],
				(Int)is_solid[idx + 2],
				(Int)is_solid[idx + 1],
				(Int)is_solid[idx]);
			mask = mm_sub(zeros, mask);

			*(QF*)&u[idx] = mm_blend(*(QF*)&u[idx], *(QF*)&zeros, *(QF*)&mask);
			*(QF*)&v[idx] = mm_blend(*(QF*)&v[idx], *(QF*)&zeros, *(QF*)&mask);
		}
	}
}

void Simulation::ResetEdges()
{
	for (int x = 1; x < nx - 1; x += 4)
	{
		*(QF*)&p[x] = *(QF*)&p[x + nx];
		*(QF*)&u[x] = mm_set(kZeroF);
		*(QF*)&v[x] = mm_set(kZeroF);

		*(QF*)&p[x + (ny - 1) * nx] = *(QF*)&p[x + (ny - 2) * nx];
		*(QF*)&u[x + (ny - 1) * nx] = mm_set(kZeroF);
		*(QF*)&v[x + (ny - 1) * nx] = mm_set(kZeroF);
	}

	for (int y = 1; y < ny - 1; ++y)
	{
		p[y * nx] = p[y * nx + nx - 2];
		u[y * nx] = u[y * nx + nx - 2];
		v[y * nx] = v[y * nx + nx - 2];

		p[y * nx + nx - 1] = p[y * nx];
		u[y * nx + nx - 1] = u[y * nx];
		v[y * nx + nx - 1] = v[y * nx];
	}
}


///////////////////////////////////////////////////////////////////////////////
// Solving Navier-Stokes equation
///////////////////////////////////////////////////////////////////////////////

void Simulation::UpdateVelocities()
{
#pragma omp parallel for
	for (int y = 1; y < ny - 1; ++y)
	{
		for (int x = 1; x < nx - 1; x += 4)
		{
			const int idx = y * nx + x;

			if (*(int*)&is_solid[idx] == 0x01010101)
				continue;

			const QF u_left = *(QF*)&un[idx - 1];
			const QF u_right = *(QF*)&un[idx + 1];
			const QF u_down = *(QF*)&un[idx - nx];
			const QF u_up = *(QF*)&un[idx + nx];

			const QF v_left = *(QF*)&vn[idx - 1];
			const QF v_right = *(QF*)&vn[idx + 1];
			const QF v_down = *(QF*)&vn[idx - nx];
			const QF v_up = *(QF*)&vn[idx + nx];

			QF p_centre = *(QF*)&pn[idx];
			QF p_left = *(QF*)&pn[idx - 1];
			QF p_right = *(QF*)&pn[idx + 1];
			QF p_down = *(QF*)&pn[idx - nx];
			QF p_up = *(QF*)&pn[idx + nx];

			QI mask_left = mm_set(
				(Int)is_solid[idx - 1 + 3],
				(Int)is_solid[idx - 1 + 2],
				(Int)is_solid[idx - 1 + 1],
				(Int)is_solid[idx - 1]);
			mask_left = mm_sub(zeros, mask_left);

			QI mask_right = mm_set(
				(Int)is_solid[idx + 1 + 3],
				(Int)is_solid[idx + 1 + 2],
				(Int)is_solid[idx + 1 + 1],
				(Int)is_solid[idx + 1]);
			mask_right = mm_sub(zeros, mask_right);

			QI mask_down = mm_set(
				(Int)is_solid[idx - nx + 3],
				(Int)is_solid[idx - nx + 2],
				(Int)is_solid[idx - nx + 1],
				(Int)is_solid[idx - nx]);
			mask_down = mm_sub(zeros, mask_down);

			QI mask_up = mm_set(
				(Int)is_solid[idx + nx + 3],
				(Int)is_solid[idx + nx + 2],
				(Int)is_solid[idx + nx + 1],
				(Int)is_solid[idx + nx]);
			mask_up = mm_sub(zeros, mask_up);

			p_left = mm_blend(p_left, p_centre, *(QF*)&mask_left);
			p_right = mm_blend(p_right, p_centre, *(QF*)&mask_right);
			p_down = mm_blend(p_down, p_centre, *(QF*)&mask_down);
			p_up = mm_blend(p_up, p_centre, *(QF*)&mask_up);

			// x velocity
			*(QF*)&u[idx] = *(QF*)&un[idx]
				// convection
				- *(QF*)&un[idx] * (dt_qf / dx_qf) * (*(QF*)&un[idx] - u_left)
				- *(QF*)&vn[idx] * (dt_qf / dy_qf) * (*(QF*)&un[idx] - u_down)
				// diffusion
				+ (viscosity_qf * dt_qf / dx2_qf) * (u_right - kTwoQF * *(QF*)&un[idx] + u_right)
				+ (viscosity_qf * dt_qf / dy2_qf) * (u_up - kTwoQF * *(QF*)&un[idx] + u_up)
				// add forces
				+ force_u_qf / (dx2_qf * density_qf) * dt_qf;
				// pressure
				//- dt_qf / (density_qf * kTwoQF * dx_qf) * (p_right - p_left);

			// y velocity
			*(QF*)&v[idx] = *(QF*)&vn[idx]
				// convection
				- *(QF*)&un[idx] * (dt_qf / dx_qf) * (*(QF*)&vn[idx] - v_left)
				- *(QF*)&vn[idx] * (dt_qf / dy_qf) * (*(QF*)&vn[idx] - v_down)
				// diffusion
				+ (viscosity_qf * dt_qf / dx2_qf) * (v_right - kTwoQF * *(QF*)&vn[idx] + v_right)
				+ (viscosity_qf * dt_qf / dy2_qf) * (v_up - kTwoQF * *(QF*)&vn[idx] + v_up)
				// add forces
				+ force_v_qf / (dx2_qf * density_qf) * dt_qf;
				// pressure
				//- dt_qf / (density_qf * kTwoQF * dy_qf) * (p_up - p_down);
		}
	}
}

void Simulation::SolveForPressure()
{
	// iteratively solve for pressure using the poisson equation
	for (int i = 0; i < niter_jacobi; ++i)
	{
		/////////////////////////////////////////////////////////////////////////////
#pragma omp parallel for
		for (int y = 1; y < ny - 1; ++y)
		{
			for (int x = 1; x < nx - 1; x += 4)
			{
				const int idx = y * nx + x;

				if (*(int*)&is_solid[idx] == 0x01010101)
					continue;

				QF p_centre = *(QF*)&pn[idx];
				QF p_left = *(QF*)&pn[idx - 1];
				QF p_right = *(QF*)&pn[idx + 1];
				QF p_down = *(QF*)&pn[idx - nx];
				QF p_up = *(QF*)&pn[idx + nx];

				QI mask_left = mm_set(
					(Int)is_solid[idx - 1 + 3],
					(Int)is_solid[idx - 1 + 2],
					(Int)is_solid[idx - 1 + 1],
					(Int)is_solid[idx - 1]);
				mask_left = mm_sub(zeros, mask_left);

				QI mask_right = mm_set(
					(Int)is_solid[idx + 1 + 3],
					(Int)is_solid[idx + 1 + 2],
					(Int)is_solid[idx + 1 + 1],
					(Int)is_solid[idx + 1]);
				mask_right = mm_sub(zeros, mask_right);

				QI mask_down = mm_set(
					(Int)is_solid[idx - nx + 3],
					(Int)is_solid[idx - nx + 2],
					(Int)is_solid[idx - nx + 1],
					(Int)is_solid[idx - nx]);
				mask_down = mm_sub(zeros, mask_down);

				QI mask_up = mm_set(
					(Int)is_solid[idx + nx + 3],
					(Int)is_solid[idx + nx + 2],
					(Int)is_solid[idx + nx + 1],
					(Int)is_solid[idx + nx]);
				mask_up = mm_sub(zeros, mask_up);

				p_left = mm_blend(p_left, p_centre, *(QF*)&mask_left);
				p_right = mm_blend(p_right, p_centre, *(QF*)&mask_right);
				p_down = mm_blend(p_down, p_centre, *(QF*)&mask_down);
				p_up = mm_blend(p_up, p_centre, *(QF*)&mask_up);

				const QF u_left = *(QF*)&un[idx - 1];
				const QF u_right = *(QF*)&un[idx + 1];

				const QF v_down = *(QF*)&vn[idx - nx];
				const QF v_up = *(QF*)&vn[idx + nx];

				const QF alpha = dx2_qf * mm_set(-kOneF);
				const QF beta = mm_set(4.0f);
				const QF div_w = (u_right - u_left) / (kTwoQF * dx_qf) + (v_up - v_down) / (kTwoQF * dy_qf);

				// jacobi iteration
				*(QF*)&p[idx] = (p_left + p_right + p_down + p_up + alpha * div_w) / beta;
				//((p_right + p_left) * dy2_qf + (p_up + p_down) * dx2_qf) / (kTwoQF * (dx2_qf + dy2_qf)) - *(QF*)&b_[idx];
			}
		}

		// update old pressure
		memcpy(pn, p, nbf);
	}
}

void Simulation::SubtractPressureGradient()
{
#pragma omp parallel for
	for (int y = 1; y < ny - 1; ++y)
	{
		for (int x = 1; x < nx - 1; x += 4)
		{
			const int idx = y * nx + x;

			if (*(int*)&is_solid[idx] == 0x01010101)
				continue;

			QF p_centre = *(QF*)&pn[idx];
			QF p_left = *(QF*)&pn[idx - 1];
			QF p_right = *(QF*)&pn[idx + 1];
			QF p_down = *(QF*)&pn[idx - nx];
			QF p_up = *(QF*)&pn[idx + nx];

			QI mask_left = mm_set(
				(Int)is_solid[idx - 1 + 3],
				(Int)is_solid[idx - 1 + 2],
				(Int)is_solid[idx - 1 + 1],
				(Int)is_solid[idx - 1]);
			mask_left = mm_sub(zeros, mask_left);

			QI mask_right = mm_set(
				(Int)is_solid[idx + 1 + 3],
				(Int)is_solid[idx + 1 + 2],
				(Int)is_solid[idx + 1 + 1],
				(Int)is_solid[idx + 1]);
			mask_right = mm_sub(zeros, mask_right);

			QI mask_down = mm_set(
				(Int)is_solid[idx - nx + 3],
				(Int)is_solid[idx - nx + 2],
				(Int)is_solid[idx - nx + 1],
				(Int)is_solid[idx - nx]);
			mask_down = mm_sub(zeros, mask_down);

			QI mask_up = mm_set(
				(Int)is_solid[idx + nx + 3],
				(Int)is_solid[idx + nx + 2],
				(Int)is_solid[idx + nx + 1],
				(Int)is_solid[idx + nx]);
			mask_up = mm_sub(zeros, mask_up);

			p_left = mm_blend(p_left, p_centre, *(QF*)&mask_left);
			p_right = mm_blend(p_right, p_centre, *(QF*)&mask_right);
			p_down = mm_blend(p_down, p_centre, *(QF*)&mask_down);
			p_up = mm_blend(p_up, p_centre, *(QF*)&mask_up);

			// x velocity
			*(QF*)&u[idx] = *(QF*)&u[idx] - (p_right - p_left) / (kTwoQF * dx_qf);

		// y velocity
			*(QF*)&v[idx] = *(QF*)&v[idx] - (p_up - p_down) / (kTwoQF * dy_qf);
		}
	}
}