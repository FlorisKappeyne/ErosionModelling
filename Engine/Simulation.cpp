#include "Simulation.h"
#include <algorithm>
#include <memory>
#include <Windows.h>


///////////////////////////////////////////////////////////////////////////////
// Simulation
///////////////////////////////////////////////////////////////////////////////

Simulation::Simulation(Graphics& gfx, Float viscosity, Float density,
	Float ds, Float delta_time)
	:
	viscosity_(viscosity),
	density_(density),
	force_u_(0.0f),
	force_v_(0.0f),
	nx(gfx.ScreenHeight),
	ny(gfx.ScreenHeight),
	nc(nx * ny),
	dx(ds),
	dy(ds),
	dx2(dx * dx),
	dy2(dy * dy),
	dt(delta_time),
	time_passed(kZeroF),
	viscosity_qf(mm_set(viscosity)),
	density_qf(mm_set(density)),
	force_u_qf(mm_set(force_u_)),
	force_v_qf(mm_set(force_v_)),
	dx_qf(mm_set(dx)),
	dy_qf(mm_set(dy)),
	dx2_qf(mm_set(dx * dx)),
	dy2_qf(mm_set(dy * dy)),
	dt_qf(mm_set(dt)),
	ones(mm_set(Int(-1))),
	zeros(mm_set(Int(0))),
	kTwoQF(mm_set(kTwoF)),
	kOneQF(mm_set(kOneF)),
	kZeroQF(mm_set(kZeroF)),
	min_mag(kZeroF),
	max_mag(kZeroF),
	min_p(kZeroF),
	max_p(kZeroF),
	drawing_vars_initialized(false),
	gfx_(gfx)
{
	p = (Float*)_aligned_malloc(nc * sizeof(Float), 64);
	pn = (Float*)_aligned_malloc(nc * sizeof(Float), 64);
	u = (Float*)_aligned_malloc((nx - 1) * ny * sizeof(Float), 64);
	un = (Float*)_aligned_malloc((nx - 1) * ny * sizeof(Float), 64);
	v = (Float*)_aligned_malloc((ny - 1) * nx * sizeof(Float), 64);
	vn = (Float*)_aligned_malloc((ny - 1) * nx * sizeof(Float), 64);
	is_solid = (bool*)_aligned_malloc(nc * sizeof(bool), 64);

	memset(p, 0, nc * sizeof(Float));
	memset(pn, 0, nc * sizeof(Float));
	memset(u, 0, (nx - 1) * ny * sizeof(Float));
	memset(un, 0, (nx - 1) * ny * sizeof(Float));
	memset(v, 0, (ny - 1) * nx * sizeof(Float));
	memset(vn, 0, (ny - 1) * nx * sizeof(Float));
	memset(is_solid, 0, nc * sizeof(bool));

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
	_aligned_free(is_solid);
}

void Simulation::Step()
{
	// reset boundary conditions
	ResetBoundaryConditions();

	// update old_cells
	memcpy(un, u, (nx - 1) * ny * sizeof(Float));
	memcpy(vn, v, (ny - 1) * nx * sizeof(Float));
	memcpy(pn, p, nc * sizeof(Float));

	//////////////////////////////////////////////////////////////////////////////
	// Update velocities
	UpdateVelocities();

	// reset boundary conditions
	ResetBoundaryConditions();

	//////////////////////////////////////////////////////////////////////////////
	// solve the Poisson Pressure equation using the iterative jacobi method
	SolveForPressure();

	//////////////////////////////////////////////////////////////////////////////
	// subtract the pressure gradient
	SubtractPressureGradient();

	///////////////////////////////////////////////////////////////////////////
	// update delta time for the next frame
	time_passed += dt;
//	Float dt_n = dt;
//	Float max_speed = kZeroF;
//	for (int y = 0; y < ny - 1; ++y)
//	{
//		for (int x = 0; x < nx - 1; ++x)
//		{
//			int idx = y * nx + x;
//			max_speed = Max(u[IndexU(x, y)], max_speed);
//			max_speed = Max(v[IndexV(x, y)], max_speed);
//		}
//	}
//	dt = dx / max_speed;
//	dt = Min(dt_n * kTwoF, dt); // make sure the dt doesn't grow too much
//	dt_qf = mm_set(dt);
}

void Simulation::Draw()
{
	// plot magnitude of u
	if (drawing_vars_initialized == false)
	{
		min_mag = Vec2(u[0], v[0]).Magnitude();
		max_mag = min_mag;

		min_p = p[0];
		max_p = min_p;

		drawing_vars_initialized = true;
	}

	for (int y = 1; y < ny - 1; ++y)
	{
		for (int x = 1; x < nx - 1; ++x)
		{
			Float avg_u = (u[IndexU(x - 1, y)] + u[IndexU(x, y)]) / kTwoF;
			Float avg_v = (v[IndexV(x, y - 1)] + v[IndexV(x, y)]) / kTwoF;
			Float mag = Vec2(avg_u, avg_v).Magnitude();

			min_mag = std::min(min_mag, mag);
			max_mag = std::max(max_mag, mag);
		}
	}

	Float inv_delta_mag = 1 / (max_mag - min_mag);
	for (int y = 1; y < ny - 1; ++y)
	{
		for (int x = 1; x < nx - 1; ++x)
		{
			int idx = y * nx + x;
			if (is_solid[idx])
			{
				gfx_.PutPixel(x, ny - y - 1, Colors::Green * 0.3f);
				continue;
			}

			Float mag = Vec2(
				(u[IndexU(x - 1, y)] + u[IndexU(x, y)]) / kTwoF,
				(v[IndexV(x, y - 1)] + v[IndexV(x, y)]) / kTwoF
			).Magnitude();

			Color res =
				(Cell::mc1 * (max_mag - mag) * inv_delta_mag) +
					Cell::mc2 * ((mag - min_mag) * inv_delta_mag);
			gfx_.PutPixel(x, ny - y - 1, res); // left top
		}
	}

	// plot pressure
	for (int y = 1; y < ny - 1; ++y)
	{
		for (int x = 1; x < nx - 1; ++x)
		{
			Float pressure = p[IndexP(x, y)];
			min_p = std::min(min_p, pressure);
			max_p = std::max(max_p, pressure);
		}
	}
	Float inv_delta_p = 1 / (max_p - min_p);

	for (int y = 0; y < ny; ++y)
	{
		for (int x = 0; x < nx; ++x)
		{
			const int idx = IndexP(x, y);
			if (is_solid[idx])
			{
				gfx_.PutPixel(x + nx, ny - y - 1, Colors::Gray * 0.8f);
				continue;
			}
			Float pressure = p[idx];
			Color res =
				(Cell::pc1 * (max_p - pressure) * inv_delta_p) +
				Cell::pc2 * ((pressure - min_p) * inv_delta_p);
			gfx_.PutPixel(x + nx, ny - y - 1, res); // right top
		}
	}

	OutputDebugStringA(("Min p = " + std::to_string(min_p) + ", max p = " + std::to_string(max_p) + "\n").c_str());
	OutputDebugStringA(("Min vel = " + std::to_string(min_mag) + ", max vel = " + std::to_string(max_mag) + "\n").c_str());
	OutputDebugStringA(("time passed = " + std::to_string(time_passed) + ", dt = " + std::to_string(dt) + "\n").c_str());
	OutputDebugStringA("\n");
}


///////////////////////////////////////////////////////////////////////////////
// Boundary control
///////////////////////////////////////////////////////////////////////////////

void Simulation::InitField()
{
	//for (int y = 1; y < ny - 1; ++y)
	//{
	//	Float height_fraction = Float(y) / Float(ny);
	//	for (int x = 1; x < nx - 1; ++x)
	//	{
	//		// control variables
	//		const Float periods = kOneF;
	//		const Float width = 0.2f;
	//		const Float offset = 0.1f;
	//		const int num_sines = 15;
	//
	//		int idx = y * nx + x;
	//		is_solid[idx] = true;
	//
	//		for (int i = 0; i < num_sines; ++i)
	//		{
	//			Float fraction = Float(x + i - num_sines / 2) / nx;
	//			Float height = Sin(fraction * kTwoPi * periods);
	//			height *= height;
	//			height *= kOneF - kTwoF * width;
	//			height += width;
	//
	//			if (int((height + offset) * ny) >= y && int((height - offset) * ny) <= y)
	//			{
	//				is_solid[idx] = false;
	//				u[idx] = kZeroF;
	//				v[idx] = kZeroF;
	//			}
	//		}
	//	}
	//}

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
			const int idx = IndexP(x, y);
			const int idx_u = IndexU(x, y);
			const int idx_v = IndexV(x, y);

			if (*(int*)&is_solid[idx] == 0)
				continue;

			QI mask = mm_set(
				(Int)is_solid[idx + 3],
				(Int)is_solid[idx + 2],
				(Int)is_solid[idx + 1],
				(Int)is_solid[idx]);
			mask = mm_sub(zeros, mask);

			*(QF*)&u[idx_u - 1] = mm_blend(*(QF*)&u[idx_u - 1], *(QF*)&zeros, *(QF*)&mask);
			*(QF*)&u[idx_u] = mm_blend(*(QF*)&u[idx_u], *(QF*)&zeros, *(QF*)&mask);
			*(QF*)&v[idx_v - nx] = mm_blend(*(QF*)&v[idx - nx], *(QF*)&zeros, *(QF*)&mask);
			*(QF*)&v[idx_v] = mm_blend(*(QF*)&v[idx], *(QF*)&zeros, *(QF*)&mask);
		}
	}
}

void Simulation::ResetEdges()
{
	for (int x = 1; x < nx - 1; x += 4)
	{
		*(QF*)&p[x] = kZeroQF;
		*(QF*)&u[x] = kOneQF;
		*(QF*)&v[x] = kZeroQF;

		*(QF*)&p[IndexP(x, ny - 1)] = *(QF*)&p[IndexP(x, ny - 2)];
		*(QF*)&u[IndexU(x, ny - 1)] = kZeroQF;
		*(QF*)&v[IndexV(x, ny - 2)] = kZeroQF;
	}

	for (int y = 1; y < ny - 1; ++y)
	{
		p[IndexP(0, y)] = p[IndexP(1, y)];
		u[IndexU(0, y)] = kZeroF;
		v[IndexV(0, y)] = kZeroF;

		p[IndexP(nx - 1, y)] = p[IndexP(nx - 2, y)];
		u[IndexU(nx - 2, y)] = kZeroF;
		v[IndexV(nx - 1, y)] = kZeroF;
	}
}


///////////////////////////////////////////////////////////////////////////////
// Solving Navier-Stokes equation
///////////////////////////////////////////////////////////////////////////////

void Simulation::UpdateVelocities()
{
//#pragma omp parallel for
	for (int y = 1; y < ny - 1; ++y)
	{
		for (int x = 1; x < nx - 1; x += 4)
		{
			const int idx = IndexP(x, y);
			const int idx_u = IndexU(x, y);
			const int idx_v = IndexV(x, y);

			if (*(int32*)&is_solid[idx] == 0x01010101)
				continue;

			const QF u_center = *(QF*)&un[idx_u];
			const QF u_left = *(QF*)&un[idx_u - 1];
			const QF u_right = *(QF*)&un[idx_u + 1];
			const QF u_down = *(QF*)&un[idx_u - nx + 1];
			const QF u_up = *(QF*)&un[idx_u + nx - 1];
			const QF avg_v = Average4(
				*(QF*)&vn[idx_v - 1],
				*(QF*)&vn[idx_v],
				*(QF*)&vn[idx_v + nx - 1],
				*(QF*)&vn[idx_v + nx]);

			const QF v_center = *(QF*)&vn[idx_v];
			const QF v_left = *(QF*)&vn[idx_v - 1];
			const QF v_right = *(QF*)&vn[idx_v + 1];
			const QF v_down = *(QF*)&vn[idx_v - nx];
			const QF v_up = *(QF*)&vn[idx_v + nx];
			const QF avg_u = Average4(
				*(QF*)&un[idx_u - 1],
				*(QF*)&un[idx_u],
				*(QF*)&un[idx_u + nx - 1],
				*(QF*)&un[idx_u + nx]);

			// precalculate the gradients
			Vec2T<QF> u_gradient = Gradient(u_left, u_right, u_down, u_up);
			Vec2T<QF> v_gradient = Gradient(v_left, v_right, v_down, v_up);

			// precalculate the laplacian
			QF u_laplacian = Laplacian(u_center, u_left, u_right, u_down, u_up);
			QF v_laplacian = Laplacian(v_center, v_left, v_right, v_down, v_up);

			// advection
			QF u_advec_qf =
				u_center * dt_qf * u_gradient.x +
				avg_v * dt_qf * u_gradient.y;
			QF v_advec_qf =
				avg_u * dt_qf * v_gradient.x +
				v_center * dt_qf * v_gradient.y;

			// diffusion
			QF u_diff_qf = viscosity_qf * dt_qf * u_laplacian;
			QF v_diff_qf = viscosity_qf * dt_qf * v_laplacian;

			// acceleration through constant force
			QF u_acc_qf = force_u_qf / (dx_qf * dy_qf * density_qf) * dt_qf;
			QF v_acc_qf = force_v_qf / (dx_qf * dy_qf * density_qf) * dt_qf;

			// update the velocities with the calculated terms
			*(QF*)&u[idx_u] = *(QF*)&un[idx_u] - u_advec_qf + u_diff_qf + u_acc_qf;
			*(QF*)&v[idx_v] = *(QF*)&vn[idx_v] - v_advec_qf + v_diff_qf + v_acc_qf;
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
				const int idx = IndexP(x, y);
				const int idx_u = IndexU(x, y);
				const int idx_v = IndexV(x, y);

				if (*(int*)&is_solid[idx] == 0x01010101)
					continue;

				QF p_center = *(QF*)&pn[idx];
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

				p_left = mm_blend(p_left, p_center, *(QF*)&mask_left);
				p_right = mm_blend(p_right, p_center, *(QF*)&mask_right);
				p_down = mm_blend(p_down, p_center, *(QF*)&mask_down);
				p_up = mm_blend(p_up, p_center, *(QF*)&mask_up);

				const QF u_left = *(QF*)&u[idx_u - 1];
				const QF u_right = *(QF*)&u[idx_u];

				const QF v_down = *(QF*)&v[idx_v - nx];
				const QF v_up = *(QF*)&v[idx_v];

				const QF alpha = dx2_qf * dy2_qf * mm_set(-kOneF);
				const QF beta = kTwoQF * (dx2_qf + dy2_qf);
				const QF div_w = Divergence(u_left, u_right, v_down, v_up);

				// jacobi iteration
				*(QF*)&p[idx] = ((p_left + p_right) * dy2_qf + (p_down + p_up) * dx2_qf + alpha * div_w) / beta;
			}
		}

		// update old pressure
		memcpy(pn, p, nc * sizeof(Float));
	}
}

void Simulation::SubtractPressureGradient()
{
//#pragma omp parallel for
	for (int y = 1; y < ny - 1; ++y)
	{
		for (int x = 1; x < nx - 1; x += 4)
		{
			const int idx = IndexP(x, y);
			const int idx_u = IndexU(x, y);
			const int idx_v = IndexV(x, y);

			if (*(int*)&is_solid[idx] == 0x01010101)
				continue;

			QF p_center = *(QF*)&pn[idx];
			QF p_right = *(QF*)&pn[idx + 1];
			QF p_up = *(QF*)&pn[idx + nx];

			QI mask_right = mm_set(
				(Int)is_solid[idx + 1 + 3],
				(Int)is_solid[idx + 1 + 2],
				(Int)is_solid[idx + 1 + 1],
				(Int)is_solid[idx + 1]);
			mask_right = mm_sub(zeros, mask_right);

			QI mask_up = mm_set(
				(Int)is_solid[idx + nx + 3],
				(Int)is_solid[idx + nx + 2],
				(Int)is_solid[idx + nx + 1],
				(Int)is_solid[idx + nx]);
			mask_up = mm_sub(zeros, mask_up);

			p_right = mm_blend(p_right, p_center, *(QF*)&mask_right);
			p_up = mm_blend(p_up, p_center, *(QF*)&mask_up);

			// pressure is always non-staggered
			Vec2T<QF> gradient = GradientStaggered(p_center, p_right, p_center, p_up);

			// x velocity
			*(QF*)&u[idx_u] = *(QF*)&u[idx_u] - gradient.x;

			// y velocity
			*(QF*)&v[idx_v] = *(QF*)&v[idx_v] - gradient.y;
		}
	}
}