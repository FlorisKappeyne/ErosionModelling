#include "Simulation.h"
#include <algorithm>
#include <memory>
#include <Windows.h>

/*
todo
 - boundary conditions checken
 - double check all indices
*/


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
	gfx_(gfx)
{
	p = (Float*)_aligned_malloc(nc * sizeof(Float), 64);
	pn = (Float*)_aligned_malloc(nc * sizeof(Float), 64);
	u = (Float*)_aligned_malloc((nx + 1) * ny * sizeof(Float), 64);
	un = (Float*)_aligned_malloc((nx + 1) * ny * sizeof(Float), 64);
	v = (Float*)_aligned_malloc((ny + 1) * nx * sizeof(Float), 64);
	vn = (Float*)_aligned_malloc((ny + 1) * nx * sizeof(Float), 64);
	is_solid = (bool*)_aligned_malloc(nc * sizeof(bool), 64);

	memset(p, 0, nc * sizeof(Float));
	memset(pn, 0, nc * sizeof(Float));
	memset(u, 0, (nx + 1) * ny * sizeof(Float));
	memset(un, 0, (nx + 1) * ny * sizeof(Float));
	memset(v, 0, (ny + 1) * nx * sizeof(Float));
	memset(vn, 0, (ny + 1) * nx * sizeof(Float));
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
	memcpy(un, u, (nx + 1) * ny * sizeof(Float));
	memcpy(vn, v, (ny + 1) * nx * sizeof(Float));
	memcpy(pn, p, nc * sizeof(Float));

	//////////////////////////////////////////////////////////////////////////////
	// Update velocities
	UpdateVelocities();

	//ResetBoundaryConditions();

	//////////////////////////////////////////////////////////////////////////////
	// solve the Poisson Pressure equation using the iterative jacobi method
	SolveForPressure();

	//////////////////////////////////////////////////////////////////////////////
	// subtract the pressure gradient
	SubtractPressureGradient();

	///////////////////////////////////////////////////////////////////////////
	// update delta time for the next frame
	time_passed += dt;
	//Float dt_n = dt;
	//Float max_speed = kZeroF;
	//for (int y = 1; y < ny - 1; ++y)
	//{
	//	for (int x = 1; x < nx - 1; ++x)
	//	{
	//		int idx = y * nx + x;
	//		max_speed = Max(u[idx], max_speed);
	//		max_speed = Max(v[idx], max_speed);
	//	}
	//}
	//dt = dx / max_speed;
	//dt *= 0.1f;
	//dt = Min(dt_n * kTwoF, dt); // make sure the dt doesn't grow too much
	//dt_qf = mm_set(dt);
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

	for (int y = 0; y < ny; ++y)
	{
		for (int x = 0; x < nx; ++x)
		{
			int idx = y * nx + x;
			if (is_solid[idx])
			{
				gfx_.PutPixel(x, ny - y - 1, Colors::Green * 0.3f);
				continue;
			}
			Float inv_delta = 1 / (max_mag - min_mag);
			Float mag = Vec2(u[idx], v[idx]).Magnitude();
			Color res = (Cell::mc1 * ((max_mag - mag) * inv_delta) + Cell::mc2 * ((mag - min_mag) * inv_delta));
			gfx_.PutPixel(x, ny - y - 1, res); // left top
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

	for (int y = 0; y < ny; ++y)
	{
		for (int x = 0; x < nx; ++x)
		{
			int idx = y * nx + x;
			if (is_solid[idx])
			{
				gfx_.PutPixel(x + nx, ny - y - 1, Colors::Gray * 0.8f);
				continue;
			}
			Float inv_delta = 1 / (max_p - min_p);
			Float pressure = p[idx];
			Color res = (Cell::pc1 * (max_p - pressure) * inv_delta) + Cell::pc2 * ((pressure - min_p) * inv_delta);
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
			int idx = y * nx + x;
			int idx_u = idx + y;

			if (*(int*)&is_solid[idx] == 0)
				continue;

			QI mask = mm_set(
				(Int)is_solid[idx + 3],
				(Int)is_solid[idx + 2],
				(Int)is_solid[idx + 1],
				(Int)is_solid[idx]);
			mask = mm_sub(zeros, mask);

			*(QF*)&u[idx_u] = mm_blend(*(QF*)&u[idx_u], *(QF*)&zeros, *(QF*)&mask);
			*(QF*)&u[idx_u + 1] = mm_blend(*(QF*)&u[idx_u + 1], *(QF*)&zeros, *(QF*)&mask);
			*(QF*)&v[idx] = mm_blend(*(QF*)&v[idx], *(QF*)&zeros, *(QF*)&mask);
			*(QF*)&v[idx + nx] = mm_blend(*(QF*)&v[idx + nx], *(QF*)&zeros, *(QF*)&mask);
		}
	}
}

void Simulation::ResetEdges()
{
	for (int x = 1; x < nx - 1; x += 4)
	{
		*(QF*)&p[x] = mm_set(kZeroF);
		*(QF*)&u[x] = mm_set(kOneF);
		*(QF*)&v[x] = mm_set(kZeroF);

		*(QF*)&p[x + (ny - 1) * nx] = *(QF*)&p[x + (ny - 1) * nx];
		*(QF*)&u[x + (ny - 1) * nx] = mm_set(kZeroF);
		*(QF*)&v[x + (ny - 1) * nx] = mm_set(kZeroF);
	}
	// handle speeds at right edge that are missed by simd loop because there's an extra column
	u[nx] = kZeroF;
	u[nx + (ny - 1) * (nx + 1)] = kZeroF;

	// not ny - 1 because there's an extra row
	for (int y = 1; y < ny; ++y)
	{
		p[y * nx] = p[y * nx + nx - 2];
		u[y * nx] = kZeroF;
		v[y * nx] = kZeroF;

		p[y * nx + nx - 1] = p[y * nx + 1];
		u[y * nx + nx - 1] = kZeroF;
		v[y * nx + nx - 1] = kZeroF;
	}
	v[ny * nx - 1] = kZeroF;
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
			const int idx_u = idx + y;

			if (*(int32*)&is_solid[idx] == 0x01010101)
				continue;

			const QF u_center = *(QF*)&un[idx_u];
			const QF u_left = *(QF*)&un[idx_u - 1];
			const QF u_right = *(QF*)&un[idx_u + 1];
			const QF u_down = *(QF*)&un[idx_u - nx - 1];
			const QF u_up = *(QF*)&un[idx_u + nx + 1];
			const QF avg_v = Average4(
				*(QF*)&vn[idx - 1],
				*(QF*)&vn[idx],
				*(QF*)&vn[idx + nx - 1],
				*(QF*)&vn[idx + nx]);

			const QF v_center = *(QF*)&vn[idx];
			const QF v_left = *(QF*)&vn[idx - 1];
			const QF v_right = *(QF*)&vn[idx + 1];
			const QF v_down = *(QF*)&vn[idx - nx];
			const QF v_up = *(QF*)&vn[idx + nx];
			const QF avg_u = Average4(
				*(QF*)&un[idx_u],
				*(QF*)&un[idx_u + 1],
				*(QF*)&un[idx_u - nx - 1],
				*(QF*)&un[idx_u - nx]);

			// precalculate the gradients
			Vec2T<QF> u_gradient = Gradient(u_left, u_right, u_down, u_up);
			Vec2T<QF> v_gradient = Gradient(v_left, v_right, v_down, v_up);

			// precalculate the laplacian
			QF u_laplacian = Laplacian(u_center, u_left, u_right, u_down, u_up);
			QF v_laplacian = Laplacian(v_center, v_left, v_right, v_down, v_up);

			// convection
			QF u_convec_qf =
				u_center * dt_qf * u_gradient.x +
				avg_v * dt_qf * u_gradient.y;
			QF v_convec_qf =
				avg_u * dt_qf * v_gradient.x +
				v_center * dt_qf * v_gradient.y;

			// diffusion
			QF u_diff_qf = viscosity_qf * dt_qf * u_laplacian;
			QF v_diff_qf = viscosity_qf * dt_qf * v_laplacian;

			// acceleration through constant force
			QF u_acc_qf = force_u_qf / (dx_qf * dy_qf * density_qf) * dt_qf;
			QF v_acc_qf = force_v_qf / (dx_qf * dy_qf * density_qf) * dt_qf;

			// update the velocities with the calculated terms
			*(QF*)&u[idx] = *(QF*)&un[idx] - u_convec_qf + u_diff_qf + u_acc_qf;
			*(QF*)&v[idx] = *(QF*)&vn[idx] - v_convec_qf + v_diff_qf + v_acc_qf;
		}
		// for the last column of u values
		const int idx_u = y * (nx + 1) + nx - 1;
		Float u_center = un[idx_u];
		Float u_left = un[idx_u - 1];
		Float u_right = un[idx_u + 1];
		Float u_down = un[idx_u - nx - 1];
		Float u_up = un[idx_u + nx + 1];

		Float avg_v = Average4(
			vn[y * nx - 2],
			vn[y * nx - 1],
			vn[y * nx + nx - 2],
			vn[y * nx + nx - 1]);
		Float u_conv = 
			u_center * dt * (u_right - u_left) / (kTwoF * dx) + 
			avg_v * dt * (u_up - u_down) / (kTwoF * dy);
		Float u_diff = viscosity_ * dt * Laplacian(u_center, u_left, u_right, u_down, u_up);
		Float u_acc = force_u_ / (dx * dy * density_) * dt;
		u[idx_u] = un[idx_u] - u_conv + u_diff + u_acc;
	}

	// for the last row of v values
	for (int x = 1; x < nx - 1; ++x)
	{
		// for the last column of u values
		const int idx = (ny - 2) * nx + x;
		const int idx_u = idx + ny - 2;

		Float v_center = vn[idx];
		Float v_left = vn[idx - 1];
		Float v_right = vn[idx + 1];
		Float v_down = vn[idx - nx];
		Float v_up = vn[idx + nx];

		Float avg_u = Average4(
			vn[idx_u],
			vn[idx_u + 1],
			vn[idx_u - nx - 1],
			vn[idx_u - nx]);

		Float v_conv =
			avg_u * dt * (v_right - v_left) / (kTwoF * dx) +
			v_center * dt * (v_up - v_down) / (kTwoF * dy);
		Float v_diff = viscosity_ * dt * Laplacian(v_center, v_left, v_right, v_down, v_up);
		Float v_acc = force_v_ / (dx * dy * density_) * dt;
		u[idx_u] = un[idx_u] - v_conv + v_diff + v_acc;
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
				const int idx_u = idx + y;

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

				const QF u_left = *(QF*)&u[idx_u];
				const QF u_right = *(QF*)&u[idx_u + 1];

				const QF v_down = *(QF*)&v[idx];
				const QF v_up = *(QF*)&v[idx + nx];

				const QF alpha = dx2_qf * dy2_qf * mm_set(-kOneF);
				const QF beta = mm_set(2.0) * dx2_qf * dy2_qf;
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
#pragma omp parallel for
	for (int y = 1; y < ny - 1; ++y)
	{
		for (int x = 1; x < nx - 1; x += 4)
		{
			const int idx = y * nx + x;
			const int idx_u = idx + y;

			if (*(int*)&is_solid[idx] == 0x01010101)
				continue;

			QF p_centre = *(QF*)&pn[idx];
			QF p_left = *(QF*)&pn[idx - 1];
			QF p_down = *(QF*)&pn[idx - nx];

			QI mask_left = mm_set(
				(Int)is_solid[idx - 1 + 3],
				(Int)is_solid[idx - 1 + 2],
				(Int)is_solid[idx - 1 + 1],
				(Int)is_solid[idx - 1]);
			mask_left = mm_sub(zeros, mask_left);

			QI mask_down = mm_set(
				(Int)is_solid[idx - nx + 3],
				(Int)is_solid[idx - nx + 2],
				(Int)is_solid[idx - nx + 1],
				(Int)is_solid[idx - nx]);
			mask_down = mm_sub(zeros, mask_down);

			p_left = mm_blend(p_left, p_centre, *(QF*)&mask_left);\
			p_down = mm_blend(p_down, p_centre, *(QF*)&mask_down);

			// pressure is always non-staggered
			Vec2T<QF> gradient = GradientStaggered(p_left, p_centre, p_down, p_centre);

			// x velocity
			*(QF*)&u[idx_u] = *(QF*)&u[idx_u] - gradient.x;

			// y velocity
			*(QF*)&v[idx] = *(QF*)&v[idx] - gradient.y;
		}
		// for the last column of u values
		Float p_centre = p[y * nx + nx - 1];
		Float p_left = is_solid[y * nx + nx - 1] ? p[y * nx + nx - 1] : p[y * nx + nx - 2];
		u[y * (nx + 1) + nx - 1] = u[y * (nx + 1) + nx - 1] - (p_centre - p_left) / dx;
	}
	
	// for the last row of v values
	for (int x = 1; x < nx - 1; ++x)
	{
		// for the last column of u values
		Float p_centre = p[(ny - 1) * nx + x];
		Float p_down = p[(ny - 2) * nx + x];
		p_down = is_solid[(ny - 1) * nx + x] ? p_centre : p_down;
		v[(ny - 1) * nx + x] = v[ny * nx + x] - (p_centre - p_down) / dy;
	}
}