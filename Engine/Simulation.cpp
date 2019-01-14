#include "Simulation.h"
#include <algorithm>
#include <memory>
#include <string>
#include <Windows.h>

Simulation::Simulation(Graphics& gfx, Float viscosity, Float density,
	const Vec2I& dim, Float dt)
	:
	viscosity_(viscosity),
	density_(density),
	dim_(dim),
	nx(gfx.ScreenWidth),
	ny(gfx.ScreenHeight),
	nc(nx * ny),
	dx(dim_.x / (Float)nx),
	dy(dim_.y / (Float)ny),
	dt_(dt),
	gfx_(gfx)
{
	cur_cells_ = new Cell[nc];
	old_cells_ = new Cell[nc];
	b_ = new Float[nc];

	memset(b_, 0, nc * sizeof(Float));

	// initialize the field
	InitField();
}

Simulation::~Simulation()
{
	delete[] cur_cells_;
	delete[] old_cells_;
	delete[] b_;
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
#pragma omp parallel for
		for (int i = 0; i < nx * ny; ++i)
		{
			old_cells_[i] = cur_cells_[i];
		}

		/////////////////////////////////////////////////////////////////////////////
		// precalculate const values in pressure calculation, store in a buffer
#pragma omp parallel for
		for (int y = 1; y < ny - 1; ++y)
		{
			for (int x = 1; x < nx - 1; ++x)
			{
				const int idx = y * nx + x;

				if (cur_cells_[idx].is_solid_ == true)
				{
					continue;
				}

				const Float u_left = old_cells_[idx - 1].u_.x;
				const Float u_right = old_cells_[idx + 1].u_.x;
				const Float u_down = old_cells_[idx - nx].u_.x;
				const Float u_up = old_cells_[idx + nx].u_.x;

				const Float v_left = old_cells_[idx - 1].u_.y;
				const Float v_right = old_cells_[idx + 1].u_.y;
				const Float v_down = old_cells_[idx - nx].u_.y;
				const Float v_up = old_cells_[idx + nx].u_.y;

				b_[idx] = (density_ * dx2 * dy2) / (kTwoF * (dx2 + dy2))
					* ((kOneF / dt_) * ((u_right - u_left) / (kTwoF * dx) + (v_up - v_down) / (kTwoF * dy))
						- ((u_right - u_left) / (kTwoF * dx)) * ((u_right - u_left) / (kTwoF * dx))
						- kTwoF * (((u_up - u_down) / (kTwoF * dy)) * ((v_right - v_left) / (kTwoF * dx)))
						- ((v_up - v_down) / (kTwoF * dy)) * ((v_up - v_down) / (kTwoF * dy)));
			}
		}

		// iteratively solve for pressure using the poisson equation
		Float l1norm = l1norm_target + kOneF;
		while (l1norm > l1norm_target)
		{
			// update old_cells_
#pragma omp parallel for
			for (int i = 0; i < nc; ++i)
			{
				old_cells_[i].p_ = cur_cells_[i].p_;
			}

			/////////////////////////////////////////////////////////////////////////////
#pragma omp parallel for
			for (int y = 1; y < ny - 1; ++y)
			{
				for (int x = 1; x < nx - 1; ++x)
				{
					const int idx = y * nx + x;

					if (cur_cells_[idx].is_solid_ == true)
					{
						continue;
					}

					const Float p_left = old_cells_[idx - 1].p_;
					const Float p_right = old_cells_[idx + 1].p_;
					const Float p_down = old_cells_[idx - nx].p_;
					const Float p_up = old_cells_[idx + nx].p_;

					cur_cells_[idx].p_ = ((p_right + p_left) * dy2 + (p_up + p_down) * dx2) / (kTwoF * (dx2 + dy2)) -b_[idx];
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

				const Float u = old_cells_[idx].u_.x;
				const Float u_left = old_cells_[idx - 1].u_.x;
				const Float u_right = old_cells_[idx + 1].u_.x;
				const Float u_down = old_cells_[idx - nx].u_.x;
				const Float u_up = old_cells_[idx + nx].u_.x;

				const Float v = old_cells_[idx].u_.y;
				const Float v_left = old_cells_[idx - 1].u_.y;
				const Float v_right = old_cells_[idx + 1].u_.y;
				const Float v_down = old_cells_[idx - nx].u_.y;
				const Float v_up = old_cells_[idx + nx].u_.y;

				const Float p_left = cur_cells_[idx - 1].is_solid_ ? cur_cells_[idx].p_ : cur_cells_[idx - 1].p_;
				const Float p_right = cur_cells_[idx + 1].is_solid_? cur_cells_[idx].p_ : cur_cells_[idx + 1].p_;
				const Float p_down = cur_cells_[idx - nx].is_solid_ ? cur_cells_[idx].p_ : cur_cells_[idx - nx].p_;
				const Float p_up = cur_cells_[idx + nx].is_solid_ ? cur_cells_[idx].p_ : cur_cells_[idx + nx].p_;

				// x velocity
				cur_cells_[idx].u_.x = u
					// convection
					- u * (dt_ / dx) * (u - u_left)
					- v * (dt_ / dy) * (u - u_down)
					// diffusion
					+ (viscosity_ * dt_ / dx2) * (u_right - kTwoF * u + u_left)
					+ (viscosity_ * dt_ / dy2) * (u_up - kTwoF * u + u_down)
					// pressure
					- dt_ / (density_ * kTwoF * dx) * (p_right - p_left);

				// y velocity
				cur_cells_[idx].u_.y = v
					// convection
					- u * (dt_ / dx) * (v - v_left)
					- v * (dt_ / dy) * (v - v_down)
					// diffusion
					+ (viscosity_ * dt_ / dx2) * (v_right - kTwoF * v + v_left)
					+ (viscosity_ * dt_ / dy2) * (v_up - kTwoF * v + v_down)
					// pressure
					- dt_ / (density_ * kTwoF * dy) * (p_up - p_down);
			}
		}
	}
}

void Simulation::Draw()
{
	// plot magnitude of u
	Float min_mag = cur_cells_[0].u_.Magnitude();
	Float max_mag = cur_cells_[0].u_.Magnitude();
	for (int i = 0; i < nc; ++i)
	{
			min_mag = std::min(min_mag, cur_cells_[i].u_.Magnitude());
			max_mag = std::max(max_mag, cur_cells_[i].u_.Magnitude());
	}

	for (int y = 0; y < ny; y += 2)
	{
		for (int x = 0; x < nx; x += 2)
		{
			int idx = y * nx + x;
			Float inv_delta = 1 / (max_mag - min_mag);
			Float mag[4] = {
			 cur_cells_[idx].u_.Magnitude(),
			 cur_cells_[idx + 1].u_.Magnitude(),
			 cur_cells_[idx + nx].u_.Magnitude(),
			 cur_cells_[idx + nx + 1].u_.Magnitude()
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
	Float min_p = cur_cells_[0].p_;
	Float max_p = cur_cells_[0].p_;
	for (int i = 0; i < nc; ++i)
	{
			min_p = std::min(min_p, cur_cells_[i].p_);
			max_p = std::max(max_p, cur_cells_[i].p_);
	}

	for (int y = 0; y < ny; y += 2)
	{
		for (int x = 0; x < nx; x += 2)
		{
			int idx = y * nx + x;
			Float inv_delta = 1 / (max_p - min_p);
			Float p[4] = {
			 cur_cells_[idx].p_,
			 cur_cells_[idx + 1].p_,
			 cur_cells_[idx + nx].p_,
			 cur_cells_[idx + nx + 1].p_
			};
			Color res =
				(Cell::pc1 * (max_p - p[0]) * inv_delta) + Cell::pc2 * ((p[0] - min_p) * inv_delta) * 0.25f +
				(Cell::pc1 * (max_p - p[1]) * inv_delta) + Cell::pc2 * ((p[1] - min_p) * inv_delta) * 0.25f +
				(Cell::pc1 * (max_p - p[2]) * inv_delta) + Cell::pc2 * ((p[2] - min_p) * inv_delta) * 0.25f +
				(Cell::pc1 * (max_p - p[3]) * inv_delta) + Cell::pc2 * ((p[3] - min_p) * inv_delta) * 0.25f;
			gfx_.PutPixel(nx / 2 + x / 2, ny - y / 2 - 1, res); // right top
		}
	}

	// plot b
	Float min_b = b_[0];
	Float max_b = b_[0];

	for (int i = 0; i < nc; ++i)
	{
			min_b = std::min(min_b, b_[i]);
			max_b = std::max(max_b, b_[i]);
	}

	for (int y = 0; y < ny; y += 2)
	{
		for (int x = 0; x < nx; x += 2)
		{
			int idx = y * nx + x;
			Float inv_delta = 1 / (max_b - min_b);
			Float b[4] = {
			 b_[idx],
			 b_[idx + 1],
			 b_[idx + nx],
			 b_[idx + nx + 1]
			};
			Color res =
				(Cell::bc1 * ((max_b - b[0]) * inv_delta) + Cell::bc2 * ((b[0] - min_b) * inv_delta)) * 0.25f +
				(Cell::bc1 * ((max_b - b[1]) * inv_delta) + Cell::bc2 * ((b[1] - min_b) * inv_delta)) * 0.25f +
				(Cell::bc1 * ((max_b - b[2]) * inv_delta) + Cell::bc2 * ((b[2] - min_b) * inv_delta)) * 0.25f +
				(Cell::bc1 * ((max_b - b[3]) * inv_delta) + Cell::bc2 * ((b[3] - min_b) * inv_delta)) * 0.25f;
			gfx_.PutPixel(x / 2, ny - ny / 2 - y / 2 - 1, res); // right top
		}
	}

	OutputDebugStringA(("Min p = " + std::to_string(min_p) + ", max p = " + std::to_string(max_p) + "\n").c_str());
	OutputDebugStringA(("Min b = " + std::to_string(min_b) + ", max b = " + std::to_string(max_b) + "\n").c_str());
	OutputDebugStringA(("Min vel = " + std::to_string(min_mag) + ", max vel = " + std::to_string(max_mag) + "\n").c_str());
	OutputDebugStringA("\n");
}

void Simulation::InitField()
{
	for (int x = 0; x < nx; ++x)
	{
		cur_cells_[x].p_ = cur_cells_[x + nx].p_;
		cur_cells_[x].u_ = Vec2(kZeroF);
		cur_cells_[x].is_solid_ = true;
		cur_cells_[x + (ny - 1) * nx].p_ = cur_cells_[x + (ny - 2) * nx].p_;
		cur_cells_[x + (ny - 1) * nx].u_ = Vec2(kZeroF);
		cur_cells_[x + (ny - 1) * nx].is_solid_ = true;
	}
	for (int y = 1; y < ny - 1; ++y)
	{
		cur_cells_[y * nx].p_ = const_pressure;
		cur_cells_[y * nx].u_ = cur_cells_[y * nx + 1].u_;
		cur_cells_[y * nx].is_solid_ = true;
		cur_cells_[y * nx + nx - 1].u_ = cur_cells_[y * nx + nx - 2].u_;
		cur_cells_[y * nx + nx - 1].p_ = -const_pressure;
		cur_cells_[y * nx + nx - 1].is_solid_ = true;
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
			Cell& c = cur_cells_[y * nx + x];
			c.u_ = Vec2();
			c.is_solid_ = true;
		}
	}

	ResetBoundaryConditions();
}

void Simulation::ResetBoundaryConditions()
{
	for (int x = 0; x < nx; ++x)
	{
		cur_cells_[x].p_ = cur_cells_[x + nx].p_;
		cur_cells_[x + (ny - 1) * nx].p_ = cur_cells_[x + (ny - 2) * nx].p_;
	}

	for (int y = 1; y < ny - 1; ++y)
	{
		cur_cells_[y * nx].p_ = const_pressure;
		cur_cells_[y * nx].u_ = cur_cells_[y * nx + 1].u_;
		cur_cells_[y * nx + nx - 1].u_ = cur_cells_[y * nx + nx - 2].u_;
		cur_cells_[y * nx + nx - 1].p_ = -const_pressure;
	}

#pragma omp parallel for
	for (int y = 1; y < ny - 1; ++y)
	{
		for (int x = 1; x < nx - 1; ++x)
		{
			int idx = y * nx + x;
			Cell& mid = cur_cells_[idx];

			if (mid.is_solid_ == false)
			{
				continue;
			}

			int count = 0;
			Float avg_pressure = kZeroF;

			// left
			const Cell& left = cur_cells_[idx - 1];
			count += (int)left.is_solid_;
			avg_pressure += left.is_solid_ ? kZeroF : left.p_;

			// right
			const Cell& right = cur_cells_[idx + 1];
			count += (int)right.is_solid_;
			avg_pressure += right.is_solid_ ? kZeroF : right.p_;

			// down
			const Cell& down = cur_cells_[idx - nx];
			count += (int)down.is_solid_;
			avg_pressure += down.is_solid_ ? kZeroF : down.p_;

			// up
			const Cell& up = cur_cells_[idx + nx];
			count += (int)up.is_solid_;
			avg_pressure += up.is_solid_ ? kZeroF : up.p_;

			if (count == 0)
			{
				continue;
			}

			mid.u_ = Vec2();
			mid.p_ = avg_pressure / (Float)count;
		}
	}
}

Float Simulation::GetL1Norm()
{
	Float sum_diff = 0.0f;
	Float sum_old = 0.0f;

	const int loop_count = nx * ny;

#pragma omp parallel for default(shared) reduction(+:sum_old, sum_diff)
	for (int i = 0; i < loop_count; ++i)
	{
		if (cur_cells_[i].is_solid_ == true)
		{
			continue;
		}

		const Float old = Abs(old_cells_[i].p_);
		const Float cur = Abs(cur_cells_[i].p_);
		sum_diff += cur - old;
		sum_old += old;
	}

	return sum_diff / sum_old;
}
