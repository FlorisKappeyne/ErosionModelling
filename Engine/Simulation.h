#pragma once
#include "Cell.h"
#include "Graphics.h"

class Simulation
{
public:
	Simulation(Graphics& gfx, Float viscosity, Float density, 
		Float ds, Float delta_time);
	~Simulation();

	void Step();
	void Draw();

private:

	void InitField();
	void ResetBoundaryConditions();
	void ResetEdges();

	void SolveForPressure();
	void UpdateVelocities();
	void SubtractPressureGradient();

	inline QF Average2(const QF& a, const QF& b)
	{
		return (a + b) / kTwoQF;
	}
	inline QF Average4(const QF& a, const QF& b, const QF& c, const QF& d)
	{
		return (a + b + c + d) / mm_set(Float(4.0));
	}
	inline Float Average4(const Float& a, const Float& b, const Float& c, const Float& d)
	{
		return (a + b + c + d) / Float(4.0);
	}

	// utility functions
	inline Vec2T<QF> Gradient(const QF& left, const QF& right, const QF& down, const QF& up)
	{
		return GradientNonstaggered(left, right, down, up);
	}
	inline Vec2T<QF> GradientNonstaggered(const QF& left, const QF& right, const QF& down, const QF& up)
	{
		return  Vec2T<QF>((right - left) / dx_qf * kTwoQF, (up - down) / dy_qf * kTwoQF);
	}
	inline Vec2T<QF> GradientStaggered(const QF& left, const QF& right, const QF& down, const QF& up)
	{
		return  Vec2T<QF>((right - left) / dx_qf, (up - down) / dy_qf);
	}

	inline QF Divergence(const QF& u_left, const QF& u_right, const QF& v_down, const QF& v_up)
	{
		return DivergenceStaggered(u_left, u_right, v_down, v_up);
	}
	inline QF DivergenceNonstaggered(const QF& u_left, const QF& u_right, const QF& v_down, const QF& v_up)
	{
		return (u_right - u_left) / (dx_qf * kTwoQF) + (v_up - v_down) / (dy_qf * kTwoQF);
	}
	inline QF DivergenceStaggered(const QF& u_left, const QF& u_right, const QF& v_down, const QF& v_up)
	{
		return (u_right - u_left) / dx_qf + (v_up - v_down) / dy_qf;
	}

	inline QF Laplacian(const QF& center, const QF& left, const QF& right, const QF& down, const QF& up)
	{
		return LaplacianNonstaggered(center, left, right, down, up);
	}
	inline QF LaplacianNonstaggered(const QF& center, const QF& left, const QF& right, const QF& down, const QF& up)
	{
		return (left + right - kTwoQF * center) / dx2_qf + (down + up - kTwoQF * center) / dy2_qf;
	}
	inline QF LaplacianStaggered(const QF& center, const QF& left, const QF& right, const QF& down, const QF& up)
	{
		return ((left + right - kTwoQF * center) / dx2_qf + (down + up - kTwoQF * center) / dy2_qf) * mm_set(Float(4.0));
	}
	inline Float Laplacian(const Float& center, const Float& left, const Float& right, const Float& down, const Float& up)
	{
		return (left + right - kTwoF * center) / dx2 + (down + up - kTwoF * center) / dy2;
	}

private:
	Float* p, *pn;
	Float* u, *un;
	Float* v, *vn;
	bool* is_solid;

	const Float viscosity_;
	const Float density_;
	const Float force_u_;
	const Float force_v_;
	const Vec2I dim_;
	const int nx, ny;
	const int nc; // number of cells
	const Float dx, dy;
	const Float dx2, dy2;
	Float dt;
	Float time_passed;

	const QF viscosity_qf;
	const QF density_qf;
	const QF force_u_qf;
	const QF force_v_qf;
	const QF dx_qf, dy_qf;
	const QF dx2_qf, dy2_qf;
	QF dt_qf;

	QI ones = mm_set(Int(-1));
	QI zeros;
	QF kTwoQF;
	Graphics& gfx_;

	static constexpr Float const_pressure = 10.0f;
	static constexpr int niter_jacobi = 80;
};