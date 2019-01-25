#pragma once
#include "Cell.h"
#include "Graphics.h"
#include "Params.h"
#include <string>

class Simulation
{
public:
	Simulation(Graphics& gfx, Params& params);
	~Simulation();

	void Step();
	void Draw();

	Float GetTimePassed()
	{
		return time_passed;
	}

private:
	// field control
	void InitField(const std::string& file_name);
	void ResetBoundaryConditions();
	void ResetEdges();

	// simulation
	void SolveForPressure();
	void UpdateVelocities();
	void SubtractPressureGradient();
	void UpdateErosionProcess();
	void UpdateDeltaTime();
	void CalculateShearStress();
	void ErodeGeometry(Vec2I pos);
	void Sedimentate(int n_cells);

private:
	// buffers
	Float* p, *pn;
	Float* u, *un;
	Float* v, *vn;
	Float* s;
	bool* is_solid;

	// I/O
	std::string file_name_input;
	std::string file_name_output;

	// fluid params
	const Float viscosity;
	const Float density;

	// time control params
	Float dt;
	Float time_passed;
	Float time_until_erosion;
	Float convergence_sim_seconds;

	// simulation params
	const Vec2 dim;
	const int nx, ny;
	const int nc; // number of cells
	const Float dx, dy;
	const Float dx2, dy2;
	const Float force_u;
	const Float force_v;
	const Float lid_speed;
	const Float inlet_velocity;
	const Float outlet_pressure;
	int erosion_radius;
	int niter_jacobi;
	bool do_cavity_flow;

	// quadfloat precalculations
	const QF viscosity_qf;
	const QF density_qf;
	const QF force_u_qf;
	const QF force_v_qf;
	const QF dx_qf, dy_qf;
	const QF dx2_qf, dy2_qf;
	QF dt_qf;

	const QI ones;
	const QI zeros;
	const QF kTwoQF;
	const QF kOneQF;
	const QF kZeroQF;

	// graphics
	Float min_mag;
	Float max_mag;
	Float min_p;
	Float max_p;
	bool drawing_vars_initialized;
	Graphics& gfx;

	static constexpr bool visualize_stress_rt = true;

private:
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

	inline int IndexP(int x, int y)
	{
		return y * nx + x;
	}
	inline int IndexU(int x, int y)
	{
		return y * (nx - 1) + x;
	}
	inline int IndexV(int x, int y)
	{
		return y * nx + x;
	}
};