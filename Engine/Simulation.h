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
	void AddForces();
	void Convect();
	void Diffuse();
	void SubtractPressureGradient();

private:
	Float* p, *pn;
	Float* u, *un;
	Float* v, *vn;
	Float* b_;
	bool* is_solid;

	const Float viscosity_;
	const Float density_;
	const Float force_u_;
	const Float force_v_;
	const Vec2I dim_;
	const int nx, ny;
	const int nc; // number of cells
	const int nbf; // number of bytes for Float
	const int nbb; // number of bytes for bool
	const Float dx, dy;
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