#pragma once
#include "Cell.h"
#include "Graphics.h"

class Simulation
{
public:
	Simulation(Graphics& gfx, Float viscosity, Float density, 
		const Vec2I& dim, Float dt);
	~Simulation();

	void Step();
	void Draw();

private:
	void InitField();
	void ResetBoundaryConditions();
	Float GetL1Norm();

private:
	Float* p, *pn;
	Float* u, *un;
	Float* v, *vn;
	Float* b_;
	bool* is_solid;

	const Float viscosity_;
	const Float density_;
	const Vec2I dim_;
	const int nx, ny;
	const int nc; // number of cells
	const int nbf; // number of bytes for Float
	const int nbb; // number of bytes for bool
	const Float dx, dy;
	const Float dt_;

	Graphics& gfx_;

	static constexpr Float l1norm_target = 1e-3f;
	static constexpr Float const_pressure = 10.0f;
};