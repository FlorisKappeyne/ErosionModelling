#pragma once
#include "Typedefs.h"
#include <string>

struct Params
{
public:
	// I/O
	std::string file_name_input;
	std::string file_name_output;

	// fluid params
	Float viscosity;
	Float density;

	// time control
	Float steps_per_second;
	Float init_time;
	Float erosion_step_time; 

	// simulation params
	int nx;
	int ny;
	Float force_u;
	Float force_v;
	Float field_size_x;
	Float field_size_y;
	Float lid_speed;
	Float inlet_velocity;
	Float outlet_pressure;
	Float erosion_percentile;
	int niter_jacobi;
};

Params* LoadParams(const std::string& file_name, int& n_params);
