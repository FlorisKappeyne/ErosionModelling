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
	int erosion_radius;
	int niter_jacobi;
};

// todo:
/*
better output
param sweep functionality
output of variables
some way to read and process the datadump
analyze numerical precision, try to optimize for it.
*/