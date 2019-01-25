#include "Params.h"

Params * LoadParams(const std::string & file_name, int& n_params)
{
	Params* res = new Params[1];

	res->file_name_input = "test.png";
	res->file_name_output = "output.buf";
	res->viscosity = 0.089f;
	res->density = 997.0f;
	res->steps_per_second = 1500.0f;
	res->init_time = 5.0f;
	res->erosion_step_time = 1.0f;
	res->nx = 256;
	res->ny = 256;
	res->force_u = 0.0f;
	res->force_v = 0.0f;
	res->field_size_x = 16.0f;
	res->field_size_y = 16.0f;
	res->lid_speed = 1.0f;
	res->inlet_velocity = 1.0f;
	res->outlet_pressure = 0.0f;
	res->erosion_radius = 4;
	res->niter_jacobi = 160;

	n_params = 1;
	return res;
}
