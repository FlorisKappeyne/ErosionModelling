#include "Results.h"
#include <fstream>

Results::Results(Params& params)
	:
	params(params)
{
	int nx = params.nx + 2, ny = params.ny + 2;
	nb_p = nx * ny * sizeof(Float);
	nb_u = (nx - 1) * ny * sizeof(Float);
	nb_v = (ny - 1) * nx * sizeof(Float);
	nb_s = nx * ny * sizeof(Float);
	nb_is_solid = nx * ny * sizeof(Int);
}

Results::~Results()
{
	for (int i = 0; i < buf_p.size(); ++i)
	{
		free(buf_p[i]);
		free(buf_u[i]);
		free(buf_v[i]);
		free(buf_s[i]);
		free(buf_is_solid[i]);
	}
}

void Results::AddSnapshot(Float* u, Float* v, Float* p, Float* s, Int* is_solid, Float time_stamp)
{
	// allocate memory
	Float* p_copy = (Float*)malloc(nb_p);
	Float* u_copy = (Float*)malloc(nb_u);
	Float* v_copy = (Float*)malloc(nb_v);
	Float* s_copy = (Float*)malloc(nb_s);
	Int* is_solid_copy = (Int*)malloc(nb_is_solid);

	// copy over the data from the simulatin
	memcpy(p_copy, p, nb_p);
	memcpy(u_copy, u, nb_u);
	memcpy(v_copy, v, nb_v);
	memcpy(s_copy, s, nb_s);
	memcpy(is_solid_copy, is_solid, nb_is_solid);

	// store that data long term
	buf_p.push_back(p_copy);
	buf_u.push_back(u_copy);
	buf_v.push_back(v_copy);
	buf_s.push_back(s_copy);
	buf_is_solid.push_back(is_solid_copy);
	time_stamps.push_back(time_stamp);
}

void Results::OuputToFile(const std::string& file_name)
{
	std::ofstream out_file = std::ofstream(file_name, std::ios::binary);

	/*
	buffer file contents:

	sizeof(Float)
	sizeof(int)
	sizeof(bool)
	0
	size of params.file_name_input
	params.file_name_input
	size of params.file_name_output
	params.file_name_output 
	params.viscosity
	...
	params.niter_jacobi
	num_snapshots
	all pressure buffers (nb_p * buf_p.size())
	all u buffers (nb_u * buf_p.size())
	...
	all is_solid buffers (nb_is_solid * buf_if_solid.size())
	all timestamps
	*/

	int size_of_Float = sizeof(Float);
	int size_of_int = sizeof(int);
	int size_of_bool = sizeof(bool);

	out_file.write((char*)&size_of_Float, 1);
	out_file.write((char*)&size_of_int, 1);
	out_file.write((char*)&size_of_bool, 1);
	out_file.write("\0", 1); // 4 byte alignment

	int size_of_input = params.file_name_input.size();
	out_file.write((char*)&size_of_input, sizeof(int));
	out_file.write(params.file_name_input.c_str(), size_of_input);

	int size_of_output = params.file_name_output.size();
	out_file.write((char*)&size_of_output, sizeof(int));
	out_file.write(params.file_name_output.c_str(), size_of_output);

	out_file.write((char*)&params.viscosity, sizeof(Float));
	out_file.write((char*)&params.density, sizeof(Float));
	out_file.write((char*)&params.steps_per_second, sizeof(Float));
	out_file.write((char*)&params.init_time, sizeof(Float));
	out_file.write((char*)&params.erosion_step_time, sizeof(Float));
	out_file.write((char*)&params.nx, sizeof(int));
	out_file.write((char*)&params.ny, sizeof(int));
	out_file.write((char*)&params.force_u, sizeof(Float));
	out_file.write((char*)&params.force_v, sizeof(Float));
	out_file.write((char*)&params.field_size_x, sizeof(Float));
	out_file.write((char*)&params.field_size_y, sizeof(Float));
	out_file.write((char*)&params.lid_speed, sizeof(Float));
	out_file.write((char*)&params.inlet_velocity, sizeof(Float));
	out_file.write((char*)&params.outlet_pressure, sizeof(Float));
	out_file.write((char*)&params.erosion_percentile, sizeof(Float));
	out_file.write((char*)&params.niter_jacobi, sizeof(int));

	int num_snapshots = buf_p.size();
	out_file.write((char*)&num_snapshots, sizeof(int));

	for (int i = 0; i < buf_p.size(); ++i)
	{
		out_file.write((char*)buf_p[i], nb_p);
	}

	for (int i = 0; i < buf_u.size(); ++i)
	{
		out_file.write((char*)buf_u[i], nb_u);
	}

	for (int i = 0; i < buf_v.size(); ++i)
	{
		out_file.write((char*)buf_v[i], nb_v);
	}

	for (int i = 0; i < buf_s.size(); ++i)
	{
		out_file.write((char*)buf_s[i], nb_s);
	}

	for (int i = 0; i < buf_is_solid.size(); ++i)
	{
		out_file.write((char*)buf_is_solid[i], nb_is_solid);
	}

	out_file.write((char*)time_stamps.data(), time_stamps.size() * sizeof(Float));

	// close and flush are called in the destructor, so its exception safe.
}
