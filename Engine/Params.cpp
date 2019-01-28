#include "Params.h"
#include <fstream>
#include <vector>

void GetAll(const std::string& line, std::vector<std::string>& output)
{

}
void GetAll(const std::string& line, std::vector<Float>& output)
{

}
void GetAll(const std::string& line, std::vector<int>& output)
{

}

Params * LoadParams(const std::string & file_name, int& n_params)
{
	std::string line;
	std::ifstream readFile("params.txt");
	int permutations = 1;

	std::vector<std::string> file_name_input;
	std::vector<std::string> file_name_output;
	std::vector<Float> viscosity;
	std::vector<Float> density;
	std::vector<Float> steps_per_second;
	std::vector<Float> init_time;
	std::vector<Float> erosion_step_timef;
	std::vector<int> nx;
	std::vector<int> ny;
	std::vector<Float> force_u;
	std::vector<Float> force_v;
	std::vector<Float> field_size_x;
	std::vector<Float> field_size_y;
	std::vector<Float> lid_speed;
	std::vector<Float> inlet_velocity;
	std::vector<Float> outlet_pressure;
	std::vector<Float> erosion_percentile;
	std::vector<int> niter_jacobi;

	while (std::getline(readFile, line, '='))
	{
		if (line == "file_name_input")
		{
			//GetAll(std::getline(readFile, line), file_name_input);
		}
		else if (line == "file_name_output")
		{

		}
		else if (line == "viscosity")
		{

		}
		else if (line == "density")
		{

		}
		else if (line == "steps_per_second")
		{

		}
		else if (line == "init_time")
		{

		}
		else if (line == "erosion_step_time")
		{

		}
		else if (line == "nx")
		{

		}
		else if (line == "ny")
		{

		}
		else if (line == "force_u")
		{

		}
		else if (line == "force_v")
		{

		}
		else if (line == "field_size_x")
		{

		}
		else if (line == "field_size_y")
		{

		}
		else if (line == "lid_speed")
		{

		}
		else if (line == "inlet_velocity")
		{

		}
		else if (line == "outlet_pressure")
		{

		}
		else if (line == "erosion_percentile")
		{

		}
		else if (line == "niter_jacobi")
		{

		}
	}

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
