#include "Params.h"
#include <fstream>
#include <vector>

void GetAll(const std::string& line, std::vector<std::string>& output)
{
	std::string::size_type new_start = 0;
	std::string::size_type pos = line.find(';');
	do
	{
		output.push_back(line.substr(new_start, pos - new_start));
		new_start = pos + 1;
		pos = line.find(';', new_start);
	} while (pos != std::string::npos);
}
void GetAll(const std::string& line, std::vector<Float>& output)
{
	std::string::size_type new_start = 0;
	std::string::size_type pos = line.find(';');
	do
	{
		// get the substring from the line and convert it to a Float using std::stod
		std::string sub_string = line.substr(new_start, pos - new_start);
		output.push_back((Float)std::stod(sub_string));
		new_start = pos + 1;
		pos = line.find(';', new_start);
	} while (pos != std::string::npos);
}
void GetAll(const std::string& line, std::vector<int>& output)
{
	std::string::size_type new_start = 0;
	std::string::size_type pos = line.find(';');
	do
	{
		// get the substring from the line and convert it to a Float using std::stod
		std::string sub_string = line.substr(new_start, pos - new_start);
		output.push_back((int)std::stoi(sub_string));
		new_start = pos + 1;
		pos = line.find(';', new_start);
	} while (pos != std::string::npos);
}

Params* LoadParams(const std::string & file_name, int& n_params)
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
	std::vector<Float> erosion_step_time;
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
			std::getline(readFile, line);
			GetAll(line, file_name_input);
			permutations *= (int)file_name_input.size();
		}
		else if (line == "file_name_output")
		{
			std::getline(readFile, line);
			GetAll(line, file_name_output);
			permutations *= (int)file_name_output.size();
		}
		else if (line == "viscosity")
		{
			std::getline(readFile, line);
			GetAll(line, viscosity);
			permutations *= (int)viscosity.size();
		}
		else if (line == "density")
		{
			std::getline(readFile, line);
			GetAll(line, density);
			permutations *= (int)density.size();
		}
		else if (line == "steps_per_second")
		{
			std::getline(readFile, line);
			GetAll(line, steps_per_second);
			permutations *= (int)steps_per_second.size();
		}
		else if (line == "init_time")
		{
			std::getline(readFile, line);
			GetAll(line, init_time);
			permutations *= (int)init_time.size();
		}
		else if (line == "erosion_step_time")
		{
			std::getline(readFile, line);
			GetAll(line, erosion_step_time);
			permutations *= (int)erosion_step_time.size();
		}
		else if (line == "nx")
		{
			std::getline(readFile, line);
			GetAll(line, nx);
			permutations *= (int)nx.size();
		}
		else if (line == "ny")
		{
			std::getline(readFile, line);
			GetAll(line, ny);
			permutations *= (int)ny.size();
		}
		else if (line == "force_u")
		{
			std::getline(readFile, line);
			GetAll(line, force_u);
			permutations *= (int)force_u.size();
		}
		else if (line == "force_v")
		{
			std::getline(readFile, line);
			GetAll(line, force_v);
			permutations *= (int)force_v.size();
		}
		else if (line == "field_size_x")
		{
			std::getline(readFile, line);
			GetAll(line, field_size_x);
			permutations *= (int)field_size_x.size();
		}
		else if (line == "field_size_y")
		{
			std::getline(readFile, line);
			GetAll(line, field_size_y);
			permutations *= (int)field_size_y.size();
		}
		else if (line == "lid_speed")
		{
			std::getline(readFile, line);
			GetAll(line, lid_speed);
			permutations *= (int)lid_speed.size();
		}
		else if (line == "inlet_velocity")
		{
			std::getline(readFile, line);
			GetAll(line, inlet_velocity);
			permutations *= (int)inlet_velocity.size();
		}
		else if (line == "outlet_pressure")
		{
			std::getline(readFile, line);
			GetAll(line, outlet_pressure);
			permutations *= (int)outlet_pressure.size();
		}
		else if (line == "erosion_percentile")
		{
			std::getline(readFile, line);
			GetAll(line, erosion_percentile);
			permutations *= (int)erosion_percentile.size();
		}
		else if (line == "niter_jacobi")
		{
			std::getline(readFile, line);
			GetAll(line, niter_jacobi);
			permutations *= (int)niter_jacobi.size();
		}
	}

	Params* res;

	// if something went wrong
	if (permutations == 0)
	{
		n_params = 1;
		res = new Params[n_params];

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
		res->erosion_percentile = 0.05f;
		res->niter_jacobi = 160;
	}
	else
	{
		n_params = permutations;
		res = new Params[n_params];

		int idx = 0;
		for (int a = 0; a < file_name_input.size(); ++a)
		{
			for (int b = 0; b < file_name_output.size(); ++b)
			{
				for (int c = 0; c < viscosity.size(); ++c)
				{
					for (int d = 0; d < density.size(); ++d)
					{
						for (int e = 0; e < steps_per_second.size(); ++e)
						{
							for (int f = 0; f < init_time.size(); ++f)
							{
								for (int g = 0; g < erosion_step_time.size(); ++g)
								{
									for (int h = 0; h < nx.size(); ++h)
									{
										for (int i = 0; i < ny.size(); ++i)
										{
											for (int j = 0; j < force_u.size(); ++j)
											{
												for (int k = 0; k < force_v.size(); ++k)
												{
													for (int l = 0; l < field_size_x.size(); ++l)
													{
														for (int m = 0; m < field_size_y.size(); ++m)
														{
															for (int n = 0; n < lid_speed.size(); ++n)
															{
																for (int o = 0; o < inlet_velocity.size(); ++o)
																{
																	for (int p = 0; p < outlet_pressure.size(); ++p)
																	{
																		for (int q = 0; q < erosion_percentile.size(); ++q)
																		{
																			for (int r = 0; r < niter_jacobi.size(); ++r, ++idx)
																			{
																				res[idx].file_name_input = file_name_input[a];
																				res[idx].file_name_output = file_name_output[b];
																				res[idx].viscosity = viscosity[c];
																				res[idx].density = density[d];
																				res[idx].steps_per_second = steps_per_second[e];
																				res[idx].init_time = init_time[f];
																				res[idx].erosion_step_time = erosion_step_time[g];
																				res[idx].nx = nx[h];
																				res[idx].ny = ny[i];
																				res[idx].force_u = force_u[j];
																				res[idx].force_v = force_v[k];
																				res[idx].field_size_x = field_size_x[l];
																				res[idx].field_size_y = field_size_y[m];
																				res[idx].lid_speed = lid_speed[n];
																				res[idx].inlet_velocity = inlet_velocity[o];
																				res[idx].outlet_pressure = outlet_pressure[p];
																				res[idx].erosion_percentile = erosion_percentile[q];
																				res[idx].niter_jacobi = niter_jacobi[r];
																			}
																		}
																	}
																}
															}
														}
													}
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
	return res;
}
