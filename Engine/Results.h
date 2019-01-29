#pragma once
#include "Params.h"
#include <vector>

class Results
{
public:
	Results(Params& params);
	~Results();

	void AddSnapshot(Float* u, Float* v, Float* p, Float* s, Int* is_solid, Float time_stamp);
	void OuputToFile(const std::string& file_name);

public:
	Params params;
	std::vector<Float*> buf_p;
	std::vector<Float*> buf_u;
	std::vector<Float*> buf_v;
	std::vector<Float*> buf_s; // shear stress
	std::vector<Int*> buf_is_solid;
	std::vector<Float> time_stamps;

	int nb_p, nb_u, nb_v, nb_s, nb_is_solid;
};
