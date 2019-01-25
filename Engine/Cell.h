#pragma once
#include "Vec2.h"
#include "Colors.h"

class Cell
{
public:
	Cell();
	Cell(const Vec2& u, Float p, bool is_solid);

public:
	Vec2 u_;
	Float p_;
	bool is_solid_;

	// colors to interpolate between based on magnitude
	static constexpr Color mc1 = Colors::Black;
	static constexpr Color mc2 = Colors::Blue;

	// colors to interpolate between based on pressure
	static constexpr Color pc1 = Colors::Black;
	static constexpr Color pc2 = Colors::Green;

	// colors to interpolate between speed
	static constexpr Color uc1 = Colors::Black;
	static constexpr Color uc2 = Colors::Red;
	static constexpr Color vc1 = Colors::Black;
	static constexpr Color vc2 = Colors::Green;
};