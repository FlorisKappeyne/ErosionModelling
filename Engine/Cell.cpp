#include "Cell.h"

Cell::Cell()
	:
	u_(Vec2(kZeroF)),
	p_(kZeroF),
	is_solid_(false)
{
}

Cell::Cell(const Vec2& u, Float p, bool is_solid)
	:
	u_(u),
	p_(p),
	is_solid_(is_solid)
{
}