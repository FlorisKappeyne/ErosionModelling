/******************************************************************************************
 *	Chili DirectX Framework Version 16.07.20											  *
 *	Game.cpp																			  *
 *	Copyright 2016 PlanetChili.net <http://www.planetchili.net>							  *
 *																						  *
 *	This file is part of The Chili DirectX Framework.									  *
 *																						  *
 *	The Chili DirectX Framework is free software: you can redistribute it and/or modify	  *
 *	it under the terms of the GNU General Public License as published by				  *
 *	the Free Software Foundation, either version 3 of the License, or					  *
 *	(at your option) any later version.													  *
 *																						  *
 *	The Chili DirectX Framework is distributed in the hope that it will be useful,		  *
 *	but WITHOUT ANY WARRANTY; without even the implied warranty of						  *
 *	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the						  *
 *	GNU General Public License for more details.										  *
 *																						  *
 *	You should have received a copy of the GNU General Public License					  *
 *	along with The Chili DirectX Framework.  If not, see <http://www.gnu.org/licenses/>.  *
 ******************************************************************************************/
#include "MainWindow.h"
#include "Game.h"
#include "Vec2.h"
#include "MathUtilities.h"

Game::Game(MainWindow& wnd, Params* parameters, int params_count)
	:
	sim_done(false),
	wnd(wnd),
	gfx(wnd),
	params(parameters),
	n_params(params_count),
	param_iter(0),
	sim(new Simulation(gfx, params[0])),
	steps_per_second(params[0].steps_per_second),
	steps_per_frame((int)Max(steps_per_second / (Float)60.0f, kOneF))
{
}

void Game::Go()
{
	gfx.BeginFrame();
	UpdateModel();
	ComposeFrame();
	gfx.EndFrame();
}

void Game::UpdateModel()
{
	if (sim->GetTimePassed() > 300.0f)
	{
		param_iter++;
		if (param_iter == n_params)
		{
			sim_done = true;
			delete sim;
			return;
		}
		delete sim;
		sim = new Simulation(gfx, params[param_iter]);
		steps_per_second = params[param_iter].steps_per_second;
		steps_per_frame = (int)Max(steps_per_second / (Float)60.0f, kOneF);
	}

	for (int i = 0; i < steps_per_frame; ++i)
	{
		sim->Step();
	}
	//sim.Step();
}

void Game::ComposeFrame()
{
	sim->Draw();
}
