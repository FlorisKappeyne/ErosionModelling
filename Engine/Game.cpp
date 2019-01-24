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

Game::Game(MainWindow& wnd)
	:
	wnd(wnd),
	gfx(wnd),
	sim_(gfx, 8.9e-4f, 997.0f, 16.0f / (Float)gfx.ScreenHeight, kOneF / (Float)steps_per_second),
	steps_per_frame(Max(steps_per_second / 60, 1)),
	stable_sim()
{
	stable_sim.init();
	stable_sim.reset();
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
	//for (int i = 0; i < steps_per_frame * time_multiplier; ++i)
	//{
	//	sim_.Step();
	//}

	sim_.Step();
	//static int steps = 0;
	//if (steps < 100)
	//{
	//	stable_sim.setVX0(50, 50, 0.5f * (kOneF / 60.0f));
	//	stable_sim.setVY0(50, 50, 0.5f * (kOneF / 60.0f));
	//	stable_sim.addSource();
	//}
	//	steps++;
	//
	//stable_sim.vortConfinement();
	//stable_sim.animVel();
	//stable_sim.animDen();
}

void Game::ComposeFrame()
{
	sim_.Draw();
	//// velocity
	//Float *px = stable_sim.getPX();
	//Float *py = stable_sim.getPY();
	//Float *vx = stable_sim.getVX();
	//Float *vy = stable_sim.getVY();
	//
	//// plot magnitude of u
	//Float min_mag = Vec2(vx[0], vy[0]).Magnitude();
	//Float max_mag = min_mag;
	//
	//for (int y = 0, i = 0; y < 128; ++y)
	//{
	//	for (int x = 0; x < 128; ++x, ++i)
	//	{
	//		min_mag = std::min(min_mag, Vec2(vx[i], vy[i]).Magnitude());
	//		max_mag = std::max(max_mag, Vec2(vx[i], vy[i]).Magnitude());
	//	}
	//}
	//
	//for (int y = 0; y < 128; ++y)
	//{
	//	for (int x = 0; x < 128; ++x)
	//	{
	//		int idx = y * 128 + x;
	//
	//		Float inv_delta = 1 / (max_mag - min_mag);
	//		Float mag = Vec2(vx[idx], vy[idx]).Magnitude();
	//		Color res = (Cell::mc1 * ((max_mag - mag) * inv_delta) + Cell::mc2 * ((mag - min_mag) * inv_delta));
	//		gfx.PutPixel(x, gfx.ScreenHeight - y - 1, res); // left top
	//	}
	//}
	//
	//// density
	//Float x;
	//Float y;
	//Float min_d = stable_sim.getDens(1, 1);
	//Float max_d = min_d;
	//
	//for (int i = 1; i <= 128 - 2; i++)
	//{
	//	for (int j = 1; j <= 128 - 2; j++)
	//	{
	//		Float d = stable_sim.getDens(i, j);
	//		min_d = Min(d, min_d);
	//		max_d = Max(d, max_d);
	//	}
	//}
	//
	//for (int i = 1; i <= 128 - 2; i++)
	//{
	//	x = (float)i;
	//	for (int j = 1; j <= 128 - 2; j++)
	//	{
	//		y = (float)j;
	//
	//		int idx = y * 128 + x;
	//		Float inv_delta = 1 / (max_d - min_d);
	//		Float d = stable_sim.getDens(i, j);
	//		Color res = (Cell::pc1 * (max_d - d) * inv_delta) + Cell::pc2 * ((d - min_d) * inv_delta);
	//		gfx.PutPixel(x + 128, gfx.ScreenHeight - y - 1, res); // right top
	//
	//	}
	//}
	//
	//OutputDebugStringA(("Min d = " + std::to_string(min_d) + ", max d = " + std::to_string(max_d) + "\n").c_str());
	//OutputDebugStringA(("Min vel = " + std::to_string(min_mag) + ", max vel = " + std::to_string(max_mag) + "\n").c_str());
	//OutputDebugStringA("\n");
}
