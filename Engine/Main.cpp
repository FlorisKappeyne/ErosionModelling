/******************************************************************************************
*	Chili DirectX Framework Version 16.07.20											  *
*	Main.cpp																			  *
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
#include "ChiliException.h"
#include <algorithm>

int WINAPI wWinMain( HINSTANCE hInst,HINSTANCE,LPWSTR pArgs,INT )
{

	Params params;
	params.file_name_input = "test.png";
	params.file_name_output = "output.buf";
	params.viscosity = 0.089f;
	params.density = 997.0f;
	params.steps_per_second = 3000.0f;
	params.init_time = 5.0f;
	params.erosion_step_time = 1.0f;
	params.nx = 256;
	params.ny = 256;
	params.force_u = 0.0f;
	params.force_v = 0.0f;
	params.field_size_x = 16.0f;
	params.field_size_y = 16.0f;
	params.lid_speed = 1.0f;
	params.inlet_velocity = 1.0f;
	params.outlet_pressure = 0.0f;
	params.erosion_radius = 4;
	params.niter_jacobi = 80;

	try
	{
		MainWindow wnd( hInst,pArgs, params.nx, params.ny );		
		try
		{
			Game theGame( wnd, params );
			while( wnd.ProcessMessage() )
			{
				theGame.Go();
			}
		}
		catch( const ChiliException& e )
		{
			const std::wstring eMsg = e.GetFullMessage() + 
				L"\n\nException caught at Windows message loop.";
			wnd.ShowMessageBox( e.GetExceptionType(),eMsg,MB_ICONERROR );
		}
		catch( const std::exception& e )
		{
			// need to convert std::exception what() string from narrow to wide string
			const std::string whatStr( e.what() );
			const std::wstring eMsg = std::wstring( whatStr.begin(),whatStr.end() ) + 
				L"\n\nException caught at Windows message loop.";
			wnd.ShowMessageBox( L"Unhandled STL Exception",eMsg,MB_ICONERROR );
		}
		catch( ... )
		{
			wnd.ShowMessageBox( L"Unhandled Non-STL Exception",
				L"\n\nException caught at Windows message loop.",MB_ICONERROR );
		}
	}
	catch( const ChiliException& e )
	{
		const std::wstring eMsg = e.GetFullMessage() +
			L"\n\nException caught at main window creation.";
		MessageBox( nullptr,eMsg.c_str(),e.GetExceptionType().c_str(),MB_ICONERROR );
	}
	catch( const std::exception& e )
	{
		// need to convert std::exception what() string from narrow to wide string
		const std::string whatStr( e.what() );
		const std::wstring eMsg = std::wstring( whatStr.begin(),whatStr.end() ) +
			L"\n\nException caught at main window creation.";
		MessageBox( nullptr,eMsg.c_str(),L"Unhandled STL Exception",MB_ICONERROR );
	}
	catch( ... )
	{
		MessageBox( nullptr,L"\n\nException caught at main window creation.",
			L"Unhandled Non-STL Exception",MB_ICONERROR );
	}

	return 0;
}