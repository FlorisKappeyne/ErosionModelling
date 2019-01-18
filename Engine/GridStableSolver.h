/** File:    GridStableSolver.h
 ** Author:  Dongli Zhang
 ** Contact: dongli.zhang0129@gmail.com
 **
 ** Copyright (C) Dongli Zhang 2013
 **
 ** This program is free software;  you can redistribute it and/or modify
 ** it under the terms of the GNU General Public License as published by
 ** the Free Software Foundation; either version 2 of the License, or
 ** (at your option) any later version.
 **
 ** This program is distributed in the hope that it will be useful,
 ** but WITHOUT ANY WARRANTY;  without even the implied warranty of
 ** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See
 ** the GNU General Public License for more details.
 **
 ** You should have received a copy of the GNU General Public License
 ** along with this program;  if not, write to the Free Software 
 ** Foundation, 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
 */
#pragma once
#include "Typedefs.h"

class StableSolver
{
public:
    StableSolver();
    ~StableSolver();
    void init();
    void reset();
    void cleanBuffer();
    void start(){ running=1; }
    void stop(){ running=0; }
    int isRunning(){ return running; }

    //animation
    void setBoundary(Float *value, int flag);
    void projection();
    void advection(Float *value, Float *value0, Float *u, Float *v, int flag);
    void diffusion(Float *value, Float *value0, Float rate, int flag);
    void vortConfinement();
    void addSource();
    void animVel();
    void animDen();

    //getter
    int getRowSize(){ return rowSize; }
    int getColSize(){ return colSize; }
    int getTotSize(){ return totSize; }
    Float getH(){ return h; }
    Float getSimSizeX(){ return simSizeX; }
    Float getSimSizeY(){ return simSizeY; }
    Float* getVX(){ return vx; }
    Float* getVY(){ return vy; }
    Float* getD(){ return d; }
    Float* getPX(){ return px; }
    Float* getPY(){ return py; }
    Float getDens(int i, int j){ return (d[cIdx(i-1, j-1)]+d[cIdx(i, j-1)]+d[cIdx(i-1, j)]+d[cIdx(i, j)])/4.0f; }

    //setter
    void setVX0(int i, int j, Float value){ vx0[cIdx(i, j)]=value; }
    void setVY0(int i, int j, Float value){ vy0[cIdx(i, j)]=value; }
    void setD0(int i, int j, Float value){ d0[cIdx(i, j)]=value; }

private:
    int cIdx(int i, int j){ return j*rowSize+i; }

private:
    int rowSize;
    int colSize;
    int totSize;
    Float h;
    Float simSizeX;
    Float simSizeY;
    Float minX;
    Float maxX;
    Float minY;
    Float maxY;

    //params
    int running;
    Float visc;
    Float diff;
    Float vorticity;
    Float timeStep;

    Float *vx;
    Float *vy;
    Float *vx0;
    Float *vy0;
    Float *d;
    Float *d0;
    Float *px;
    Float *py;
    Float *div;
    Float *p;

    //vorticity confinement
    Float *vort;
    Float *absVort;
    Float *gradVortX;
    Float *gradVortY;
    Float *lenGrad;
    Float *vcfx;
    Float *vcfy;
};
