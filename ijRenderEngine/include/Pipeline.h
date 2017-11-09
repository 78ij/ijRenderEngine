/*
*                  The main pipeline header of the engine.
*
*                       Created by 78ij in 2017.11
*/

#ifndef PIPELINE_H
#define PIPELINE_H

#include "Bases.h"

void Line(IJVector,IJVector,IJColor *);
IJPatch *VertexShaderStage1(IJWorld);
IJPatch *VertexShaderStage2(IJWorld, IJPatch *);
IJPatch *RasterizationStage1(IJWorld, IJPatch *);
IJPatch *RasterizationStage2(IJWorld, IJPatch *);
void FreePatch(IJPatch *, IJWorld);

#endif