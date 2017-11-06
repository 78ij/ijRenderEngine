

#ifndef PIPELINE_H
#define PIPELINE_H

#include "Bases.h"
void Line(IJVector, IJVector);

IJPatch *VertexShaderStage1(IJWorld );
IJPatch *VertexShaderStage2(IJWorld, IJPatch *);
IJPatch *RasterizationStage1(IJWorld, IJPatch *);
void FreePatch(IJPatch *, IJWorld);
#endif