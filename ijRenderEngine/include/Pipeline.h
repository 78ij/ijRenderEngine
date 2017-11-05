

#ifndef PIPELINE_H
#define PIPELINE_H

#include "Bases.h"

IJPatch *VertexShaderStage1(IJWorld world);
IJPatch *VertexShaderStage2(IJWorld world, IJPatch *data);
IJPatch *RasterizationStage1(IJWorld world, IJPatch *data);

#endif