// SPDX-FileCopyrightText: Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
// SPDX-License-Identifier: BSD-3-Clause
#include "vtkFixedPointVolumeRayCastMIPHelper.h"

#include "vtkCommand.h"
#include "vtkDataArray.h"
#include "vtkFixedPointRayCastImage.h"
#include "vtkFixedPointVolumeRayCastMapper.h"
#include "vtkImageData.h"
#include "vtkObjectFactory.h"
#include "vtkRectilinearGrid.h"
#include "vtkRenderWindow.h"
#include "vtkVolume.h"
#include "vtkVolumeProperty.h"

#include <cmath>

VTK_ABI_NAMESPACE_BEGIN
vtkStandardNewMacro(vtkFixedPointVolumeRayCastMIPHelper);

// Construct a new vtkFixedPointVolumeRayCastMIPHelper with default values
vtkFixedPointVolumeRayCastMIPHelper::vtkFixedPointVolumeRayCastMIPHelper() = default;

// Destruct a vtkFixedPointVolumeRayCastMIPHelper - clean up any memory used
vtkFixedPointVolumeRayCastMIPHelper::~vtkFixedPointVolumeRayCastMIPHelper() = default;

// This method is called when the interpolation type is nearest neighbor and
// the data contains one component. In the inner loop we will compute the
// maximum value (in native type). After we have a maximum value for the ray
// we will convert it to unsigned short using the scale/shift, then use this
// index to lookup the final color/opacity.
template <class T>
void vtkFixedPointMIPHelperGenerateImageOneNN(T* data, int threadID, int threadCount,
  vtkFixedPointVolumeRayCastMapper* mapper, vtkVolume* vtkNotUsed(vol))
{
  VTKKWRCHelper_InitializationAndLoopStartNN();
  VTKKWRCHelper_InitializeMIPOneNN();
  VTKKWRCHelper_SpaceLeapSetup();

  if (cropping)
  {
    int maxValueDefined = 0;
    unsigned short maxIdx = 0;

    for (k = 0; k < numSteps; k++)
    {
      if (k)
      {
        mapper->FixedPointIncrement(pos, dir);
      }

      VTKKWRCHelper_MIPSpaceLeapCheck(maxIdx, maxValueDefined, mapper->GetFlipMIPComparison());

      if (!mapper->CheckIfCropped(pos))
      {
        mapper->ShiftVectorDown(pos, spos);
        dptr = data + spos[0] * inc[0] + spos[1] * inc[1] + spos[2] * inc[2];
        if (!maxValueDefined ||
          ((mapper->GetFlipMIPComparison() && *dptr < maxValue) ||
            (!mapper->GetFlipMIPComparison() && *dptr > maxValue)))
        {
          maxValue = *dptr;
          maxIdx = static_cast<unsigned short>((maxValue + shift[0]) * scale[0]);
          maxValueDefined = 1;
        }
      }
    }

    if (maxValueDefined)
    {
      VTKKWRCHelper_LookupColorMax(colorTable[0], scalarOpacityTable[0], maxIdx, imagePtr);
    }
    else
    {
      imagePtr[0] = imagePtr[1] = imagePtr[2] = imagePtr[3] = 0;
    }
  }
  else
  {
    unsigned short maxIdx = static_cast<unsigned short>((maxValue + shift[0]) * scale[0]);

    for (k = 0; k < numSteps; k++)
    {
      if (k)
      {
        mapper->FixedPointIncrement(pos, dir);
      }

      VTKKWRCHelper_MIPSpaceLeapCheck(maxIdx, 1, mapper->GetFlipMIPComparison());

      mapper->ShiftVectorDown(pos, spos);
      dptr = data + spos[0] * inc[0] + spos[1] * inc[1] + spos[2] * inc[2];
      if (mapper->GetFlipMIPComparison())
      {
        maxValue = (*dptr < maxValue) ? (*dptr) : (maxValue);
      }
      else
      {
        maxValue = (*dptr > maxValue) ? (*dptr) : (maxValue);
      }

      maxIdx = static_cast<unsigned short>((maxValue + shift[0]) * scale[0]);
    }

    VTKKWRCHelper_LookupColorMax(colorTable[0], scalarOpacityTable[0], maxIdx, imagePtr);
  }

  VTKKWRCHelper_IncrementAndLoopEnd();
}

// This method is called when the interpolation type is nearest neighbor and
// the data has two or four dependent components. If it is four, they must be
// unsigned char components.  Compute max of last components in native type,
// then use first component to look up a color (2 component data) or first three
// as the color directly (four component data). Lookup alpha off the last component.
template <class T>
void vtkFixedPointMIPHelperGenerateImageDependentNN(T* data, int threadID, int threadCount,
  vtkFixedPointVolumeRayCastMapper* mapper, vtkVolume* vtkNotUsed(vol))
{
  VTKKWRCHelper_InitializationAndLoopStartNN();
  VTKKWRCHelper_InitializeMIPMultiNN();
  VTKKWRCHelper_SpaceLeapSetup();

  int maxValueDefined = 0;
  unsigned short maxIdxS = 0;

  for (k = 0; k < numSteps; k++)
  {
    if (k)
    {
      mapper->FixedPointIncrement(pos, dir);
    }

    VTKKWRCHelper_MIPSpaceLeapCheck(maxIdxS, maxValueDefined, mapper->GetFlipMIPComparison());
    VTKKWRCHelper_CroppingCheckNN(pos);

    mapper->ShiftVectorDown(pos, spos);
    dptr = data + spos[0] * inc[0] + spos[1] * inc[1] + spos[2] * inc[2];
    if (!maxValueDefined ||
      ((mapper->GetFlipMIPComparison() && *(dptr + components - 1) < maxValue[components - 1]) ||
        (!mapper->GetFlipMIPComparison() && *(dptr + components - 1) > maxValue[components - 1])))
    {
      for (c = 0; c < components; c++)
      {
        maxValue[c] = *(dptr + c);
      }
      maxIdxS = static_cast<unsigned short>(
        (maxValue[components - 1] + shift[components - 1]) * scale[components - 1]);
      maxValueDefined = 1;
    }
  }

  if (maxValueDefined)
  {
    unsigned short maxIdx[4] = { 0, 0, 0, 0 };
    if (components == 2)
    {
      maxIdx[0] = static_cast<unsigned short>((maxValue[0] + shift[0]) * scale[0]);
      maxIdx[1] = static_cast<unsigned short>((maxValue[1] + shift[1]) * scale[1]);
    }
    else
    {
      maxIdx[0] = static_cast<unsigned short>(maxValue[0]);
      maxIdx[1] = static_cast<unsigned short>(maxValue[1]);
      maxIdx[2] = static_cast<unsigned short>(maxValue[2]);
      maxIdx[3] = static_cast<unsigned short>((maxValue[3] + shift[3]) * scale[3]);
    }

    for (c = 0; c < components; c++)
    {
    }
    VTKKWRCHelper_LookupDependentColorUS(
      colorTable[0], scalarOpacityTable[0], maxIdx, components, imagePtr);
  }
  else
  {
    imagePtr[0] = imagePtr[1] = imagePtr[2] = imagePtr[3] = 0;
  }

  VTKKWRCHelper_IncrementAndLoopEnd();
}

// This method is called when the interpolation type is nearest neighbor and
// the data has more than one independent components. We compute the max of
// each component along the ray in native type, then use the scale/shift to
// convert this into an unsigned short index value. We use the index values
// to lookup the color/opacity per component, then use the component weights to
// blend these into one final color.
template <class T>
void vtkFixedPointMIPHelperGenerateImageIndependentNN(
  T* data, int threadID, int threadCount, vtkFixedPointVolumeRayCastMapper* mapper, vtkVolume* vol)
{
  VTKKWRCHelper_InitializeWeights();
  VTKKWRCHelper_InitializationAndLoopStartNN();
  VTKKWRCHelper_InitializeMIPMultiNN();
  VTKKWRCHelper_SpaceLeapSetupMulti();

  int maxValueDefined = 0;
  unsigned short maxIdx[4] = { 0, 0, 0, 0 };

  for (k = 0; k < numSteps; k++)
  {
    if (k)
    {
      mapper->FixedPointIncrement(pos, dir);
    }
    VTKKWRCHelper_CroppingCheckNN(pos);
    VTKKWRCHelper_MIPSpaceLeapPopulateMulti(maxIdx, mapper->GetFlipMIPComparison())

      mapper->ShiftVectorDown(pos, spos);
    dptr = data + spos[0] * inc[0] + spos[1] * inc[1] + spos[2] * inc[2];

    if (!maxValueDefined)
    {
      for (c = 0; c < components; c++)
      {
        maxValue[c] = *(dptr + c);
        maxIdx[c] = static_cast<unsigned short>((maxValue[c] + shift[c]) * scale[c]);
      }
      maxValueDefined = 1;
    }
    else
    {
      for (c = 0; c < components; c++)
      {
        if (VTKKWRCHelper_MIPSpaceLeapCheckMulti(c, mapper->GetFlipMIPComparison()) &&
          ((mapper->GetFlipMIPComparison() && *(dptr + c) < maxValue[c]) ||
            (!mapper->GetFlipMIPComparison() && *(dptr + c) > maxValue[c])))
        {
          maxValue[c] = *(dptr + c);
          maxIdx[c] = static_cast<unsigned short>((maxValue[c] + shift[c]) * scale[c]);
        }
      }
    }
  }

  imagePtr[0] = imagePtr[1] = imagePtr[2] = imagePtr[3] = 0;
  if (maxValueDefined)
  {
    VTKKWRCHelper_LookupAndCombineIndependentColorsMax(
      colorTable, scalarOpacityTable, maxIdx, weights, components, imagePtr);
  }

  VTKKWRCHelper_IncrementAndLoopEnd();
}

// This method is called when the interpolation type is linear, the
// data contains one component and scale = 1.0 and shift = 0.0. This is
// the simple case were we do not need to apply scale/shift in the
// inner loop. In the inner loop we compute the eight cell vertex values
// (if we have changed cells). We compute our weights within the cell
// according to our fractional position within the cell, and apply trilinear
// interpolation to compute the index. We find the maximum index along
// the ray, and then use this to look up a final color.
template <class T>
void vtkFixedPointMIPHelperGenerateImageOneSimpleTrilin(T* dataPtr, int threadID, int threadCount,
  vtkFixedPointVolumeRayCastMapper* mapper, vtkVolume* vtkNotUsed(vol))
{
  VTKKWRCHelper_InitializationAndLoopStartTrilin();
  VTKKWRCHelper_InitializeMIPOneTrilin();
  VTKKWRCHelper_SpaceLeapSetup();

  int maxValueDefined = 0;
  unsigned short maxIdx = 0;
  unsigned int maxScalar = 0;

  for (k = 0; k < numSteps; k++)
  {
    if (k)
    {
      mapper->FixedPointIncrement(pos, dir);
    }

    VTKKWRCHelper_MIPSpaceLeapCheck(maxIdx, maxValueDefined, mapper->GetFlipMIPComparison());
    VTKKWRCHelper_CroppingCheckTrilin(pos);

    mapper->ShiftVectorDown(pos, spos);
    if (spos[0] != oldSPos[0] || spos[1] != oldSPos[1] || spos[2] != oldSPos[2])
    {
      oldSPos[0] = spos[0];
      oldSPos[1] = spos[1];
      oldSPos[2] = spos[2];

      dptr = dataPtr + spos[0] * inc[0] + spos[1] * inc[1] + spos[2] * inc[2];
      VTKKWRCHelper_GetCellScalarValuesSimple(dptr);
      if (mapper->GetFlipMIPComparison())
      {
        maxScalar = (A < B) ? (A) : (B);
        maxScalar = (C < maxScalar) ? (C) : (maxScalar);
        maxScalar = (D < maxScalar) ? (D) : (maxScalar);
        maxScalar = (E < maxScalar) ? (E) : (maxScalar);
        maxScalar = (F < maxScalar) ? (F) : (maxScalar);
        maxScalar = (G < maxScalar) ? (G) : (maxScalar);
        maxScalar = (H < maxScalar) ? (H) : (maxScalar);
      }
      else
      {
        maxScalar = (A > B) ? (A) : (B);
        maxScalar = (C > maxScalar) ? (C) : (maxScalar);
        maxScalar = (D > maxScalar) ? (D) : (maxScalar);
        maxScalar = (E > maxScalar) ? (E) : (maxScalar);
        maxScalar = (F > maxScalar) ? (F) : (maxScalar);
        maxScalar = (G > maxScalar) ? (G) : (maxScalar);
        maxScalar = (H > maxScalar) ? (H) : (maxScalar);
      }
    }

    if (!maxValueDefined ||
      ((mapper->GetFlipMIPComparison() && maxScalar < static_cast<unsigned int>(maxValue)) ||
        (!mapper->GetFlipMIPComparison() && maxScalar > static_cast<unsigned int>(maxValue))))
    {
      VTKKWRCHelper_ComputeWeights(pos);
      VTKKWRCHelper_InterpolateScalar(val);

      if (!maxValueDefined ||
        ((mapper->GetFlipMIPComparison() && val < maxValue) ||
          (!mapper->GetFlipMIPComparison() && val > maxValue)))
      {
        maxValue = val;
        maxIdx = static_cast<unsigned short>(maxValue);
        maxValueDefined = 1;
      }
    }
  }

  if (maxValueDefined)
  {
    VTKKWRCHelper_LookupColorMax(colorTable[0], scalarOpacityTable[0], maxIdx, imagePtr);
  }
  else
  {
    imagePtr[0] = imagePtr[1] = imagePtr[2] = imagePtr[3] = 0;
  }

  VTKKWRCHelper_IncrementAndLoopEnd();
}

// This method is called when the interpolation type is linear, the
// data contains one component and scale != 1.0 or shift != 0.0. This
// means that we need to apply scale/shift in the inner loop to compute
// an unsigned short index value. In the inner loop we compute the eight cell
// vertex values (as unsigned short indices, if we have changed cells). We
// compute our weights within the cell according to our fractional position
// within the cell, and apply trilinear interpolation to compute the index.
// We find the maximum index along the ray, and then use this to look up a
// final color.
template <class T>
void vtkFixedPointMIPHelperGenerateImageOneTrilin(T* dataPtr, int threadID, int threadCount,
  vtkFixedPointVolumeRayCastMapper* mapper, vtkVolume* vtkNotUsed(vol))
{
  VTKKWRCHelper_InitializationAndLoopStartTrilin();
  VTKKWRCHelper_InitializeMIPOneTrilin();
  VTKKWRCHelper_SpaceLeapSetup();

  int maxValueDefined = 0;
  unsigned short maxIdx = 0;
  for (k = 0; k < numSteps; k++)
  {
    if (k)
    {
      mapper->FixedPointIncrement(pos, dir);
    }

    VTKKWRCHelper_CroppingCheckTrilin(pos);
    VTKKWRCHelper_MIPSpaceLeapCheck(maxIdx, maxValueDefined, mapper->GetFlipMIPComparison());

    mapper->ShiftVectorDown(pos, spos);
    if (spos[0] != oldSPos[0] || spos[1] != oldSPos[1] || spos[2] != oldSPos[2])
    {
      oldSPos[0] = spos[0];
      oldSPos[1] = spos[1];
      oldSPos[2] = spos[2];

      dptr = dataPtr + spos[0] * inc[0] + spos[1] * inc[1] + spos[2] * inc[2];
      VTKKWRCHelper_GetCellScalarValues(dptr, scale[0], shift[0]);
    }

    VTKKWRCHelper_ComputeWeights(pos);
    VTKKWRCHelper_InterpolateScalar(val);

    if (!maxValueDefined ||
      ((mapper->GetFlipMIPComparison() && val < maxValue) ||
        (!mapper->GetFlipMIPComparison() && val > maxValue)))
    {
      maxValue = val;
      maxIdx = static_cast<unsigned short>(maxValue);
      maxValueDefined = 1;
    }
  }

  if (maxValueDefined)
  {
    VTKKWRCHelper_LookupColorMax(colorTable[0], scalarOpacityTable[0], maxIdx, imagePtr);
  }
  else
  {
    imagePtr[0] = imagePtr[1] = imagePtr[2] = imagePtr[3] = 0;
  }

  VTKKWRCHelper_IncrementAndLoopEnd();
}

// This method is used when the interpolation type is linear, the data has
// two or four components and the components are not considered independent.
// For four component d>>(VTKKW_FP_SHIFT - 8));ata, the data must be unsigned char in type. In the
// inner loop we get the data value for the eight cell corners (if we have
// changed cells) for all components as unsigned shorts (we use the
// scale/shift to ensure the correct range). We compute our weights within
// the cell according to our fractional position within the cell, and apply
// trilinear interpolation to compute the index values. For two component data,
// We use the first index to lookup a color and the second to look up an opacity
// for this sample. For four component data we use the first three components
// directly as a color, then we look up the opacity using the fourth component.
// We then composite this into the color computed so far along the ray, and
// check if we can terminate at this point (if the accumulated opacity is
// higher than some threshold).
template <class T>
void vtkFixedPointMIPHelperGenerateImageDependentTrilin(T* dataPtr, int threadID, int threadCount,
  vtkFixedPointVolumeRayCastMapper* mapper, vtkVolume* vtkNotUsed(vol))
{
  VTKKWRCHelper_InitializationAndLoopStartTrilin();
  VTKKWRCHelper_InitializeMIPMultiTrilin();
  VTKKWRCHelper_SpaceLeapSetup();

  int maxValueDefined = 0;
  unsigned short maxIdx = 0;
  for (k = 0; k < numSteps; k++)
  {
    if (k)
    {
      mapper->FixedPointIncrement(pos, dir);
    }

    VTKKWRCHelper_CroppingCheckTrilin(pos);
    VTKKWRCHelper_MIPSpaceLeapCheck(maxIdx, maxValueDefined, mapper->GetFlipMIPComparison());

    mapper->ShiftVectorDown(pos, spos);
    if (spos[0] != oldSPos[0] || spos[1] != oldSPos[1] || spos[2] != oldSPos[2])
    {
      oldSPos[0] = spos[0];
      oldSPos[1] = spos[1];
      oldSPos[2] = spos[2];

      if (components == 2)
      {
        for (c = 0; c < components; c++)
        {
          dptr = dataPtr + spos[0] * inc[0] + spos[1] * inc[1] + spos[2] * inc[2] + c;
          VTKKWRCHelper_GetCellComponentScalarValues(dptr, c, scale[c], shift[c]);
        }
      }
      else
      {
        for (c = 0; c < 3; c++)
        {
          dptr = dataPtr + spos[0] * inc[0] + spos[1] * inc[1] + spos[2] * inc[2] + c;
          VTKKWRCHelper_GetCellComponentRawScalarValues(dptr, c);
        }
        dptr = dataPtr + spos[0] * inc[0] + spos[1] * inc[1] + spos[2] * inc[2] + c;
        VTKKWRCHelper_GetCellComponentScalarValues(dptr, 3, scale[3], shift[3]);
      }
    }

    VTKKWRCHelper_ComputeWeights(pos);
    VTKKWRCHelper_InterpolateScalarComponent(val, c, components);

    if (!maxValueDefined ||
      ((mapper->GetFlipMIPComparison() && (val[components - 1] < maxValue[components - 1])) ||
        (!mapper->GetFlipMIPComparison() && (val[components - 1] > maxValue[components - 1]))))
    {
      for (c = 0; c < components; c++)
      {
        maxValue[c] = val[c];
      }
      maxIdx = static_cast<unsigned short>(
        (maxValue[components - 1] + shift[components - 1]) * scale[components - 1]);
      maxValueDefined = 1;
    }
  }

  if (maxValueDefined)
  {
    VTKKWRCHelper_LookupDependentColorUS(
      colorTable[0], scalarOpacityTable[0], maxValue, components, imagePtr);
  }
  else
  {
    imagePtr[0] = imagePtr[1] = imagePtr[2] = imagePtr[3] = 0;
  }

  VTKKWRCHelper_IncrementAndLoopEnd();
}

// This method is used when the interpolation type is linear, the data has
// more than one component and the components are considered independent. In
// the inner loop we get the data value for the eight cell corners (if we have
// changed cells) for all components as an unsigned shorts (we have to use the
// scale/shift to ensure that we obtained unsigned short indices) We compute
// our weights within the cell according to our fractional position within the
// cell, and apply trilinear interpolation to compute a value for each
// component. We do this for each sample along the ray to find a maximum value
// per component, then we look up a color/opacity for each component and blend
// them according to the component weights.
template <class T>
void vtkFixedPointMIPHelperGenerateImageIndependentTrilin(T* dataPtr, int threadID, int threadCount,
  vtkFixedPointVolumeRayCastMapper* mapper, vtkVolume* vol)
{
  VTKKWRCHelper_InitializeWeights();
  VTKKWRCHelper_InitializationAndLoopStartTrilin();
  VTKKWRCHelper_InitializeMIPMultiTrilin();

  int maxValueDefined = 0;
  for (k = 0; k < numSteps; k++)
  {
    if (k)
    {
      mapper->FixedPointIncrement(pos, dir);
    }

    VTKKWRCHelper_CroppingCheckTrilin(pos);

    mapper->ShiftVectorDown(pos, spos);
    if (spos[0] != oldSPos[0] || spos[1] != oldSPos[1] || spos[2] != oldSPos[2])
    {
      oldSPos[0] = spos[0];
      oldSPos[1] = spos[1];
      oldSPos[2] = spos[2];

      for (c = 0; c < components; c++)
      {
        dptr = dataPtr + spos[0] * inc[0] + spos[1] * inc[1] + spos[2] * inc[2] + c;
        VTKKWRCHelper_GetCellComponentScalarValues(dptr, c, scale[c], shift[c]);
      }
    }

    VTKKWRCHelper_ComputeWeights(pos);
    VTKKWRCHelper_InterpolateScalarComponent(val, c, components);

    if (!maxValueDefined)
    {
      for (c = 0; c < components; c++)
      {
        maxValue[c] = val[c];
      }
      maxValueDefined = 1;
    }
    else
    {
      for (c = 0; c < components; c++)
      {
        if ((mapper->GetFlipMIPComparison() && val[c] < maxValue[c]) ||
          (!mapper->GetFlipMIPComparison() && val[c] > maxValue[c]))
        {
          maxValue[c] = val[c];
        }
      }
    }
  }

  imagePtr[0] = imagePtr[1] = imagePtr[2] = imagePtr[3] = 0;
  if (maxValueDefined)
  {
    VTKKWRCHelper_LookupAndCombineIndependentColorsMax(
      colorTable, scalarOpacityTable, maxValue, weights, components, imagePtr);
  }

  VTKKWRCHelper_IncrementAndLoopEnd();
}

void vtkFixedPointVolumeRayCastMIPHelper::GenerateImage(
  int threadID, int threadCount, vtkVolume* vol, vtkFixedPointVolumeRayCastMapper* mapper)
{
  void* dataPtr = mapper->GetCurrentScalars()->GetVoidPointer(0);
  int scalarType = mapper->GetCurrentScalars()->GetDataType();

  // Nearest Neighbor interpolate
  if (mapper->ShouldUseNearestNeighborInterpolation(vol))
  {
    // One component data
    if (mapper->GetCurrentScalars()->GetNumberOfComponents() == 1)
    {
      switch (scalarType)
      {
        vtkTemplateMacro(vtkFixedPointMIPHelperGenerateImageOneNN(
          static_cast<VTK_TT*>(dataPtr), threadID, threadCount, mapper, vol));
      }
    }
    // More that one independent components
    else if (vol->GetProperty()->GetIndependentComponents())
    {
      switch (scalarType)
      {
        vtkTemplateMacro(vtkFixedPointMIPHelperGenerateImageIndependentNN(
          static_cast<VTK_TT*>(dataPtr), threadID, threadCount, mapper, vol));
      }
    }
    // Dependent (color) components
    else
    {
      switch (scalarType)
      {
        vtkTemplateMacro(vtkFixedPointMIPHelperGenerateImageDependentNN(
          static_cast<VTK_TT*>(dataPtr), threadID, threadCount, mapper, vol));
      }
    }
  }
  // Trilinear Interpolation
  else
  {
    // One component
    if (mapper->GetCurrentScalars()->GetNumberOfComponents() == 1)
    {
      // Scale == 1.0 and shift == 0.0 - simple case (faster)
      if (mapper->GetTableScale()[0] == 1.0 && mapper->GetTableShift()[0] == 0.0)
      {
        switch (scalarType)
        {
          vtkTemplateMacro(vtkFixedPointMIPHelperGenerateImageOneSimpleTrilin(
            static_cast<VTK_TT*>(dataPtr), threadID, threadCount, mapper, vol));
        }
      }
      // Scale != 1.0 or shift != 0.0 - must apply scale/shift in inner loop
      else
      {
        switch (scalarType)
        {
          vtkTemplateMacro(vtkFixedPointMIPHelperGenerateImageOneTrilin(
            static_cast<VTK_TT*>(dataPtr), threadID, threadCount, mapper, vol));
        }
      }
    }
    // Independent components (more than one)
    else if (vol->GetProperty()->GetIndependentComponents())
    {
      switch (scalarType)
      {
        vtkTemplateMacro(vtkFixedPointMIPHelperGenerateImageIndependentTrilin(
          static_cast<VTK_TT*>(dataPtr), threadID, threadCount, mapper, vol));
      }
    }
    // Dependent components
    else
    {
      switch (scalarType)
      {
        vtkTemplateMacro(vtkFixedPointMIPHelperGenerateImageDependentTrilin(
          static_cast<VTK_TT*>(dataPtr), threadID, threadCount, mapper, vol));
      }
    }
  }
}

// Print method for vtkFixedPointVolumeRayCastMIPHelper
void vtkFixedPointVolumeRayCastMIPHelper::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}

//-----------------------------------------------------------------------------//

vtkStandardNewMacro(vqFixedPointVolumeRayCastAIPHelper);

// Construct a new VQFixedPointVolumeRayCastAIPHelper with default values
vqFixedPointVolumeRayCastAIPHelper::vqFixedPointVolumeRayCastAIPHelper() {}

// Destruct a VQFixedPointVolumeRayCastAIPHelper - clean up any memory used
vqFixedPointVolumeRayCastAIPHelper ::~vqFixedPointVolumeRayCastAIPHelper() {}

// This method is called when the interpolation type is nearest neighbor and
// the data contains one component. In the inner loop we will compute the
// maximum value (in native type). After we have a maximum value for the ray
// we will convert it to unsigned short using the scale/shift, then use this
// index to lookup the final color/opacity.
template <class T>
void vqFixedPointAIPHelperGenerateImageOneNN(T* data, int threadID, int threadCount,
  vtkFixedPointVolumeRayCastMapper* mapper, vtkVolume* vtkNotUsed(vol))
{
  VTKKWRCHelper_InitializationAndLoopStartNN();
  VTKKWRCHelper_InitializeMIPOneNN();
  VTKKWRCHelper_SpaceLeapSetup();

  if (cropping)
  {
    unsigned long long aveSum = 0;
    unsigned long long numSamples = 1;
    unsigned short aveIdx = 0;

    for (k = 0; k < numSteps; k++)
    {
      if (k)
      {
        mapper->FixedPointIncrement(pos, dir);
      }

      if (!mapper->CheckIfCropped(pos))
      {
        mapper->ShiftVectorDown(pos, spos);
        dptr = data + spos[0] * inc[0] + spos[1] * inc[1] + spos[2] * inc[2];
        aveSum += *dptr;
        numSamples++;
      }
    }

    if (numSamples != 1)
    {
      aveIdx = aveSum / numSamples;
      VTKKWRCHelper_LookupColorMax(colorTable[0], scalarOpacityTable[0], aveIdx, imagePtr);
    }
    else
    {
      imagePtr[0] = imagePtr[1] = imagePtr[2] = imagePtr[3] = 0;
    }
  }
  else
  {
    unsigned long long aveSum = 0;
    unsigned long long numSamples = 1;
    unsigned short aveIdx = 0;

    for (k = 0; k < numSteps; k++)
    {
      if (k)
      {
        mapper->FixedPointIncrement(pos, dir);
      }

      mapper->ShiftVectorDown(pos, spos);
      dptr = data + spos[0] * inc[0] + spos[1] * inc[1] + spos[2] * inc[2];
      aveSum += static_cast<unsigned short>((*dptr + shift[0]) * scale[0]);
      numSamples++;
    }

    aveIdx = aveSum / numSamples;
    VTKKWRCHelper_LookupColorMax(colorTable[0], scalarOpacityTable[0], aveIdx, imagePtr);
  }

  VTKKWRCHelper_IncrementAndLoopEnd();
}

// This method is called when the interpolation type is nearest neighbor and
// the data has two or four dependent components. If it is four, they must be
// unsigned char components.  Compute max of last components in native type,
// then use first component to look up a color (2 component data) or first three
// as the color directly (four component data). Lookup alpha off the last component.
template <class T>
void vqFixedPointAIPHelperGenerateImageDependentNN(T* data, int threadID, int threadCount,
  vtkFixedPointVolumeRayCastMapper* mapper, vtkVolume* vtkNotUsed(vol))
{
  VTKKWRCHelper_InitializationAndLoopStartNN();
  VTKKWRCHelper_InitializeMIPMultiNN();
  VTKKWRCHelper_SpaceLeapSetup();

  unsigned long long aveSum[4] = { 0, 0, 0, 0 };
  unsigned long long numSamples = 1;
  unsigned short aveIdx[4] = { 0, 0, 0, 0 };

  for (k = 0; k < numSteps; k++)
  {
    if (k)
    {
      mapper->FixedPointIncrement(pos, dir);
    }

    VTKKWRCHelper_CroppingCheckNN(pos);

    mapper->ShiftVectorDown(pos, spos);
    dptr = data + spos[0] * inc[0] + spos[1] * inc[1] + spos[2] * inc[2];
    for (c = 0; c < components; c++)
    {
      aveSum[c] += *(dptr + c);
    }
    numSamples++;
  }

  if (numSamples != 1)
  {
    for (c = 0; c < components; c++)
    {
      aveIdx[c] = (aveSum[c] / numSamples);
    }
    aveIdx[components - 1] =
      (aveIdx[components - 1] + shift[components - 1]) * scale[components - 1];

    VTKKWRCHelper_LookupDependentColorUS(
      colorTable[0], scalarOpacityTable[0], aveIdx, components, imagePtr);
  }
  else
  {
    imagePtr[0] = imagePtr[1] = imagePtr[2] = imagePtr[3] = 0;
  }

  VTKKWRCHelper_IncrementAndLoopEnd();
}

// This method is called when the interpolation type is nearest neighbor and
// the data has more than one independent components. We compute the max of
// each component along the ray in native type, then use the scale/shift to
// convert this into an unsigned short index value. We use the index values
// to lookup the color/opacity per component, then use the component weights to
// blend these into one final color.
template <class T>
void vqFixedPointAIPHelperGenerateImageIndependentNN(
  T* data, int threadID, int threadCount, vtkFixedPointVolumeRayCastMapper* mapper, vtkVolume* vol)
{
  VTKKWRCHelper_InitializeWeights();
  VTKKWRCHelper_InitializationAndLoopStartNN();
  VTKKWRCHelper_InitializeMIPMultiNN();
  VTKKWRCHelper_SpaceLeapSetupMulti();

  unsigned long long aveSum[4] = { 0, 0, 0, 0 };
  unsigned long long numSamples = 1;
  unsigned short aveIdx[4] = { 0, 0, 0, 0 };
  for (k = 0; k < numSteps; k++)
  {
    if (k)
    {
      mapper->FixedPointIncrement(pos, dir);
    }
    VTKKWRCHelper_CroppingCheckNN(pos);

    mapper->ShiftVectorDown(pos, spos);
    dptr = data + spos[0] * inc[0] + spos[1] * inc[1] + spos[2] * inc[2];

    for (c = 0; c < components; c++)
    {
      unsigned int value = *(dptr + c);
      aveSum[c] += static_cast<unsigned short>((value + shift[c]) * scale[c]);
    }
    numSamples++;
  }

  imagePtr[0] = imagePtr[1] = imagePtr[2] = imagePtr[3] = 0;
  if (numSamples != 1)
  {
    for (c = 0; c < components; c++)
    {
      aveIdx[c] = aveSum[c] / numSamples;
    }
    VTKKWRCHelper_LookupAndCombineIndependentColorsMax(
      colorTable, scalarOpacityTable, aveIdx, weights, components, imagePtr);
  }

  VTKKWRCHelper_IncrementAndLoopEnd();
}

// This method is called when the interpolation type is linear, the
// data contains one component and scale = 1.0 and shift = 0.0. This is
// the simple case were we do not need to apply scale/shift in the
// inner loop. In the inner loop we compute the eight cell vertex values
// (if we have changed cells). We compute our weights within the cell
// according to our fractional position within the cell, and apply trilinear
// interpolation to compute the index. We find the maximum index along
// the ray, and then use this to look up a final color.
template <class T>
void vqFixedPointAIPHelperGenerateImageOneSimpleTrilin(T* dataPtr, int threadID, int threadCount,
  vtkFixedPointVolumeRayCastMapper* mapper, vtkVolume* vtkNotUsed(vol))
{
  VTKKWRCHelper_InitializationAndLoopStartTrilin();
  VTKKWRCHelper_InitializeMIPOneTrilin();
  VTKKWRCHelper_SpaceLeapSetup();

  unsigned long long numSamples = 1;
  unsigned short aveIdx = 0;
  unsigned int aveSum = 0;

  for (k = 0; k < numSteps; k++)
  {
    if (k)
    {
      mapper->FixedPointIncrement(pos, dir);
    }

    VTKKWRCHelper_CroppingCheckTrilin(pos);

    mapper->ShiftVectorDown(pos, spos);
    if (spos[0] != oldSPos[0] || spos[1] != oldSPos[1] || spos[2] != oldSPos[2])
    {
      oldSPos[0] = spos[0];
      oldSPos[1] = spos[1];
      oldSPos[2] = spos[2];

      dptr = dataPtr + spos[0] * inc[0] + spos[1] * inc[1] + spos[2] * inc[2];
      val = static_cast<unsigned short>(*dptr);
      VTKKWRCHelper_GetCellScalarValuesSimple(dptr);
    }

    aveSum += val;
    numSamples++;
  }

  if (numSamples != 1)
  {
    aveIdx = aveSum / numSamples;
    VTKKWRCHelper_LookupColorMax(colorTable[0], scalarOpacityTable[0], aveIdx, imagePtr);
  }
  else
  {
    imagePtr[0] = imagePtr[1] = imagePtr[2] = imagePtr[3] = 0;
  }

  VTKKWRCHelper_IncrementAndLoopEnd();
}

// This method is called when the interpolation type is linear, the
// data contains one component and scale != 1.0 or shift != 0.0. This
// means that we need to apply scale/shift in the inner loop to compute
// an unsigned short index value. In the inner loop we compute the eight cell
// vertex values (as unsigned short indices, if we have changed cells). We
// compute our weights within the cell according to our fractional position
// within the cell, and apply trilinear interpolation to compute the index.
// We find the maximum index along the ray, and then use this to look up a
// final color.
template <class T>
void vqFixedPointAIPHelperGenerateImageOneTrilin(T* dataPtr, int threadID, int threadCount,
  vtkFixedPointVolumeRayCastMapper* mapper, vtkVolume* vtkNotUsed(vol))
{
  VTKKWRCHelper_InitializationAndLoopStartTrilin();
  VTKKWRCHelper_InitializeMIPOneTrilin();
  VTKKWRCHelper_SpaceLeapSetup();

  unsigned long long sumValue = 0;
  unsigned short aveId = 0;

  unsigned long long numSamples = 1;
  for (k = 0; k < numSteps; k++)
  {
    if (k)
    {
      mapper->FixedPointIncrement(pos, dir);
    }

    VTKKWRCHelper_CroppingCheckTrilin(pos);

    mapper->ShiftVectorDown(pos, spos);
    if (spos[0] != oldSPos[0] || spos[1] != oldSPos[1] || spos[2] != oldSPos[2])
    {
      oldSPos[0] = spos[0];
      oldSPos[1] = spos[1];
      oldSPos[2] = spos[2];

      dptr = dataPtr + spos[0] * inc[0] + spos[1] * inc[1] + spos[2] * inc[2];
      VTKKWRCHelper_GetCellScalarValues(dptr, scale[0], shift[0]);
    }

    VTKKWRCHelper_ComputeWeights(pos);
    VTKKWRCHelper_InterpolateScalar(val);

    sumValue += val;
    numSamples++;
  }

  if (numSamples != 1)
  {
    aveId = sumValue / numSamples;
    VTKKWRCHelper_LookupColorMax(colorTable[0], scalarOpacityTable[0], aveId, imagePtr);
  }
  else
  {
    imagePtr[0] = imagePtr[1] = imagePtr[2] = imagePtr[3] = 0;
  }

  VTKKWRCHelper_IncrementAndLoopEnd();
}

// This method is used when the interpolation type is linear, the data has
// two or four components and the components are not considered independent.
// For four component d>>(VTKKW_FP_SHIFT - 8));ata, the data must be unsigned char in type. In the
// inner loop we get the data value for the eight cell corners (if we have
// changed cells) for all components as unsigned shorts (we use the
// scale/shift to ensure the correct range). We compute our weights within
// the cell according to our fractional position within the cell, and apply
// trilinear interpolation to compute the index values. For two component data,
// We use the first index to lookup a color and the second to look up an opacity
// for this sample. For four component data we use the first three components
// directly as a color, then we look up the opacity using the fourth component.
// We then composite this into the color computed so far along the ray, and
// check if we can terminate at this point (if the accumulated opacity is
// higher than some threshold).
template <class T>
void vqFixedPointAIPHelperGenerateImageDependentTrilin(T* dataPtr, int threadID, int threadCount,
  vtkFixedPointVolumeRayCastMapper* mapper, vtkVolume* vtkNotUsed(vol))
{
  VTKKWRCHelper_InitializationAndLoopStartTrilin();
  VTKKWRCHelper_InitializeMIPMultiTrilin();
  VTKKWRCHelper_SpaceLeapSetup();

  unsigned long long sumValue[4] = { 0, 0, 0, 0 };
  unsigned short aveId[4] = { 0, 0, 0, 0 };
  unsigned long long numSamples = 1;
  for (k = 0; k < numSteps; k++)
  {
    if (k)
    {
      mapper->FixedPointIncrement(pos, dir);
    }

    VTKKWRCHelper_CroppingCheckTrilin(pos);

    mapper->ShiftVectorDown(pos, spos);
    if (spos[0] != oldSPos[0] || spos[1] != oldSPos[1] || spos[2] != oldSPos[2])
    {
      oldSPos[0] = spos[0];
      oldSPos[1] = spos[1];
      oldSPos[2] = spos[2];

      if (components == 2)
      {
        for (c = 0; c < components; c++)
        {
          dptr = dataPtr + spos[0] * inc[0] + spos[1] * inc[1] + spos[2] * inc[2] + c;
          VTKKWRCHelper_GetCellComponentScalarValues(dptr, c, scale[c], shift[c]);
        }
      }
      else
      {
        for (c = 0; c < 3; c++)
        {
          dptr = dataPtr + spos[0] * inc[0] + spos[1] * inc[1] + spos[2] * inc[2] + c;
          VTKKWRCHelper_GetCellComponentRawScalarValues(dptr, c);
        }
        dptr = dataPtr + spos[0] * inc[0] + spos[1] * inc[1] + spos[2] * inc[2] + c;
        VTKKWRCHelper_GetCellComponentScalarValues(dptr, 3, scale[3], shift[3]);
      }
    }

    VTKKWRCHelper_ComputeWeights(pos);
    VTKKWRCHelper_InterpolateScalarComponent(val, c, components);

    for (c = 0; c < components; c++)
    {
      sumValue[c] += val[c];
    }
    numSamples++;
  }

  if (numSamples != 1)
  {
    for (c = 0; c < components; c++)
    {
      aveId[c] = sumValue[c] / numSamples;
    }
    VTKKWRCHelper_LookupDependentColorUS(
      colorTable[0], scalarOpacityTable[0], aveId, components, imagePtr);
  }
  else
  {
    imagePtr[0] = imagePtr[1] = imagePtr[2] = imagePtr[3] = 0;
  }

  VTKKWRCHelper_IncrementAndLoopEnd();
}

// This method is used when the interpolation type is linear, the data has
// more than one component and the components are considered independent. In
// the inner loop we get the data value for the eight cell corners (if we have
// changed cells) for all components as an unsigned shorts (we have to use the
// scale/shift to ensure that we obtained unsigned short indices) We compute
// our weights within the cell according to our fractional position within the
// cell, and apply trilinear interpolation to compute a value for each
// component. We do this for each sample along the ray to find a maximum value
// per component, then we look up a color/opacity for each component and blend
// them according to the component weights.
template <class T>
void vqFixedPointAIPHelperGenerateImageIndependentTrilin(T* dataPtr, int threadID, int threadCount,
  vtkFixedPointVolumeRayCastMapper* mapper, vtkVolume* vol)
{
  VTKKWRCHelper_InitializeWeights();
  VTKKWRCHelper_InitializationAndLoopStartTrilin();
  VTKKWRCHelper_InitializeMIPMultiTrilin();

  unsigned long long sumValue[4] = { 0, 0, 0, 0 };
  unsigned short aveId[4] = { 0, 0, 0, 0 };
  unsigned long long numSamples = 1;
  for (k = 0; k < numSteps; k++)
  {
    if (k)
    {
      mapper->FixedPointIncrement(pos, dir);
    }

    VTKKWRCHelper_CroppingCheckTrilin(pos);

    mapper->ShiftVectorDown(pos, spos);
    if (spos[0] != oldSPos[0] || spos[1] != oldSPos[1] || spos[2] != oldSPos[2])
    {
      oldSPos[0] = spos[0];
      oldSPos[1] = spos[1];
      oldSPos[2] = spos[2];

      for (c = 0; c < components; c++)
      {
        dptr = dataPtr + spos[0] * inc[0] + spos[1] * inc[1] + spos[2] * inc[2] + c;
        VTKKWRCHelper_GetCellComponentScalarValues(dptr, c, scale[c], shift[c]);
      }
    }

    VTKKWRCHelper_ComputeWeights(pos);
    VTKKWRCHelper_InterpolateScalarComponent(val, c, components);

    for (c = 0; c < components; c++)
    {
      sumValue[c] += val[c];
    }
    numSamples++;
  }

  imagePtr[0] = imagePtr[1] = imagePtr[2] = imagePtr[3] = 0;
  if (numSamples != 1)
  {
    for (c = 0; c < components; c++)
    {
      aveId[c] = sumValue[c] / numSamples;
    }
    VTKKWRCHelper_LookupAndCombineIndependentColorsMax(
      colorTable, scalarOpacityTable, aveId, weights, components, imagePtr);
  }

  VTKKWRCHelper_IncrementAndLoopEnd();
}

void vqFixedPointVolumeRayCastAIPHelper::GenerateImage(
  int threadID, int threadCount, vtkVolume* vol, vtkFixedPointVolumeRayCastMapper* mapper)
{
  void* dataPtr = mapper->GetCurrentScalars()->GetVoidPointer(0);
  int scalarType = mapper->GetCurrentScalars()->GetDataType();

  // Nearest Neighbor interpolate
  if (mapper->ShouldUseNearestNeighborInterpolation(vol))
  {
    // One component data
    if (mapper->GetCurrentScalars()->GetNumberOfComponents() == 1)
    {
      switch (scalarType)
      {
        vtkTemplateMacro(vqFixedPointAIPHelperGenerateImageOneNN(
          static_cast<VTK_TT*>(dataPtr), threadID, threadCount, mapper, vol));
      }
    }
    // More that one independent components
    else if (vol->GetProperty()->GetIndependentComponents())
    {
      switch (scalarType)
      {
        vtkTemplateMacro(vqFixedPointAIPHelperGenerateImageIndependentNN(
          static_cast<VTK_TT*>(dataPtr), threadID, threadCount, mapper, vol));
      }
    }
    // Dependent (color) components
    else
    {
      switch (scalarType)
      {
        vtkTemplateMacro(vqFixedPointAIPHelperGenerateImageDependentNN(
          static_cast<VTK_TT*>(dataPtr), threadID, threadCount, mapper, vol));
      }
    }
  }
  // Trilinear Interpolation
  else
  {
    // One component
    if (mapper->GetCurrentScalars()->GetNumberOfComponents() == 1)
    {
      // Scale == 1.0 and shift == 0.0 - simple case (faster)
      if (mapper->GetTableScale()[0] == 1.0 && mapper->GetTableShift()[0] == 0.0)
      {
        switch (scalarType)
        {
          vtkTemplateMacro(vqFixedPointAIPHelperGenerateImageOneSimpleTrilin(
            static_cast<VTK_TT*>(dataPtr), threadID, threadCount, mapper, vol));
        }
      }
      // Scale != 1.0 or shift != 0.0 - must apply scale/shift in inner loop
      else
      {
        switch (scalarType)
        {
          vtkTemplateMacro(vqFixedPointAIPHelperGenerateImageOneTrilin(
            static_cast<VTK_TT*>(dataPtr), threadID, threadCount, mapper, vol));
        }
      }
    }
    // Indepedent components (more than one)
    else if (vol->GetProperty()->GetIndependentComponents())
    {
      switch (scalarType)
      {
        vtkTemplateMacro(vqFixedPointAIPHelperGenerateImageIndependentTrilin(
          static_cast<VTK_TT*>(dataPtr), threadID, threadCount, mapper, vol));
      }
    }
    // Dependent components
    else
    {
      switch (scalarType)
      {
        vtkTemplateMacro(vqFixedPointAIPHelperGenerateImageDependentTrilin(
          static_cast<VTK_TT*>(dataPtr), threadID, threadCount, mapper, vol));
      }
    }
  }
}

// Print method for vqFixedPointVolumeRayCastAIPPHelper
void vqFixedPointVolumeRayCastAIPHelper::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}

//----------------------------------------------------------------------------------------------------
using namespace vqRender;
vtkStandardNewMacro(vqFixedPointVolumeRayCastISOHelper);

// Construct a new VQFixedPointVolumeRayCastAIPHelper with default values
vqFixedPointVolumeRayCastISOHelper::vqFixedPointVolumeRayCastISOHelper() {}

// Destruct a VQFixedPointVolumeRayCastAIPHelper - clean up any memory used
vqFixedPointVolumeRayCastISOHelper ::~vqFixedPointVolumeRayCastISOHelper() {}

#define ComputeVolumeIntensity(Idx, Pos, VAL)                                                      \
  T* imgPtr = dataPtr + Idx[0] * inc[0] + Idx[1] * inc[1] + Idx[2] * inc[2];                       \
  VTKKWRCHelper_GetCellScalarValues(imgPtr, scale[0], shift[0]);                                   \
  VTKKWRCHelper_ComputeWeights(Pos);                                                               \
  VTKKWRCHelper_InterpolateScalar(VAL);

#define ComputeSmoothVolmeIntensity(Idx, Pos, Smoothness, VAL)                                     \
  if (Smoothness == 0)                                                                             \
  {                                                                                                \
    ComputeVolumeIntensity(Idx, Pos, VAL);                                                         \
  }                                                                                                \
  else                                                                                             \
  {                                                                                                \
    unsigned long long volSum = 0;                                                                 \
    unsigned int cnt = 0;                                                                          \
    for (int i = -Smoothness; i <= Smoothness; i++)                                                \
    {                                                                                              \
      for (int j = -Smoothness; j <= Smoothness; j++)                                              \
      {                                                                                            \
        for (int k = -Smoothness; k <= Smoothness; k++)                                            \
        {                                                                                          \
          unsigned int idx1[3] = { Idx[0] + i, Idx[1] + j, Idx[2] + k };                           \
          unsigned int pos1[3] = { Pos[0] + i * 0x8000, Pos[1] + j * 0x8000,                       \
            Pos[2] + k * 0x8000 };                                                                 \
          unsigned short intensity3 = 0;                                                           \
          {                                                                                        \
            ComputeVolumeIntensity(idx2, pos2, intensity3);                                        \
          }                                                                                        \
          cnt++;                                                                                   \
          volSum += intensity3;                                                                    \
        }                                                                                          \
      }                                                                                            \
    }                                                                                              \
    VAL = volSum / cnt;                                                                            \
  }

#define VTKVRCHelper_ComputeGradientNorm(Event, Smoothness)                                        \
  unsigned int sdPos[3];                                                                           \
  mapper->ShiftVectorDown(Event.dPos, sdPos);                                                      \
  int direction[3] = { 2, 2, 2 };                                                                  \
  if (sdPos[0] <= (direction[0] + 1) || sdPos[1] <= (direction[1] + 1) ||                          \
    sdPos[2] <= (direction[2] + 1) || sdPos[0] >= (dim[0] - (direction[0] + 1)) ||                 \
    sdPos[1] >= (dim[1] - (direction[0] + 1)) || sdPos[2] >= dim[2] - (direction[0] + 1))          \
  {                                                                                                \
    Event.norm = vec3(0.0f);                                                                       \
  }                                                                                                \
  else                                                                                             \
  {                                                                                                \
    int dir2[3] = { -1, 1, 1 };                                                                    \
    for (int i = 0; i < 3; i++)                                                                    \
    {                                                                                              \
      unsigned int idx1[3] = { sdPos[0], sdPos[1], sdPos[2] };                                     \
      unsigned int pos1[3] = { Event.dPos[0], Event.dPos[1], Event.dPos[2] };                      \
      unsigned int idx2[3] = { sdPos[0], sdPos[1], sdPos[2] };                                     \
      unsigned int pos2[3] = { Event.dPos[0], Event.dPos[1], Event.dPos[2] };                      \
      unsigned short intensity1 = 1, intensity2 = 0;                                               \
      idx1[i] = idx1[i] - dir2[i] * direction[i];                                                  \
      pos1[i] = pos1[i] - dir2[i] * direction[i] * 0x8000;                                         \
      {                                                                                            \
        ComputeSmoothVolmeIntensity(idx1, pos1, Smoothness, intensity1);                           \
      }                                                                                            \
      idx2[i] = idx2[i] + dir2[i] * direction[i];                                                  \
      pos2[i] = pos2[i] + dir2[i] * direction[i] * 0x8000;                                         \
      {                                                                                            \
        ComputeSmoothVolmeIntensity(idx2, pos2, Smoothness, intensity2);                           \
      }                                                                                            \
      event.norm[i] = intensity2 - intensity1;                                                     \
    }                                                                                              \
    Event.norm.normalize();                                                                        \
  }

#define VTKKWRCHelper_InitShader(Event)                                                            \
  shader.reflectRay = Event.refRay;                                                                \
  unsigned short emission[4] = { 0, 0, 0, 0 };                                                     \
  emission[3] = SHRT_MAX;                                                                          \
  emission[0] = static_cast<unsigned short>(                                                       \
    (colorTable[0][3 * Event.intensity[0]] * emission[3] + 0x7fff) >> (VTKKW_FP_SHIFT));           \
  emission[1] = static_cast<unsigned short>(                                                       \
    (colorTable[0][3 * Event.intensity[0] + 1] * emission[3] + 0x7fff) >> (VTKKW_FP_SHIFT));       \
  emission[2] = static_cast<unsigned short>(                                                       \
    (colorTable[0][3 * Event.intensity[0] + 2] * emission[3] + 0x7fff) >> (VTKKW_FP_SHIFT));       \
  shader.emissionColor = vec4(static_cast<float>(emission[0]), static_cast<float>(emission[1]),    \
    static_cast<float>(emission[2]), static_cast<float>(emission[3]));                             \
  shader.diffuseColor = shader.emissionColor;                                                      \
  shader.specularColor = vec4(static_cast<float>(USHRT_MAX), static_cast<float>(USHRT_MAX),        \
    static_cast<float>(USHRT_MAX), shader.emissionColor[3]);                                       \
  shader.ior = 0.5f;                                                                               \
  shader.glossiness = glossinessCoef;

#define VTKKWRCHelper_InitLight()                                                                  \
  light.type = 1;                                                                                  \
  light.surfaceP = vec3(static_cast<float>(vqCamPos[0]), static_cast<float>(vqCamPos[1]),          \
    static_cast<float>(vqCamPos[2]));                                                              \
  light.surfaceN = vec3(event.inRayDir);

#define VTKKWRCHelper_uniformSampleOneLight(Event)                                                 \
  VTKKWRCHelper_InitShader(Event);                                                                 \
  VTKKWRCHelper_InitLight();                                                                       \
  vec3 scatterPointToLight = light.surfaceP - Event.nPos;                                          \
  scatterPointToLight.normalize();                                                                 \
  vec3 normU = Event.norm.cross(Event.refRay);                                                     \
  normU.normalize();                                                                               \
  vec3 normV = Event.norm.cross(normU);                                                            \
  vec3 localRefRay(                                                                                \
    Event.refRay.dot(normU), Event.refRay.dot(normV), Event.refRay.dot(Event.norm));               \
  vec3 localScatterPointtoLight(scatterPointToLight.dot(normU), scatterPointToLight.dot(normV),    \
    scatterPointToLight.dot(Event.norm));                                                          \
  float localReflectAngle = std::abs(localScatterPointtoLight[2]);                                 \
  vec4 diffuseColor = shader.diffuseColor * localReflectAngle;                                     \
  float lambertianpdf = localReflectAngle * localReflectAngle * INV_PI_F;                          \
  float cosTheta = std::abs(localRefRay.dot(localScatterPointtoLight));                            \
  float microfacepdf =                                                                             \
    (shader.glossiness) * pow(cosTheta, shader.glossiness) / (2.0 * PI_F * 4.0 * cosTheta);        \
  float lightpdf = 1.0 / std::abs(scatterPointToLight.dot(light.surfaceN) * 4.0);                  \
  float weightpdf = lightpdf * lightpdf / (lightpdf * lightpdf + lambertianpdf * lambertianpdf);   \
  weightpdf = (1.0 - weightpdf);                                                                   \
  vec4 specularColor = shader.specularColor * (microfacepdf + microfacepdf) *                      \
    std::abs(scatterPointToLight.dot(light.surfaceN)) * weightpdf;                                 \
  imagePtr[0] += Event.lightEnergy *                                                               \
    (ambientCoef * shader.emissionColor[0] + diffuseCoef * diffuseColor[0] +                       \
      specularCoef * specularColor[0]);                                                            \
  imagePtr[1] += Event.lightEnergy *                                                               \
    (ambientCoef * shader.emissionColor[1] + diffuseCoef * diffuseColor[1] +                       \
      specularCoef * specularColor[1]);                                                            \
  imagePtr[2] += Event.lightEnergy *                                                               \
    (ambientCoef * shader.emissionColor[2] + diffuseCoef * diffuseColor[2] +                       \
      specularCoef * specularColor[2]);                                                            \
  imagePtr[3] = shader.emissionColor[3];

#define VTKKWRCHelper_SampleVolume(MaxStep)                                                        \
  event.event = 0.0f;                                                                              \
  event.valid = 1.0f;                                                                              \
  float energy = 20;                                                                               \
  float remainEnergy = 0.0f;                                                                       \
  float sigmaT = 0.0f;                                                                             \
  float densityScale = 20.0f;                                                                      \
  unsigned int cnt = 0;                                                                            \
  unsigned int zoomStepRate = 1;                                                                   \
  unsigned int zoomDir[3] = { 0, 0, 0 };                                                           \
  vec3 initPos(0.0f), endPos(0.0f);                                                                \
  while (cnt < MaxStep && remainEnergy < energy)                                                   \
  {                                                                                                \
    if (cnt)                                                                                       \
    {                                                                                              \
      zoomDir[0] = ((dir[0] & 0x7fffffff)) / zoomStepRate + (dir[0] & 0x80000000);                 \
      zoomDir[1] = ((dir[1] & 0x7fffffff)) / zoomStepRate + (dir[1] & 0x80000000);                 \
      zoomDir[2] = ((dir[2] & 0x7fffffff)) / zoomStepRate + (dir[2] & 0x80000000);                 \
      mapper->FixedPointIncrement(pos, zoomDir);                                                   \
      endPos[0] = pos[0];                                                                          \
      endPos[1] = pos[1];                                                                          \
      endPos[2] = pos[2];                                                                          \
    }                                                                                              \
    else                                                                                           \
    {                                                                                              \
      initPos[0] = pos[0];                                                                         \
      initPos[1] = pos[1];                                                                         \
      initPos[2] = pos[2];                                                                         \
    }                                                                                              \
    VTKKWRCHelper_CroppingCheckTrilin(pos);                                                        \
    mapper->ShiftVectorDown(pos, spos);                                                            \
    if (spos[0] != oldSPos[0] || spos[1] != oldSPos[1] || spos[2] != oldSPos[2])                   \
    {                                                                                              \
      oldSPos[0] = spos[0];                                                                        \
      oldSPos[1] = spos[1];                                                                        \
      oldSPos[2] = spos[2];                                                                        \
      dptr = dataPtr + spos[0] * inc[0] + spos[1] * inc[1] + spos[2] * inc[2];                     \
      VTKKWRCHelper_GetCellScalarValues(dptr, scale[0], shift[0]);                                 \
    }                                                                                              \
    VTKKWRCHelper_ComputeWeights(pos);                                                             \
    VTKKWRCHelper_InterpolateScalar(val);                                                          \
    unsigned short opacity = scalarOpacityTable[0][val];                                           \
    sigmaT = static_cast<float>(opacity) / SHRT_MAX * densityScale;                                \
    remainEnergy += sigmaT;                                                                        \
    if (sigmaT > 0.01)                                                                             \
    {                                                                                              \
      event.intensity[0] = val;                                                                    \
      event.dPos[0] = pos[0];                                                                      \
      event.dPos[1] = pos[1];                                                                      \
      event.dPos[2] = pos[2];                                                                      \
      event.nPos[0] = spos[0] / dim[0];                                                            \
      event.nPos[1] = spos[1] / dim[1];                                                            \
      event.nPos[2] = spos[2] / dim[2];                                                            \
      {                                                                                            \
        VTKVRCHelper_ComputeGradientNorm(event, 0);                                                \
      }                                                                                            \
      vec3 rayDir(0.0f);                                                                           \
      rayDir = initPos - endPos;                                                                   \
      rayDir.normalize();                                                                          \
      rayDir[0] = -rayDir[0];                                                                      \
      event.inRayDir = rayDir;                                                                     \
      event.refRay = (rayDir - vec3(2.0f) * (vec3(rayDir.dot(event.norm)) * event.norm));          \
      event.refRay.normalize();                                                                    \
      event.lightEnergy = sigmaT / energy;                                                         \
      {                                                                                            \
        VTKKWRCHelper_uniformSampleOneLight(event);                                                \
      }                                                                                            \
    }                                                                                              \
    cnt++;                                                                                         \
  }

#define VTKKWRCHelper_SampleRay(MaxStep) VTKKWRCHelper_SampleVolume(MaxStep);

template <class T>
void vqFixedPointISOHelperGenerateImageOneTrilin(T* dataPtr, int threadID, int threadCount,
  vtkFixedPointVolumeRayCastMapper* mapper, vtkVolume* vol)
{
  // VQFixedPointVolumeRayCastMapper* vqMapper =
  // dynamic_cast<VQFixedPointVolumeRayCastMapper*>(mapper);

  vtkVolumeProperty* property = vol->GetProperty();
  double ambientCoef = property->GetAmbient();
  double diffuseCoef = property->GetDiffuse();
  double specularCoef = property->GetSpecular();
  double glossinessCoef = property->GetSpecularPower();
  double* vqCamPos = mapper->vqGetWorldCameraPosition();
  VTKKWRCHelper_InitializationAndLoopStartTrilin();
  VTKKWRCHelper_InitializeMIPOneTrilin();
  VTKKWRCHelper_SpaceLeapSetup();
  shaderInfo shader;
  lightInfo light;
  scatterEvent event;
  imagePtr[0] = imagePtr[1] = imagePtr[2] = imagePtr[3] = 0;
  VTKKWRCHelper_SampleVolume(numSteps);
  if (!event.valid && event.event == 0)
  {
    // imagePtr[0] = imagePtr[1] = imagePtr[2] = imagePtr[3] = 0;
  }
  if (event.valid && event.event == 1)
  {
    /*light color here*/
  }
  if (event.valid && event.event == 2)
  {
    /*shader for object eg. roi*/
  }

  VTKKWRCHelper_IncrementAndLoopEnd();
}

#define ComputeVolumeIntensityIndependent(Idx, Pos, Val)                                           \
  int c = 0;                                                                                       \
  for (c = 0; c < components; c++)                                                                 \
  {                                                                                                \
    T* imgPtr = dataPtr + Idx[0] * inc[0] + Idx[1] * inc[1] + Idx[2] * inc[2] + c;                 \
    VTKKWRCHelper_GetCellComponentScalarValues(imgPtr, c, scale[c], shift[c]);                     \
  }                                                                                                \
  VTKKWRCHelper_ComputeWeights(Pos);                                                               \
  VTKKWRCHelper_InterpolateScalarComponent(Val, c, components);

#define VTKVRCHelper_ComputeGradientNormIndependent(Event, Smoothness)                             \
  unsigned int sdPos[3];                                                                           \
  mapper->ShiftVectorDown(Event.dPos, sdPos);                                                      \
  int direction[3] = { 2, 2, 2 };                                                                  \
  if (sdPos[0] <= (direction[0] + 1) || sdPos[1] <= (direction[1] + 1) ||                          \
    sdPos[2] <= (direction[2] + 1) || sdPos[0] >= (dim[0] - (direction[0] + 1)) ||                 \
    sdPos[1] >= (dim[1] - (direction[0] + 1)) || sdPos[2] >= dim[2] - (direction[0] + 1))          \
  {                                                                                                \
    Event.norm = vec3(0.0f);                                                                       \
  }                                                                                                \
  else                                                                                             \
  {                                                                                                \
    int dir2[3] = { -1, 1, 1 };                                                                    \
    for (int i = 0; i < 3; i++)                                                                    \
    {                                                                                              \
      unsigned int idx1[3] = { sdPos[0], sdPos[1], sdPos[2] };                                     \
      unsigned int pos1[3] = { Event.dPos[0], Event.dPos[1], Event.dPos[2] };                      \
      unsigned int idx2[3] = { sdPos[0], sdPos[1], sdPos[2] };                                     \
      unsigned int pos2[3] = { Event.dPos[0], Event.dPos[1], Event.dPos[2] };                      \
      unsigned short intensity1[3], intensity2[3];                                                 \
      idx1[i] = idx1[i] - dir2[i] * direction[i];                                                  \
      pos1[i] = pos1[i] - dir2[i] * direction[i] * 0x8000;                                         \
      {                                                                                            \
        ComputeVolumeIntensityIndependent(idx1, pos1, intensity1);                                 \
      }                                                                                            \
      idx2[i] = idx2[i] + dir2[i] * direction[i];                                                  \
      pos2[i] = pos2[i] + dir2[i] * direction[i] * 0x8000;                                         \
      {                                                                                            \
        ComputeVolumeIntensityIndependent(idx2, pos2, intensity2);                                 \
      }                                                                                            \
      event.norm[i] = abslength(intensity2) - abslength(intensity1);                               \
    }                                                                                              \
    Event.norm.normalize();                                                                        \
  }

#define VTKKWRCHelper_SampleVolumeIndependent(MaxStep)                                             \
  event.event = 0.0f;                                                                              \
  event.valid = 1.0f;                                                                              \
  float energy = 20.0f;                                                                            \
  float remainEnergy = 0.0f;                                                                       \
  float sigmaT = 0.0f;                                                                             \
  float densityScale = 20.0f;                                                                      \
  unsigned int cnt = 0;                                                                            \
  unsigned int zoomStepRate = 1;                                                                   \
  unsigned int zoomDir[3] = { 0, 0, 0 };                                                           \
  vec3 initPos(0.0f), endPos(0.0f);                                                                \
  while (cnt < MaxStep && remainEnergy < energy)                                                   \
  {                                                                                                \
    if (cnt)                                                                                       \
    {                                                                                              \
      zoomDir[0] = ((dir[0] & 0x7fffffff)) / zoomStepRate + (dir[0] & 0x80000000);                 \
      zoomDir[1] = ((dir[1] & 0x7fffffff)) / zoomStepRate + (dir[1] & 0x80000000);                 \
      zoomDir[2] = ((dir[2] & 0x7fffffff)) / zoomStepRate + (dir[2] & 0x80000000);                 \
      mapper->FixedPointIncrement(pos, zoomDir);                                                   \
      endPos[0] = pos[0];                                                                          \
      endPos[1] = pos[1];                                                                          \
      endPos[2] = pos[2];                                                                          \
    }                                                                                              \
    else                                                                                           \
    {                                                                                              \
      initPos[0] = pos[0];                                                                         \
      initPos[1] = pos[1];                                                                         \
      initPos[2] = pos[2];                                                                         \
    }                                                                                              \
    VTKKWRCHelper_CroppingCheckTrilin(pos);                                                        \
    mapper->ShiftVectorDown(pos, spos);                                                            \
    unsigned int opacity = 0;                                                                      \
    if (spos[0] != oldSPos[0] || spos[1] != oldSPos[1] || spos[2] != oldSPos[2])                   \
    {                                                                                              \
      oldSPos[0] = spos[0];                                                                        \
      oldSPos[1] = spos[1];                                                                        \
      oldSPos[2] = spos[2];                                                                        \
      for (c = 0; c < components; c++)                                                             \
      {                                                                                            \
        dptr = dataPtr + spos[0] * inc[0] + spos[1] * inc[1] + spos[2] * inc[2] + c;               \
        VTKKWRCHelper_GetCellComponentScalarValues(dptr, c, scale[c], shift[c]);                   \
      }                                                                                            \
    }                                                                                              \
    VTKKWRCHelper_ComputeWeights(pos);                                                             \
    VTKKWRCHelper_InterpolateScalarComponent(val, c, components);                                  \
    for (c = 0; c < components; c++)                                                               \
    {                                                                                              \
      opacity += scalarOpacityTable[c][val[c]];                                                    \
    }                                                                                              \
    opacity /= components;                                                                         \
    sigmaT = float(opacity) / SHRT_MAX * densityScale;                                             \
    /*zoomStepRate =std::max(1, std::min(10, int(std::round(opacity / 10000 * densityScale ))));   \
     */                                                                                            \
    remainEnergy += sigmaT;                                                                        \
    cnt++;                                                                                         \
  }                                                                                                \
  for (c = 0; c < components; c++)                                                                 \
  {                                                                                                \
    event.intensity[c] = val[c];                                                                   \
  }                                                                                                \
  event.dPos[0] = pos[0];                                                                          \
  event.dPos[1] = pos[1];                                                                          \
  event.dPos[2] = pos[2];                                                                          \
  event.nPos[0] = float(spos[0]) / float(dim[0]);                                                  \
  event.nPos[1] = float(spos[1]) / float(dim[1]);                                                  \
  event.nPos[2] = float(spos[2]) / float(dim[2]);                                                  \
  {                                                                                                \
    VTKVRCHelper_ComputeGradientNormIndependent(event, 0);                                         \
  }                                                                                                \
  vec3 rayDir(0.0f);                                                                               \
  rayDir = initPos - endPos;                                                                       \
  rayDir.normalize();                                                                              \
  rayDir[0] = -rayDir[0];                                                                          \
  event.inRayDir = rayDir;                                                                         \
  event.refRay = (rayDir - vec3(2.0f) * (vec3(rayDir.dot(event.norm)) * event.norm));              \
  event.refRay.normalize();

#define VTKKWRCHelper_InitShaderIndependent(Event)                                                 \
  shader.reflectRay = Event.refRay;                                                                \
  unsigned short emission[4] = { 0, 0, 0, 0 };                                                     \
  emission[3] = USHRT_MAX;                                                                         \
  emission[0] = static_cast<unsigned short>(                                                       \
    (colorTable[0][3 * Event.intensity[0]] * emission[3] + 0x7fff) >> (VTKKW_FP_SHIFT));           \
  emission[1] = static_cast<unsigned short>(                                                       \
    (colorTable[1][3 * Event.intensity[1] + 1] * emission[3] + 0x7fff) >> (VTKKW_FP_SHIFT));       \
  emission[2] = static_cast<unsigned short>(                                                       \
    (colorTable[2][3 * Event.intensity[2] + 2] * emission[3] + 0x7fff) >> (VTKKW_FP_SHIFT));       \
  shader.emissionColor =                                                                           \
    vec4(float(emission[0]), float(emission[1]), float(emission[2]), float(emission[3]));          \
  shader.diffuseColor = shader.emissionColor;                                                      \
  shader.specularColor = vec4(static_cast<float>(USHRT_MAX), static_cast<float>(USHRT_MAX),        \
    static_cast<float>(USHRT_MAX), shader.emissionColor[3]);                                       \
  shader.ior = 0.5;                                                                                \
  shader.glossiness = glossinessCoef;

#define VTKKWRCHelper_uniformSampleOneLightIndependent(Event)                                      \
  VTKKWRCHelper_InitShaderIndependent(Event);                                                      \
  VTKKWRCHelper_InitLight();                                                                       \
  vec3 scatterPointToLight = light.surfaceP - Event.nPos;                                          \
  scatterPointToLight.normalize();                                                                 \
  vec3 normU = Event.norm.cross(Event.refRay);                                                     \
  normU.normalize();                                                                               \
  vec3 normV = Event.norm.cross(normU);                                                            \
  vec3 localRefRay(                                                                                \
    Event.refRay.dot(normU), Event.refRay.dot(normV), Event.refRay.dot(Event.norm));               \
  vec3 localScatterPointtoLight(scatterPointToLight.dot(normU), scatterPointToLight.dot(normV),    \
    scatterPointToLight.dot(Event.norm));                                                          \
  float localReflectAngle = std::abs(localScatterPointtoLight[2]);                                 \
  vec4 diffuseColor = shader.diffuseColor * localReflectAngle;                                     \
  float lambertianpdf = localReflectAngle * localReflectAngle * INV_PI_F;                          \
  float cosTheta = std::abs(localRefRay.dot(localScatterPointtoLight));                            \
  float microfacepdf =                                                                             \
    (shader.glossiness) * pow(cosTheta, shader.glossiness) / (2.0f * PI_F * 4.0f * cosTheta);      \
  float lightpdf = 1.0f / std::abs(scatterPointToLight.dot(light.surfaceN) * 4.0f);                \
  float weightpdf = lightpdf * lightpdf / (lightpdf * lightpdf + lambertianpdf * lambertianpdf);   \
  weightpdf = (1.0f - weightpdf);                                                                  \
  vec4 specularColor = shader.specularColor * (microfacepdf + microfacepdf) *                      \
    std::abs(scatterPointToLight.dot(light.surfaceN)) * weightpdf;                                 \
  imagePtr[0] = ambientCoef * shader.emissionColor[0] + diffuseCoef * diffuseColor[0] +            \
    specularCoef * specularColor[0];                                                               \
  imagePtr[1] = ambientCoef * shader.emissionColor[1] + diffuseCoef * diffuseColor[1] +            \
    specularCoef * specularColor[1];                                                               \
  imagePtr[2] = ambientCoef * shader.emissionColor[2] + diffuseCoef * diffuseColor[2] +            \
    specularCoef * specularColor[2];

template <class T>
void vqFixedPointISOHelperGenerateImageIndependentTrilin(T* dataPtr, int threadID, int threadCount,
  vtkFixedPointVolumeRayCastMapper* mapper, vtkVolume* vol)
{
  double* vqCamPos = mapper->vqGetWorldCameraPosition();
  vtkVolumeProperty* property = vol->GetProperty();
  double ambientCoef = property->GetAmbient();
  double diffuseCoef = property->GetDiffuse();
  double specularCoef = property->GetSpecular();
  double glossinessCoef = property->GetSpecularPower();
  VTKKWRCHelper_InitializeWeights();
  VTKKWRCHelper_InitializationAndLoopStartTrilin();
  VTKKWRCHelper_InitializeMIPMultiTrilin();

  shaderInfo shader;
  lightInfo light;
  scatterEvent event;
  VTKKWRCHelper_SampleVolumeIndependent(numSteps);
  imagePtr[0] = imagePtr[1] = imagePtr[2] = imagePtr[3] = 0;
  if (event.valid && event.event == 0)
  {
    VTKKWRCHelper_uniformSampleOneLightIndependent(event);
  }
  if (event.valid && event.event == 1)
  {
    /*light color here*/
  }
  if (event.valid && event.event == 2)
  {
    /*shader for object eg. roi*/
  }

  VTKKWRCHelper_IncrementAndLoopEnd();
}

void vqFixedPointVolumeRayCastISOHelper::GenerateImage(
  int threadID, int threadCount, vtkVolume* vol, vtkFixedPointVolumeRayCastMapper* mapper)
{
  void* dataPtr = mapper->GetCurrentScalars()->GetVoidPointer(0);
  int scalarType = mapper->GetCurrentScalars()->GetDataType();

  // Nearest Neighbor interpolate
  if (mapper->ShouldUseNearestNeighborInterpolation(vol))
  {
    // Not Implement
  }
  // Trilinear Interpolation
  else
  {
    // One component
    if (mapper->GetCurrentScalars()->GetNumberOfComponents() == 1)
    {
      switch (scalarType)
      {
        vtkTemplateMacro(vqFixedPointISOHelperGenerateImageOneTrilin(
          static_cast<VTK_TT*>(dataPtr), threadID, threadCount, mapper, vol));
      }
    }
    // Indepedent components (more than one)
    else if (vol->GetProperty()->GetIndependentComponents())
    {
      switch (scalarType)
      {
        vtkTemplateMacro(vqFixedPointISOHelperGenerateImageIndependentTrilin(
          static_cast<VTK_TT*>(dataPtr), threadID, threadCount, mapper, vol));
      }
    }
    // Dependent components
    else
    {
      // not implemented
    }
  }
}

// Print method for vqFixedPointVolumeRayCastISOHelper
void vqFixedPointVolumeRayCastISOHelper::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}
VTK_ABI_NAMESPACE_END