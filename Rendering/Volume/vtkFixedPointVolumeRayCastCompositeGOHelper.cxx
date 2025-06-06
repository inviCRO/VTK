/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkFixedPointVolumeRayCastCompositeGOHelper.cxx
  Language:  C++

  Copyright (c) 1993-2002 Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkFixedPointVolumeRayCastCompositeGOHelper.h"

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

vtkStandardNewMacro(vtkFixedPointVolumeRayCastCompositeGOHelper);

// Construct a new vtkFixedPointVolumeRayCastCompositeGOHelper with default values
vtkFixedPointVolumeRayCastCompositeGOHelper::vtkFixedPointVolumeRayCastCompositeGOHelper()
{
  m_compositeMethod = RegularComposite;
  m_invert = false;
  m_weight = 0.1;
  m_threshold = 0.5;
  m_transPeriod = -2;
  for (int i = 0; i < 3; i++)
  {
    m_channelWeight[i] = 1.0;
    m_channelWeight[2] = 3.0;
  }
}

// Destruct a vtkFixedPointVolumeRayCastCompositeGOHelper - clean up any memory used
vtkFixedPointVolumeRayCastCompositeGOHelper::~vtkFixedPointVolumeRayCastCompositeGOHelper() =
  default;

// This method is used when the interpolation type is nearest neighbor and
// the data has one component and scale == 1.0 and shift == 0.0. In the inner
// loop we get the data value as an unsigned short, and use this index to
// lookup a color and opacity for this sample. We then composite this into
// the color computed so far along the ray, and check if we can terminate at
// this point (if the accumulated opacity is higher than some threshold).
// Finally we move on to the next sample along the ray.
template <class T>
void vtkFixedPointCompositeGOHelperGenerateImageOneSimpleNN(
  T* data, int threadID, int threadCount, vtkFixedPointVolumeRayCastMapper* mapper, vtkVolume* vol)
{
  VTKKWRCHelper_InitializationAndLoopStartGONN();
  VTKKWRCHelper_InitializeCompositeOneNN();
  VTKKWRCHelper_InitializeCompositeGONN();
  VTKKWRCHelper_SpaceLeapSetup();

  for (k = 0; k < numSteps; k++)
  {
    if (k)
    {
      VTKKWRCHelper_MoveToNextSampleGONN();
    }

    VTKKWRCHelper_SpaceLeapCheck();
    VTKKWRCHelper_CroppingCheckNN(pos);

    unsigned short val = static_cast<unsigned short>(((*dptr)));
    unsigned char mag = *magPtr;

    VTKKWRCHelper_LookupColorGOUS(
      colorTable[0], scalarOpacityTable[0], gradientOpacityTable[0], val, mag, tmp);

    if (tmp[3])
    {
      VTKKWRCHelper_CompositeColorAndCheckEarlyTermination(color, tmp, remainingOpacity);
    }
  }

  VTKKWRCHelper_SetPixelColor(imagePtr, color, remainingOpacity);
  VTKKWRCHelper_IncrementAndLoopEnd();
}

// This method is used when the interpolation type is nearest neighbor and
// the data has one component. In the inner loop we get the data value as
// an unsigned short using the scale/shift, and use this index to lookup
// a color and opacity for this sample. We then composite this into the
// color computed so far along the ray, and check if we can terminate at
// this point (if the accumulated opacity is higher than some threshold).
// Finally we move on to the next sample along the ray.
template <class T>
void vtkFixedPointCompositeGOHelperGenerateImageOneNN(
  T* data, int threadID, int threadCount, vtkFixedPointVolumeRayCastMapper* mapper, vtkVolume* vol)
{
  VTKKWRCHelper_InitializationAndLoopStartGONN();
  VTKKWRCHelper_InitializeCompositeOneNN();
  VTKKWRCHelper_InitializeCompositeGONN();
  VTKKWRCHelper_SpaceLeapSetup();

  for (k = 0; k < numSteps; k++)
  {
    if (k)
    {
      VTKKWRCHelper_MoveToNextSampleGONN();
    }

    VTKKWRCHelper_SpaceLeapCheck();
    VTKKWRCHelper_CroppingCheckNN(pos);

    unsigned short val = static_cast<unsigned short>(((*dptr) + shift[0]) * scale[0]);
    unsigned char mag = *magPtr;

    VTKKWRCHelper_LookupColorGOUS(
      colorTable[0], scalarOpacityTable[0], gradientOpacityTable[0], val, mag, tmp);

    if (tmp[3])
    {
      VTKKWRCHelper_CompositeColorAndCheckEarlyTermination(color, tmp, remainingOpacity);
    }
  }

  VTKKWRCHelper_SetPixelColor(imagePtr, color, remainingOpacity);
  VTKKWRCHelper_IncrementAndLoopEnd();
}

// This method is used when the interpolation type is nearest neighbor and
// the data has two components which are not considered independent. In the
// inner loop we compute the two unsigned short index values from the data
// values (using the scale/shift). We use the first index to lookup a color,
// and we use the second index to look up the opacity. We then composite
// the color into the color computed so far along this ray, and check to
// see if we can terminate here (if the opacity accumulated exceed some
// threshold). Finally we move to the next sample along the ray.
template <class T>
void vtkFixedPointCompositeGOHelperGenerateImageTwoDependentNN(
  T* data, int threadID, int threadCount, vtkFixedPointVolumeRayCastMapper* mapper, vtkVolume* vol)
{
  VTKKWRCHelper_InitializationAndLoopStartGONN();
  VTKKWRCHelper_InitializeCompositeOneNN();
  VTKKWRCHelper_InitializeCompositeGONN();
  VTKKWRCHelper_SpaceLeapSetup();

  for (k = 0; k < numSteps; k++)
  {
    if (k)
    {
      VTKKWRCHelper_MoveToNextSampleGONN();
    }

    VTKKWRCHelper_SpaceLeapCheck();
    VTKKWRCHelper_CroppingCheckNN(pos);

    unsigned short val[2];
    val[0] = static_cast<unsigned short>(((*(dptr)) + shift[0]) * scale[0]);
    val[1] = static_cast<unsigned short>(((*(dptr + 1)) + shift[1]) * scale[1]);
    unsigned char mag = *magPtr;

    tmp[3] =
      (scalarOpacityTable[0][val[1]] * gradientOpacityTable[0][mag] + 0x3fff) >> (VTKKW_FP_SHIFT);
    if (!tmp[3])
    {
      continue;
    }

    tmp[0] = static_cast<unsigned short>(
      (colorTable[0][3 * val[0]] * tmp[3] + 0x7fff) >> (VTKKW_FP_SHIFT));
    tmp[1] = static_cast<unsigned short>(
      (colorTable[0][3 * val[0] + 1] * tmp[3] + 0x7fff) >> (VTKKW_FP_SHIFT));
    tmp[2] = static_cast<unsigned short>(
      (colorTable[0][3 * val[0] + 2] * tmp[3] + 0x7fff) >> (VTKKW_FP_SHIFT));

    VTKKWRCHelper_CompositeColorAndCheckEarlyTermination(color, tmp, remainingOpacity);
  }

  VTKKWRCHelper_SetPixelColor(imagePtr, color, remainingOpacity);
  VTKKWRCHelper_IncrementAndLoopEnd();
}

// This method is used when the interpolation type is nearest neighbor and
// the data has four components which are not considered independent . This
// means that the first three components directly represent color, and this
// data must be of unsigned char type. In the inner loop we directly access
// the four data values (no scale/shift is needed). The first three are the
// color of this sample and the fourth is used to look up an opacity in the
// scalar opacity transfer function. We then composite this color into the
// color we have accumulated so far along the ray, and check if we can
// terminate here (if our accumulated opacity has exceed some threshold).
// Finally we move onto the next sample along the ray.
template <class T>
void vtkFixedPointCompositeGOHelperGenerateImageFourDependentNN(
  T* data, int threadID, int threadCount, vtkFixedPointVolumeRayCastMapper* mapper, vtkVolume* vol)
{
  VTKKWRCHelper_InitializationAndLoopStartGONN();
  VTKKWRCHelper_InitializeCompositeOneNN();
  VTKKWRCHelper_InitializeCompositeGONN();
  VTKKWRCHelper_SpaceLeapSetup();

  for (k = 0; k < numSteps; k++)
  {
    if (k)
    {
      VTKKWRCHelper_MoveToNextSampleGONN();
    }

    VTKKWRCHelper_SpaceLeapCheck();
    VTKKWRCHelper_CroppingCheckNN(pos);

    unsigned short val[4];
    val[0] = *(dptr);
    val[1] = *(dptr + 1);
    val[2] = *(dptr + 2);
    val[3] = static_cast<unsigned short>(((*(dptr + 3)) + shift[3]) * scale[3]);

    unsigned char mag = *magPtr;

    tmp[3] =
      (scalarOpacityTable[0][val[3]] * gradientOpacityTable[0][mag] + 0x3fff) >> (VTKKW_FP_SHIFT);
    if (!tmp[3])
    {
      continue;
    }

    tmp[0] = (val[0] * tmp[3] + 0x7f) >> (8);
    tmp[1] = (val[1] * tmp[3] + 0x7f) >> (8);
    tmp[2] = (val[2] * tmp[3] + 0x7f) >> (8);

    VTKKWRCHelper_CompositeColorAndCheckEarlyTermination(color, tmp, remainingOpacity);
  }

  VTKKWRCHelper_SetPixelColor(imagePtr, color, remainingOpacity);
  VTKKWRCHelper_IncrementAndLoopEnd();
}

// This method is used when the interpolation type is nearest neighbor and
// the data has more than one component and the components are considered to
// be independent. In the inner loop we access each component value, using
// the scale/shift to turn the data value into an unsigned short index. We
// then lookup the color/opacity for each component and combine them according
// to the weighting value for each component. We composite this resulting
// color into the color already accumulated for this ray, and we check
// whether we can terminate here (if the accumulated opacity exceeds some
// threshold). Finally we increment to the next sample on the ray.
template <class T>
void vtkFixedPointCompositeGOHelperGenerateImageIndependentNN(
  T* data, int threadID, int threadCount, vtkFixedPointVolumeRayCastMapper* mapper, vtkVolume* vol)
{
  VTKKWRCHelper_InitializeWeights();
  VTKKWRCHelper_InitializationAndLoopStartGONN();
  VTKKWRCHelper_InitializeCompositeMultiNN();
  VTKKWRCHelper_InitializeCompositeGONN();

  for (k = 0; k < numSteps; k++)
  {
    if (k)
    {
      VTKKWRCHelper_MoveToNextSampleGONN();
    }

    VTKKWRCHelper_CroppingCheckNN(pos);

    unsigned char mag[4] = { 1, 1, 1, 1 };
    for (c = 0; c < components; c++)
    {
      val[c] = static_cast<unsigned short>(((*(dptr + c)) + shift[c]) * scale[c]);
      mag[c] = static_cast<unsigned short>(*(magPtr + c));
    }

    VTKKWRCHelper_LookupAndCombineIndependentColorsGOUS(
      colorTable, scalarOpacityTable, gradientOpacityTable, val, mag, weights, components, tmp);

    if (tmp[3])
    {
      VTKKWRCHelper_CompositeColorAndCheckEarlyTermination(color, tmp, remainingOpacity);
    }
  }

  VTKKWRCHelper_SetPixelColor(imagePtr, color, remainingOpacity);
  VTKKWRCHelper_IncrementAndLoopEnd();
}

// This method is used when the interpolation type is linear and the data
// has one component and scale = 1.0 and shift = 0.0. In the inner loop we
// get the data value for the eight cell corners (if we have changed cells)
// as an unsigned short (the range must be right and we don't need the
// scale/shift). We compute our weights within the cell according to our
// fractional position within the cell, apply trilinear interpolation to
// compute the index, and use this index to lookup a color and opacity for
// this sample. We then composite this into the color computed so far along
// the ray, and check if we can terminate at this point (if the accumulated
// opacity is higher than some threshold). Finally we move on to the next
// sample along the ray.
template <class T>
void vtkFixedPointCompositeGOHelperGenerateImageOneSimpleTrilin(T* data, int threadID,
  int threadCount, vtkFixedPointVolumeRayCastMapper* mapper, vtkVolume* vol, vtkFixedPointVolumeRayCastCompositeGOHelper::VQStruct1* VQ1 = nullptr)
{
  VTKKWRCHelper_InitializationAndLoopStartGOTrilin();
  VTKKWRCHelper_InitializeCompositeOneTrilin();
  VTKKWRCHelper_InitializeCompositeOneGOTrilin();
  VTKKWRCHelper_SpaceLeapSetup();

  int needToSampleGO = 0;
  for (k = 0; k < numSteps; k++)
  {
    if (k)
    {
      mapper->FixedPointIncrement(pos, dir);
    }

    VTKKWRCHelper_SpaceLeapCheck();
    VTKKWRCHelper_CroppingCheckTrilin(pos);

    mapper->ShiftVectorDown(pos, spos);
    if (spos[0] != oldSPos[0] || spos[1] != oldSPos[1] || spos[2] != oldSPos[2])
    {
      oldSPos[0] = spos[0];
      oldSPos[1] = spos[1];
      oldSPos[2] = spos[2];

      dptr = data + spos[0] * inc[0] + spos[1] * inc[1] + spos[2] * inc[2];
      VTKKWRCHelper_GetCellScalarValuesSimple(dptr);
      magPtrABCD = gradientMag[spos[2]] + spos[0] * mInc[0] + spos[1] * mInc[1];
      magPtrEFGH = gradientMag[spos[2] + 1] + spos[0] * mInc[0] + spos[1] * mInc[1];
      needToSampleGO = 1;
    }

    VTKKWRCHelper_ComputeWeights(pos);
    VTKKWRCHelper_InterpolateScalar(val);

    tmp[3] = scalarOpacityTable[0][val];
    if (!tmp[3])
    {
      continue;
    }

    if (needToSampleGO)
    {
      VTKKWRCHelper_GetCellMagnitudeValues(magPtrABCD, magPtrEFGH);
      needToSampleGO = 0;
    }

    VTKKWRCHelper_InterpolateMagnitude(mag);

    if (VQ1)
    {
      double d_val = static_cast<double>(val) / 255.0;
      double d_mag = static_cast<double>(mag) / 255.0;
      double weighttmp = ((VQ1->weight * d_val) + (0x7fff - VQ1->weight) * d_mag);
      double otmp = (weighttmp - VQ1->threshold);
      double op = 1.0 / (1.0 + exp(VQ1->transPeriod * otmp));
      unsigned short ops = op * 0x7fff;
      tmp[3] = ops > 0x7fff ? 0x7fff : ops;
    }
    else
    {
      tmp[3] = (tmp[3] * gradientOpacityTable[0][mag] + 0x7fff) >> VTKKW_FP_SHIFT;
    }

    if (!tmp[3])
    {
      continue;
    }

    tmp[0] =
      static_cast<unsigned short>((colorTable[0][3 * val] * tmp[3] + 0x7fff) >> (VTKKW_FP_SHIFT));
    tmp[1] = static_cast<unsigned short>(
      (colorTable[0][3 * val + 1] * tmp[3] + 0x7fff) >> (VTKKW_FP_SHIFT));
    tmp[2] = static_cast<unsigned short>(
      (colorTable[0][3 * val + 2] * tmp[3] + 0x7fff) >> (VTKKW_FP_SHIFT));

    VTKKWRCHelper_CompositeColorAndCheckEarlyTermination(color, tmp, remainingOpacity);
  }

  VTKKWRCHelper_SetPixelColor(imagePtr, color, remainingOpacity);
  VTKKWRCHelper_IncrementAndLoopEnd();
}

// This method is used when the interpolation type is linear and the data
// has one component and scale != 1.0 or shift != 0.0. In the inner loop we
// get the data value for the eight cell corners (if we have changed cells)
// as an unsigned short (we use the scale/shift to ensure the correct range).
// We compute our weights within the cell according to our fractional position
// within the cell, apply trilinear interpolation to compute the index, and use
// this index to lookup a color and opacity for this sample. We then composite
// this into the color computed so far along the ray, and check if we can
// terminate at this point (if the accumulated opacity is higher than some
// threshold). Finally we move on to the next sample along the ray.
template <class T>
void vtkFixedPointCompositeGOHelperGenerateImageOneTrilin(T* data, int threadID, int threadCount,
  vtkFixedPointVolumeRayCastMapper* mapper, vtkVolume* vol,
  vtkFixedPointVolumeRayCastCompositeGOHelper::VQStruct1* VQ1 = nullptr)
{
  VTKKWRCHelper_InitializationAndLoopStartGOTrilin();
  VTKKWRCHelper_InitializeCompositeOneTrilin();
  VTKKWRCHelper_InitializeCompositeOneGOTrilin();
  VTKKWRCHelper_SpaceLeapSetup();

  int needToSampleGO = 0;
  for (k = 0; k < numSteps; k++)
  {
    if (k)
    {
      mapper->FixedPointIncrement(pos, dir);
    }

    VTKKWRCHelper_SpaceLeapCheck();
    VTKKWRCHelper_CroppingCheckTrilin(pos);

    mapper->ShiftVectorDown(pos, spos);
    if (spos[0] != oldSPos[0] || spos[1] != oldSPos[1] || spos[2] != oldSPos[2])
    {
      oldSPos[0] = spos[0];
      oldSPos[1] = spos[1];
      oldSPos[2] = spos[2];

      dptr = data + spos[0] * inc[0] + spos[1] * inc[1] + spos[2] * inc[2];
      VTKKWRCHelper_GetCellScalarValues(dptr, scale[0], shift[0]);
      magPtrABCD = gradientMag[spos[2]] + spos[0] * mInc[0] + spos[1] * mInc[1];
      magPtrEFGH = gradientMag[spos[2] + 1] + spos[0] * mInc[0] + spos[1] * mInc[1];
      needToSampleGO = 1;
    }

    VTKKWRCHelper_ComputeWeights(pos);
    VTKKWRCHelper_InterpolateScalar(val);

    tmp[3] = scalarOpacityTable[0][val];
    if (!tmp[3])
    {
      continue;
    }

    if (needToSampleGO)
    {
      VTKKWRCHelper_GetCellMagnitudeValues(magPtrABCD, magPtrEFGH);
      needToSampleGO = 0;
    }
    VTKKWRCHelper_InterpolateMagnitude(mag);

    if (VQ1)
    {
      double d_val = static_cast<double>(val) / 255.0;
      double d_mag = static_cast<double>(mag) / 255.0;

      double weighttmp = ((VQ1->weight * d_val) + (1 - VQ1->weight) * d_mag);
      double otmp = (weighttmp - VQ1->threshold);
      double op = 1.0 / (1.0 + exp(VQ1->transPeriod * otmp));
      unsigned short ops = op * 0x7fff;
      tmp[3] = ops > 0x7fff ? 0x7fff : ops;
    }
    else
    {
      tmp[3] = (tmp[3] * gradientOpacityTable[0][mag] + 0x7fff) >> VTKKW_FP_SHIFT;
    }

    if (!tmp[3])
    {
      continue;
    }

    tmp[0] =
      static_cast<unsigned short>((colorTable[0][3 * val] * tmp[3] + 0x7fff) >> (VTKKW_FP_SHIFT));
    tmp[1] = static_cast<unsigned short>(
      (colorTable[0][3 * val + 1] * tmp[3] + 0x7fff) >> (VTKKW_FP_SHIFT));
    tmp[2] = static_cast<unsigned short>(
      (colorTable[0][3 * val + 2] * tmp[3] + 0x7fff) >> (VTKKW_FP_SHIFT));

    VTKKWRCHelper_CompositeColorAndCheckEarlyTermination(color, tmp, remainingOpacity);
  }

  VTKKWRCHelper_SetPixelColor(imagePtr, color, remainingOpacity);
  VTKKWRCHelper_IncrementAndLoopEnd();
}

// This method is used when the interpolation type is linear, the data has
// two components and the components are not considered independent. In the
// inner loop we get the data value for the eight cell corners (if we have
// changed cells) for both components as an unsigned shorts (we use the
// scale/shift to ensure the correct range). We compute our weights within
// the cell according to our fractional position within the cell, and apply
// trilinear interpolation to compute the two index value. We use the first
// index to lookup a color and the second to look up an opacity for this sample.
// We then composite this into the color computed so far along the ray, and
// check if we can terminate at this point (if the accumulated opacity is
// higher than some threshold). Finally we move on to the next sample along
// the ray.
template <class T>
void vtkFixedPointCompositeGOHelperGenerateImageTwoDependentTrilin(
  T* data, int threadID, int threadCount, vtkFixedPointVolumeRayCastMapper* mapper, vtkVolume* vol)
{
  VTKKWRCHelper_InitializationAndLoopStartGOTrilin();
  VTKKWRCHelper_InitializeCompositeMultiTrilin();
  VTKKWRCHelper_InitializeCompositeOneGOTrilin();
  VTKKWRCHelper_SpaceLeapSetup();

  int needToSampleGO = 0;
  for (k = 0; k < numSteps; k++)
  {
    if (k)
    {
      mapper->FixedPointIncrement(pos, dir);
    }

    VTKKWRCHelper_SpaceLeapCheck();
    VTKKWRCHelper_CroppingCheckTrilin(pos);

    mapper->ShiftVectorDown(pos, spos);
    if (spos[0] != oldSPos[0] || spos[1] != oldSPos[1] || spos[2] != oldSPos[2])
    {
      oldSPos[0] = spos[0];
      oldSPos[1] = spos[1];
      oldSPos[2] = spos[2];

      dptr = data + spos[0] * inc[0] + spos[1] * inc[1] + spos[2] * inc[2];
      VTKKWRCHelper_GetCellComponentScalarValues(dptr, 0, scale[0], shift[0]);

      dptr++;
      VTKKWRCHelper_GetCellComponentScalarValues(dptr, 1, scale[1], shift[1]);

      magPtrABCD = gradientMag[spos[2]] + spos[0] * mInc[0] + spos[1] * mInc[1];
      magPtrEFGH = gradientMag[spos[2] + 1] + spos[0] * mInc[0] + spos[1] * mInc[1];
      needToSampleGO = 1;
    }

    VTKKWRCHelper_ComputeWeights(pos);
    VTKKWRCHelper_InterpolateScalarComponent(val, c, 2);

    tmp[3] = scalarOpacityTable[0][val[1]];
    if (!tmp[3])
    {
      continue;
    }

    if (needToSampleGO)
    {
      VTKKWRCHelper_GetCellMagnitudeValues(magPtrABCD, magPtrEFGH);
      needToSampleGO = 0;
    }

    VTKKWRCHelper_InterpolateMagnitude(mag);
    tmp[3] = (tmp[3] * gradientOpacityTable[0][mag] + 0x7fff) >> VTKKW_FP_SHIFT;
    if (!tmp[3])
    {
      continue;
    }

    tmp[0] = static_cast<unsigned short>(
      (colorTable[0][3 * val[0]] * tmp[3] + 0x7fff) >> (VTKKW_FP_SHIFT));
    tmp[1] = static_cast<unsigned short>(
      (colorTable[0][3 * val[0] + 1] * tmp[3] + 0x7fff) >> (VTKKW_FP_SHIFT));
    tmp[2] = static_cast<unsigned short>(
      (colorTable[0][3 * val[0] + 2] * tmp[3] + 0x7fff) >> (VTKKW_FP_SHIFT));

    VTKKWRCHelper_CompositeColorAndCheckEarlyTermination(color, tmp, remainingOpacity);
  }

  VTKKWRCHelper_SetPixelColor(imagePtr, color, remainingOpacity);
  VTKKWRCHelper_IncrementAndLoopEnd();
}

// This method is used when the interpolation type is linear, the data has
// four components and the components are not considered independent. In the
// inner loop we get the data value for the eight cell corners (if we have
// changed cells) for all components as an unsigned shorts (we don't have to
// use the scale/shift because only unsigned char data is supported for four
// component data when the components are not independent). We compute our
// weights within the cell according to our fractional position within the cell,
// and apply trilinear interpolation to compute a value for each component. We
// use the first three directly as the color of the sample, and the fourth is
// used to look up an opacity for this sample. We then composite this into the
// color computed so far along the ray, and check if we can terminate at this
// point (if the accumulated opacity is higher than some threshold). Finally we
// move on to the next sample along the ray.
template <class T>
void vtkFixedPointCompositeGOHelperGenerateImageFourDependentTrilin(
  T* data, int threadID, int threadCount, vtkFixedPointVolumeRayCastMapper* mapper, vtkVolume* vol)
{
  VTKKWRCHelper_InitializationAndLoopStartGOTrilin();
  VTKKWRCHelper_InitializeCompositeMultiTrilin();
  VTKKWRCHelper_InitializeCompositeOneGOTrilin();
  VTKKWRCHelper_SpaceLeapSetup();

  int needToSampleGO = 0;
  for (k = 0; k < numSteps; k++)
  {
    if (k)
    {
      mapper->FixedPointIncrement(pos, dir);
    }

    VTKKWRCHelper_SpaceLeapCheck();
    VTKKWRCHelper_CroppingCheckTrilin(pos);

    mapper->ShiftVectorDown(pos, spos);
    if (spos[0] != oldSPos[0] || spos[1] != oldSPos[1] || spos[2] != oldSPos[2])
    {
      oldSPos[0] = spos[0];
      oldSPos[1] = spos[1];
      oldSPos[2] = spos[2];

      dptr = data + spos[0] * inc[0] + spos[1] * inc[1] + spos[2] * inc[2];
      VTKKWRCHelper_GetCellComponentRawScalarValues(dptr, 0);

      dptr++;
      VTKKWRCHelper_GetCellComponentRawScalarValues(dptr, 1);

      dptr++;
      VTKKWRCHelper_GetCellComponentRawScalarValues(dptr, 2);

      dptr++;
      VTKKWRCHelper_GetCellComponentScalarValues(dptr, 3, scale[3], shift[3]);

      magPtrABCD = gradientMag[spos[2]] + spos[0] * mInc[0] + spos[1] * mInc[1];
      magPtrEFGH = gradientMag[spos[2] + 1] + spos[0] * mInc[0] + spos[1] * mInc[1];
      needToSampleGO = 1;
    }

    VTKKWRCHelper_ComputeWeights(pos);
    VTKKWRCHelper_InterpolateScalarComponent(val, c, 4);

    tmp[3] = scalarOpacityTable[0][val[3]];
    if (!tmp[3])
    {
      continue;
    }

    if (needToSampleGO)
    {
      VTKKWRCHelper_GetCellMagnitudeValues(magPtrABCD, magPtrEFGH);
      needToSampleGO = 0;
    }

    VTKKWRCHelper_InterpolateMagnitude(mag);
    tmp[3] = (tmp[3] * gradientOpacityTable[0][mag] + 0x7fff) >> VTKKW_FP_SHIFT;
    if (!tmp[3])
    {
      continue;
    }

    tmp[0] = (val[0] * tmp[3] + 0x7f) >> 8;
    tmp[1] = (val[1] * tmp[3] + 0x7f) >> 8;
    tmp[2] = (val[2] * tmp[3] + 0x7f) >> 8;

    VTKKWRCHelper_CompositeColorAndCheckEarlyTermination(color, tmp, remainingOpacity);
  }

  VTKKWRCHelper_SetPixelColor(imagePtr, color, remainingOpacity);
  VTKKWRCHelper_IncrementAndLoopEnd();
}

// This method is used when the interpolation type is linear, the data has
// more than one component and the components are considered independent. In
// the inner loop we get the data value for the eight cell corners (if we have
// changed cells) for all components as an unsigned shorts (we have to use the
// scale/shift to ensure that we obtained unsigned short indices) We compute our
// weights within the cell according to our fractional position within the cell,
// and apply trilinear interpolation to compute a value for each component. We
// look up a color/opacity for each component and blend them according to the
// component weights. We then composite this resulting color into the
// color computed so far along the ray, and check if we can terminate at this
// point (if the accumulated opacity is higher than some threshold). Finally we
// move on to the next sample along the ray.
template <class T>
void vtkFixedPointCompositeGOHelperGenerateImageIndependentTrilin(T* data, int threadID,
  int threadCount, vtkFixedPointVolumeRayCastMapper* mapper, vtkVolume* vol,
  vtkFixedPointVolumeRayCastCompositeGOHelper::VQStruct2* VQ2 = nullptr)
{
  VTKKWRCHelper_InitializeWeights();
  VTKKWRCHelper_InitializationAndLoopStartGOTrilin();
  VTKKWRCHelper_InitializeCompositeMultiTrilin();
  VTKKWRCHelper_InitializeCompositeMultiGOTrilin();

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

      dptr = data + spos[0] * inc[0] + spos[1] * inc[1] + spos[2] * inc[2];
      VTKKWRCHelper_GetCellComponentScalarValues(dptr, 0, scale[0], shift[0]);

      dptr++;
      VTKKWRCHelper_GetCellComponentScalarValues(dptr, 1, scale[1], shift[1]);

      if (components > 2)
      {
        dptr++;
        VTKKWRCHelper_GetCellComponentScalarValues(dptr, 2, scale[2], shift[2]);
        if (components > 3)
        {
          dptr++;
          VTKKWRCHelper_GetCellComponentScalarValues(dptr, 3, scale[3], shift[3]);
        }
      }

      magPtrABCD = gradientMag[spos[2]] + spos[0] * mInc[0] + spos[1] * mInc[1];
      magPtrEFGH = gradientMag[spos[2] + 1] + spos[0] * mInc[0] + spos[1] * mInc[1];
      VTKKWRCHelper_GetCellComponentMagnitudeValues(magPtrABCD, magPtrEFGH, 0);

      magPtrABCD++;
      magPtrEFGH++;
      VTKKWRCHelper_GetCellComponentMagnitudeValues(magPtrABCD, magPtrEFGH, 1);

      if (components > 2)
      {
        magPtrABCD++;
        magPtrEFGH++;
        VTKKWRCHelper_GetCellComponentMagnitudeValues(magPtrABCD, magPtrEFGH, 2);
        if (components > 3)
        {
          magPtrABCD++;
          magPtrEFGH++;
          VTKKWRCHelper_GetCellComponentMagnitudeValues(magPtrABCD, magPtrEFGH, 3);
        }
      }
    }

    VTKKWRCHelper_ComputeWeights(pos);
    VTKKWRCHelper_InterpolateScalarComponent(val, c, components);
    VTKKWRCHelper_InterpolateMagnitudeComponent(mag, c, components);

    if (VQ2)
    {
      double op = 0.0;
      if (!VQ2->colorProjection)
      {
        double gradientSum = 0.0;
        for (int _idx = 0; _idx < components; _idx++)
        {
          gradientSum += static_cast<double>(mag[_idx]) / 100.0;
        }
        gradientSum /= static_cast<double>(components);
        if (gradientSum < 0.0001)
          continue;

        double scalarSum = 0.0;
        double scalarSum2 = 0.0;
        for (int _idx = 0; _idx < components; _idx++)
        {
          scalarSum += static_cast<double>(val[_idx]) / 255.0 * VQ2->channelWeight[_idx];
          scalarSum2 += static_cast<double>(val[_idx]) / 255.0;
        }
        scalarSum /= static_cast<double>(components);
        scalarSum2 /= static_cast<double>(components);

        double pow = ((VQ2->weight * scalarSum + (1 - VQ2->weight) * gradientSum) - VQ2->threshold);
        op = 1.0 / (1.0 + exp(VQ2->transPeriod * pow));

        if (VQ2->invert)
        {
          op *= (1 - scalarSum2);
        }
        else
        {
          op *= scalarSum2;
        }
      }
      else
      {
        double distance = 0.0;
        for (int idx = 0; idx < components; idx++)
        {
          distance += pow(val[idx] - VQ2->channelWeight[idx], 2);
        }
        distance = sqrt(distance);
        op = distance / (1.732 * 255) * VQ2->transPeriod;
        if (op > 1)
          op = 1;
        if (op < 0)
          op = 0;
        op = 1 - op;
      }

      unsigned short ops = op * 32767;
      ops = ops > 0x7fff ? 0x7fff : ops;
      if (!ops)
      {
        continue;
      }
      unsigned int _tmp[4] = { 0, 0, 0, 0 };

      {
        for (int _idx = 0; _idx < components; _idx++)
        {
          _tmp[0] += static_cast<unsigned short>(
            ((colorTable[_idx][3 * val[_idx]]) * ops + 0x7fff) >> (VTKKW_FP_SHIFT));
          _tmp[1] += static_cast<unsigned short>(
            ((colorTable[_idx][3 * val[_idx] + 1]) * ops + 0x7fff) >> (VTKKW_FP_SHIFT));
          _tmp[2] += static_cast<unsigned short>(
            ((colorTable[_idx][3 * val[_idx] + 2]) * ops + 0x7fff) >> (VTKKW_FP_SHIFT));
        }
      }
      _tmp[3] = ops;
      tmp[0] = (_tmp[0] > 32767) ? (32767) : (_tmp[0]);
      tmp[1] = (_tmp[1] > 32767) ? (32767) : (_tmp[1]);
      tmp[2] = (_tmp[2] > 32767) ? (32767) : (_tmp[2]);
      tmp[3] = (_tmp[3] > 32767) ? (32767) : (_tmp[3]);
    }
    else
    {
      VTKKWRCHelper_LookupAndCombineIndependentColorsGOUS(
        colorTable, scalarOpacityTable, gradientOpacityTable, val, mag, weights, components, tmp);
    }

    VTKKWRCHelper_CompositeColorAndCheckEarlyTermination(color, tmp, remainingOpacity);
  }

  if (VQ2 && VQ2->colorProjection)
  {
    for (int com = 0; com < components; com++)
    {
      color[com] = static_cast<unsigned int>(static_cast<double>(color[com]) * VQ2->weight);
    }
  }

  VTKKWRCHelper_SetPixelColor(imagePtr, color, remainingOpacity);
  VTKKWRCHelper_IncrementAndLoopEnd();
}

void vtkFixedPointVolumeRayCastCompositeGOHelper::GenerateImage(
  int threadID, int threadCount, vtkVolume* vol, vtkFixedPointVolumeRayCastMapper* mapper)
{
  void* data = mapper->GetCurrentScalars()->GetVoidPointer(0);
  int scalarType = mapper->GetCurrentScalars()->GetDataType();

  // Nearest Neighbor interpolate
  if (mapper->ShouldUseNearestNeighborInterpolation(vol))
  {
    // One component data
    if (mapper->GetCurrentScalars()->GetNumberOfComponents() == 1)
    {
      // Scale == 1.0 and shift == 0.0 - simple case (faster)
      if (mapper->GetTableScale()[0] == 1.0 && mapper->GetTableShift()[0] == 0.0)
      {
        switch (scalarType)
        {
          vtkTemplateMacro(vtkFixedPointCompositeGOHelperGenerateImageOneSimpleNN(
            static_cast<VTK_TT*>(data), threadID, threadCount, mapper, vol));
        }
      }
      else
      {
        switch (scalarType)
        {
          vtkTemplateMacro(vtkFixedPointCompositeGOHelperGenerateImageOneNN(
            static_cast<VTK_TT*>(data), threadID, threadCount, mapper, vol));
        }
      }
    }
    // More that one independent components
    else if (vol->GetProperty()->GetIndependentComponents())
    {
      switch (scalarType)
      {
        vtkTemplateMacro(vtkFixedPointCompositeGOHelperGenerateImageIndependentNN(
          static_cast<VTK_TT*>(data), threadID, threadCount, mapper, vol));
      }
    }
    // Dependent (color) components
    else
    {
      // Two components - the first specifies color (through a lookup table) and
      // the second specified opacity (through a lookup table)
      if (mapper->GetCurrentScalars()->GetNumberOfComponents() == 2)
      {
        switch (scalarType)
        {
          vtkTemplateMacro(vtkFixedPointCompositeGOHelperGenerateImageTwoDependentNN(
            static_cast<VTK_TT*>(data), threadID, threadCount, mapper, vol));
        }
      }
      // Four components - they must be unsigned char, the first three directly
      // specify color and the fourth specifies opacity (through a lookup table)
      else
      {
        if (scalarType == VTK_UNSIGNED_CHAR)
        {
          vtkFixedPointCompositeGOHelperGenerateImageFourDependentNN(
            static_cast<unsigned char*>(data), threadID, threadCount, mapper, vol);
        }
        else
        {
          vtkErrorMacro("Four component dependent must be unsigned char!");
        }
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
        if (m_compositeMethod == RegularComposite)
        {
            switch (scalarType)
            {
              vtkTemplateMacro(vtkFixedPointCompositeGOHelperGenerateImageOneSimpleTrilin(
                static_cast<VTK_TT*>(data), threadID, threadCount, mapper, vol));
            }
        }
        else
        {
          vtkFixedPointVolumeRayCastCompositeGOHelper::VQStruct1 temp(
            m_weight, m_threshold, m_transPeriod);

          switch (scalarType)
          {
            vtkTemplateMacro(vtkFixedPointCompositeGOHelperGenerateImageOneSimpleTrilin(
              static_cast<VTK_TT*>(data), threadID, threadCount, mapper, vol, &temp));
          }
        }
      }
      // Scale != 1.0 or shift != 0.0 - must apply scale/shift in inner loop
      else
      {
        if (m_compositeMethod == RegularComposite)
        {
            switch (scalarType)
            {
              vtkTemplateMacro(vtkFixedPointCompositeGOHelperGenerateImageOneTrilin(
                static_cast<VTK_TT*>(data), threadID, threadCount, mapper, vol));
            }
        }
        else
        {
          vtkFixedPointVolumeRayCastCompositeGOHelper::VQStruct1 temp(
            m_weight, m_threshold, m_transPeriod);

          switch (scalarType)
          {
            vtkTemplateMacro(vtkFixedPointCompositeGOHelperGenerateImageOneTrilin(
              static_cast<VTK_TT*>(data), threadID, threadCount, mapper, vol, &temp));
          }
        }
      }
    }
    // Independent components (more than one)
    else if (vol->GetProperty()->GetIndependentComponents())
    {
      if (m_compositeMethod == ColorProjection)
      {
        vtkFixedPointVolumeRayCastCompositeGOHelper::VQStruct2 temp(
          m_weight, m_threshold, m_transPeriod, false, m_channelWeight, true);

        switch (scalarType)
        {
          vtkTemplateMacro(vtkFixedPointCompositeGOHelperGenerateImageIndependentTrilin(
            static_cast<VTK_TT*>(data), threadID, threadCount, mapper, vol, &temp));
        }
      }
      else if (m_compositeMethod == FeatureDetection)
      {
        vtkFixedPointVolumeRayCastCompositeGOHelper::VQStruct2 temp(
          m_weight, m_threshold, m_transPeriod, !m_invert, m_channelWeight, false);

        switch (scalarType)
        {
          vtkTemplateMacro(vtkFixedPointCompositeGOHelperGenerateImageIndependentTrilin(
            static_cast<VTK_TT*>(data), threadID, threadCount, mapper, vol, &temp));
        }
      }
      else
      {
        switch (scalarType)
        {
          vtkTemplateMacro(vtkFixedPointCompositeGOHelperGenerateImageIndependentTrilin(
            static_cast<VTK_TT*>(data), threadID, threadCount, mapper, vol));
        }
      }
    }
    // Dependent components
    else
    {
      // Two components - the first specifies color (through a lookup table)
      // and the second specified opacity (through a lookup table)
      if (mapper->GetCurrentScalars()->GetNumberOfComponents() == 2)
      {
        switch (scalarType)
        {
          vtkTemplateMacro(vtkFixedPointCompositeGOHelperGenerateImageTwoDependentTrilin(
            static_cast<VTK_TT*>(data), threadID, threadCount, mapper, vol));
        }
      }
      // Four components - they must be unsigned char, the first three directly
      // specify color and the fourth specifies opacity (through a lookup
      // table)
      else
      {
        if (scalarType == VTK_UNSIGNED_CHAR)
        {
          vtkFixedPointCompositeGOHelperGenerateImageFourDependentTrilin(
            static_cast<unsigned char*>(data), threadID, threadCount, mapper, vol);
        }
        else
        {
          vtkErrorMacro("Four component dependent must be unsigned char!");
        }
      }
    }
  }
}

// Print method for vtkFixedPointVolumeRayCastCompositeGOHelper
void vtkFixedPointVolumeRayCastCompositeGOHelper::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}
