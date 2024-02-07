/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkFixedPointVolumeRayCastCompositeGOHelper.h
  Language:  C++

  Copyright (c) 1993-2002 Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

/**
 * @class   vtkFixedPointVolumeRayCastCompositeGOHelper
 * @brief   A helper that generates composite images for the volume ray cast mapper
 *
 * This is one of the helper classes for the vtkFixedPointVolumeRayCastMapper.
 * It will generate composite images using an alpha blending operation.
 * This class should not be used directly, it is a helper class for
 * the mapper and has no user-level API.
 *
 * @sa
 * vtkFixedPointVolumeRayCastMapper
*/

#ifndef vtkFixedPointVolumeRayCastCompositeGOHelper_h
#define vtkFixedPointVolumeRayCastCompositeGOHelper_h

#include "vtkRenderingVolumeModule.h" // For export macro
#include "vtkFixedPointVolumeRayCastHelper.h"

class vtkFixedPointVolumeRayCastMapper;
class vtkVolume;

class VTKRENDERINGVOLUME_EXPORT vtkFixedPointVolumeRayCastCompositeGOHelper : public vtkFixedPointVolumeRayCastHelper
{
public:
  static vtkFixedPointVolumeRayCastCompositeGOHelper *New();
  vtkTypeMacro(vtkFixedPointVolumeRayCastCompositeGOHelper,vtkFixedPointVolumeRayCastHelper);
  void PrintSelf( ostream& os, vtkIndent indent ) VTK_OVERRIDE;

  void  GenerateImage( int threadID,
                               int threadCount,
                               vtkVolume *vol,
                               vtkFixedPointVolumeRayCastMapper *mapper) VTK_OVERRIDE;
  struct VQStruct1
  {
      double weight = 0.0;
      double threshold = 0.0;
      double transPeriod = 0.0;
  };
  struct VQStruct2
  {
      double weight = 0.0;
      double threshold = 0.0;
      double transPeriod = 0.0;
      bool invert = false;
      double* channelWeight = nullptr;
      bool colorProjection = false;
  };

  enum CompositeMethod {
      RegularComposite = 0,
      FeatureDetection = 1,
      ColorProjection = 2
  };

  void setFeatureDetectionWeight(double weight) { m_weight = weight; }
  void setFeatureDetectionThreshold(double threshold) { m_threshold = threshold; }
  void setFeatureDetectionTransitionPeriod(double transPeriod) { m_transPeriod = transPeriod; }
  void setFeatureDetectionEnabled(bool toEnable) {
      if (toEnable) {
          m_compositeMethod = FeatureDetection;
      }
      else {
          if (m_compositeMethod == FeatureDetection) {
              m_compositeMethod = RegularComposite;
          }
      }
  }
  void setColorDetectionWeight(double* weight) { for (int i = 0; i < 3; i++) { m_channelWeight[i] = weight[i]; } }
  void setColorProjectionEnabled(bool toEnable) {
      if (toEnable) {
          m_compositeMethod = ColorProjection;
      }
      else {
          if (m_compositeMethod == ColorProjection) {
              m_compositeMethod = RegularComposite;
          }
      }
  }
  void setCompositeOpacityInverted(bool toEnable) { m_invert = toEnable; }

protected:
  vtkFixedPointVolumeRayCastCompositeGOHelper();
  ~vtkFixedPointVolumeRayCastCompositeGOHelper() VTK_OVERRIDE;

private:
  vtkFixedPointVolumeRayCastCompositeGOHelper(const vtkFixedPointVolumeRayCastCompositeGOHelper&) VTK_DELETE_FUNCTION;
  void operator=(const vtkFixedPointVolumeRayCastCompositeGOHelper&) VTK_DELETE_FUNCTION;

  CompositeMethod m_compositeMethod;
  bool m_invert;
  double m_weight;
  double m_threshold;
  double m_transPeriod;
  double m_channelWeight[3];
};


#endif


