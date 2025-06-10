// SPDX-FileCopyrightText: Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
// SPDX-License-Identifier: BSD-3-Clause

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

#include "vtkFixedPointVolumeRayCastHelper.h"
#include "vtkRenderingVolumeModule.h" // For export macro

VTK_ABI_NAMESPACE_BEGIN
class vtkFixedPointVolumeRayCastMapper;
class vtkVolume;

class VTKRENDERINGVOLUME_EXPORT vtkFixedPointVolumeRayCastCompositeGOHelper
  : public vtkFixedPointVolumeRayCastHelper
{
public:
  static vtkFixedPointVolumeRayCastCompositeGOHelper* New();
  vtkTypeMacro(vtkFixedPointVolumeRayCastCompositeGOHelper, vtkFixedPointVolumeRayCastHelper);
  void PrintSelf(ostream& os, vtkIndent indent) override;

  void GenerateImage(int threadID, int threadCount, vtkVolume* vol,
    vtkFixedPointVolumeRayCastMapper* mapper) override;

  struct VQStruct1
  {
    double weight = 0.0;
    double threshold = 0.0;
    double transPeriod = 0.0;

    VQStruct1(double w, double t, double tp)
      : weight(w)
      , threshold(t)
      , transPeriod(tp)
    {
    }
  };
  struct VQStruct2
  {
    double weight = 0.0;
    double threshold = 0.0;
    double transPeriod = 0.0;
    bool invert = false;
    double* channelWeight = nullptr;
    bool colorProjection = false;

    VQStruct2(double w, double t, double tp, bool inv, double* cw, bool cp)
      : weight(w)
      , threshold(t)
      , transPeriod(tp)
      , invert(inv)
      , channelWeight(cw)
      , colorProjection(cp)
    {
    }
  };

  enum CompositeMethod
  {
    RegularComposite = 0,
    FeatureDetection = 1,
    ColorProjection = 2
  };

  void setFeatureDetectionWeight(double weight) { m_weight = weight; }
  void setFeatureDetectionThreshold(double threshold) { m_threshold = threshold; }
  void setFeatureDetectionTransitionPeriod(double transPeriod) { m_transPeriod = transPeriod; }
  void setFeatureDetectionEnabled(bool toEnable)
  {
    if (toEnable)
    {
      m_compositeMethod = FeatureDetection;
    }
    else
    {
      if (m_compositeMethod == FeatureDetection)
      {
        m_compositeMethod = RegularComposite;
      }
    }
  }
  void setColorDetectionWeight(double* weight)
  {
    for (int i = 0; i < 3; i++)
    {
      m_channelWeight[i] = weight[i];
    }
  }
  void setColorProjectionEnabled(bool toEnable)
  {
    if (toEnable)
    {
      m_compositeMethod = ColorProjection;
    }
    else
    {
      if (m_compositeMethod == ColorProjection)
      {
        m_compositeMethod = RegularComposite;
      }
    }
  }
  void setCompositeOpacityInverted(bool toEnable) { m_invert = toEnable; }

protected:
  vtkFixedPointVolumeRayCastCompositeGOHelper();
  ~vtkFixedPointVolumeRayCastCompositeGOHelper() override;

private:
  vtkFixedPointVolumeRayCastCompositeGOHelper(
    const vtkFixedPointVolumeRayCastCompositeGOHelper&) = delete;
  void operator=(const vtkFixedPointVolumeRayCastCompositeGOHelper&) = delete;

  CompositeMethod m_compositeMethod;
  bool m_invert;
  double m_weight;
  double m_threshold;
  double m_transPeriod;
  double m_channelWeight[3];
};

VTK_ABI_NAMESPACE_END
#endif
