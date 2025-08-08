// SPDX-FileCopyrightText: Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
// SPDX-License-Identifier: BSD-3-Clause
/**
 * @class   vtkFixedPointVolumeRayCastMIPHelper
 * @brief   A helper that generates MIP images for the volume ray cast mapper
 *
 * This is one of the helper classes for the vtkFixedPointVolumeRayCastMapper.
 * It will generate maximum intensity images.
 * This class should not be used directly, it is a helper class for
 * the mapper and has no user-level API.
 *
 * @sa
 * vtkFixedPointVolumeRayCastMapper
 */

#ifndef vtkFixedPointVolumeRayCastMIPHelper_h
#define vtkFixedPointVolumeRayCastMIPHelper_h

#include "vtkFixedPointVolumeRayCastHelper.h"
#include "vtkRenderingVolumeModule.h" // For export macro

#include <array>
#include "vtkMath.h"
#include "vtkVector.h"

VTK_ABI_NAMESPACE_BEGIN
class vtkFixedPointVolumeRayCastMapper;
class vtkVolume;

class VTKRENDERINGVOLUME_EXPORT vtkFixedPointVolumeRayCastMIPHelper
  : public vtkFixedPointVolumeRayCastHelper
{
public:
  static vtkFixedPointVolumeRayCastMIPHelper* New();
  vtkTypeMacro(vtkFixedPointVolumeRayCastMIPHelper, vtkFixedPointVolumeRayCastHelper);
  void PrintSelf(ostream& os, vtkIndent indent) override;

  void GenerateImage(int threadID, int threadCount, vtkVolume* vol,
    vtkFixedPointVolumeRayCastMapper* mapper) override;

protected:
  vtkFixedPointVolumeRayCastMIPHelper();
  ~vtkFixedPointVolumeRayCastMIPHelper() override;

private:
  vtkFixedPointVolumeRayCastMIPHelper(const vtkFixedPointVolumeRayCastMIPHelper&) = delete;
  void operator=(const vtkFixedPointVolumeRayCastMIPHelper&) = delete;
};

class VTKRENDERINGVOLUME_EXPORT vqFixedPointVolumeRayCastAIPHelper
  : public vtkFixedPointVolumeRayCastHelper
{
public:
  static vqFixedPointVolumeRayCastAIPHelper* New();
  vtkTypeMacro(vqFixedPointVolumeRayCastAIPHelper, vtkFixedPointVolumeRayCastHelper);
  void PrintSelf(ostream& os, vtkIndent indent) override;

  void GenerateImage(int threadID, int threadCount, vtkVolume* vol,
    vtkFixedPointVolumeRayCastMapper* mapper) override;

protected:
  vqFixedPointVolumeRayCastAIPHelper();
  ~vqFixedPointVolumeRayCastAIPHelper() override;

private:
  vqFixedPointVolumeRayCastAIPHelper(const vqFixedPointVolumeRayCastAIPHelper&) = delete;
  void operator=(const vqFixedPointVolumeRayCastAIPHelper&) = delete;
};

namespace vqRender
{
template <typename T, size_t size>
class vector
{
public:
  vector()
  {
    for (int i = 0; i < size; i++)
    {
      m_data[i] = 0;
    }
  }
  vector(T s)
  {
    for (int i = 0; i < size; i++)
    {
      m_data[i] = s;
    }
  }

  template <typename... Args>
  vector(Args const&... args)
  {
    const size_t packSize = sizeof...(args);
    static_assert(packSize == size, "pack size must equal to vector size");
    T x[] = { args... };
    for (int i = 0; i < packSize; i++)
    {
      m_data[i] = x[i];
    }
  }
  T operator[](const size_t& i) const { return m_data[i]; }

  T& operator[](const size_t& i) { return m_data[i]; }

  vector operator+(const vector& v) const
  {
    vector result;
    for (int i = 0; i < size; i++)
    {
      result[i] = this->m_data[i] + v[i];
    }
    return result;
  }

  vector operator-(const vector& v) const
  {
    vector result;
    for (int i = 0; i < size; i++)
    {
      result[i] = this->m_data[i] - v[i];
    }
    return result;
  }

  vector operator*(const vector& v) const
  {
    vector result;
    for (int i = 0; i < size; i++)
    {
      result[i] = this->m_data[i] * v[i];
    }
    return result;
  }

  T min()
  {
    T min = m_data[0];
    for (int i = 1; i < size; i++)
    {
      if (m_data[i] < min)
      {
        min = m_data[i];
      }
    }
  }

  void clamp(const T& min, const T& max)
  {
    for (int i = 0; i < size; i++)
    {
      this->m_data[i] = max(min, min(m_data[i], max));
    }
  }

  float length() const
  {
    return sqrt(m_data[0] * m_data[0] + m_data[1] * m_data[1] + m_data[2] * m_data[2]);
  }

  void scale(T f)
  {
    for (int i = 0; i < size; i++)
    {
      this->m_data[i] = this->m_data[i] * f;
    }
  }

  float dot(const vector& v) const
  {
    return (m_data[0] * v[0] + m_data[1] * v[1] + m_data[2] * v[2]);
  }

  vector cross(const vector& v) const
  {
    return vector((m_data[1] * v[2] - m_data[2] * v[1]), (m_data[2] * v[0] - m_data[0] * v[2]),
      (m_data[0] * v[1] - m_data[1] * v[0]));
  }

  float theta(const vector& v) const
  {
    float costheta = this->dot(v) / (this->length() * v.length());
    float theta = acos(costheta) * 180 / 3.1415;
    if (theta > 180)
    {
      theta -= 180;
    }
    return theta;
  }

  void normalize()
  {
    float len = length();
    if (len <= 0.0)
    {
      for (int i = 0; i < size; i++)
      {
        this->m_data[i] = 0;
      }
    }
    else
    {
    }
    for (int i = 0; i < size; i++)
    {
      this->m_data[i] = this->m_data[i] / len;
    }
  }

  void reverse()
  {
    for (int i = 0; i < size; i++)
    {
      this->m_data[i] = -this->m_data[i];
    }
  }

private:
  T m_data[size];
};

template <typename T, size_t size>
float abslength(const T (&arr)[size])
{
  float len = 0;
  for (int i = 0; i < size; i++)
  {
    len += arr[i] * arr[i];
  }
  return sqrt(len);
}
}

typedef vqRender::vector<float, 3> vec3;
typedef vqRender::vector<float, 4> vec4;
typedef vqRender::vector<int, 3> vec3i;
typedef vqRender::vector<unsigned int, 3> vec3ui;
typedef vqRender::vector<unsigned short, 4> vec4us;

const float INV_PI_F = 0.31830988618379067154;
const float INV_TWO_PI_F = 0.15915494309189533577;
const float PI_F = 3.141592654;

struct shaderInfo
{
  vec3 reflectRay;
  vec4 emissionColor;
  vec4 diffuseColor;
  vec4 specularColor;
  float ior;
  float glossiness;
};

struct lightInfo
{
  vec3 surfaceN;
  vec3 surfaceP;
  vec3 pos;
  int type; // 0 for point and 1 for plane
};

struct scatterEvent
{
  vec3 minT;
  vec3 nPos;
  unsigned int dPos[3];
  vec3 norm;
  unsigned short intensity[3];
  float normWeight;
  vec3 inRayDir;
  vec3 refRay;
  float lightEnergy;
  int event; // 0 for volume, 1 for light, 2 for other object
  int valid;
};

class VTKRENDERINGVOLUME_EXPORT vqFixedPointVolumeRayCastISOHelper
  : public vtkFixedPointVolumeRayCastHelper
{
public:
  static vqFixedPointVolumeRayCastISOHelper* New();
  vtkTypeMacro(vqFixedPointVolumeRayCastISOHelper, vtkFixedPointVolumeRayCastHelper);
  void PrintSelf(ostream& os, vtkIndent indent) override;

  void GenerateImage(int threadID, int threadCount, vtkVolume* vol,
    vtkFixedPointVolumeRayCastMapper* mapper) override;

protected:
  vqFixedPointVolumeRayCastISOHelper();
  ~vqFixedPointVolumeRayCastISOHelper() override;

private:
  vqFixedPointVolumeRayCastISOHelper(const vqFixedPointVolumeRayCastISOHelper&) = delete;
  void operator=(const vqFixedPointVolumeRayCastISOHelper&) = delete;
};

VTK_ABI_NAMESPACE_END
#endif
