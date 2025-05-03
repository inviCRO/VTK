#ifndef VTKOPENGLGUPVOLUMERENDERINGHELPER_H
#define VTKOPENGLGUPVOLUMERENDERINGHELPER_H

/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkOpenGLGPUVolumeRayCastMapper.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

// Include compiled shader code
#include <vtkBoundingBox.h>
#include <vtkCamera.h>
#include <vtkCellArray.h>
#include <vtkClipConvexPolyData.h>
#include <vtkColorTransferFunction.h>
#include <vtkCommand.h>
#include <vtkDataArray.h>
#include <vtkDensifyPolyData.h>
#include <vtkFloatArray.h>
#include <vtk_glew.h>
#include <vtkImageData.h>
#include <vtkLight.h>
#include <vtkLightCollection.h>
#include <vtkMath.h>
#include <vtkMatrix4x4.h>
#include <vtkNew.h>
#include <vtkObjectFactory.h>
#include <vtkOpenGLError.h>
#include <vtkOpenGLShaderCache.h>
#include <vtkOpenGLRenderWindow.h>
#include <vtkPerlinNoise.h>
#include <vtkPlaneCollection.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkShader.h>
#include <vtkShaderProgram.h>
#include <vtkSmartPointer.h>
#include <vtkTessellatedBoxSource.h>
#include <vtkTextureObject.h>
#include <vtkTimerLog.h>
#include <vtkTransform.h>
#include <vtkUnsignedCharArray.h>
#include <vtkUnsignedIntArray.h>
#include <vtkVolumeMask.h>
#include <vtkVolumeProperty.h>
#include <vtkWeakPointer.h>
#include <vtkSphericalDirectionEncoder.h>


#define VOLUME_GRADIENT


#ifdef VOLUME_GRADIENT
template <class T>
void vtkFixedPointVolumeRayCastMapperComputeGradients(T *dataPtr,
    int dim[3],
    double spacing[3],
    int components,
    int independent,
    double scalarRange[4][2],
    unsigned short **gradientNormal,
    unsigned char  **gradientMagnitude,
    vtkDirectionEncoder *directionEncoder)
{
    int                 x, y, z, c;
    vtkIdType           yinc, zinc;
    int                 x_start, x_limit;
    int                 y_start, y_limit;
    int                 z_start, z_limit;
    T                   *dptr, *cdptr;
    float               n[3], t;
    float               gvalue = 0;
    int                 xlow, xhigh;
    double              aspect[3];
    vtkIdType           xstep, ystep, zstep;
    float               scale[4];
    unsigned short      *dirPtr, *cdirPtr;
    unsigned char       *magPtr, *cmagPtr;

    double avgSpacing = (spacing[0] + spacing[1] + spacing[2]) / 3.0;

    // adjust the aspect
    aspect[0] = spacing[0] * 2.0 / avgSpacing;
    aspect[1] = spacing[1] * 2.0 / avgSpacing;
    aspect[2] = spacing[2] * 2.0 / avgSpacing;

    // compute the increments
    yinc = static_cast<vtkIdType>(dim[0]);
    zinc = yinc*static_cast<vtkIdType>(dim[1]);

    // Compute steps through the volume in x, y, and z
    xstep = components;
    ystep = components*yinc;
    zstep = components*zinc;

    if (!independent)
    {
        if (scalarRange[components - 1][1] - scalarRange[components - 1][0])
        {
            scale[0] = 255.0 / (0.25*(scalarRange[components - 1][1] - scalarRange[components - 1][0]));
        }
        else
        {
            scale[0] = 0.0;
        }
    }
    else
    {
        for (c = 0; c < components; c++)
        {
            if (scalarRange[c][1] - scalarRange[c][0])
            {
                scale[c] = 255.0 / (0.25*(scalarRange[c][1] - scalarRange[c][0]));
            }
            else
            {
                scale[c] = 1.0;
            }
        }
    }

    x_start = 0;
    x_limit = dim[0];
    y_start = 0;
    y_limit = dim[1];
    z_start = 0;
    z_limit = dim[2];

    int increment = (independent) ? (components) : (1);

    float tolerance[4];
    for (c = 0; c < components; c++)
    {
        tolerance[c] = .00001 * (scalarRange[c][1] - scalarRange[c][0]);
    }

    // Loop through all the data and compute the encoded normal and
    // gradient magnitude for each scalar location
    for (z = z_start; z < z_limit; z++)
    {
        unsigned short *gradientDirPtr = gradientNormal[z];
        unsigned char *gradientMagPtr = gradientMagnitude[z];

        for (y = y_start; y < y_limit; y++)
        {
            xlow = x_start;
            xhigh = x_limit;

            dptr = dataPtr + components*(z * zinc + y * yinc + xlow);

            dirPtr = gradientDirPtr + (y * yinc + xlow)*increment;
            magPtr = gradientMagPtr + (y * yinc + xlow)*increment;

            for (x = xlow; x < xhigh; x++)
            {
                for (c = 0; (independent && c < components) || c == 0; c++)
                {
                    cdptr = dptr + ((independent) ? (c) : (components - 1));
                    cdirPtr = dirPtr + ((independent) ? (c) : (0));
                    cmagPtr = magPtr + ((independent) ? (c) : (0));

                    // Allow up to 3 tries to find the gadient - looking out at a distance of
                    // 1, 2, and 3 units.
                    int foundGradient = 0;
                    for (int d = 1; d <= 3 && !foundGradient; d++)
                    {
                        // Use a central difference method if possible,
                        // otherwise use a forward or backward difference if
                        // we are on the edge
                        // Compute the X component
                        if (x < d)
                        {
                            n[0] = 2.0*((float)*(cdptr)-(float)*(cdptr + d*xstep));
                        }
                        else if (x >= dim[0] - d)
                        {
                            n[0] = 2.0*((float)*(cdptr - d*xstep) - (float)*(cdptr));
                        }
                        else
                        {
                            n[0] = (float)*(cdptr - d*xstep) - (float)*(cdptr + d*xstep);
                        }

                        // Compute the Y component
                        if (y < d)
                        {
                            n[1] = 2.0*((float)*(cdptr)-(float)*(cdptr + d*ystep));
                        }
                        else if (y >= dim[1] - d)
                        {
                            n[1] = 2.0*((float)*(cdptr - d*ystep) - (float)*(cdptr));
                        }
                        else
                        {
                            n[1] = (float)*(cdptr - d*ystep) - (float)*(cdptr + d*ystep);
                        }

                        // Compute the Z component
                        if (z < d)
                        {
                            n[2] = 2.0*((float)*(cdptr)-(float)*(cdptr + d*zstep));
                        }
                        else if (z >= dim[2] - d)
                        {
                            n[2] = 2.0*((float)*(cdptr - d*zstep) - (float)*(cdptr));
                        }
                        else
                        {
                            n[2] = (float)*(cdptr - d*zstep) - (float)*(cdptr + d*zstep);
                        }

                        // Take care of the aspect ratio of the data
                        // Scaling in the vtkVolume is isotropic, so this is the
                        // only place we have to worry about non-isotropic scaling.
                        n[0] /= d*aspect[0];
                        n[1] /= d*aspect[1];
                        n[2] /= d*aspect[2];

                        // Compute the gradient magnitude
                        t = sqrt((double)(n[0] * n[0] +
                            n[1] * n[1] +
                            n[2] * n[2]));


                        // Encode this into an 8 bit value
                        gvalue = t * scale[c];

                        if (d > 1)
                        {
                            gvalue = 0;
                        }

                        gvalue = (gvalue<0.0) ? (0.0) : (gvalue);
                        gvalue = (gvalue>255.0) ? (255.0) : (gvalue);

                        // Normalize the gradient direction
                        if (t > tolerance[c])
                        {
                            n[0] /= t;
                            n[1] /= t;
                            n[2] /= t;
                            foundGradient = 1;
                        }
                        else
                        {
                            n[0] = n[1] = n[2] = 0.0;
                        }
                    }


                    *cmagPtr = static_cast<unsigned char>(gvalue + 0.5);
                    *cdirPtr = directionEncoder->GetEncodedDirection(n);
                }

                dptr += components;
                dirPtr += increment;
                magPtr += increment;
            }
        }
        if (z % 8 == 7)
        {
            double args[1];
            args[0] =
                static_cast<float>(z - z_start) /
                static_cast<float>(z_limit - z_start - 1);
        }
    }
}

#endif


#endif //VTKOPENGLGUPVOLUMERENDERINGHELPER_H
