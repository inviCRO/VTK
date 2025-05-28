//VTK::System::Dec

/*=========================================================================

  Program:   Visualization Toolkit
  Module:    raycasterfs.glsl

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

//////////////////////////////////////////////////////////////////////////////
///
/// Inputs
///
//////////////////////////////////////////////////////////////////////////////

/// 3D texture coordinates form vertex shader
varying vec3 ip_textureCoords;
varying vec3 ip_vertexPos;

//////////////////////////////////////////////////////////////////////////////
///
/// Outputs
///
//////////////////////////////////////////////////////////////////////////////

vec4 g_fragColor = vec4(0.0);

//////////////////////////////////////////////////////////////////////////////
///
/// Uniforms, attributes, and globals
///
//////////////////////////////////////////////////////////////////////////////
vec3 g_dataPos;
vec3 g_dirStep;
vec4 g_srcColor;
vec4 g_eyePosObj;
bool g_exit;
bool g_skip;
float g_currentT;
float g_terminatePointMax;

uniform vec4 in_volume_scale;
uniform vec4 in_volume_bias;

//VTK::Output::Dec

      
// Volume dataset      
uniform sampler3D in_volume;      
uniform int in_noOfComponents;      
uniform int in_independentComponents;      
      
uniform sampler2D in_noiseSampler;      
#ifndef GL_ES      
uniform sampler2D in_depthSampler;      
#endif      
      
// Camera position      
uniform vec3 in_cameraPos;      
      
// view and model matrices      
uniform mat4 in_volumeMatrix;      
uniform mat4 in_inverseVolumeMatrix;      
uniform mat4 in_projectionMatrix;      
uniform mat4 in_inverseProjectionMatrix;      
uniform mat4 in_modelViewMatrix;      
uniform mat4 in_inverseModelViewMatrix;      
uniform mat4 in_textureDatasetMatrix;      
uniform mat4 in_inverseTextureDatasetMatrix;      
varying mat4 ip_inverseTextureDataAdjusted;      
uniform vec3 in_texMin;      
uniform vec3 in_texMax;      
uniform mat4 in_textureToEye;      
      
// Ray step size      
uniform vec3 in_cellStep;      
uniform vec2 in_scalarsRange[4];      
uniform vec3 in_cellSpacing;      
      
// Sample distance      
uniform float in_sampleDistance;      
      
// Scales      
uniform vec3 in_cellScale;      
uniform vec2 in_windowLowerLeftCorner;      
uniform vec2 in_inverseOriginalWindowSize;      
uniform vec2 in_inverseWindowSize;      
uniform vec3 in_textureExtentsMax;      
uniform vec3 in_textureExtentsMin;      
      
// Material and lighting      
uniform vec3 in_diffuse[4];      
uniform vec3 in_ambient[4];      
uniform vec3 in_specular[4];      
uniform float in_shininess[4];      
      
// Others      
uniform bool in_cellFlag;      
uniform bool in_useJittering;      
vec3 g_rayJitter = vec3(0.0);      
uniform bool in_clampDepthToBackface;      
      
uniform vec2 in_averageIPRange;

//VQ ADDED
uniform mat4 in_originMatrix;
uniform mat4 in_InverseOriginMatrix;
uniform mat4 in_volumeUserMatrix;
uniform mat4 in_volumeInvUserMatrix;
uniform mat4 in_flipMatrix;
uniform mat4 in_textureOriginMatrix;
mat4 invTextureOriginMatrix = inverse(in_textureOriginMatrix);

  vec4 l_bb_min; 
  vec4 l_bb_max; 
  uniform float fuseCoef; 
  uniform int m_lightonly;
  uniform int m_inverted;
  uniform float mc_weight;
  uniform float mc_threshold;
  uniform float mc_transPeriod;
  uniform vec3 mc_channelWeight;
  uniform sampler3D in_gradient;        
uniform bool in_twoSidedLighting;        
vec3 g_xvec;        
vec3 g_yvec;        
vec3 g_zvec;        
vec3 g_aspect;        
vec3 g_cellSpacing;        
float g_avgSpacing;        
uniform vec3 in_lightAmbientColor[1];        
uniform vec3 in_lightDiffuseColor[1];        
uniform vec3 in_lightSpecularColor[1];        
vec4 g_lightPosObj;        
vec3 g_ldir;        
vec3 g_vdir;        
vec3 g_h;        
uniform vec4 in_componentWeight;

      
 const float g_opacityThreshold = 1.0 - 1.0 / 255.0;


  uniform float in_croppingPlanes[6];

      
 int clippingPlanesSize;      
 vec3 objRayDir;      
 mat4 textureToObjMat;

        
 bool l_firstValue;        
 vec4 l_sampledValue;

//VTK::BinaryMask::Dec

//VTK::CompositeMask::Dec


 uniform sampler2D in_opacityTransferFunc;
 uniform sampler2D in_opacityTransferFunc1;
 uniform sampler2D in_opacityTransferFunc2;          
float computeOpacity(vec4 scalar, int component)          
  {            
  if (component == 0)            
    {            
    return texture2D(in_opacityTransferFunc,vec2(scalar[0],0)).r;            
    }            
  if (component == 1)            
    {            
    return texture2D(in_opacityTransferFunc1,vec2(scalar[1],0)).r;            
    }            
  if (component == 2)            
    {            
    return texture2D(in_opacityTransferFunc2,vec2(scalar[2],0)).r;            
    }
  }


 uniform sampler2D in_gradientTransferFunc;
 uniform sampler2D in_gradientTransferFunc1;
 uniform sampler2D in_gradientTransferFunc2;        
float computeGradientOpacity(vec4 grad, int component)        
  {          
  if (component == 0)          
    {          
    return texture2D(in_gradientTransferFunc, vec2(grad.w, 0.0)).r;          
    }          
  if (component == 1)          
    {          
    return texture2D(in_gradientTransferFunc1, vec2(grad.w, 0.0)).r;          
    }          
  if (component == 2)          
    {          
    return texture2D(in_gradientTransferFunc2, vec2(grad.w, 0.0)).r;          
    }        
  }        
// c is short for component        
vec4 computeGradient(int c)        
  {        
  // Approximate Nabla(F) derivatives with central differences.        
  vec3 g1; // F_front        
  vec3 g2; // F_back        
  g1.x = texture3D(in_volume, vec3(g_dataPos + g_xvec))[c];        
  g1.y = texture3D(in_volume, vec3(g_dataPos + g_yvec))[c];        
  g1.z = texture3D(in_volume, vec3(g_dataPos + g_zvec))[c];        
  g2.x = texture3D(in_volume, vec3(g_dataPos - g_xvec))[c];        
  g2.y = texture3D(in_volume, vec3(g_dataPos - g_yvec))[c];        
  g2.z = texture3D(in_volume, vec3(g_dataPos - g_zvec))[c];        
        
  // Apply scale and bias to the fetched values.        
  g1 = g1 * in_volume_scale[c] + in_volume_bias[c];        
  g2 = g2 * in_volume_scale[c] + in_volume_bias[c];        
        
  // Central differences: (F_front - F_back) / 2h        
  // This version of computeGradient() is only used for lighting        
  // calculations (only direction matters), hence the difference is        
  // not scaled by 2h and a dummy gradient mag is returned (-1.).        
  return vec4((g1 - g2), -1.0);        
  }

      
vec4 computeLighting(vec4 color, int component)      
  {      
  vec4 finalColor = vec4(0.0);        
  // Compute gradient function only once        
  vec4 gradient = computeGradient(component);          
  vec3 diffuse = vec3(0.0);          
  vec3 specular = vec3(0.0);          
  vec3 normal = gradient.xyz / in_cellSpacing;          
  float normalLength = length(normal);          
  if (normalLength > 0.0)          
    {          
    normal = normalize(normal);          
    }          
  else          
    {          
    normal = vec3(0.0, 0.0, 0.0);          
    }          
   float nDotL = dot(normal, g_ldir);          
   float nDotH = dot(normal, g_h);          
   if (nDotL < 0.0 && in_twoSidedLighting)          
     {          
     nDotL = -nDotL;          
     }          
   if (nDotH < 0.0 && in_twoSidedLighting)          
     {          
     nDotH = -nDotH;          
     }          
   if (nDotL > 0.0)          
     {          
     diffuse = nDotL * in_diffuse[component] *          
               in_lightDiffuseColor[0] * color.rgb;          
     }          
    specular = pow(nDotH, in_shininess[component]) *          
                 in_specular[component] *          
                 in_lightSpecularColor[0];          
  // For the headlight, ignore the light's ambient color          
  // for now as it is causing the old mapper tests to fail          
  finalColor.xyz = in_ambient[component] * color.rgb +          
                   diffuse + specular;      
  finalColor.a = color.a;      
  return finalColor;      
  }


 uniform sampler2D in_colorTransferFunc;
 uniform sampler2D in_colorTransferFunc1;
 uniform sampler2D in_colorTransferFunc2;          
vec4 computeColor(vec4 scalar, float opacity, int component)          
  {            
  if (component == 0)            
    {            
    return computeLighting(vec4(texture2D(            
      in_colorTransferFunc, vec2(            
      scalar[0],0.0)).xyz,            
      opacity),0);            
    }            
  if (component == 1)            
    {            
    return computeLighting(vec4(texture2D(            
      in_colorTransferFunc1, vec2(            
      scalar[1],0.0)).xyz,            
      opacity),1);            
    }            
  if (component == 2)            
    {            
    return computeLighting(vec4(texture2D(            
      in_colorTransferFunc2, vec2(            
      scalar[2],0.0)).xyz,            
      opacity),2);            
    }
  }

        
vec3 computeRayDirection()        
  {        
  return normalize(ip_vertexPos.xyz - g_eyePosObj.xyz);        
  }

//VTK::Picking::Dec

//VTK::RenderToImage::Dec

//VTK::DepthPeeling::Dec

/// We support only 8 clipping planes for now
/// The first value is the size of the data array for clipping
/// planes (origin, normal)
uniform float in_clippingPlanes[49];
uniform float in_scale;
uniform float in_bias;

//////////////////////////////////////////////////////////////////////////////
///
/// Helper functions
///
//////////////////////////////////////////////////////////////////////////////

/**
 * Transform window coordinate to NDC.
 */
vec4 WindowToNDC(const float xCoord, const float yCoord, const float zCoord)
{
  vec4 NDCCoord = vec4(0.0, 0.0, 0.0, 1.0);

  NDCCoord.x = (xCoord - in_windowLowerLeftCorner.x) * 2.0 *
    in_inverseWindowSize.x - 1.0;
  NDCCoord.y = (yCoord - in_windowLowerLeftCorner.y) * 2.0 *
    in_inverseWindowSize.y - 1.0;
  NDCCoord.z = (2.0 * zCoord - (gl_DepthRange.near + gl_DepthRange.far)) /
    gl_DepthRange.diff;

  return NDCCoord;
}

/**
 * Transform NDC coordinate to window coordinates.
 */
vec4 NDCToWindow(const float xNDC, const float yNDC, const float zNDC)
{
  vec4 WinCoord = vec4(0.0, 0.0, 0.0, 1.0);

  WinCoord.x = (xNDC + 1.f) / (2.f * in_inverseWindowSize.x) +
    in_windowLowerLeftCorner.x;
  WinCoord.y = (yNDC + 1.f) / (2.f * in_inverseWindowSize.y) +
    in_windowLowerLeftCorner.y;
  WinCoord.z = (zNDC * gl_DepthRange.diff +
    (gl_DepthRange.near + gl_DepthRange.far)) / 2.f;

  return WinCoord;
}

//////////////////////////////////////////////////////////////////////////////
///
/// Ray-casting
///
//////////////////////////////////////////////////////////////////////////////

/**
 * Global initialization. This method should only be called once per shader
 * invocation regardless of whether castRay() is called several times (e.g.
 * vtkDualDepthPeelingPass). Any castRay() specific initialization should be
 * placed within that function.
 */
void initializeRayCast()
{
    /// Initialize g_fragColor (output) to 0
    g_fragColor = vec4(0.0);
    g_dirStep = vec3(0.0);
    g_srcColor = vec4(0.0);
    g_exit = false;

        
  bool l_adjustTextureExtents =  !in_cellFlag;        
  // Get the 3D texture coordinates for lookup into the in_volume dataset        
  //g_dataPos = ip_textureCoords.xyz;
 vec4 texremapped = vec4(ip_textureCoords, 1.0); 
 texremapped = invTextureOriginMatrix * in_flipMatrix * in_textureOriginMatrix * texremapped; 
 g_dataPos = texremapped.xyz;         
        
  // Eye position in dataset space        
  //g_eyePosObj = (in_inverseVolumeMatrix * vec4(in_cameraPos, 1.0));
  g_eyePosObj = ( in_originMatrix  * in_volumeInvUserMatrix * in_InverseOriginMatrix* vec4(in_cameraPos, 1.0));        
  if (g_eyePosObj.w != 0.0)        
    {        
    g_eyePosObj.x /= g_eyePosObj.w;        
    g_eyePosObj.y /= g_eyePosObj.w;        
    g_eyePosObj.z /= g_eyePosObj.w;        
    g_eyePosObj.w = 1.0;        
    }        
        
  // Getting the ray marching direction (in dataset space);        
  vec3 rayDir = computeRayDirection();        
        
  // Multiply the raymarching direction with the step size to get the        
  // sub-step size we need to take at each raymarching step        
  //g_dirStep = (ip_inverseTextureDataAdjusted *        
  //            vec4(rayDir, 0.0)).xyz * in_sampleDistance;
  g_dirStep = (in_inverseTextureDatasetMatrix * invTextureOriginMatrix * in_flipMatrix * in_textureOriginMatrix * vec4(rayDir, 0.0)).xyz * in_sampleDistance;        
        
  // 2D Texture fragment coordinates [0,1] from fragment coordinates.        
  // The frame buffer texture has the size of the plain buffer but         
  // we use a fraction of it. The texture coordinate is less than 1 if        
  // the reduction factor is less than 1.        
  // Device coordinates are between -1 and 1. We need texture        
  // coordinates between 0 and 1. The in_noiseSampler and in_depthSampler        
  // buffers have the original size buffer.        
  vec2 fragTexCoord = (gl_FragCoord.xy - in_windowLowerLeftCorner) *        
                      in_inverseWindowSize;        
        
  if (in_useJittering)        
  {        
    float jitterValue = texture2D(in_noiseSampler, fragTexCoord).x;        
    g_rayJitter = g_dirStep * jitterValue;        
    g_dataPos += g_rayJitter;        
  }        
  else        
  {        
    g_dataPos += g_dirStep;        
  }        
        
  // Flag to deternmine if voxel should be considered for the rendering        
  g_skip = false;          
  // Light position in dataset space          
  g_lightPosObj = (in_inverseVolumeMatrix *          
                      vec4(in_cameraPos, 1.0));          
  if (g_lightPosObj.w != 0.0)          
    {          
    g_lightPosObj.x /= g_lightPosObj.w;          
    g_lightPosObj.y /= g_lightPosObj.w;          
    g_lightPosObj.z /= g_lightPosObj.w;          
    g_lightPosObj.w = 1.0;          
    }          
  g_ldir = normalize(g_lightPosObj.xyz - ip_vertexPos);          
  g_vdir = normalize(g_eyePosObj.xyz - ip_vertexPos);          
  g_h = normalize(g_ldir + g_vdir);        
  g_xvec = vec3(in_cellStep[0], 0.0, 0.0);        
  g_yvec = vec3(0.0, in_cellStep[1], 0.0);        
  g_zvec = vec3(0.0, 0.0, in_cellStep[2]);        
  g_cellSpacing = vec3(in_cellSpacing[0],        
                       in_cellSpacing[1],        
                       in_cellSpacing[2]);        
  g_avgSpacing = (g_cellSpacing[0] +        
                  g_cellSpacing[1] +        
                  g_cellSpacing[2])/3.0;        
  // Adjust the aspect        
  g_aspect.x = g_cellSpacing[0] * 2.0 / g_avgSpacing;        
  g_aspect.y = g_cellSpacing[1] * 2.0 / g_avgSpacing;        
  g_aspect.z = g_cellSpacing[2] * 2.0 / g_avgSpacing;

        
  // Flag to indicate if the raymarch loop should terminate       
  bool stop = false;      
      
  g_terminatePointMax = 0.0;      
      
#ifdef GL_ES      
  vec4 l_depthValue = vec4(1.0,1.0,1.0,1.0);      
#else      
  vec4 l_depthValue = texture2D(in_depthSampler, fragTexCoord);      
#endif      
  // Depth test      
  if(gl_FragCoord.z >= l_depthValue.x)      
    {      
    discard;      
    }      
      
  // color buffer or max scalar buffer have a reduced size.      
  fragTexCoord = (gl_FragCoord.xy - in_windowLowerLeftCorner) *      
                 in_inverseOriginalWindowSize;      
      
  // Compute max number of iterations it will take before we hit      
  // the termination point      
      
  // Abscissa of the point on the depth buffer along the ray.      
  // point in texture coordinates      
  vec4 terminatePoint = WindowToNDC(gl_FragCoord.x, gl_FragCoord.y, l_depthValue.x);      
      
  // From normalized device coordinates to eye coordinates.      
  // in_projectionMatrix is inversed because of way VT      
  // From eye coordinates to texture coordinates      
  terminatePoint = ip_inverseTextureDataAdjusted *
  in_originMatrix * 
  in_volumeInvUserMatrix * 
  in_InverseOriginMatrix *       
                   //in_inverseVolumeMatrix *      
                   in_inverseModelViewMatrix *      
                   in_inverseProjectionMatrix *      
                   terminatePoint;      
  terminatePoint /= terminatePoint.w;      

  g_terminatePointMax = length(terminatePoint.xyz - ip_textureCoords.xyz) /      
  //g_terminatePointMax = length(terminatePoint.xyz - g_dataPos.xyz) /      
                        length(g_dirStep);      
  g_currentT = 0.0;

  
  l_bb_min = vec4(in_croppingPlanes[0], in_croppingPlanes[2], in_croppingPlanes[4], 1.0);
  l_bb_min = in_inverseTextureDatasetMatrix * l_bb_min;
  if (l_bb_min.w != 0.0)
  {
         l_bb_min /= l_bb_min.w;
  }
  l_bb_max = vec4(in_croppingPlanes[1], in_croppingPlanes[3], in_croppingPlanes[5], 1.0);
  l_bb_max = in_inverseTextureDatasetMatrix * l_bb_max;
  if (l_bb_max.w != 0.0)
  {
         l_bb_max /= l_bb_max.w;
  }

//VTK::Clipping::Init

//VTK::RenderToImage::Init

//VTK::DepthPass::Init
}

/**
 * March along the ray direction sampling the volume texture.  This function
 * takes a start and end point as arguments but it is up to the specific render
 * pass implementation to use these values (e.g. vtkDualDepthPeelingPass). The
 * mapper does not use these values by default, instead it uses the number of
 * steps defined by g_terminatePointMax.
 */
vec4 castRay(const float zStart, const float zEnd)
{
  //VTK::DepthPeeling::Ray::Init

  //VTK::DepthPeeling::Ray::PathCheck

          
  l_firstValue = true;        
  l_sampledValue = vec4(0.0);        
  float remainOpacity = 1.0;                    
  l_sampledValue.w = 0.5;            

  /// For all samples along the ray
  while (!g_exit)
  {
          
    g_skip = false;

    
  g_skip = dot(sign(g_dataPos - l_bb_min.xyz), sign(l_bb_max.xyz - g_dataPos)) < 3.0;

    //VTK::Clipping::Impl

    //VTK::BinaryMask::Impl

    //VTK::CompositeMask::Impl

          
    if (!g_skip)      
      {      
      vec4 scalar = texture3D(in_volume, g_dataPos);        
      scalar = scalar*in_volume_scale + in_volume_bias;
          vec4 gradient = texture3D(in_gradient, g_dataPos); 
          vec4 gradientOpacity = vec4(0.0);
          vec4 scalarOpacity = vec4(0.0);
          if (l_firstValue)
          {
              gradient = vec4(0.0);
          }
          float totalAlpha = 0.0;
          vec4 tmpAlpha = vec4(0.0);
          vec4 tmpGradient = vec4(0.0);
          for(int i = 0; i < in_noOfComponents; i++){
              scalarOpacity[i] = computeOpacity(scalar,i);
              tmpGradient = vec4(gradient[i]);
              gradientOpacity[i] = computeGradientOpacity(tmpGradient,i);
          }
          for (int i = 0; i < in_noOfComponents; ++i)
          {
               tmpAlpha[i] = scalarOpacity[i] * gradientOpacity[i];
               totalAlpha += scalarOpacity[i] * gradientOpacity[i];
          }
          if(totalAlpha > 0 && !g_skip) { 
              for (int i = 0; i < in_noOfComponents; ++i)
              {
                  vec4 tmpcolor = computeColor(scalar, 1.0, i);
                  l_sampledValue[0] += tmpcolor[0] * tmpAlpha[i];
                  l_sampledValue[1] += tmpcolor[1] * tmpAlpha[i];
                  l_sampledValue[2] += tmpcolor[2] * tmpAlpha[i];
                  l_sampledValue[3] += tmpAlpha[i] * tmpAlpha[i] / totalAlpha;
              }
          }
          remainOpacity *= 1 - l_sampledValue[3];
          if (remainOpacity < 0.01)
          {
              break;
          }
          if (l_firstValue)
          {
              l_firstValue = false;
          }
                        
      }

    //VTK::RenderToImage::Impl

    //VTK::DepthPass::Impl

    /// Advance ray
    g_dataPos += g_dirStep;

          
    if(any(greaterThan(g_dataPos, in_texMax)) ||      
      any(lessThan(g_dataPos, in_texMin)))      
      {      
      break;      
      }      
      
    // Early ray termination      
    // if the currently composited colour alpha is already fully saturated      
    // we terminated the loop or if we have hit an obstacle in the      
    // direction of they ray (using depth buffer) we terminate as well.      
    if((g_fragColor.a > g_opacityThreshold) ||       
       g_currentT >= g_terminatePointMax)      
      {      
      break;      
      }      
    ++g_currentT;
  }

  
  g_srcColor = vec4(0);
   g_srcColor.a = 1.0;
  for (int i = 0; i < in_noOfComponents; ++i)
    {
    vec4 tmp = computeColor(l_sampledValue, computeOpacity(l_sampledValue, i), i);
    g_srcColor[0] += tmp[0] * tmp[3] * in_componentWeight[i];
    g_srcColor[1] += tmp[1] * tmp[3] * in_componentWeight[i];
    g_srcColor[2] += tmp[2] * tmp[3] * in_componentWeight[i];
    g_srcColor[3] += tmp[3] * in_componentWeight[i];
    }
  g_fragColor = g_srcColor;

  return g_fragColor;
}

/**
 * Finalize specific modes and set output data.
 */
void finalizeRayCast()
{
  //VTK::Base::Exit

  //VTK::Terminate::Exit

  //VTK::Cropping::Exit

  //VTK::Clipping::Exit

  //VTK::Picking::Exit

    g_fragColor.r = g_fragColor.r * in_scale + in_bias * g_fragColor.a;
    g_fragColor.g = g_fragColor.g * in_scale + in_bias * g_fragColor.a;
    g_fragColor.b = g_fragColor.b * in_scale + in_bias * g_fragColor.a;
    if(m_lightonly == 0) {
        g_fragColor.a = fuseCoef;
    } else {
        g_fragColor.a = (g_fragColor.r + g_fragColor.g + g_fragColor.b)/3;
    }

    gl_FragData[0] = g_fragColor;

  //VTK::RenderToImage::Exit

  //VTK::DepthPass::Exit
}

//////////////////////////////////////////////////////////////////////////////
///
/// Main
///
//////////////////////////////////////////////////////////////////////////////
void main()
{
      
  initializeRayCast();    
  castRay(-1.0, -1.0);    
  finalizeRayCast();
}

