//VTK::System::Dec

varying vec3 ip_textureCoords;
varying vec3 ip_vertexPos;

vec4 g_fragColor = vec4(0.0);

vec3 g_dataPos;
vec3 g_dirStep;
float g_rayStepLength;
vec4 g_srcColor;
vec4 g_eyePosObj;
bool g_exit;
bool g_skip = false;
uniform vec4 in_volume_scale;
uniform vec4 in_volume_bias;

//VTK::Output::Dec
// Volume dataset
uniform sampler3D in_volume;
uniform int in_noOfComponents;

uniform int in_independentComponents;

uniform sampler2D in_noiseSampler;
uniform sampler2D in_depthSampler;

// Camera position
uniform vec3 in_cameraPos;

// view and model matrices
uniform mat4 in_projectionMatrix;
uniform mat4 in_inverseProjectionMatrix;
uniform mat4 in_modelViewMatrix;
uniform mat4 in_inverseModelViewMatrix;
uniform mat4 in_textureDatasetMatrix;
uniform mat4 in_inverseTextureDatasetMatrix;
uniform mat4 in_textureOriginMatrix;
uniform mat4 in_inverseVolumeMatrix;
uniform mat4 in_volumeMatrix;
mat4 invTextureOriginMatrix = inverse(in_textureOriginMatrix);

// Ray step size
uniform vec3 in_cellStep;
uniform vec2 in_scalarsRange[4];
uniform vec3 in_cellSpacing;
uniform float in_croppingPlanes[6];

// Sample distance
uniform float in_sampleDistance;

// Scales
uniform vec3 in_cellScale;
uniform vec2 in_windowLowerLeftCorner;
uniform vec2 in_inverseOriginalWindowSize;
uniform vec2 in_inverseWindowSize;
uniform vec3 in_textureExtentsMax;
uniform vec3 in_textureExtentsMin;

//ray and light angleuniform int in_twoSidedLighting;
vec3 g_xvec;
vec3 g_yvec;
vec3 g_zvec;

// Material and lighting
uniform vec3 in_diffuse[4];
uniform vec3 in_ambient[4];
uniform vec3 in_specular[4];
uniform float in_shininess[4];
uniform vec3 in_lightAmbientColor[1];
uniform vec3 in_lightDiffuseColor[1];
uniform vec3 in_lightSpecularColor[1];
vec4 g_lightPosObj;
vec3 g_ldir;
vec3 g_vdir;
vec3 g_h;
vec3 g_aspect;

vec3 computeRayDirection() {
    return normalize(ip_vertexPos.xyz - g_eyePosObj.xyz);
}

//VTK::Picking::Dec

uniform float in_clippingPlanes[49];
uniform float in_scale;
uniform float in_bias;


/*===================================VQ===============================*/
// cropping box coordinate in texture domain
//g_skip = false;
vec4 l_bb_min = vec4(in_croppingPlanes[0], in_croppingPlanes[2], in_croppingPlanes[4], 1.0);
// Maximum texture access coordinate
vec4 l_bb_max = vec4(in_croppingPlanes[1], in_croppingPlanes[3], in_croppingPlanes[5], 1.0);
/*===================================VQ===============================*/
//VQ:variables to support vtkvolume usermatrix
uniform mat4 in_originMatrix;
uniform mat4 in_InverseOriginMatrix;
uniform mat4 in_volumeUserMatrix;
uniform mat4 in_volumeInvUserMatrix;
uniform mat4 in_flipMatrix;

uniform sampler2D in_opacityTransferFunc;
float computeOpacity(vec4 scalar) {
    return texture2D(in_opacityTransferFunc, vec2(scalar.w, 0)).r;
}

uniform sampler3D in_gradient;

uniform sampler2D in_gradientTransferFunc;
float computeGradientOpacity(vec4 grad) {
    return texture2D(in_gradientTransferFunc, vec2(grad.w,0)).r;
}

uniform sampler2D in_colorTransferFunc;
vec4 computeColor(vec3 pos, float opacity) {
    vec4 scalar = texture3D(in_volume, pos);
    scalar = scalar * in_volume_scale + in_volume_bias;
    scalar = vec4(scalar.r, scalar.r, scalar.r, scalar.r);
    return vec4(texture2D(in_colorTransferFunc, vec2(scalar.w,0)).xyz, opacity);
}

/*===================================VQ===============================*/
//VQ:tissue textures, will need to do string replace dynamically
//uniform sampler1D in_custom1DTransferFunc0;
//uniform sampler2D in_custom2DTransferFunc0;

/*===================================VQ===============================*/
//VQ:define scatter event, this is the location where light get absored by tissue
struct scatterEvent {
    vec3 minT;
    vec3 pos;
    vec3 norm;
    float normWeight;
    float lightWeight;
    vec3 refRay;
    int event; //0 for volume, 1 for light, 2 for other object
    int valid;
};

float random (vec2 st) {
    return fract(sin(dot(st.xy, vec2(12.9898,78.233)))* 43758.5453123);
}

/*===================================VQ===============================*/
//VQ:get smoothed volume intensity with mean filter
vec4 ComputeSmoothVolumeIntensity(vec3 dataPos, int smoothness) {
    vec4 intensity = vec4(1.0);
    if (smoothness == 0) {
        intensity = texture3D(in_volume, dataPos);
    } else {
        vec4 sum = vec4(0.0);
        float cnt = 0.0;
        for (int i = -smoothness; i <= smoothness; i++) {
            for (int j = -smoothness; j <= smoothness; j++) {
                for (int k = -smoothness; k <= smoothness; k++) {
                    vec3 newPos = vec3(dataPos.x + float(i) * in_cellStep[0], dataPos.y + float(j) * in_cellStep[1], dataPos.z + float(k) * in_cellStep[2]);
                    intensity = texture3D(in_volume, newPos);
                    sum += intensity;
                    cnt += 1.0;
                }
            }
        }
        intensity = sum / cnt;
    }
    intensity = intensity * in_volume_scale + in_volume_bias;
    return intensity;
}

/*===================================VQ===============================*/
//VQ:compute surface gradient
vec3 ComputeGradientNormal(vec3 dataPos, int smoothness) {
    vec4 gradientNormal = vec4(0.0);
    if (dataPos.x <= 0.0 ||
        dataPos.x >= 1.0 ||
        dataPos.y <= 0.0 ||
        dataPos.y >= 1.0 ||
        dataPos.z <= 0.0 ||
        dataPos.z >= 1.0) {
        return gradientNormal.rgb;
    }

    vec3 direction = vec3(-1.0, 1.0, 1.0);
    float stepscale = 2.0;
    for (int i = 0; i < 3; i++) {
        vec3 g1dataPos = dataPos;
        g1dataPos[i] = g1dataPos[i] - stepscale * direction[i] * in_cellStep[i];
        vec3 g2dataPos = dataPos;
        g2dataPos[i] = g2dataPos[i] + stepscale * direction[i] * in_cellStep[i];
        vec4 g1 = ComputeSmoothVolumeIntensity(g1dataPos, smoothness);
        vec4 g2 = ComputeSmoothVolumeIntensity(g2dataPos, smoothness);
        g1 = g1 * in_volume_scale + in_volume_bias;
        g2 = g2 * in_volume_scale + in_volume_bias;
        gradientNormal[i] = length(g2.rgb) - length(g1.rgb);
    }
    return gradientNormal.rgb;
}

/*===================================VQ===============================*/
//VQ:compute ray intercept with light, now unused
scatterEvent intersectLights(vec3 ray) {
    scatterEvent lightEvent;
    lightEvent.event = 1;
    lightEvent.valid = 0;

    return lightEvent;
}

/*===================================VQ===============================*/
//VQ:compute ray intercept with object like roi, now unused
scatterEvent intersectObject(vec3 pos, vec3 rayDir, int totalSamplePoint) {
    scatterEvent objectEvent;
    objectEvent.event = 2;
    if (totalSamplePoint > 0) {
        objectEvent.valid = 1;
        objectEvent.pos = pos + rayDir * float(totalSamplePoint);
    } else {
        objectEvent.valid = 0;
    }
    objectEvent.valid = 0;
    return objectEvent;
}

/*===================================VQ===============================*/
//VQ:struct that store shader color
struct shaderInfo {
    vec3 reflectRay;
    vec4 emissionColor;
    vec4 diffuseColor;
    vec4 specularColor;
    float ior;
    float glossiness;
};

/*===================================VQ===============================*/
//VQ:light information
struct lightInfo {
    vec3 surfaceN;
    vec3 surfaceP;
    vec3 pos;
    int type; //0 for point and 1 for plane
};

const float INV_PI_F = 0.31830988618379067154;
const float INV_TWO_PI_F = 0.15915494309189533577;
const float PI_F = 3.141592654;

struct normInfo {
    vec3 ldir;
    vec3 vdir;
    vec3 h;
};

/*===================================VQ===============================*/
//VQ:light information.lightPos is in volume space and eyepos is in world space
normInfo computeLightCameraNormal(vec3 lightpos) {

    // Light position in dataset space
    vec4 t_lightPosObj = (in_inverseVolumeMatrix * vec4(lightpos, 1.0));
    if (t_lightPosObj.w != 0.0) {
        t_lightPosObj.x /= t_lightPosObj.w;
        t_lightPosObj.y /= t_lightPosObj.w;
        t_lightPosObj.z /= t_lightPosObj.w;
        t_lightPosObj.w = 1.0;
    }

    normInfo g_info;
    g_info.ldir = normalize(t_lightPosObj.xyz - ip_vertexPos);
    g_info.vdir = normalize(g_eyePosObj.xyz - ip_vertexPos);
    g_info.h = normalize(g_info.ldir + g_info.vdir);
    return g_info;
}

void computeLightNormal() {
    // Light position in dataset space
    g_lightPosObj = (in_inverseVolumeMatrix * vec4(in_cameraPos, 1.0));
    if (g_lightPosObj.w != 0.0) {
        g_lightPosObj.x /= g_lightPosObj.w;
        g_lightPosObj.y /= g_lightPosObj.w;
        g_lightPosObj.z /= g_lightPosObj.w;
        g_lightPosObj.w = 1.0;
    }
    g_ldir = normalize(g_lightPosObj.xyz - ip_vertexPos);
    g_vdir = normalize(g_eyePosObj.xyz - ip_vertexPos);
    g_h = normalize(g_ldir + g_vdir);
    return;
}

/*===================================VQ===============================*/
//VQ:get shader information, it map volume intensity to final color with transfer funciton or texture
shaderInfo getShader(scatterEvent event) {
    shaderInfo shader;
    shader.reflectRay = event.refRay;
    shader.emissionColor = computeColor(event.pos, 1.0);
    float mapping = length(cross(event.norm + 0.5, event.pos - 0.5));
    float radius = 6.0;
    //vec2 coorMapping = vec2(sqrt(radius * 2.0 / (radius + event.pos.z)) * event.pos.x, sqrt(radius * 2.0 / (radius + event.pos.z)) * event.pos.y);
    //vec4 diffusecol = texture2D(in_custom2DTransferFunc0, coorMapping);
    shader.diffuseColor = shader.emissionColor;
    shader.specularColor = vec4(0.9, 0.9, 0.9, shader.emissionColor.a);
    shader.ior = 0.5;
    shader.glossiness = in_shininess[0];
    return shader;
}

/*===================================VQ===============================*/
//VQ:get light info, now only one spot light from camera
lightInfo getLightInfo() {
    lightInfo planeLight;
    planeLight.type = 1;
    planeLight.surfaceP = g_eyePosObj.xyz;
    planeLight.surfaceN = g_vdir;
    return planeLight;
}

/*===================================VQ===============================*/
//VQ: shader, compute three lights
vec4 uniformSampleOneLight(scatterEvent event) {
    //init shader
    shaderInfo shader = getShader(event);

    //init light
    lightInfo planeLight = getLightInfo();
    vec3 scatterPointtoLight = normalize(planeLight.surfaceP - event.pos);

    //BRDF information
    vec3 normN = normalize(event.norm);
    vec3 normU = normalize(cross(event.norm, event.refRay));
    vec3 normV = normalize(cross(event.norm, normU));
    vec3 localRefRay = vec3(dot(event.refRay, normU), dot(event.refRay, normV), dot(event.refRay, normN));
    vec3 localScatterPointtoLight = vec3(dot(scatterPointtoLight, normU), dot(scatterPointtoLight, normV), dot(scatterPointtoLight, normN));

    //diffuse effect
    float reflactfactor = abs(localScatterPointtoLight.z);
    vec4 diffuseColor = shader.diffuseColor * (reflactfactor);
    float cosThetaO = abs(localRefRay.z);
    float cosThetaI = abs(localScatterPointtoLight.z);
    vec4 reflectColor = diffuseColor;

    //specular effect, the method is from exposure render
    float lambertianpdf = pow(cosThetaI, 2.0) * INV_PI_F;
    float bsfpdf = lambertianpdf;
    if (cosThetaO * cosThetaI > 0.0) {
        vec3 wh = normalize(localRefRay + localScatterPointtoLight);
        float cosTheta = abs(dot(localRefRay, localScatterPointtoLight));
        float microfacepdf = (shader.glossiness) * pow(cosTheta, shader.glossiness) / (2.0 * PI_F * 4.0 * cosTheta);
        bsfpdf += microfacepdf;
    }
    vec4 specularColor = vec4(0.0);
    if (bsfpdf > 0.0) {
        float lightpdf = normalize(distance(planeLight.surfaceP, event.pos)) / abs(dot(planeLight.surfaceN, scatterPointtoLight) * 4.0);
        float weightpdf = lightpdf * lightpdf / (lightpdf * lightpdf + bsfpdf * bsfpdf);
        weightpdf = 1.0 - weightpdf;
        specularColor = bsfpdf * abs(dot(scatterPointtoLight, planeLight.surfaceN)) * weightpdf * shader.specularColor;
    }
    vec4 finalColor = in_ambient[0][0] * shader.emissionColor + in_diffuse[0][0] * diffuseColor + in_specular[0][0] * specularColor;

    //tone mapping
    bool toneMapping = false;
    if (toneMapping) {
        float exposure = 0.25;
        float invexposure = 1.0 / exposure;
        finalColor = 1.0 - exp(-finalColor * invexposure);
    }
    finalColor = clamp(finalColor, vec4(0.0), vec4(1.0));
    return finalColor;
}

/*===================================VQ===============================*/
//VQ:compute the location where light get absored by tissue, now it only store one scatter point
vec4 sampleVolume(vec3 pos, vec3 rayDir, vec3 dirStep, float totalObjectSamplePoint) {


    float energy = 20.0;
    float remainEnergy = 0.0;
    float sigmaT = 0.0;
    float densityScale = 20.0;
    bool stop = false;
    float normalizedIntensity = 1.0;

    const vec3 tex_min = vec3(0.0);
    const vec3 tex_max = vec3(1.0);
    vec4 color = vec4(0.0);

    float l_currentT = 0.0;
    while (remainEnergy < energy && !stop && l_currentT <= totalObjectSamplePoint) {
        g_skip = dot(sign(pos - l_bb_min.xyz), sign(l_bb_max.xyz - pos)) < 3.0;
        vec4 intensity = texture3D(in_volume, pos);
        intensity = intensity * in_volume_scale + in_volume_bias;
        intensity.rgba = intensity.rrrr;
        float opacity = computeOpacity(intensity);
        normalizedIntensity = intensity.r;
        sigmaT = opacity * densityScale;
        if(sigmaT > 0.01 && !g_skip) {
            remainEnergy += sigmaT;
            scatterEvent volumeEvent;
            volumeEvent.event = 0;
            volumeEvent.valid = 1;
            vec3 norm = ComputeGradientNormal(pos, 0);
            norm = norm.xyz / in_cellSpacing;
            volumeEvent.pos = pos;
            volumeEvent.norm = normalize(norm);
            volumeEvent.normWeight = normalizedIntensity;
            volumeEvent.refRay = normalize(rayDir - 2.0 * (dot(rayDir, volumeEvent.norm)) * volumeEvent.norm);
            float weight = sigmaT / energy;
            vec4 scatterColor = uniformSampleOneLight(volumeEvent);
            color += weight * scatterColor;
            l_currentT += 1.0;
        }
        float zoomStepRate = clamp(densityScale * 10.0 * opacity, 1.0, 20.0);
        pos += dirStep / zoomStepRate;
        stop = dot(sign(pos - tex_min), sign(tex_max - pos)) < 3.0;
    }
    return color;
}

/*===================================VQ===============================*/
//VQ:test ray intercept with either volume, light, or object, now only volume
vec4 sampleRay(vec3 pos, vec3 rayDir, vec3 dirStep, float totalObjectSamplePoint) {

    //scatterEvent event[3];
    vec4 color = sampleVolume(pos, rayDir, dirStep, totalObjectSamplePoint);
    //event[1] = intersectLights(pos);
    //event[2] = intersectObject(pos, rayDir, totalObjectSamplePoint);
    return color;
}

/*===================================VQ===============================*/
//VQ: select scatter position
vec4 multipleScattering(vec3 pos, vec3 rayDir, vec3 dirStep, float totalObjectSamplePoint) {
    vec4 color = sampleRay(pos, rayDir, dirStep, totalObjectSamplePoint);
    color.a = 1.0;
    return color;
}
/*******************************isosurface, vtk**********************************/

void main() {
    /// Initialize g_fragColor (output) to 0
    g_fragColor = vec4(0.0);
    g_dirStep = vec3(0.0);
    g_srcColor = vec4(0.1);
    g_exit = false;

/*compute bounding box minmax coordinate in texture domain*/
    l_bb_min = in_inverseTextureDatasetMatrix * l_bb_min;
    if (l_bb_min.w != 0.0)
    {
       l_bb_min.x /= l_bb_min.w;
       l_bb_min.y /= l_bb_min.w;
       l_bb_min.z /= l_bb_min.w;
       l_bb_min.w = 1.0;
    }
    l_bb_max = in_inverseTextureDatasetMatrix * l_bb_max;
    if (l_bb_max.w != 0.0)
    {
       l_bb_max.x /= l_bb_max.w;
       l_bb_max.y /= l_bb_max.w;
       l_bb_max.z /= l_bb_max.w;
       l_bb_max.w = 1.0;
    }
    // Get the 3D texture coordinates for lookup into the in_volume dataset
    vec4 texremapped = vec4(ip_textureCoords, 1.0);
    texremapped = invTextureOriginMatrix * in_flipMatrix * in_textureOriginMatrix * texremapped;
    g_dataPos = texremapped.xyz;

    // Eye position in object space
    g_eyePosObj = (in_InverseOriginMatrix * in_volumeInvUserMatrix * in_originMatrix * vec4(in_cameraPos, 1.0));
    if (g_eyePosObj.w != 0.0) {
        g_eyePosObj.x /= g_eyePosObj.w;
        g_eyePosObj.y /= g_eyePosObj.w;
        g_eyePosObj.z /= g_eyePosObj.w;
        g_eyePosObj.w = 1.0;
    }

    // Getting the ray marching direction (in object space);
    vec3 rayDir = computeRayDirection();

    g_dirStep = (in_inverseTextureDatasetMatrix * invTextureOriginMatrix * in_flipMatrix * in_textureOriginMatrix * vec4(rayDir, 0.0)).xyz * in_sampleDistance;
    g_rayStepLength = length(g_dirStep);

    g_dataPos += g_dirStep * (texture2D(in_noiseSampler, g_dataPos.xy).x);

    bool stop = false;

    vec2 fragTexCoord = (gl_FragCoord.xy - in_windowLowerLeftCorner) * in_inverseWindowSize;
    vec4 l_depthValue = texture2D(in_depthSampler, fragTexCoord);
    float l_terminatePointMax = 0.0;

    // Depth test
    if (gl_FragCoord.z >= l_depthValue.x) {
        discard;
    }

    vec4 terminatePoint;
    terminatePoint.x = (gl_FragCoord.x - in_windowLowerLeftCorner.x) * 2.0 * in_inverseWindowSize.x - 1.0;
    terminatePoint.y = (gl_FragCoord.y - in_windowLowerLeftCorner.y) * 2.0 * in_inverseWindowSize.y - 1.0;
    terminatePoint.z = (2.0 * l_depthValue.x - (gl_DepthRange.near + gl_DepthRange.far)) / gl_DepthRange.diff;
    terminatePoint.w = 1.0;

    terminatePoint = in_inverseTextureDatasetMatrix *
        in_originMatrix * in_volumeInvUserMatrix * in_InverseOriginMatrix * in_inverseModelViewMatrix * in_inverseProjectionMatrix * terminatePoint;
    terminatePoint /= terminatePoint.w;
    l_terminatePointMax = length(terminatePoint.xyz - ip_textureCoords.xyz) / g_rayStepLength;

    computeLightNormal();

    // For all samples along the ray
    g_srcColor = multipleScattering(g_dataPos, rayDir, g_dirStep, l_terminatePointMax);

    g_fragColor.rgb = g_srcColor.rgb * g_srcColor.a;
    g_fragColor.a = g_srcColor.a;

    g_fragColor.r = g_fragColor.r * in_scale + in_bias * g_fragColor.a;
    g_fragColor.g = g_fragColor.g * in_scale + in_bias * g_fragColor.a;
    g_fragColor.b = g_fragColor.b * in_scale + in_bias * g_fragColor.a;

    gl_FragData[0] = g_fragColor;
}
