//VTK::System::Dec

varying vec3 ip_textureCoords;
varying vec3 ip_vertexPos;

vec4 g_fragColor = vec4(0.0);

vec3 g_dataPos;
vec3 g_dirStep;
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
mat4 invTextureOriginMatrix = inverse(in_textureOriginMatrix);

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

// scale and shift
uniform float in_clippingPlanes[49];
uniform float in_scale;
uniform float in_bias;
uniform float in_croppingPlanes[6];

/*===================================VQ===============================*/
// cropping box coordinate in texture domain
//g_skip = false;
vec4 l_bb_min = vec4(in_croppingPlanes[0], in_croppingPlanes[2], in_croppingPlanes[4], 1.0);
// Maximum texture access coordinate
vec4 l_bb_max = vec4(in_croppingPlanes[1], in_croppingPlanes[3], in_croppingPlanes[5], 1.0);

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
uniform vec4 in_componentWeight;

uniform float fuseCoef;
uniform int m_lightonly;
uniform int m_inverted;
uniform float mc_weight;
uniform float mc_threshold;
uniform float mc_transPeriod;
uniform vec3 mc_channelWeight;
uniform sampler2D in_opacityTransferFunc;
uniform sampler2D in_opacityTransferFunc1;
uniform sampler2D in_opacityTransferFunc2;

/*===================================VQ===============================*/
//VQ:variables to support vtkvolume usermatrix
uniform mat4 in_originMatrix;
uniform mat4 in_InverseOriginMatrix;
uniform mat4 in_volumeUserMatrix;
uniform mat4 in_volumeInvUserMatrix;
uniform mat4 in_flipMatrix;

vec4 computeOpacity(vec4 scalar) {
    vec4 opacity = vec4(0.0);
    opacity.r = texture2D(in_opacityTransferFunc, vec2(scalar.r, 0)).r;
    opacity.g = texture2D(in_opacityTransferFunc1, vec2(scalar.g, 0)).r;
    opacity.b = texture2D(in_opacityTransferFunc2, vec2(scalar.b, 0)).r;
    return opacity;
}

float computeOpacitySingle(vec4 scalar, int component) {
    if (component == 0) {
        return texture2D(in_opacityTransferFunc, vec2(scalar[0], 0)).r;
    }
    if (component == 1) {
        return texture2D(in_opacityTransferFunc1, vec2(scalar[1], 0)).r;
    }
    if (component == 2) {
        return texture2D(in_opacityTransferFunc2, vec2(scalar[2], 0)).r;
    }
}

//VQ::RGBRayCastingMethod::GradientDefinition

uniform sampler2D in_colorTransferFunc;
uniform sampler2D in_colorTransferFunc1;
uniform sampler2D in_colorTransferFunc2;

vec4 computeColor(vec4 scalar, vec4 opacity) {
    vec4 color = vec4(1.0);
    color.r = texture2D(in_colorTransferFunc, vec2(scalar.r,0)).r * opacity.r;
    color.g = texture2D(in_colorTransferFunc, vec2(scalar.g,0)).r * opacity.g;
    color.b = texture2D(in_colorTransferFunc, vec2(scalar.b,0)).r * opacity.b;
    return color;
}

vec3 computeRayDirection() {
    return normalize(ip_vertexPos.xyz - g_eyePosObj.xyz);
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

    vec3 direction = vec3(-1.0, 1.0, 1.0); //vec3(-sign(in_cameraPos.x), sign(in_cameraPos.z), sign(in_cameraPos.y));
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

/**************************************************************/

const float INV_PI_F = 0.31830988618379067154;
const float INV_TWO_PI_F = 0.15915494309189533577;
const float PI_F = 3.141592654;

/*===================================VQ===============================*/
//VQ:define scatter event, this is the location where light get absored by tissue
struct scatterEvent {
    vec3 minT;
    vec3 pos;
    vec3 norm;
    float normWeight;
    vec3 refRay;
    int event; //0 for volume, 1 for light, 2 for other object
    int valid;
};

/*===================================VQ===============================*/
//VQ:struct that store shader color
struct shaderInfo {
    vec3 norm;
    vec3 reflectRay;
    vec4 diffuseColor;
    vec4 specularColor;
    vec4 emissionColor;
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

/*===================================VQ===============================*/
//VQ:light information.lightPos is in volume space and eyepos is in world space
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
};

/*===================================VQ===============================*/
//VQ:get shader information, it map volume intensity to final color with transfer funciton or texture
shaderInfo getShader(scatterEvent event) {
    //define shader information
    shaderInfo shader;
    vec4 intensity = texture3D(in_volume, event.pos);
    intensity = intensity * in_volume_scale + in_volume_bias;
    shader.emissionColor = computeColor(intensity, vec4(1.0));
    shader.norm = event.norm;
    shader.reflectRay = event.refRay;
    shader.diffuseColor = shader.emissionColor;
    shader.specularColor = vec4(0.9, 0.9, 0.9, shader.emissionColor.a);
    shader.ior = 0.5;
    shader.glossiness = in_shininess[0];
    return shader;
}

/*===================================VQ===============================*/
//VQ:get light info, now only one spot light from camera
lightInfo getLight() {
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
    lightInfo planeLight = getLight();
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

    //specular effect
    float cosThetaO = abs(localRefRay.z);
    float cosThetaI = abs(localScatterPointtoLight.z);
    vec4 reflectColor = diffuseColor;
    float bsfpdf = 0.0;
    if (cosThetaO * cosThetaI > 0.0) {
        float lambertianpdf = pow(cosThetaI, 2.0) * INV_PI_F;
        bsfpdf += lambertianpdf;
    }
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
    vec4 finalColor = vec4(in_ambient[0], 1.0) * shader.emissionColor + vec4(in_diffuse[0], 1.0) * diffuseColor + vec4(in_specular[0], 1.0) * specularColor;

    //tone mapping
    if (false) {
        float exposure = 0.25;
        float invexposure = 1.0 / exposure;
        finalColor = 1.0 - exp(-finalColor * invexposure);
    }
    finalColor = clamp(finalColor, vec4(0.0), vec4(1.0));
    return finalColor;
}

/*===================================VQ===============================*/
//VQ:compute the location where light get absored by tissue, now it only store one scatter point
scatterEvent sampleVolume(vec3 pos, vec3 rayDir, vec3 dirStep) {
    scatterEvent volumeEvent;
    volumeEvent.event = 0;
    volumeEvent.valid = 1;

    float densityScale = 10.0;
    float energy = 20.0;
    float remainEnergy = 0.0;
    float sigmaT = 0.0;
    int cnt = 0;
    bool stop = false;
    vec4 normalizedIntensity = vec4(1.0);

    const vec3 tex_min = vec3(0);
    const vec3 tex_max = vec3(1);

    while (remainEnergy < energy && !stop) {
        g_skip = dot(sign(pos - l_bb_min.xyz), sign(l_bb_max.xyz - pos)) < 3.0;
        vec4 intensity = texture3D(in_volume, pos);
        intensity = intensity * in_volume_scale + in_volume_bias;
        vec4 opacity = computeOpacity(intensity);
        if(!g_skip) {normalizedIntensity = intensity;
        sigmaT = length(opacity) * densityScale;
        remainEnergy += sigmaT;}
        float zoomStepRate = clamp(in_sampleDistance * densityScale * 10.0 * length(opacity), 1.0, 20.0);
        pos += dirStep / zoomStepRate;
        stop = dot(sign(pos - tex_min), sign(tex_max - pos)) < 3.0;
        cnt++;
    }
    vec3 norm = ComputeGradientNormal(pos, 2);
    norm = norm.xyz / in_cellSpacing;
    volumeEvent.pos = pos;
    volumeEvent.norm = normalize(norm);
    volumeEvent.normWeight = length(normalizedIntensity);
    volumeEvent.refRay = normalize(rayDir - 2.0 * (dot(rayDir, volumeEvent.norm)) * volumeEvent.norm);
    return volumeEvent;
}

scatterEvent intersectLights(vec3 ray) {
    scatterEvent lightEvent;
    lightEvent.event = 1;
    lightEvent.valid = 0;

    return lightEvent;
}

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
//VQ:test ray intercept with either volume, light, or object, now only volume
scatterEvent sampleRay(vec3 pos, vec3 rayDir, vec3 dirStep, int totalObjectSamplePoint) {

    scatterEvent event[3];
    event[0] = sampleVolume(pos, rayDir, dirStep);
    event[1] = intersectLights(pos);
    event[2] = intersectObject(pos, rayDir, totalObjectSamplePoint);

    int idx = -1;
    vec3 intersect = vec3(1.0);
    for (int i = 0; i < 3; i++) {
        if (event[i].valid == 1 && (length(event[i].pos) < length(intersect))) {
            idx = i;
        }
    }
    if (idx < 3 && idx >= 0) {
        return event[idx];
    } else {
        scatterEvent invalidEvent;
        invalidEvent.valid = 0;
        return invalidEvent;
    }
}

/*===================================VQ===============================*/
//VQ: select scatter position
vec4 singleScattering(vec3 pos, vec3 rayDir, vec3 dirStep, int totalObjectSamplePoint) {
    scatterEvent event = sampleRay(pos, rayDir, dirStep, totalObjectSamplePoint);

    vec4 color = vec4(0.0);
    if (event.valid == 1 && event.event == 0) { // volume
        color += uniformSampleOneLight(event);
    }
    if (event.valid == 1 && event.event == 1) { // light
        //color += vec4(0.9);
    }
    if (event.valid == 1 && event.event == 2) { // object
        //color += computeLightShadingColor(event.pos, 1);
    }

    color.a = 1.0;
    return color;
}

/**************************************************************/

void main() {
    /// Initialize g_fragColor (output) to 0
    g_fragColor = vec4(0.0);
    g_dirStep = vec3(0.0);
    g_srcColor = vec4(0.0);
    g_exit = false;

/*compute bounding box minmax coordinate in texture domain*/    l_bb_min = in_inverseTextureDatasetMatrix * l_bb_min;
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

    g_dataPos += g_dirStep * (texture2D(in_noiseSampler, g_dataPos.xy).x);

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
    l_terminatePointMax = length(terminatePoint.xyz - ip_textureCoords.xyz) / length(g_dirStep);

    // We get data between 0.0 - 1.0 range
    computeLightNormal();

    // For all samples along the ray
    g_srcColor = singleScattering(g_dataPos, rayDir, g_dirStep, int(l_terminatePointMax));

    g_fragColor = g_srcColor;

    g_fragColor.r = g_fragColor.r * in_scale + in_bias * g_fragColor.a;
    g_fragColor.g = g_fragColor.g * in_scale + in_bias * g_fragColor.a;
    g_fragColor.b = g_fragColor.b * in_scale + in_bias * g_fragColor.a;
    gl_FragData[0] = g_fragColor;

}
