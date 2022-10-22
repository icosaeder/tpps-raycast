/**
 * raycast.glsl
 * 
 * Raycasting-Based Renderer of Three-Parameter Point Set Solid Models
 * This code is suitable for running at ShaderToy: https://www.shadertoy.com/
 * 
 * Written by Konstantin Ryabinin under terms of MIT license.
 *
 * The MIT License (MIT)
 * Copyright (c) 2018 Konstantin Ryabinin
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy of this software 
 * and associated documentation files (the "Software"), to deal in the Software without restriction, 
 * including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, 
 * and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, 
 * subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in all copies or substantial 
 * portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT 
 * LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE 
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

#define MAX_DIST 100.0
#define PI 3.1415926
#define EPSILON 1.0e-3
#define MAX_NEWTON_STEP 100
#define U0 0
#define U1 1
#define V0 2
#define V1 3
#define W0 4
#define W1 5

vec4 uvwt(vec3 val, int sel)
{
    switch (sel)
    {
        case U0: return vec4(0.0, val.x, val.y, val.z);
        case U1: return vec4(1.0, val.x, val.y, val.z);
        case V0: return vec4(val.x, 0.0, val.y, val.z);
        case V1: return vec4(val.x, 1.0, val.y, val.z);
        case W0: return vec4(val.x, val.y, 0.0, val.z);
        case W1: return vec4(val.x, val.y, 1.0, val.z);
        default: return vec4(0.0);
    }
}

vec3 F(vec3 val, vec3 pP, vec3 pQ, int sel)
{
    vec4 param = uvwt(val, sel);

    // HINT: this is to be changed for particular solid model.
    const vec3 pA = vec3(0, 0, 0);
    const vec3 pB = vec3(3, 5, 0);
    const vec3 pC = vec3(8, 3, 0);
    const vec3 pD = vec3(4, 2, 7);
    float u = param.s;
    float v = param.t;
    float w = param.p;
    float t = param.q;
    float nu = 1.0 - u;
    float nv = 1.0 - v;
    float nw = 1.0 - w;
    return pA * u * v * nw + pB * nv * nw + pC * nu * v * nw + pD * w - pP - (pQ - pP) * t;
}

vec3 center()
{
    // HINT: this is to be changed for particular solid model.
    const vec3 pA = vec3(0, 0, 0);
    const vec3 pB = vec3(3, 5, 0);
    const vec3 pC = vec3(8, 3, 0);
    const vec3 pD = vec3(4, 2, 7);
    
    return (pA + pB + pC + pD) / 4.0;
}

float orbit()
{
    // HINT: this is to be changed for particular solid model.
    const vec3 pA = vec3(0, 0, 0);
    const vec3 pB = vec3(3, 5, 0);
    const vec3 pC = vec3(8, 3, 0);
    const vec3 pD = vec3(4, 2, 7);

    vec3 maxPt = max(max(max(pA, pB), pC), pD);
    return max(max(maxPt.x, maxPt.y), maxPt.z) * 2.0;
}

mat3 jacobian(vec3 val, vec3 pP, vec3 pQ, int sel)
{
    const float h = 2.0 * EPSILON;
    const vec3 dx = vec3(EPSILON, 0.0, 0.0);
    const vec3 dy = vec3(0.0, EPSILON, 0.0);
    const vec3 dz = vec3(0.0, 0.0, EPSILON);
    vec3 dxL = val - dx;
    vec3 dxR = val + dx;
    vec3 dyL = val - dy;
    vec3 dyR = val + dy;
    vec3 dzL = val - dz;
    vec3 dzR = val + dz;
    vec3 dFdx = (F(dxR, pP, pQ, sel) - F(dxL, pP, pQ, sel)) / h;
    vec3 dFdy = (F(dyR, pP, pQ, sel) - F(dyL, pP, pQ, sel)) / h;
    vec3 dFdz = (F(dzR, pP, pQ, sel) - F(dzL, pP, pQ, sel)) / h;
    return mat3(dFdx, dFdy, dFdz);
}

vec3 solveNewton(vec3 pP, vec3 pQ, int sel)
{
    vec3 result = vec3(0.0);
    float norm = 1.0;
    for (int i = 0; i < MAX_NEWTON_STEP && norm > EPSILON; ++i)
    {
        mat3 mJ = jacobian(result, pP, pQ, sel);
        if (abs(determinant(mJ)) < EPSILON)
        {
            if (i == 0)
                mJ = mat3(1.0);
            else
                return vec3(2.0, 2.0, 2.0);
        }
        vec3 newResult = result - inverse(mJ) * F(result, pP, pQ, sel);
        norm = distance(newResult, result);
        result = newResult;
    }
    return result;
}

float calcDistance(vec3 rO, vec3 rD, float maxDist)
{
    vec3 pP = rO;
    vec3 pQ = rO + rD * maxDist;
    float t = 2.0;
    for (int i = 0, j = 0; i < 6; ++i)
    {
        vec3 result = solveNewton(pP, pQ, i);
        if (result.x > 0.0 && result.x < 1.0 &&
            result.y > 0.0 && result.y < 1.0 &&
            result.z > 0.0 && result.z < 1.0)
        {
            t = min(t, result.z);
            ++j;
        }
        if (j == 2)
            break;
    }
    return t * maxDist;
}

vec3 calcNormal(vec3 rO, vec3 rD, float dist, float maxDist, vec3 i, vec3 j, vec3 k)
{
    return normalize(vec3(calcDistance(rO + i * EPSILON, rD, maxDist),
                          calcDistance(rO + j * EPSILON, rD, maxDist),
                          calcDistance(rO + k * EPSILON, rD, maxDist)) -
                     vec3(dist));
}

void mainImage(out vec4 fragColor, in vec2 fragCoord)
{
    vec2 uv = (fragCoord - iResolution.xy * 0.5) / iResolution.y;
    vec2 mouse = iMouse.xy / iResolution.xy * PI * 2.0;
  
    vec3 camPos = vec3(sin(mouse.x) * sin(mouse.y), cos(mouse.y), cos(mouse.x) * sin(mouse.y)) * orbit();
    vec3 camTarget = center();
  
    vec3 camZ = normalize(camTarget - camPos);
    vec3 camX = normalize(cross(vec3(0.0, 1.0, 0.0), camZ));
    vec3 camY = normalize(cross(camZ, camX));

    vec3 color = vec3(0);
    vec3 rO = camPos;
    vec3 rD = normalize(camX * uv.x + camY * uv.y + camZ);

    float d = calcDistance(rO, rD, MAX_DIST);

    if (d > MAX_DIST)
        color = vec3(0.8);
    else
    {
        vec3 normal = calcNormal(rO, rD, d, MAX_DIST, camX, camY, camZ);
        vec3 lightPosition = vec3(0.0, 3.0, 0.0);
        vec3 surfacePoint = vec3(uv, 1.0) * d;
        vec3 lightDirection = normalize(lightPosition - surfacePoint);
        float diffuse = max(dot(normal, lightDirection), 0.1);
        color = pow(vec3(205.0/255.0, 127.0/255.0, 50.0/255.0) * diffuse, vec3(0.45));
    }

    fragColor = vec4(color, 1.0);
}
