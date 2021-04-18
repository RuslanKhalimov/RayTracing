#include "Light.h"


Light::Light(const vec3& origin, double intensity, const SpectralValues& kd)
    : origin_(origin), intensity_(intensity), kd_(kd) {}


SpectralValues Light::luminanceFromPoint(const vec3& lightPoint,
                                         const vec3& hitPoint,
                                         const vec3& N,
                                         const std::vector<Triangle>& triangles,
                                         int hitTriangleId) const {
  assert(hitTriangleId != 1);

  vec3 lightDirection = (lightPoint - hitPoint).normalize();
  double dist = (lightPoint - hitPoint).length();
  double cosTheta = lightDirection.dot(N);
  if (cosTheta <= 0) { // back side
    return SpectralValues(0);
  }

  Ray shadowRay(hitPoint, lightDirection);
  for (int id = 0; id < triangles.size(); ++id) {
    if (id == hitTriangleId) {
      continue;
    }
    if (triangles[id].hitTest(shadowRay, dist)) { // shadow
      return SpectralValues(0);
    }
  }

  double E = (intensity_ / (dist * dist)) * cosTheta;
  auto material = triangles[hitTriangleId].material;
  return material->color * material->kd * kd_ * E / M_PI;
}


PointLight::PointLight(const vec3& origin, double intensity, const SpectralValues& kd)
    : Light(origin, intensity, kd) {}


SpectralValues PointLight::calculateLuminance(const vec3& hitPoint,
                                              const vec3& N,
                                              const std::vector<Triangle>& triangles,
                                              int hitTriangleId) const {
  return luminanceFromPoint(origin_, hitPoint, N, triangles, hitTriangleId);
}


RectangleLight::RectangleLight(const vec3& origin, double intensity, const SpectralValues& kd, double xSize, double ySize)
    : Light(origin, intensity, kd), xSize_(xSize), ySize_(ySize) {}


SpectralValues RectangleLight::calculateLuminance(const vec3& hitPoint,
                                                  const vec3& N,
                                                  const std::vector<Triangle>& triangles,
                                                  int hitTriangleId) const {
  SpectralValues resultLuminance(0);
  for (int i = 0; i < SUBDIVISION; ++i) {
    for (int j = 0; j < SUBDIVISION; ++j) {
      double x = origin_.x - xSize_ / 2 + (i + 1) * xSize_ / (SUBDIVISION + 1);
      double y = origin_.y - ySize_ / 2 + (j + 1) * ySize_ / (SUBDIVISION + 1);
      vec3 lightPoint(x, y, origin_.z);
      resultLuminance += luminanceFromPoint(lightPoint, hitPoint, N, triangles, hitTriangleId) / (SUBDIVISION * SUBDIVISION);
    }
  }
  return resultLuminance;
}
