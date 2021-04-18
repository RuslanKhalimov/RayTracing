#pragma once

#include <vector>
#include <random>
#include "structs.h"

class Light {
public:
  Light(const vec3& origin, double intensity, const SpectralValues& kd);
  virtual SpectralValues calculateLuminance(const vec3& hitPoint,
                                            const vec3& N,
                                            const std::vector<Triangle>& triangles,
                                            int hitTriangleId) const = 0;

protected:
  const vec3 origin_;
  const double intensity_;
  const SpectralValues kd_;

  SpectralValues luminanceFromPoint(const vec3& lightPoint,
                                    const vec3& hitPoint,
                                    const vec3& N,
                                    const std::vector<Triangle>& triangles,
                                    int hitTriangleId) const;
};

class PointLight : public Light {
public:
  PointLight(const vec3& origin, double intensity, const SpectralValues& kd);
  SpectralValues calculateLuminance(const vec3& hitPoint,
                                    const vec3& N,
                                    const std::vector<Triangle>& triangles,
                                    int hitTriangleId) const override;
};

class RectangleLight : public Light {
public:
  RectangleLight(const vec3& origin, double intensity, const SpectralValues& kd, double xSize, double ySize);
  SpectralValues calculateLuminance(const vec3& hitPoint,
                                    const vec3& N,
                                    const std::vector<Triangle>& triangles,
                                    int hitTriangleId) const override;

private:
  const double xSize_;
  const double ySize_;
  static const int SUBDIVISION = 3;
};
