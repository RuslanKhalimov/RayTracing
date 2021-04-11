#pragma once

#include <unordered_map>
#include <cmath>

struct vec3 {
  double x, y, z;

  vec3() {}

  vec3(double x, double y, double z) : x(x), y(y), z(z) {}

  vec3 operator+(const vec3& v) const {
    return vec3(x + v.x, y + v.y, z + v.z);
  }

  vec3 operator-(const vec3& v) const {
    return vec3(x - v.x, y - v.y, z - v.z);
  }

  vec3 operator*(double d) const {
    return vec3(x * d, y * d, z * d);
  }

  vec3 operator/(double d) const {
    return vec3(x / d, y / d, z / d);
  }

  double dot(const vec3& v) const {
    return x * v.x + y * v.y + z * v.z;
  }

  vec3 cross(const vec3& v) const {
    return vec3(y * v.z - z * v.y, z * v.x - x * v.z, x * v.y - y * v.x);
  }

  vec3 normalize() const {
    double len = length();
    return vec3(x / len, y / len, z / len);
  }

  double length() const {
    return sqrt(x * x + y * y + z * z);
  }
};

struct SpectralValues {
  std::unordered_map<int, double> values;

  SpectralValues() : SpectralValues(0) {}

  explicit SpectralValues(double value) {
    values = {{400, value}, {500, value}, {600, value}, {700, value}};
  }

  SpectralValues(const std::unordered_map<int, double>& kd) : values(kd) {}

  bool isZero() const {
    return values.at(400) == 0 && values.at(500) == 0 && values.at(600) == 0 && values.at(700) == 0;
  }

  SpectralValues operator*(double x) const {
    SpectralValues res(x);
    res *= *this;
    return res;
  }

  SpectralValues operator*(const SpectralValues& x) const {
    SpectralValues res = *this;
    res *= x;
    return res;
  }

  SpectralValues operator/(double x) const {
    SpectralValues res = *this;
    res.values[400] /= x;
    res.values[500] /= x;
    res.values[600] /= x;
    res.values[700] /= x;
    return res;
  }

  void operator*=(const SpectralValues& other) {
    values[400] *= other.values.at(400);
    values[500] *= other.values.at(500);
    values[600] *= other.values.at(600);
    values[700] *= other.values.at(700);
  }

  void operator+=(const SpectralValues& other) {
    values[400] += other.values.at(400);
    values[500] += other.values.at(500);
    values[600] += other.values.at(600);
    values[700] += other.values.at(700);
  }
};

struct Ray {
  vec3 origin;
  vec3 direction;
  SpectralValues kdeff;
  SpectralValues luminance;

  Ray() {}

  Ray(const vec3& origin, const vec3& direction)
    : origin(origin), direction(direction), kdeff(SpectralValues(1)), luminance(SpectralValues(0)) {}
};

struct Triangle {
  int objectId;
  vec3 v1, v2, v3;
  SpectralValues color;

  Triangle(int objectId, const vec3& v1, const vec3& v2, const vec3& v3, const SpectralValues& color)
      : objectId(objectId), v1(v1), v2(v2), v3(v3), color(color) {}

  bool hitTest(const Ray& ray, double& t) const {
    vec3 e1 = v2 - v1;
    vec3 e2 = v3 - v1;
    vec3 pvec = ray.direction.cross(e2);
    double det = e1.dot(pvec);

    if (det < 1e-8 && det > -1e-8) return false;

    double inv_det = 1 / det;
    vec3 tvec = ray.origin - v1;
    double u = tvec.dot(pvec) * inv_det;
    if (u < 0 || u > 1) return false;

    vec3 qvec = tvec.cross(e1);
    double v = ray.direction.dot(qvec) * inv_det;
    if (v < 0 || v + u > 1) return false;

    double dist = e2.dot(qvec) * inv_det;
    if (dist < t && dist > 1e-8) {
      t = dist;
      return true;
    }
    return false;
  }

  vec3 getNormal(const vec3& direction) {
    vec3 N = (v2 - v1).cross(v3 - v1).normalize();
    if (N.dot(direction) < 0) {
      N = (v3 - v1).cross(v2 - v1).normalize();
    }
    return N;
  }
};

struct Light {
  vec3 origin;
  double intensity;
  SpectralValues kd;

  Light(const vec3 origin, const double& intensity, const SpectralValues& kd)
    : origin(origin), intensity(intensity), kd(kd) {
  }

  SpectralValues calculateLuminance(const vec3& hitPoint,
                                    const vec3& N,
                                    const std::vector<Triangle>& triangles,
                                    int hitTriangleId) const {
    vec3 lightDirection = (origin - hitPoint).normalize();
    double dist = (origin - hitPoint).length();
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

    double E = (intensity / (dist * dist)) * cosTheta;
    return triangles[hitTriangleId].color * kd * E / M_PI;
  }
};

struct Camera {
  vec3 origin;
  vec3 direction;
  int width, height;

  Camera(const vec3& origin, const vec3& target, int width, int height) : origin(origin), width(width), height(height) {
    direction = (target - origin).normalize();
  }
};
