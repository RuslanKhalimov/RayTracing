#pragma once

#include <unordered_map>
#include <cmath>
#include <cassert>
#include <memory>

struct vec3 {
  double x, y, z;

  vec3() {}

  vec3(double x, double y, double z) : x(x), y(y), z(z) {}

  vec3 operator-() const {
    return vec3(-x, -y, -z);
  }

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

  double compare(const SpectralValues& other) const {
    double result = 0;
    result += (values.at(400) - other.values.at(400)) * (values.at(400) - other.values.at(400));
    result += (values.at(500) - other.values.at(500)) * (values.at(500) - other.values.at(500));
    result += (values.at(600) - other.values.at(600)) * (values.at(600) - other.values.at(600));
    result += (values.at(700) - other.values.at(700)) * (values.at(700) - other.values.at(700));
    return sqrt(result);
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
  double ks;
  SpectralValues luminance;

  Ray() {}

  Ray(const vec3& origin, const vec3& direction, double ks = 1)
    : origin(origin), direction(direction), ks(ks), luminance(SpectralValues(0)) {}
};

struct Material {
  SpectralValues color;
  double kd;
  double ks;

  Material(const SpectralValues& color, double kd, double ks) : color(color), kd(kd), ks(ks) {
    assert(kd + ks - 1 < 1e-8);
  }
};

struct Triangle {
  int objectId;
  vec3 v1, v2, v3;
  std::shared_ptr<Material> material;

  Triangle(int objectId, const vec3& v1, const vec3& v2, const vec3& v3, const std::shared_ptr<Material>& material)
      : objectId(objectId), v1(v1), v2(v2), v3(v3), material(material) {}

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

struct Camera {
  const vec3 origin;
  const vec3 direction;
  const int width, height;
  const double viewAngle = M_PI / 2;

  Camera(const vec3& origin, const vec3& target, int width, int height)
      : origin(origin), width(width), height(height), direction((target - origin).normalize()) {}

  Ray castRay(int i, int j, int resolutionScaler = 1, int iOffset = 0, int jOffset = 0) const {
    assert(iOffset < resolutionScaler);
    assert(jOffset < resolutionScaler);

    double scaledWidth = width * resolutionScaler;
    double scaledHeight = height * resolutionScaler;
    int iScaled = (i - 1) * resolutionScaler + iOffset;
    int jScaled = (j - 1) * resolutionScaler + jOffset;

    double x = jScaled + 0.5 - scaledWidth / 2;
    double y = scaledWidth / tan(viewAngle / 2);
    double z = iScaled + 0.5 - scaledHeight / 2;
    vec3 dir = -vec3(x, y, z).normalize();

    return Ray(origin, dir);
  }
};
