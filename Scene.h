#pragma once

#include <vector>
#include <memory>
#include <string>
#include "structs.h"
#include "Light.h"

class Scene {
public:
  void readSceneFromFiles(const std::string& geometryFile,
                          const std::string& materialsFile,
                          const std::string& lightsFile,
                          const std::string& cameraFile);
  void render(const std::string& outputFileName);

private:
  struct HitPoint {
    vec3 point;
    int triangleId;
  };

  HitPoint getIntersection(const Ray& ray);
  SpectralValues getLuminance(const Ray& ray);

  void readGeometry(const std::string& fileName, const std::string& materialsFileName);
  std::vector<std::shared_ptr<Material>> readMaterials(const std::string& fileName);
  void readLights(const std::string& fileName);
  void readCamera(const std::string& fileName);

  static const int ANTIALIASING_FACTOR = 4;

  std::vector<Triangle> triangles_;
  std::vector<std::unique_ptr<Light>> lights_;
  std::unique_ptr<Camera> camera_ = nullptr;
};
