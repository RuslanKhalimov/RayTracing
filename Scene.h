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
  void readGeometry(const std::string& fileName, const std::string& materialsFileName);
  std::vector<SpectralValues> readMaterials(const std::string& fileName);
  void readLights(const std::string& fileName);
  void readCamera(const std::string& fileName);

  std::vector<Triangle> triangles_;
  std::vector<std::unique_ptr<Light>> lights_;
  std::unique_ptr<Camera> camera_ = nullptr;
};
