#include "Scene.h"
#include <fstream>
#include <sstream>
#include <limits>
#include <iostream>

enum GeometryParseState {
  INITIAL,
  VERTICES,
  TRIANGLES
};

struct HitPoint {
  vec3 point;
  int triangleId;
  bool isShadow;
};

void Scene::readSceneFromFiles(const std::string& geometryFile,
                               const std::string& materialsFile,
                               const std::string& lightsFile,
                               const std::string& cameraFile) {
  readGeometry(geometryFile, materialsFile);
  readLights(lightsFile);
  readCamera(cameraFile);
}


void Scene::render(const std::string& outputFileName) {
  int width = camera_->width;
  int height = camera_->height;
  std::vector<std::vector<HitPoint>> hitPoints(height, std::vector<HitPoint>(width));

#pragma omp parallel for shared(height, width, hitPoints) default(none)
  for (int i = 0; i < height; ++i) {
    for (int j = 0; j < width; ++j) {
      Ray ray = camera_->castRay(i, j);

      double t = std::numeric_limits<double>::max();
      int triangleId = -1;
      for (int id = 0; id < triangles_.size(); ++id) {
        if (triangles_[id].hitTest(ray, t)) {
          triangleId = id;
        }
      }

      vec3 hitPoint = ray.origin + ray.direction * t;
      bool isShadow = false;
      if (triangleId != -1) {
        isShadow = true;
        vec3 N = triangles_[triangleId].getNormal(camera_->origin - hitPoint);
        for (const auto& light : lights_) {
          if (!light->calculateLuminance(hitPoint, N, triangles_, triangleId).isZero()) {
            isShadow = false;
            break;
          }
        }
      }

      hitPoints[i][j] = {hitPoint, triangleId, isShadow};
    }
  }

  // for antialiasing cast more rays on triangle and shadow edges pixel
  std::vector<std::vector<std::vector<HitPoint>>> scaledHitPoints(height, std::vector<std::vector<HitPoint>>(width));
#pragma omp parallel for shared(height, width, hitPoints, scaledHitPoints) default(none)
  for (int i = 0; i < std::max(width, height); ++i) {
    if (i < width) {
      scaledHitPoints[0][i].push_back(hitPoints[0][i]);
      scaledHitPoints[height - 1][i].push_back(hitPoints[height - 1][i]);
    }
    if (i < height) {
      scaledHitPoints[i][0].push_back(hitPoints[i][0]);
      scaledHitPoints[i][width - 1].push_back(hitPoints[i][width - 1]);
    }
  }
#pragma omp parallel for shared(height, width, hitPoints, scaledHitPoints) default(none)
  for (int i = 1; i < height - 1; ++i) {
    for (int j = 1; j < width - 1; ++j) {
      if (hitPoints[i][j].triangleId != hitPoints[i - 1][j].triangleId ||
          hitPoints[i][j].isShadow != hitPoints[i - 1][j].isShadow ||
          hitPoints[i][j].triangleId != hitPoints[i][j - 1].triangleId ||
          hitPoints[i][j].isShadow != hitPoints[i][j - 1].isShadow ||
          hitPoints[i][j].triangleId != hitPoints[i + 1][j].triangleId ||
          hitPoints[i][j].isShadow != hitPoints[i + 1][j].isShadow ||
          hitPoints[i][j].triangleId != hitPoints[i][j + 1].triangleId ||
          hitPoints[i][j].isShadow != hitPoints[i][j + 1].isShadow) {
        for (int iOffset = 0; iOffset < ANTIALIASING_FACTOR; ++iOffset) {
          for (int jOffset = 0; jOffset < ANTIALIASING_FACTOR; ++jOffset) {
            Ray ray = camera_->castRay(i, j, ANTIALIASING_FACTOR, iOffset, jOffset);

            double t = std::numeric_limits<double>::max();
            int triangleId = -1;
            for (int id = 0; id < triangles_.size(); ++id) { // TODO check only relevant triangles
              if (triangles_[id].hitTest(ray, t)) {
                triangleId = id;
              }
            }

            // `isShadow` is not used, so don't calculate it
            scaledHitPoints[i][j].push_back({ray.origin + ray.direction * t, triangleId, false});
          }
        }
      } else {
        scaledHitPoints[i][j].push_back(hitPoints[i][j]);
      }
    }
  }

  std::vector<std::vector<SpectralValues>> outLuminance(height, std::vector<SpectralValues>(width, SpectralValues(0)));
#pragma omp parallel for shared(height, width, scaledHitPoints, outLuminance) default(none)
  for (int i = 0; i < height; ++i) {
    for (int j = 0; j < width; ++j) {
      for (const HitPoint& nextHitPoint : scaledHitPoints[i][j]) {
        vec3 hitPoint = nextHitPoint.point;
        int triangleId = nextHitPoint.triangleId;
        if (triangleId == -1) {
          continue;
        }

        vec3 N = triangles_[triangleId].getNormal(camera_->origin - hitPoint);
        for (const std::unique_ptr<Light>& light : lights_) {
          outLuminance[i][j] += light->calculateLuminance(hitPoint, N, triangles_, triangleId);
        }
      }
      outLuminance[i][j] = outLuminance[i][j] / scaledHitPoints[i][j].size();
    }
  }

  std::ofstream out(outputFileName);

  for (int waveLength : {400, 500, 600, 700}) {
    out << "wave_length " << waveLength << std::endl;
    for (int i = 0; i < height; ++i) {
      for (int j = 0; j < width; ++j) {
        out << outLuminance[i][j].values[waveLength] << " ";
      }
      out << std::endl;
    }
    out << std::endl;
  }

  out.close();
}


void Scene::readGeometry(const std::string& fileName, const std::string& materialsFileName) {
  std::ifstream in(fileName);

  std::vector<SpectralValues> colors = readMaterials(materialsFileName);

  std::string line;
  GeometryParseState state = GeometryParseState::INITIAL;
  int currentObjectId = 0;
  std::vector<vec3> vertices;
  while (getline(in, line)) {
    switch (state) {
      case GeometryParseState::INITIAL: {
        if (line.find("vertices") != std::string::npos) {
          state = GeometryParseState::VERTICES;
        }
        break;
      }
      case GeometryParseState::VERTICES: {
        if (line.find("triangles") != std::string::npos) {
          state = GeometryParseState::TRIANGLES;
          break;
        }
        std::istringstream iss(line);
        double x, y, z;
        iss >> x >> z >> y;
        vertices.push_back(vec3(x, y, z));
        break;
      }
      case GeometryParseState::TRIANGLES: {
        if (line.find("parts") != std::string::npos) {
          state = GeometryParseState::INITIAL;
          vertices.clear();
          ++currentObjectId;
          break;
        }
        std::istringstream iss(line);
        int v1, v2, v3;
        iss >> v1 >> v2 >> v3;
        triangles_.push_back(Triangle(currentObjectId, vertices[v1], vertices[v2], vertices[v3], colors[currentObjectId]));
        break;
      }
    }
  }

  in.close();
}


std::vector<SpectralValues> Scene::readMaterials(const std::string& fileName) {
  std::ifstream in(fileName);

  std::vector<SpectralValues> colors;
  std::string line;
  int currentObjectId = -1;
  std::unordered_map<int, double> kds;
  while (getline(in, line)) {
    if (line.find("id") != std::string::npos) {
      if (currentObjectId != -1) {
        colors.push_back(SpectralValues(kds));
        kds.clear();
      }
      ++currentObjectId;
      continue;
    }
    std::istringstream iss(line);
    int waveLength;
    double kd;
    iss >> waveLength >> kd;
    kds[waveLength] = kd;
  }

  in.close();

  return colors;
}


void Scene::readLights(const std::string& fileName) {
  std::ifstream in(fileName);
  vec3 origin;
  double intensity;
  std::unordered_map<int, double> kds;
  in >> origin.x >> origin.y >> origin.z;
  in >> intensity;
  std::string line;
  getline(in, line); // '\n'
  while (getline(in, line)) {
    std::istringstream iss(line);
    int waveLength;
    double kd;
    iss >> waveLength >> kd;
    kds[waveLength] = kd;
  }
  lights_.push_back(std::make_unique<PointLight>(origin, intensity, kds));
  in.close();
}


void Scene::readCamera(const std::string& fileName) {
  std::ifstream in(fileName);
  vec3 origin, target;
  int width, height;
  std::unordered_map<int, double> color;
  in >> origin.x >> origin.y >> origin.z;
  in >> target.x >> target.y >> target.z;
  in >> width >> height;
  camera_ = std::make_unique<Camera>(origin, target, width, height);
  in.close();
}
