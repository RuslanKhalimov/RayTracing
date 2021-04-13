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
  std::vector<std::vector<std::pair<vec3, int>>> hitPoints(height, std::vector<std::pair<vec3, int>>(width));

  int intersects = 0;
  int radValues = 0;

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
      hitPoints[i][j] = {ray.origin + ray.direction * t, triangleId};
    }
  }
  std::vector<std::vector<SpectralValues>> outLuminance(height, std::vector<SpectralValues>(width, SpectralValues(0)));

#pragma omp parallel for shared(height, width, hitPoints, outLuminance, intersects, radValues) default(none)
  for (int i = 0; i < height; ++i) {
    for (int j = 0; j < width; ++j) {
      vec3 hitPoint = hitPoints[i][j].first;
      int triangleId = hitPoints[i][j].second;
      if (triangleId == -1) {
        continue;
      }

      intersects++;

      vec3 N = triangles_[triangleId].getNormal(camera_->origin - hitPoint);
      for (const std::unique_ptr<Light>& light : lights_) {
        outLuminance[i][j] += light->calculateLuminance(hitPoint, N, triangles_, triangleId);
      }

      if (!outLuminance[i][j].isZero()) {
        radValues++;
      }
    }
  }

  std::cout << intersects << " intersects" << std::endl;
  std::cout << radValues << " radValues" << std::endl;

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
  lights_.push_back(std::make_unique<RectangleLight>(origin, intensity, kds, 130, 105));
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
