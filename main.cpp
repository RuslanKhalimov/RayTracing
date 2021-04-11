#include "Scene.h"

int main() {
  Scene scene;

  std::string geometryFile = "../cornel_box0/cornel_box0.shp";
  std::string materialsFile = "../cornel_box0/materials.txt";
  std::string lightsFile = "../cornel_box0/lights.txt";
  std::string cameraFile = "../cornel_box0/camera.txt";

  scene.readSceneFromFiles(geometryFile, materialsFile, lightsFile, cameraFile);
  scene.render("output.txt");

  return 0;
}
