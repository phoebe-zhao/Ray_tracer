// Winter 2019

#pragma once
#include <glm/glm.hpp>
#include "Image.hpp"

class Material {
public:
    virtual ~Material();
    /*virtual glm::vec3 get_kd() = 0;
    virtual glm::vec3 get_ks() = 0;
    virtual double get_shine() = 0;*/

protected:
  Material();
};


class Texture : public Material {
public:
  Image image;
  Texture(const std::string &filename){
    image.loadPng(filename);
  }
  glm::vec3 uv_cal(double u,double v);
};


class Bump : public Material {
public:
  Image image;
  Bump(const std::string &filename){
    image.loadPng(filename);
  }
    glm::vec3 uv_cal(double u,double v);
};
