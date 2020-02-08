// Winter 2019

#pragma once

#include <glm/glm.hpp>

#include "Material.hpp"
using namespace glm;

class PhongMaterial : public Material {
public:
  PhongMaterial(const glm::vec3& kd, const glm::vec3& ks, double shininess);
  PhongMaterial(const glm::vec3& kd, const glm::vec3& ks, double shininess,bool,float);
  virtual ~PhongMaterial();
    vec3 get_kd(){
        return m_kd;
    }
    vec3 get_ks(){
        return m_ks;
    }
    double get_shine(){
        return m_shininess;
    }
    bool m_dielectric;
    float m_ior;
private:
  glm::vec3 m_kd;
  glm::vec3 m_ks;

  double m_shininess;
};
