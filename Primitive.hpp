// Fall 2018

#pragma once

#include <iostream>
#include <fstream>
#include <glm/glm.hpp>
#include <algorithm>
#include "polyroots.hpp"
#include "A4.hpp"

using namespace glm;
class Ray;
class Intersection;
struct Triangle
{
    size_t v1;
    size_t v2;
    size_t v3;

    Triangle( size_t pv1, size_t pv2, size_t pv3 )
    : v1( pv1 )
    , v2( pv2 )
    , v3( pv3 )
    {}
};


class Primitive {
public:
  virtual ~Primitive();
    double can_intersect(Ray &my_ray, vec3 p0,vec3 p1, vec3 p2,bool& is_on);
    virtual Intersection intersect(Ray& ray)=0;
    virtual int my_shape()=0;

    //change
    Material* material = NULL;
    Material* extra_material = NULL;

    int texture_type = -1;
    int material_type = -1;

};


class Sphere : public Primitive {
public:
  virtual ~Sphere();
    Intersection intersect(Ray& ray);
    int my_shape(){
      return 0;
    }
};


class Cube : public Primitive {
public:

  virtual ~Cube();
    Intersection  intersect(Ray& ray);
    std::vector<glm::vec3> m_vertices;
    std::vector<Triangle> m_faces;
    int my_shape(){
      return 1;
    }
};


class NonhierSphere : public Primitive {
public:
  NonhierSphere(const glm::vec3& pos, double radius)
    : m_pos(pos), m_radius(radius)
  {
  }
  virtual ~NonhierSphere();
    Intersection intersect(Ray& ray);
    int my_shape(){
      return 0;
    }

private:
  glm::vec3 m_pos;
  double m_radius;
};



class Cylinder : public Primitive {
public:
  Intersection intersect(Ray& ray);
  virtual ~Cylinder();
  int my_shape(){
    return 2;
  }
};


class Cone : public Primitive {
public:
  Intersection intersect(Ray& ray);
  virtual ~Cone();
  int my_shape(){
    return 3;
  }
};




class NonhierBox : public Primitive {
public:
  NonhierBox(const glm::vec3& pos, double size)
    : m_pos(pos), m_size(size)
  {
  }

  virtual ~NonhierBox();
    Intersection intersect(Ray& ray);
    int my_shape(){
      return 1;
    }


private:
  glm::vec3 m_pos;
  double m_size;
    std::vector<glm::vec3> m_vertices;
    std::vector<Triangle> m_faces;
};
