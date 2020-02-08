// Winter 2019

#include <iostream>
#include "Primitive.hpp"
#include <limits>
#include <vector>
#include "math.h"
#include <glm/gtx/normal.hpp>
#include <glm/ext.hpp>

using namespace std;
using namespace glm;

double Primitive::can_intersect(Ray& ray, vec3 p0,vec3 p1, vec3 p2,bool& is_on){
    vec3 a = vec3(ray.origin);
    vec3 b = a + vec3(ray.direction);
    float R1 = a.x-p0.x;
    float R2 = a.y-p0.y;
    float R3 = a.z-p0.z;
    vec3 my_x = vec3(p1.x-p0.x,p2.x-p0.x,a.x-b.x);
    vec3 my_y = vec3(p1.y-p0.y,p2.y-p0.y,a.y-b.y);
    vec3 my_z = vec3(p1.z-p0.z,p2.z-p0.z,a.z-b.z);

    mat3 D = mat3(my_x,my_y,my_z);


    mat3 D1 = mat3(vec3(R1,R2,R3),vec3(p2.x-p0.x,p2.y-p0.y,p2.z-p0.z),vec3(a.x-b.x,a.y-b.y,a.z-b.z));
    //D1 = transpose(D1);
    mat3 D2 = mat3(vec3(p1.x-p0.x,p1.y-p0.y,p1.z-p0.z),vec3(R1,R2,R3),vec3(a.x-b.x,a.y-b.y,a.z-b.z));
    //D2 = transpose(D2);
    mat3 D3 = mat3(vec3(p1.x-p0.x,p1.y-p0.y,p1.z-p0.z),vec3(p2.x-p0.x,p2.y-p0.y,p2.z-p0.z),vec3(R1,R2,R3));
    //D3 = transpose(D3);

    double beta = determinant(D1)/determinant(D);
    double gamma = determinant(D2)/determinant(D);
    double t = determinant(D3)/determinant(D);

    if(beta>=0&&gamma>=0&&(beta+gamma)<=1&&t>=0.05){
        is_on = true;
    }else{
        is_on = false;
    }
    // cout << "t is " << t << endl;
    return t;

}


Intersection Sphere::intersect(Ray& ray){
  Intersection cur_sec(ray);
    //cout << "here sphere " <<  endl;
    vec3 center(0,0,0);
    vec3 origin = vec3(ray.origin);
    vec3 termination = origin + vec3(ray.direction);
    double A = dot(termination-origin,termination-origin);
    double B = dot(termination-origin,origin-center)*2;
    double C = dot(origin-center,origin-center)-1;

    double roots[2];
    size_t roots_n = quadraticRoots(A,B,C,roots);

    switch (roots_n) {
        case 0:
            cur_sec.hit = false;
            break;
        case 1:
        if(roots[0]>0.001){
            cur_sec.hit = true;
            cur_sec.t = roots[0];
          } else {
            cur_sec.hit = false;
          }
          break;
        default:
        double mini = std::min(roots[0],roots[1]);
        double maxi = std::max(roots[0],roots[1]);
            if(mini>0.001){
                cur_sec.t = mini;
                cur_sec.t_max = maxi;
                cur_sec.hit = true;
                vec3 tmp_normal_max = origin + float(cur_sec.t_max)*(termination-origin);
                cur_sec.normal_max = vec4(tmp_normal_max,0);
            } else {
                cur_sec.hit = true;
                if(maxi>0.001){
                    cur_sec.t = maxi;
                } else {
                    cur_sec.hit = false;
                }
            }
            break;
    }
   // cout << "t is " << cur_sec.t << endl;

    if(cur_sec.hit){
        vec3 tmp_normal = origin + float(cur_sec.t)*(termination-origin);
        cur_sec.normal = vec4(tmp_normal,0);
        //cout << "normal is " << glm::to_string(cur_sec.normal) << endl;
        cur_sec.center = center;
    }
    return cur_sec;
  }


//from woo4.me
Intersection Cylinder::intersect(Ray& ray){
    Intersection cur_sec(ray);
    vec3 m_pos(0,0,0);

    double temp_t = std::numeric_limits<double>::infinity();
    double cur_t = 0;
    vec3 temp_normal(0,0,0);

    vec3 E = vec3(ray.origin);
    vec3 D = vec3(ray.direction);

    double A = D.x*D.x + D.y*D.y;
    double B = 2*D.x*E.x + 2*D.y*E.y;
    double C = E.x*E.x + E.y*E.y - 1;

    double roots[2];
    size_t roots_n = quadraticRoots(A,B,C,roots);

    switch (roots_n) {
      case 0:
          cur_sec.hit = false;
          break;
      case 1:
      if(roots[0]>0.001){
          cur_sec.hit = true;
          vec3 tmp_normal = E + float(roots[0])*D;
          cur_sec.normal = vec4(tmp_normal.x,tmp_normal.y,0,0);
          cur_sec.t = roots[0];
        } else {
          cur_sec.hit = false;
        }
        break;
      default:
        double t0 = std::min(roots[0],roots[1]);
        double t1 = std::max(roots[0],roots[1]);
        float y0 = E.z + t0*D.z;
        float y1 = E.z + t1*D.z;

        if(y0<-1){
          if(y1<-1){
            cur_sec.hit = false;
          }else{
            //hit the cap
            float t = t0 + (t1-t0)*(y0+1)/(y0-y1);
            if(t<=0.001){
              cur_sec.hit = false;
              break;
            } else {
              cur_sec.t = t;
              cur_sec.hit = true;
              cur_sec.normal = vec4(0,-1,0,0);
              break;
            }
          }
        } else if(y0>=-1&&y0<=1){
          //hit the cylinder bit
          if(t0<=0.001){
            cur_sec.hit = false;
            break;
          }else{
            cur_sec.t = t0;
            vec3 tmp_normal = E + t0*D;
            cur_sec.normal = vec4(tmp_normal.x,tmp_normal.y,0,0);
            cur_sec.hit = true;
          }
        } else if(y0>1){
          if(y1>1){
            cur_sec.hit = false;
            break;
          }else{
            //hit the cap
            float t = t0 + (t1-t0)*(y0-1)/(y0-y1);
            if(t<=0.001){
              cur_sec.hit = false;
              break;
            }else{
              vec3 tmp_normal = E + t*D;
              cur_sec.t = t;
              cur_sec.hit = true;
              cur_sec.normal = vec4(0,1,0,0);
            }
          }
        }else{
          cur_sec.hit = false;
        }
    }
    cur_sec.center = m_pos;
    return cur_sec;
}


Intersection Cone::intersect(Ray& ray){
  Intersection cur_sec(ray);
  vec3 m_pos(0,0,0);
  double angle = 60;

  vec3 E = vec3(ray.origin);
  vec3 D = vec3(ray.direction);

  double A = D.x*D.x + D.y*D.y - D.z*D.z;
  double B = 2*D.x*E.x + 2*D.y*E.y - 2*D.z*E.z;
  double C = E.x*E.x + E.y*E.y - E.z*E.z;

  double roots[2];
  size_t roots_n = quadraticRoots(A,B,C,roots);

  //cout << "roots 0 is " << roots[0] << " roots 1 is " << roots[1] << endl;

  switch (roots_n) {
    case 0:
        cur_sec.hit = false;
        break;
    case 1:
    if(roots[0]>0.001){
        cur_sec.hit = true;
        vec3 tmp_normal = E + float(roots[0])*D;
        cur_sec.normal = vec4(tmp_normal.x,tmp_normal.y,-tmp_normal.z,0);
        cur_sec.t = roots[0];
      } else {
        cur_sec.hit = false;
      }
      break;
    default:
    double t0 = std::min(roots[0],roots[1]);
    double t1 = std::max(roots[0],roots[1]);
    float y0 = E.z + t0*D.z;
    float y1 = E.z + t1*D.z;

    if(y0<0){
      cur_sec.hit = false;
    } else if(y0>=0&&y0<=1){
      //hit the cylinder bit
      if(t0<=0.001){
        cur_sec.hit = false;
        break;
      }else{
        cur_sec.t = t0;
        vec3 tmp_normal = E + t0*D;
        cur_sec.normal = vec4(tmp_normal.x,tmp_normal.y,-tmp_normal.z,0);
        cur_sec.hit = true;
      }
    } else if(y0>1){
      if(y1>1){
        cur_sec.hit = false;
        break;
      }else{
        //hit the cap
        float t = t0 + (t1-t0)*(y0-1)/(y0-y1);
        if(t<=0.001){
          cur_sec.hit = false;
          break;
        }else{
          vec3 tmp_normal = E + t*D;
          cur_sec.t = t;
          cur_sec.hit = true;
          cur_sec.normal = vec4(0,0,1,0);
        }
      }
    }else{
      cur_sec.hit = false;
    }
}
return cur_sec;
}


Intersection Cube::intersect(Ray& ray){
  Intersection cur_sec(ray);
    vec3 m_pos(0,0,0);

    bool is_on;
    double temp_t = std::numeric_limits<double>::infinity();
    double temp_t_max = 0;
    vec3 temp_normal_max(0,0,0);
    double cur_t = 0;
    vec3 temp_normal(0,0,0);
    m_vertices = {
        m_pos +vec3(0,0,1),
        m_pos +vec3(0,0,0),
        m_pos +vec3(1,0,0),
        m_pos +vec3(1,0,1),
        m_pos +vec3(0,1,1),
        m_pos +vec3(0,1,0),
        m_pos +vec3(1,1,0),
        m_pos +vec3(1,1,1),
    };

    m_faces = {
        Triangle(4,5,1),
        Triangle(5,6,2),
        Triangle(6,7,3),
        Triangle(4,0,7),
        Triangle(0,1,2),
        Triangle(7,6,5),
        Triangle(0,4,1),
        Triangle(1,5,2),
        Triangle(2,6,3),
        Triangle(7,0,3),
        Triangle(3,0,2),
        Triangle(4,7,5)
    };

    for (auto triangle : m_faces) {

        cur_t = can_intersect(ray,m_vertices[triangle.v1], m_vertices[triangle.v2], m_vertices[triangle.v3],is_on);
        //cout << "cur t is " << cur_t << endl;
        if(is_on){
            cur_sec.hit = true;
            if(cur_t>temp_t_max){
              temp_t_max = cur_t;
              temp_normal_max = triangleNormal(m_vertices[triangle.v1], m_vertices[triangle.v2], m_vertices[triangle.v3]);

            }
            if(cur_t<temp_t){
                temp_t = cur_t;
                temp_normal = triangleNormal(m_vertices[triangle.v1], m_vertices[triangle.v2], m_vertices[triangle.v3]);
            }
        }
    }

    cur_sec.local_normal = vec4(temp_normal,0);
    cur_sec.local_pos = ray.origin + temp_t*ray.direction;
    cur_sec.t = temp_t;
    cur_sec.t_max = temp_t_max;
    cur_sec.normal = vec4(temp_normal,0);
    cur_sec.normal_max = vec4(temp_normal_max,0);
    return cur_sec;
}

/*Intersection Cube::intersect(Ray& ray){
  Intersection cur_sec(ray);
    vec3 m_pos(0,0,0);

    bool is_on;
    double temp_t = std::numeric_limits<double>::infinity();
    double cur_t = 0;
    vec3 temp_normal(0,0,0);
    m_vertices = {
        m_pos +vec3(0,0,1),
        m_pos +vec3(0,0,0),
        m_pos +vec3(1,0,0),
        m_pos +vec3(1,0,1),
        m_pos +vec3(0,1,1),
        m_pos +vec3(0,1,0),
        m_pos +vec3(1,1,0),
        m_pos +vec3(1,1,1),
    };

    m_faces = {
        Triangle(4,5,1),
        Triangle(5,6,2),
        Triangle(6,7,3),
        Triangle(4,0,7),
        Triangle(0,1,2),
        Triangle(7,6,5),
        Triangle(0,4,1),
        Triangle(1,5,2),
        Triangle(2,6,3),
        Triangle(7,0,3),
        Triangle(3,0,2),
        Triangle(4,7,5)
    };

    for (auto triangle : m_faces) {

        cur_t = can_intersect(ray,m_vertices[triangle.v1], m_vertices[triangle.v2], m_vertices[triangle.v3],is_on);
        //cout << "cur t is " << cur_t << endl;
        if(is_on){
            cur_sec.hit = true;
            if(cur_t<temp_t){
                temp_t = cur_t;
                temp_normal = triangleNormal(m_vertices[triangle.v1], m_vertices[triangle.v2], m_vertices[triangle.v3]);
            }
        }
    }

    cur_sec.local_normal = vec4(temp_normal,0);
    cur_sec.local_pos = ray.origin + temp_t*ray.direction;
    cur_sec.t = temp_t;
    cur_sec.normal = vec4(temp_normal,0);
    return cur_sec;
}*/


Intersection NonhierSphere::intersect(Ray& ray){
    Intersection cur_sec(ray);
    vec3 center = m_pos;
    vec3 origin = vec3(ray.origin);
    vec3 termination = origin + vec3(ray.direction);
    double A = dot(termination-origin,termination-origin);
    double B = dot(termination-origin,origin-center)*2;
    double C = dot(origin-center,origin-center)-m_radius*m_radius;

    double roots[2];
    size_t roots_n = quadraticRoots(A,B,C,roots);
    switch (roots_n) {
        case 0:
            cur_sec.hit = false;
            break;
        case 1:
        if(roots[0]>0.05){
            cur_sec.hit = true;
            cur_sec.t = roots[0];
          } else {
            cur_sec.hit = false;
          }
          break;
        default:
        double mini = std::min(roots[0],roots[1]);
        double maxi = std::max(roots[0],roots[1]);
            if(mini>0.05){
                cur_sec.t = mini;
                cur_sec.t_max = maxi;
                vec3 tmp_normal_max = origin + float(cur_sec.t_max)*(termination-origin);
                cur_sec.normal_max = vec4(tmp_normal_max,0);
                cur_sec.hit = true;
            } else {
                cur_sec.hit = true;
                if(maxi>0.05){
                    cur_sec.t = maxi;
                } else {
                    cur_sec.hit = false;
                }
            }
            break;
    }

    if(cur_sec.hit){
        vec3 tmp_normal = origin + float(cur_sec.t)*(termination-origin)-m_pos;
        cur_sec.normal = vec4(tmp_normal,0);
        //cout << "normal is " << glm::to_string(cur_sec.normal) << endl;
        cur_sec.center = center;
    }
    return cur_sec;
}





Intersection NonhierBox::intersect(Ray& ray){
  Intersection cur_sec(ray);
    bool is_on;
    double temp_t = std::numeric_limits<double>::infinity();
    double temp_t_max = 0;
    vec3 temp_normal_max(0,0,0);
    double cur_t = 0;
    vec3 temp_normal(0,0,0);
    m_vertices = {
        m_pos +vec3(0,0,m_size),
        m_pos +vec3(0,0,0),
        m_pos +vec3(m_size,0,0),
        m_pos +vec3(m_size,0,m_size),
        m_pos +vec3(0,m_size,m_size),
        m_pos +vec3(0,m_size,0),
        m_pos +vec3(m_size,m_size,0),
        m_pos +vec3(m_size,m_size,m_size),
    };

    m_faces = {
        Triangle(4,5,1),
        Triangle(5,6,2),
        Triangle(6,7,3),
        Triangle(4,0,7),
        Triangle(0,1,2),
        Triangle(7,6,5),
        Triangle(0,4,1),
        Triangle(1,5,2),
        Triangle(2,6,3),
        Triangle(7,0,3),
        Triangle(3,0,2),
        Triangle(4,7,5)
    };
    for (auto triangle : m_faces) {

        cur_t = can_intersect(ray,m_vertices[triangle.v1], m_vertices[triangle.v2], m_vertices[triangle.v3],is_on);
        //cout << "cur t is " << cur_t << endl;
        if(is_on){
            cur_sec.hit = true;
            if(cur_t>temp_t_max){
              temp_t_max = cur_t;
              temp_normal_max = triangleNormal(m_vertices[triangle.v1], m_vertices[triangle.v2], m_vertices[triangle.v3]);

            }
            if(cur_t<temp_t){
                temp_t = cur_t;
                temp_normal = triangleNormal(m_vertices[triangle.v1], m_vertices[triangle.v2], m_vertices[triangle.v3]);
            }
        }
    }
    if(cur_sec.hit){
      cur_sec.local_normal = vec4(temp_normal,0);
      cur_sec.local_pos = ray.origin + temp_t*ray.direction;
      cur_sec.t = temp_t;
      cur_sec.t_max = temp_t_max;
      cur_sec.normal = vec4(temp_normal,0);
      cur_sec.normal_max = vec4(temp_normal_max,0);
      cout << "normal is " << glm::to_string(cur_sec.normal) << endl;
    }
    return cur_sec;
}







Primitive::~Primitive()
{
}

Sphere::~Sphere()
{
}

Cube::~Cube()
{
}

Cylinder::~Cylinder(){}

Cone::~Cone(){}

NonhierSphere::~NonhierSphere()
{
}

NonhierBox::~NonhierBox()
{
}
