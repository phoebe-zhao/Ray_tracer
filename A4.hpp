// Winter 2019

#pragma once

#include <glm/glm.hpp>
#include "PhongMaterial.hpp"
#include "SceneNode.hpp"
#include "Light.hpp"
#include "Image.hpp"
#include <vector>
#include <limits>
#include <random>
#include <cstdlib>
#include <utility>

using namespace std;
using namespace glm;

class Ray;
class Intersection;

//#define REFLECTION
//#define ANTIALIASING
//#define SOFTSHADOW
//#define DEPTHOFFIELD
//#define REFRACTION
//#define GLOSSY_REFLECTION

void A4_Render(
		// What to render
		SceneNode * root,

		// Image to write to, set to a given width and height
		Image & image,

		// Viewing parameters
		const glm::vec3 & eye,
		const glm::vec3 & view,
		const glm::vec3 & up,
		double fovy,

		// Lighting parameters
		const glm::vec3 & ambient,
		const std::list<Light *> & lights
);

//vec3 Ray_color(Ray& my_ray,SceneNode *root,const list<Light*> &lights,vec4 look_from,vec3 ambient,int maxhit,glm::vec3,vec3&,bool&);
vec3 uv_cal(Image image ,double ,double );
