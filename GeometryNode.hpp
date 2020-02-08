// Winter 2019

#pragma once

#include "SceneNode.hpp"
#include "Primitive.hpp"
#include "Material.hpp"



class GeometryNode : public SceneNode {
public:
	GeometryNode( const std::string & name, Primitive *prim,
		Material *mat = nullptr );

	void setMaterial( Material *material );
	void set_texture(const std::string & file);
	void set_bump(const std::string & file);
	void set_intersection(const std::string& another_node);
	void set_difference(const std::string& another_node);
	void set_glossy(int num);
	Intersection intersect(Ray& ray);

	Material *m_material;
	Primitive *m_primitive;
	Image image;
	bool extra_image;
	Image extra_texture;
	bool has_image = false;
	int image_type = -1;

	//bool has_intersection = false;
//	bool has_difference = false;
	string my_brother;
	int cur_glossy = -1;
	//GeometryNode* my_brother_node = NULL;
	//bool intersected = false;
};
