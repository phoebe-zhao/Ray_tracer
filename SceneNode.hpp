// Fall 2018

#pragma once

#include "Material.hpp"
#include "PhongMaterial.hpp"
#include <glm/glm.hpp>

#include <list>
#include <string>
#include <iostream>

class Primitive;
class GeometryNode;

enum class NodeType {
	SceneNode,
	GeometryNode,
	JointNode
};

class Ray{
public:
    glm::vec4 origin;
    glm::vec4 direction;
};

class Intersection {
public:
    Ray my_ray;
    bool hit = false;
    double t ;
		double t_max = -1;

    PhongMaterial *my_phong;
		Primitive* m_primitive;
    glm::vec4 normal;
		glm::vec4 normal_max;
		bool extra_image = false;
		bool has_image = false;
		int texture_type = -1;
		int shape = -1;
		glm::vec3 center;
		glm::vec4 local_normal = glm::vec4(0,0,0,0);
		glm::vec4 local_pos = glm::vec4(0,0,0,0);
		std::string node_name = "";
		int cur_glossy = -1;



    Intersection(Ray &r):my_ray(r),hit(false){};

};


class SceneNode {
public:
    SceneNode(const std::string & name);

	SceneNode(const SceneNode & other);

    virtual ~SceneNode();

	int totalSceneNodes() const;

    const glm::mat4& get_transform() const;
    const glm::mat4& get_inverse() const;

    void set_transform(const glm::mat4& m);

    void add_child(SceneNode* child);

    void remove_child(SceneNode* child);

	//-- Transformations:
    void rotate(char axis, float angle);
    void scale(const glm::vec3& amount);
    void translate(const glm::vec3& amount);


	friend std::ostream & operator << (std::ostream & os, const SceneNode & node);

    // Transformations
    glm::mat4 trans;
    glm::mat4 invtrans;

    std::list<SceneNode*> children;

	NodeType m_nodeType;
	std::string m_name;
	unsigned int m_nodeId;

	virtual Intersection intersect(Ray &ray);
	SceneNode* my_brother_node = NULL;
	bool intersected = false;
	bool has_intersection = false;
	bool has_difference = false;
	bool be_subtracted = false;

private:
	// The number of SceneNode instances.
	static unsigned int nodeInstanceCount;
};
