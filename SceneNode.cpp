// Winter 2019

#include "SceneNode.hpp"

#include "cs488-framework/MathUtils.hpp"

#include <iostream>
#include <sstream>
using namespace std;

#include <glm/glm.hpp>
#include <glm/ext.hpp>
#include <glm/gtx/transform.hpp>

using namespace glm;


// Static class variable
unsigned int SceneNode::nodeInstanceCount = 0;


//---------------------------------------------------------------------------------------
SceneNode::SceneNode(const std::string& name)
  : m_name(name),
	m_nodeType(NodeType::SceneNode),
	trans(mat4()),
	invtrans(mat4()),
	m_nodeId(nodeInstanceCount++)
{

}

//---------------------------------------------------------------------------------------
// Deep copy
SceneNode::SceneNode(const SceneNode & other)
	: m_nodeType(other.m_nodeType),
	  m_name(other.m_name),
	  trans(other.trans),
	  invtrans(other.invtrans)
{
	for(SceneNode * child : other.children) {
		this->children.push_front(new SceneNode(*child));
	}
}

//---------------------------------------------------------------------------------------
SceneNode::~SceneNode() {
	for(SceneNode * child : children) {
		delete child;
	}
}

//---------------------------------------------------------------------------------------
void SceneNode::set_transform(const glm::mat4& m) {
	trans = m;
	invtrans = glm::inverse(m);
}

//---------------------------------------------------------------------------------------
const glm::mat4& SceneNode::get_transform() const {
	return trans;
}

//---------------------------------------------------------------------------------------
const glm::mat4& SceneNode::get_inverse() const {
	return invtrans;
}

//---------------------------------------------------------------------------------------
void SceneNode::add_child(SceneNode* child) {
	children.push_back(child);
}

//---------------------------------------------------------------------------------------
void SceneNode::remove_child(SceneNode* child) {
	children.remove(child);
}

//---------------------------------------------------------------------------------------
void SceneNode::rotate(char axis, float angle) {
	vec3 rot_axis;

	switch (axis) {
		case 'x':
			rot_axis = vec3(1,0,0);
			break;
		case 'y':
			rot_axis = vec3(0,1,0);
	        break;
		case 'z':
			rot_axis = vec3(0,0,1);
	        break;
		default:
			break;
	}
	mat4 rot_matrix = glm::rotate(degreesToRadians(angle), rot_axis);
	set_transform( rot_matrix * trans );
}

//---------------------------------------------------------------------------------------
void SceneNode::scale(const glm::vec3 & amount) {
	set_transform( glm::scale(amount) * trans );
}

//---------------------------------------------------------------------------------------
void SceneNode::translate(const glm::vec3& amount) {
	set_transform( glm::translate(amount) * trans );
}


//---------------------------------------------------------------------------------------
int SceneNode::totalSceneNodes() const {
	return nodeInstanceCount;
}

//---------------------------------------------------------------------------------------
std::ostream & operator << (std::ostream & os, const SceneNode & node) {

	//os << "SceneNode:[NodeType: ___, name: ____, id: ____, isSelected: ____, transform: ____"
	switch (node.m_nodeType) {
		case NodeType::SceneNode:
			os << "SceneNode";
			break;
		case NodeType::GeometryNode:
			os << "GeometryNode";
			break;
		case NodeType::JointNode:
			os << "JointNode";
			break;
	}
	os << ":[";

	os << "name:" << node.m_name << ", ";
	os << "id:" << node.m_nodeId;

	os << "]\n";
	return os;
}


Intersection SceneNode::intersect(Ray &ray){
    Ray geo_ray;
    geo_ray.origin = invtrans*ray.origin;
    geo_ray.direction = invtrans*ray.direction;
    intersected = true;

    Intersection best_intersect(ray);

    for(auto child : children){
       // cout << "see child" << endl;
        //if(child->m_nodeType==NodeType::GeometryNode){
            // GeometryNode *geo_node = dynamic_cast<GeometryNode*>(child);*/
             //skip intersected node
             //if(child->intersected)continue;
             if(child->be_subtracted)continue;

             Intersection cur_intersect =child->intersect(geo_ray);
             //if hit see if i have brother that hasn't been hitted
             if(cur_intersect.hit&&child->my_brother_node!=NULL){
               Intersection brother_intersect = child->my_brother_node->intersect(geo_ray);
               if(child->has_intersection){
                 //cout << "here" << endl;
                 if(brother_intersect.hit&&cur_intersect.t<brother_intersect.t&&cur_intersect.t_max>brother_intersect.t&&cur_intersect.t_max>0.001){
                   cur_intersect = brother_intersect;
                 } else if(brother_intersect.hit&&brother_intersect.t<cur_intersect.t&&brother_intersect.t_max>cur_intersect.t&&brother_intersect.t_max>0.001){
                   //cur_intersect = cur_intersect;
                 } else {
                   cur_intersect.hit = false;
                   //continue;
                 }
               } else {
                // cout << "here" << endl;
                //!brother_intersect.hit||
                 if(cur_intersect.t<brother_intersect.t){
                   //cout << "have samller t" << endl;
                   //cur_intersect = cur_intersect;
                 } else if(brother_intersect.hit&&brother_intersect.t_max<cur_intersect.t_max&&cur_intersect.t_max>0.001&&brother_intersect.t_max>0.001){
                   cur_intersect = brother_intersect;
                   cur_intersect.t = cur_intersect.t_max;
                  // cout << "cur t max is " <<  cur_intersect.t_max << endl;
                   //cout << "normal origin is " << to_string(cur_intersect.normal) << endl;
                   cur_intersect.normal = cur_intersect.normal_max;
                   //cout << "after  is " << to_string(cur_intersect.normal_max) << endl;
                 } else{
                   cur_intersect.hit = false;
                   //continue;
                 }
                 /*if (cur_intersect.hit && !brother_intersect.hit){

                 } else {
                   cur_intersect.hit = false;
                 }*/

               }

           }






        //cout << "child name is " << m_name << endl;
        if((cur_intersect.hit&&!best_intersect.hit)||(cur_intersect.hit&&cur_intersect.t < best_intersect.t&&cur_intersect.t>0.001)){
            best_intersect = cur_intersect;
          }



      //  }
      }

    if(best_intersect.hit){
      //cout << "name is " << m_name <<endl;
      //cout << "normal before is " << glm::to_string(best_intersect.normal) << endl;
        auto invtrans_mat3 = glm::dmat3(invtrans);
        auto trans3 = glm::dmat3(trans);
        auto vec3_best = glm::vec3(best_intersect.normal);
        best_intersect.normal = vec4( transpose(invtrans_mat3)*vec3_best,0);
        best_intersect.center = trans3*best_intersect.center;
      //  best_intersect.normal = normalize(best_intersect.normal);
        //cout << "center is " << glm::to_string(best_intersect.center) << endl;
    }

    return best_intersect;
}
