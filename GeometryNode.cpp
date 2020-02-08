// Winter 2019

#include "GeometryNode.hpp"
#include <glm/glm.hpp>
#include <glm/ext.hpp>

//---------------------------------------------------------------------------------------
GeometryNode::GeometryNode(
	const std::string & name, Primitive *prim, Material *mat )
	: SceneNode( name )
	, m_material( mat )
	, m_primitive( prim )
{
	m_nodeType = NodeType::GeometryNode;
}

void GeometryNode::setMaterial( Material *mat )
{
	// Obviously, there's a potential memory leak here.  A good solution
	// would be to use some kind of reference counting, as in the
	// C++ shared_ptr.  But I'm going to punt on that problem here.
	// Why?  Two reasons:
	// (a) In practice we expect the scene to be constructed exactly
	//     once.  There's no reason to believe that materials will be
	//     repeatedly overwritten in a GeometryNode.
	// (b) A ray tracer is a program in which you compute once, and
	//     throw away all your data.  A memory leak won't build up and
	//     crash the program.

	m_material = mat;
}


Intersection GeometryNode::intersect(Ray& ray){
    Ray geo_ray;
    geo_ray.origin = invtrans*ray.origin;
    geo_ray.direction = invtrans*ray.direction;
    Intersection best_intersect = m_primitive->intersect(geo_ray);
		best_intersect.node_name = m_name;
		//intersected = true;

    if(best_intersect.hit){
        best_intersect.my_phong = static_cast<PhongMaterial*>(m_material);
				if(has_image){
						best_intersect.has_image = true;
						best_intersect.texture_type = image_type;
						best_intersect.shape = m_primitive->my_shape();
						//change
						best_intersect.m_primitive = m_primitive;
						best_intersect.m_primitive->texture_type = image_type;

						if(extra_image){
							best_intersect.extra_image = true;
						}
						if(cur_glossy>0){
							best_intersect.cur_glossy = cur_glossy;
						}
					}
    }

    for(auto child : children){
		 	//cout << "child name is " << m_name << endl;
        Intersection cur_intersect = child->intersect(geo_ray);
        if((cur_intersect.hit&&!best_intersect.hit)||(cur_intersect.hit&&cur_intersect.t < best_intersect.t&&cur_intersect.t>0.001)){
            best_intersect = cur_intersect;
        }
    }

    if(best_intersect.hit){
				//cout << "normal before is " << glm::to_string(best_intersect.normal) << endl;
				auto invtrans3 = glm::dmat3(invtrans);
				auto trans3 = glm::dmat3(trans);
				auto vec3_best = glm::vec3(best_intersect.normal);
				best_intersect.normal = vec4( transpose(invtrans3)*vec3_best,0);
				best_intersect.center = trans3*best_intersect.center;
				//cout << "after is " << glm::to_string(best_intersect.normal) << endl;
    }
    return best_intersect;

}


void GeometryNode::set_texture(const std::string & file){
	if(m_primitive->material_type==-1){
		Material* cur_type =  new Texture(file);
		m_primitive->material = cur_type;
		m_primitive->material_type = 0;
		image_type = 0;
	} else {
		Material* cur_type = new Texture(file);
		m_primitive->extra_material = cur_type;
		extra_image = true;
	}
	has_image= true;
}

	void  GeometryNode::set_bump(const std::string & file){
		Material* cur_type = new Bump(file);
		m_primitive->material = cur_type;
		m_primitive->material_type = 1;
		has_image= true;
		image_type = 1;
	}


void GeometryNode::set_intersection(const std::string& another_node){
	has_intersection = true;
	my_brother = another_node;
}


void GeometryNode::set_difference(const std::string& another_node){
	has_difference = true;
	my_brother = another_node;
}


void GeometryNode::set_glossy(int num){
	cur_glossy = num;
}
