// Winter 2019

#include <iostream>
#include <fstream>

#include <glm/ext.hpp>

// #include "cs488-framework/ObjFileDecoder.hpp"
#include "Mesh.hpp"

Mesh::Mesh( const std::string& fname )
	: m_vertices()
	, m_faces()
{
	std::string code;
	double vx, vy, vz;
	size_t s1, s2, s3;

	std::ifstream ifs( string("Assets/") + fname.c_str() );
	while( ifs >> code ) {
		if( code == "v" ) {
			ifs >> vx >> vy >> vz;
#ifdef RENDER_BOUNDING_VOLUMES
            if(vx < min_x){
                min_x = vx;
            }
            if(vx > max_x){
                max_x = vx;
            }
            if(vy < min_y){
                min_y = vy;
            }
            if(vy > max_y){
                max_y = vy;
            }
            if(vz < min_z){
                min_z = vz;
            }
            if(vz > max_z){
                max_z = vz;
            }
#endif
			m_vertices.push_back( glm::vec3( vx, vy, vz ) );
		} else if( code == "f" ) {
			ifs >> s1 >> s2 >> s3;
			m_faces.push_back( Triangle( s1 - 1, s2 - 1, s3 - 1 ) );
		}
	}

#ifdef RENDER_BOUNDING_VOLUMES
    m_vertices.clear();
    m_faces.clear();

    vec3 m_pos(min_x,min_y,min_z);

    double x_diff = max_x-min_x;
    double y_diff = max_y-min_y;
    double z_diff = max_z-min_z;
    //vec3 m_size(x_diff,y_diff,z_diff);


    m_vertices = {
        m_pos,
        m_pos +vec3(x_diff,0,0),
        m_pos +vec3(x_diff,0,z_diff),
        m_pos +vec3(0,0,z_diff),
        m_pos +vec3(0,y_diff,0),
        m_pos +vec3(x_diff,y_diff,0),
        m_pos +vec3(x_diff,y_diff,z_diff),
        m_pos +vec3(0,y_diff,z_diff)
    };

    m_faces = {
        Triangle(0,1,2),
        Triangle(0,2,3),
        Triangle(4,5,6),
        Triangle(4,6,7),
        Triangle(6,5,1),
        Triangle(6,1,2),
        Triangle(7,4,0),
        Triangle(7,0,3),
        Triangle(4,5,1),
        Triangle(4,1,0),
        Triangle(7,6,2),
        Triangle(7,2,3)
    };



#else
    // blah
#endif
}


Intersection Mesh::intersect(Ray& ray){
	  Intersection cur_sec(ray);
    bool is_on;
    double temp_t = std::numeric_limits<double>::infinity();
    double cur_t = 0;
    vec3 temp_normal(0,0,0);

    //cout << "cur t is " << cur_t << endl;
    //cout << "m_faces size is " << m_faces.size() << endl;
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
    cur_sec.t = temp_t;
    cur_sec.normal = vec4(temp_normal,0);
		return cur_sec;
    //cout << "cur t is " << cur_sec.t << endl;
     //cout << "cur normal is " << cur_sec.normal.x << ' ' << cur_sec.normal.y <<' ' <<  cur_sec.normal.z  << endl;
}


std::ostream& operator<<(std::ostream& out, const Mesh& mesh)
{
  out << "mesh {";
  /*

  for( size_t idx = 0; idx < mesh.m_verts.size(); ++idx ) {
  	const MeshVertex& v = mesh.m_verts[idx];
  	out << glm::to_string( v.m_position );
	if( mesh.m_have_norm ) {
  	  out << " / " << glm::to_string( v.m_normal );
	}
	if( mesh.m_have_uv ) {
  	  out << " / " << glm::to_string( v.m_uv );
	}
  }

*/
  out << "}";
  return out;
}
