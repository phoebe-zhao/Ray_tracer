// Winter 2019

#include <glm/ext.hpp>
#include <math.h>
#include <algorithm>
#include "A4.hpp"
#include "GeometryNode.hpp"
#include <glm/gtx/vector_angle.hpp>
#include <thread>
#include <chrono>
#include <iomanip>
#include <cstdlib>
using namespace std;
using namespace glm;

#define RAYNUM 100
#define GLOSSYNUM 5


std::random_device rd;
std::mt19937 eng(rd()); // seed the generator
std::uniform_int_distribution<> distr(0, 100); // define the range
std::uniform_int_distribution<> distr_soft(-6, 6); // define the range
std::uniform_int_distribution<> distr_one(-1, 1);
std::uniform_real_distribution<float> distribution(0.0, 1.0);

std::default_random_engine generator;


int pixels_done = 0;
int total_pixels;

void print_mat4(mat4 matrix){
  for(int i = 0;i<4;++i){
    for(int j = 0;j<4;++j){
      cout << matrix[i][j] << ' ';
    }
    cout << endl;
  }
}


int smallest(float x, float y, float z){
  x = abs(x);
  y = abs(y);
  z = abs(z);
  if(x<=y){
    if(x<=z){
      return 0;
    }else {
      return 2;
    }
  } else {
    if(x<=z){
      return 1;
    } else if(y<=z){
      return 1;
    } else {
      return 2;
    }
  }
}


SceneNode* recursive_find(SceneNode* root, string name){
    GeometryNode* temp = NULL;
    for(auto node : root->children){
      GeometryNode *geo_node = static_cast<GeometryNode*>(node);
        if(geo_node->m_name==name){
          return geo_node;
        }
    }
    return temp;
}


void find_brother(SceneNode* root){
  for(auto node : root->children){
    GeometryNode *geo_node = static_cast<GeometryNode*>(node);
    if(geo_node->has_intersection){
      geo_node->my_brother_node=recursive_find(root,geo_node->my_brother);
      //geo_node->my_brother_node->my_brother_node = geo_node;
      geo_node->my_brother_node->be_subtracted = true;
    } else if(geo_node->has_difference){
      geo_node->my_brother_node=recursive_find(root,geo_node->my_brother);
      geo_node->my_brother_node->be_subtracted = true;
    }
  }
}


float r_cal(Image image ,double u,double v){
  vec3 col(0,0,0);
  double di = (image.width() - 1)*u;
	double dj = (image.height() - 1)*v;
  int i = int(di);
  double up = di-i;
  int j = int(dj);
  double vp = dj-j;
  //cout << "image i j 0 is " << image(i,j,0) << endl;
  double col00 = image(i,j,0);
  double col01 = image(i+1,j,0);
  double col10 = image(i,j+1,0);
  double col11 = image(i+1,j+1,0);

  float result_1 = col00+up*(col01-col00);
  float result_2 = col10+up*(col11-col10);
  float result = result_1+vp*(result_2-result_1);
  //cout << "result is " << result << endl;
  return result;
}


pair<float,float> get_derivatives(Image image ,double u,double v){
  float Bu, Bv;
  double delta = 0.002;

  float col10 = r_cal(image,u+delta,v);
  float col_10 = r_cal(image,u-delta,v);
  float col01 = r_cal(image,u,v+delta);
  float col0_1 = r_cal(image,u,v-delta);

  Bu = (col10-col_10)/(2*delta);
  Bv = (col01-col0_1)/(2*delta);

  float p1 = col10-col_10;
  float p2 = col01-col0_1;
  //cout<< "\n col10 is " <<  col10 << " col_10 is " << col_10 << endl;

  //cout << " bu is " <<   Bu << " bv is " << Bv << endl;

  if(p1==0){
    Bu = 0;
  }
  if(p2==0){
    Bv = 0;
  }
  pair<float,float> result = make_pair(Bu,Bv);

  return result;
}


vector<pair<double,double> >give_random_2d(int size){
    vector<pair<double,double> > my_arr;
    for(int i = 0;i<size;++i){
        for(int j = 0;j<size;++j){
            double random_1 = ((double) rand() / (RAND_MAX));
            double random_2 = ((double) rand() / (RAND_MAX));
            pair<double,double> my_pair = make_pair(random_1,random_2);
            my_arr.push_back(my_pair);
        }
    }
    return my_arr;
}

Intersection check_intersect(Ray &ray, SceneNode* node){
    /*Ray geo_ray;
    geo_ray.origin = node->invtrans*ray.origin;
    geo_ray.direction = node->invtrans*ray.direction;
    Intersection temp_intersect(geo_ray);
    Intersection best_intersect(geo_ray);

    if(node->m_nodeType==NodeType::GeometryNode){
        GeometryNode *geo_node = static_cast<GeometryNode*>(node);
        geo_node->m_primitive->intersect(geo_ray,best_intersect);
        best_intersect.my_phong = static_cast<PhongMaterial*>(geo_node->m_material);
        best_intersect.normal = transpose(node->invtrans)*best_intersect.normal;
    }

    for(SceneNode *child: node->children) {
        Intersection cur_intersect = check_intersect(geo_ray,child);
        if((cur_intersect.hit&&!best_intersect.hit)||(cur_intersect.t>0.001&&cur_intersect.t < best_intersect.t)){
            best_intersect = cur_intersect;
            //best_intersect.normal = transpose(node->invtrans)*best_intersect.normal;
        }
    }
    return best_intersect;*/
    return node->intersect(ray);
}




//delete later caution
/*void freshnel(vec4 direction,vec4 normal,float ior,float& kr){
  float dot_val = dot(direction,normal);
  float cosi = clamp(float(-1.0),float(1.0),float(dot_val));
  float etai = 1;
  float etat = ior;
  if(cosi>0){
    std::swap(etai,etat);
  }
  float sint = etai/etat*sqrtf(std::max(float(0), 1 - cosi * cosi));
  if(sint>=1){
    kr = 1;
  }else{
    float cost = sqrtf(std::max(float(0), 1 - sint * sint));
    cosi = fabsf(cosi);
    float Rs = ((etat * cosi) - (etai * cost)) / ((etat * cosi) + (etai * cost));
    float Rp = ((etai * cosi) - (etat * cost)) / ((etai * cosi) + (etat * cost));
    kr = (Rs * Rs + Rp * Rp) / 2;
  }
}


vec4 refract(vec4 direction,vec4 normal,float ior){
  float dot_val = dot(direction,normal);
  float cosi = clamp(float(-1.0),float(1.0),float(dot_val));
  float etai = 1;
  float etat = ior;
  vec4 n = normal;
  if(cosi<0){
    cosi = -cosi;
  }else{
    std::swap(etai,etat);
    n = -normal;
  }
  float eta = etai/etat;
  float k = 1-eta *eta *(1-cosi*cosi);
  vec4 alt = eta*direction+(eta*cosi-sqrtf(k))*n;
  vec4 standard = vec4(0,0,0,0);
  return k<0 ? standard : alt;
}*/

/*void perturbed_ray(Ray& ray,Ray& return_ray){
  float a = 10;
  double random_1 = ((double) rand() / (RAND_MAX));
  double random_2 = ((double) rand() / (RAND_MAX));
  double random_3 = ((double) rand() / (RAND_MAX));
  double random_4 = ((double) rand() / (RAND_MAX));

  vec2 e1(random_1,random_2);
  vec2 e2(random_3,random_4);

  vec2 u = -a/2 + e1*a;
  vec2 v = -a/2 + e2*a;

  vec4 w = normalize(ray.direction);

  vec4 t = w;
  //t.z = 1;
  int index = smallest(t.x,t.y,t.z);
  if(index==0){
    t.x = 1;
  } else if(index==1){
    t.y = 1;
  } else {
    t.z = 1;
  }


  vec3 U = normalize(cross(vec3(t),vec3(w)));
  //cout << "u is " << to_string(U) << endl;
  vec3 V = cross(vec3(w),vec3(U));

  vec4 mult_u = vec4(u,0,1);
  vec4 mult_v = vec4(v,0,1);
  vec4 mult_U = vec4(U,1);
  vec4 mult_V = vec4(V,1);


  //Ray another_r;
  return_ray.origin  = ray.origin;
  return_ray.direction =  ray.direction +  mult_u*mult_U + mult_v*mult_V;
  //cout << "final direction is " << to_string(return_ray.direction) << endl;
  return_ray.direction.w = 0;
  //another_r.direction = normalize(another_r.direction);
  //cout << "original direction is " << to_string(ray.direction) << endl;
  //cout << "after is " << to_string(return_ray.direction) << endl;
}*/

void perturbed_ray(Ray& ray,Ray& return_ray,int glossy){
  //double x1 = ((double) rand() / (RAND_MAX));
  //double x2 = ((double) rand() / (RAND_MAX));
  //double x1 = distr_one(eng);
  //double x2 = distr_one(eng);
  float x1 = distribution(generator);
	float x2 = distribution(generator);

  vec4 reflect_dir = normalize(ray.direction);

  //float alpha = acos(pow(1-x1,1/(0.01+1)));
  vec3 U = cross(vec3(reflect_dir),vec3(0,1,0));
  vec3 V = cross(vec3(reflect_dir),U);
  float beta = 2*M_PI*x2;
  double cos_theta = pow(x1,(double)1.0/(glossy+1));
  double sin_theta = sqrt(1-cos_theta*cos_theta);
  double cos_phi = cos(beta);
  double sin_phi = sin(beta);
  vec4 new_U = vec4(U,0);
  vec4 new_V = vec4(V,0);

  //vec3 normal = vec3(0,1,ray.direction.y/ray.direction.z);
  //mat4 first_R = glm::rotate(alpha,normal);
  //mat4 second_R = glm::rotate(beta,normal);

  //vec4 temp_dir = second_R*first_R*vec4(ray.direction);
  //cout << "temp dir is " << to_string(temp_dir) << endl;
  vec4 temp_dir = reflect_dir*cos_theta+new_U*cos_phi*sin_theta+new_V*sin_phi*sin_theta;


  return_ray.origin = ray.origin;
//  return_ray.direction = normalize(temp_dir);
  return_ray.direction = temp_dir;
}




vec3 Ray_color(Ray& my_ray,SceneNode *root,const list<Light*> &lights,vec4 look_from,vec3 ambient,int maxhit,vec3 background,vec3& my_colors,bool& show_background,bool debug){
    Intersection cur_sec(my_ray);
    Intersection* best_sec = NULL;

    cur_sec = check_intersect(my_ray, root);
    cur_sec.normal = normalize(cur_sec.normal);
    //cur_sec.center = normalize(cur_sec.center);

    //cout << "in ray color" << endl;


    vec3 col = my_colors;

    if(cur_sec.hit&&cur_sec.t>0.001){
        vec3 kd  = cur_sec.my_phong->get_kd();
        vec4 cur_sec_n = cur_sec.normal;
        vec4 P = my_ray.origin+cur_sec.t*my_ray.direction;
        double u,v = 0;

      //cout << "cur shape is " << cur_sec.shape << endl;
           if(!cur_sec.has_image){
             col += cur_sec.my_phong->get_kd()*ambient;
           } else {
             if(cur_sec.shape==0||cur_sec.shape==3){
               vec3 n = -normalize(vec3(P)-cur_sec.center);
             //double u = 0.5 - atan2(n.z, n.x)/2/M_PI;
				     //double v = 0.5 + asin(n.y)/M_PI;//

             u = atan2(n.x,n.z) / (2*M_PI)+0.5;
             v = n.y*0.5+0.5;

           } else if(cur_sec.shape==2){
                  vec3 n = normalize(vec3(P)-cur_sec.center);
                 u = atan2(n.x,n.y) / (2*M_PI)+0.5;
                 v = n.z/2+0.5;
                 //u = (M_PI+atan2(n.x,n.y))/4*M_PI+0.5;
                 //v = (1+n.z)/2+0.5;

              } else if(cur_sec.shape==1){
                vec4 cur_pos = cur_sec.local_pos;
                vec3 cur_nor = vec3(cur_sec.local_normal);

                if(distance(cur_nor,vec3(0,0,1))<0.01){
                  u = 1-cur_pos.x;
                  v = 1-cur_pos.y;
                } else if(distance(cur_nor,vec3(0,0,-1))<0.01){
                  u = cur_pos.x;
                  v = cur_pos.y;
                } else if(distance(cur_nor,vec3(0,1,0))<0.01){
                  u = cur_pos.x;
                  v = cur_pos.z;
                } else if(distance(cur_nor,vec3(0,-1,0))<0.01){
                  u = 1- cur_pos.x;
                  v = 1- cur_pos.z;
                } else if(distance(cur_nor,vec3(1,0,0))<0.01){
                  u = cur_pos.y;
                  v = cur_pos.z;
                } else {
                  u = 1-cur_pos.y;
                  v = 1-cur_pos.z;
                }
              }

              vec3 get_col(0,0,0);
              if(cur_sec.texture_type==0){
                Texture *cur_texture = dynamic_cast<Texture*>(cur_sec.m_primitive->material);
                //cout << "width is " << cur_texture->image.width() << endl;
                //cout << "hieght is " << cur_texture->image.height() << endl;
                get_col = cur_texture->uv_cal(u,v);
                kd = get_col;
              } else if(cur_sec.extra_image){
                Texture *cur_texture = dynamic_cast<Texture*>(cur_sec.m_primitive->extra_material);
                get_col = cur_texture->uv_cal(u,v);
                kd = get_col;
              }
                 //cout << "color is " << to_string(get_col) << endl;
                 col += kd*ambient;

              }


        for(Light *light: lights) {
#ifdef SOFTSHADOW
            vec3 temp_col(0,0,0);
            for(int i = 0;i<RAYNUM;++i){
              vec3 one_col(0,0,0);
                Ray new_ray;
                new_ray.origin = P;
                //vec3 new_light = ballRand(5);
                vec3 new_light_ray = light->position;
                new_light_ray.x += distr(eng);
                new_light_ray.y += distr(eng);
                new_light_ray.z += distr(eng);

                new_ray.direction = vec4(new_light_ray,1)-P;
                Intersection new_intersect = check_intersect(new_ray,root);

                if(!new_intersect.hit){
                    //cout << "t is " << cur_sec.t << endl;
                    vec3 l = vec3(normalize(light->position - vec3(P)));
                    vec3 cur_normal = vec3(cur_sec.normal);
                    float n_dot_l = dot(cur_normal, l);
                    float distance_light = glm::distance(vec4(light->position, 1), P);

                    vec3 diffuse;
                    vec3 kd  = cur_sec.my_phong->get_kd();

                    //texture mapping

                    diffuse = kd * n_dot_l;
                    if (n_dot_l > 0.0) {
                        one_col.x += kd.x * n_dot_l * light->colour.r;
                        one_col.y += kd.y * n_dot_l * light->colour.g;
                        one_col.z += kd.z * n_dot_l * light->colour.b;
                    }
                    float dot_use = dot(l,cur_normal);
                    vec3 R = 2*dot_use*cur_normal-l;
                    float temp_spec = dot(vec3(-normalize(my_ray.direction)),R);
                    temp_spec = pow(temp_spec,cur_sec.my_phong->get_shine());
                    vec3 ks  = cur_sec.my_phong->get_ks();
                    one_col.x += ks.x * temp_spec * light->colour.r;
                    one_col.y += ks.y * temp_spec * light->colour.g;
                    one_col.z += ks.z * temp_spec * light->colour.b;

                    //cout << "one color is " << glm::to_string(one_col) << endl;
                    temp_col+=one_col;
                }
              }
              temp_col/=RAYNUM;
              //cout << "temp color is " << glm::to_string(temp_col) << endl;
              col += temp_col;
#else

            Ray new_ray;
            new_ray.origin = P;
            new_ray.direction = vec4(light->position, 1)-P;
            Intersection new_intersect = check_intersect(new_ray,root);

            if(!new_intersect.hit){
                vec3 l = vec3(normalize(light->position - vec3(P)));
                vec3 cur_normal = vec3(cur_sec.normal);
                float n_dot_l = dot(cur_normal, l);
                float distance_light = glm::distance(vec4(light->position, 1), P);

                vec3 diffuse;
                if(cur_sec.has_image&&cur_sec.texture_type==1){
                  //kd  = cur_sec.my_phong->get_kd();
                  vec3 axis = vec3(0,0,1);
                  vec3 normal_3 = vec3(cur_sec.normal);
                  vec3 Ou = normalize(cross(axis,vec3(normal_3)));
                  vec3 Ov = -normalize(cross(vec3(normal_3),Ou));

                  //cout << "u is " << u << " v is " << v << endl;

                  //pair<float,float> derivatives = get_derivatives(*cur_sec.my_image,u,v);
                  //cout << "x is " << derivatives.first << " y is " << derivatives.second << endl;
                  //vec3 N_new = normalize(normal_3+ derivatives.first*cross(normal_3,Ov) - derivatives.second*cross(normal_3,Ou));
                  Bump *cur_texture = dynamic_cast<Bump*>(cur_sec.m_primitive->material);
                  //Image current_image = *cur_sec.my_image;
                  //P = normalize(P);
                  vec3 Bump_bal = cur_texture->uv_cal(u,v);
                  //cout << "Bump_ball is " << to_string(Bump_bal) << endl;
                  //cout << "original N is " << to_string(original_N) << endl;
                  vec3 N_new = normalize(normal_3+ Bump_bal);
                  vec3 R_new = normalize(2*dot(l,N_new)*N_new-l);
                  float n_dot_l_new = dot(N_new, l);
                  if (n_dot_l_new> 0.0) {
                      col.x += kd.x * n_dot_l_new * light->colour.r;
                      col.y += kd.y * n_dot_l_new * light->colour.g;
                      col.z += kd.z * n_dot_l_new * light->colour.b;
                  }
                  float temp_spec = dot(vec3(-normalize(my_ray.direction)),R_new);
                  temp_spec = pow(temp_spec,cur_sec.my_phong->get_shine());
                  vec3 ks  = cur_sec.my_phong->get_ks();
                  col.x += ks.x * temp_spec * light->colour.r;
                  col.y += ks.y * temp_spec * light->colour.g;
                  col.z += ks.z * temp_spec * light->colour.b;
                } else {
                  diffuse = kd * n_dot_l;
                  if (n_dot_l > 0.0) {
                      col.x += kd.x * n_dot_l * light->colour.r;
                      col.y += kd.y * n_dot_l * light->colour.g;
                      col.z += kd.z * n_dot_l * light->colour.b;
                  }
                  float dot_use = dot(l,cur_normal);
                  vec3 R = 2*dot_use*cur_normal-l;
                  float temp_spec = dot(vec3(-normalize(my_ray.direction)),R);
                  temp_spec = pow(temp_spec,cur_sec.my_phong->get_shine());
                  vec3 ks  = cur_sec.my_phong->get_ks();
                  col.x += ks.x * temp_spec * light->colour.r;
                  col.y += ks.y * temp_spec * light->colour.g;
                  col.z += ks.z * temp_spec * light->colour.b;
                }



                //cout << "one color is " << glm::to_string(one_col) << endl;

            }

#endif
}

#ifdef REFLECTION
        if(maxhit<5){
            Ray reflected_ray;
            reflected_ray.origin = P + 0.01*cur_sec.normal;
            reflected_ray.direction = my_ray.direction - 2*cur_sec.normal*dot(my_ray.direction,cur_sec.normal);
            vec3 ks  = cur_sec.my_phong->get_ks();
            vec3 ready_col(0,0,0);
            col += ks * Ray_color(reflected_ray,root,lights,look_from,ambient,maxhit+1,background,ready_col,show_background,false);
        }
#endif

#ifdef GLOSSY_REFLECTION
if(maxhit<5&&cur_sec.cur_glossy>0){
  //cout << "in glossy " << endl;
          Ray reflected_ray;
          reflected_ray.origin = P + 0.01*cur_sec.normal;
          reflected_ray.direction = my_ray.direction - 2*cur_sec.normal*dot(my_ray.direction,cur_sec.normal);
          vec3 reflect_col(0,0,0);
          vec3 ks  = cur_sec.my_phong->get_ks();
          for(int i = 0;i<GLOSSYNUM;++i){

            Ray perturbed_r;
            //double x1 = ((double) rand() / (RAND_MAX));
          //  double x2 = ((double) rand() / (RAND_MAX));
            perturbed_ray(reflected_ray,perturbed_r,cur_sec.cur_glossy);
          //float x1 = distribution(generator);
        	//float x2 = distribution(generator);

        /*  float alpha = acos(pow(1-x1,0.001));
          float beta = 2*M_PI*x2;

          vec3 nor = vec3(0,1,reflected_ray.direction.y/reflected_ray.direction.z);
          mat4 first_R = glm::rotate(alpha,nor);
          mat4 second_R = glm::rotate(beta,nor);

          vec4 temp_dir = second_R*first_R*vec4(reflected_ray.direction);
          //cout << "temp dir is " << to_string(temp_dir) << endl;

          perturbed_r.origin = reflected_ray.origin;
          perturbed_r.direction = normalize(temp_dir);*/
          //  glossy_rays.push_back(perturbed_r );
            vec3 ready_col(0,0,0);
            reflect_col += Ray_color(perturbed_r,root,lights,look_from,ambient,maxhit+1,background,ready_col,show_background,false);
          }
          reflect_col /= GLOSSYNUM;
          col = 0.5*col + 0.5*reflect_col;//ks*Ray_color(perturbed_r,root,lights,look_from,ambient,maxhit+1,background,ready_col,show_background,false);
        }

#endif

#ifdef REFRACTION
          if(cur_sec.my_phong->m_dielectric){
            vec3 ks = cur_sec.my_phong->get_ks();
            vec4 reflection_vec = glm::reflect( my_ray.direction,cur_sec.normal);
            vec4 bold_t;
            float eta = 1/cur_sec.my_phong->m_ior;

            float cos_theta = dot( my_ray.direction,cur_sec.normal);

            //float kr;
            //freshnel(my_ray.direction,cur_sec.normal,cur_sec.my_phong->m_ior,fr);
            bool outside;
            if(dot( my_ray.direction,cur_sec.normal)<0){
              outside = true;
              bold_t = glm::refract( my_ray.direction,cur_sec.normal,eta);
              cos_theta = -cos_theta;
              ks = vec3(0.65,0.65,0.65);
              //cout << "this case" << endl;
            } else {
              outside = false;
              bold_t = glm::refract(my_ray.direction,-cur_sec.normal,eta);
              vec4 standard = vec4(0,0,0,0);
              if(bold_t !=standard){
                cos_theta = dot(bold_t,cur_sec.normal);
                //cout << "not standard" << endl;
              }
            }
            double ebs = 1;
            vec4  bias = 0.1*cur_sec.normal;
            float R0 = pow(eta-1,2)/pow(eta+1,2);
            float R = R0 + (1-R0)*pow(1-cos_theta,5);
            Ray reflected_r;
            reflected_r.origin = outside? P + bias : P-bias;
            reflected_r.direction = reflection_vec*ebs;
            Ray refracted_r;
            refracted_r.origin = outside? P - bias : P+bias;
            refracted_r.direction = bold_t*ebs;
            //vec3 k_temp(0.5,0.5,0.5);
            vec3 ready_col(0,0,0);
            vec3 ready_an(0,0,0);
            //vec3 col_another = k_temp*R*Ray_color(reflected_r,root,lights,look_from,ambient,maxhit+1,background,my_colors,show_background,false);
            //cout << "col another is " << to_string(col_another) << endl;
            col += ks*(R*Ray_color(reflected_r,root,lights,look_from,ambient,maxhit+1,background,ready_an,show_background,false)
                       + (1-R)*Ray_color(refracted_r,root,lights,look_from,ambient,maxhit+1,background,ready_col,show_background,false));

          }

#endif
        //cout << "col r is " << col.r << " col g is " << col.g << " col b is " << col.b << endl;
        my_colors = col;
         return col;
    } else {
    //cout << "col r is " << col.r << " col g is " << col.g << " col b is " << col.b << endl;
        my_colors = background;
        show_background = true;
        return background;
    }

}


void wrap_ray_color(SceneNode *root,vec3 eye, mat4 p_temp_world, const list<Light*> &lights,vec3 ambient,int maxhit,vector<vector<vec3> >& my_colors,int cur_row,int ny,int total_x,vec3 view){
    bool debug = false;
    for(int y = 0;y<ny;++y){
        //if(y==160){
        vec3 background(0.6,0.75,0.4);

        background.z = (float)y/(ny);
        bool show_background = false;
#ifdef ANTIALIASING
        if(cur_row+1<total_x&&y+1<ny){
            Ray my_ray;
            vec3 col = vec3(0,0,0);
            vec4 look_from = vec4(eye,1);
            //vec4 p_world = p_temp_world*vec4(cur_row,y,0,1);
            //my_ray.direction = normalize(p_world-look_from);
            my_ray.origin = vec4(eye,1);
            for(int i = 0;i<8;++i){
                for(int j = 0;j<8;++j){
                    double random = ((double) rand() / (RAND_MAX));
                    //cout << "random is " << random << endl;
                    vec4 p_world = p_temp_world*vec4(cur_row+(i+random)/8,y+(j+random)/8,0,1);
                    my_ray.direction = normalize(p_world-look_from);
                    Ray_color(std::ref(my_ray),root,std::ref(lights),look_from,ambient,0,background,std::ref(my_colors[cur_row][y]),ref(show_background));


                    if(show_background)break;
                }
                if(show_background)break;
            }
            if(show_background)continue;
            my_colors[cur_row][y]  = my_colors[cur_row][y]/64;

        }
       // cout << "my color is " << my_colors[cur_row][y].r << ' ' << my_colors[cur_row][y].g << ' ' << my_colors[cur_row][y].b << endl;
#else
        Ray my_ray;
        vec3 col = vec3(0,0,0);
        vec4 look_from = vec4(eye,1);
        my_ray.origin = vec4(eye,1);
        //my_ray.direction = normalize(p_world- my_ray.origin);
        vec4 p_world = p_temp_world*vec4(cur_row,y,0,1);
        my_ray.direction = normalize(p_world-my_ray.origin);

//for depth of field
        //my_ray.direction += distr(eng);
        #ifdef DEPTHOFFIELD

                double length =770;
                vec4 focal = my_ray.origin + my_ray.direction*length;

                vec3 depth_col(0,0,0);
                for(int i = 0;i<RAYNUM;++i){
                  my_ray.origin.x += distr_soft(eng);
                  my_ray.origin.y += distr_soft(eng);
                  my_ray.direction = normalize(focal- my_ray.origin);
                  vec3 ready_col(0,0,0);
                  depth_col+=Ray_color(std::ref(my_ray),root,std::ref(lights),look_from,ambient,0,background,std::ref(ready_col),show_background,debug);
                  if(show_background)break;
                }
                //if(show_background)continue;
                 if(!show_background){
                   depth_col/=RAYNUM;
                 }
                 my_colors[cur_row][y] = depth_col;


      //  cout << "after change is " << glm::to_string(my_ray.origin) << endl;
#else
        my_ray.origin = vec4(eye,1);
        my_ray.direction = normalize(p_world- my_ray.origin);
        if(cur_row==100&&y==100)debug = true;
        Ray_color(std::ref(my_ray),root,std::ref(lights),look_from,ambient,0,background,std::ref(my_colors[cur_row][y]),show_background,debug);
        debug = false;
#endif

        //cout << "my color is " << my_colors[cur_row][y].r << ' ' << my_colors[cur_row][y].g << ' ' << my_colors[cur_row][y].b << endl;
#endif
/*#else
        Ray my_ray;
        vec3 col = vec3(0,0,0);
        vec4 look_from = vec4(eye,1);
        my_ray.origin = vec4(eye,1);

        vec4 p_world = p_temp_world*vec4(cur_row,y,0,1);
#ifdef DEPTHOFFIELD
        double focal = 300;
        //p_world =
        my_ray.origin.x += distr_one(eng)*focal;
        my_ray.origin.y += distr_one(eng)*focal;
        my_ray.direction = normalize(p_world- my_ray.origin);

      //  cout << "after change is " << glm::to_string(my_ray.origin) << endl;
#else
        my_ray.origin = vec4(eye,1);
        my_ray.direction = normalize(p_world- my_ray.origin);
#endif
        if(cur_row==100&&y==100)debug = true;
        Ray_color(std::ref(my_ray),root,std::ref(lights),look_from,ambient,0,background,std::ref(my_colors[cur_row][y]),show_background,debug);
        debug = false;

        //cout << "my color is " << my_colors[cur_row][y].r << ' ' << my_colors[cur_row][y].g << ' ' << my_colors[cur_row][y].b << endl;
#endif*/



pixels_done++;
std::cout << "\r" << std::fixed << std::setprecision(2) << ((float(pixels_done) / total_pixels) * 100) << "% done";

}
        //pixels_done++;
        //std::cout << "\r" << std::fixed << std::setprecision(2) << ((float(pixels_done) / total_pixels) * 100) << "% done";
    //}
}



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
) {

  // Fill in raytracing code here

  std::cout << "Calling A4_Render(\n" <<
		  "\t" << *root <<
          "\t" << "Image(width:" << image.width() << ", height:" << image.height() << ")\n"
          "\t" << "eye:  " << glm::to_string(eye) << std::endl <<
		  "\t" << "view: " << glm::to_string(view) << std::endl <<
		  "\t" << "up:   " << glm::to_string(up) << std::endl <<
		  "\t" << "fovy: " << fovy << std::endl <<
          "\t" << "ambient: " << glm::to_string(ambient) << std::endl <<
		  "\t" << "lights{" << std::endl;

	for(const Light * light : lights) {
		std::cout << "\t\t" <<  *light << std::endl;
	}
	std::cout << "\t}" << std::endl;
	std:: cout <<")" << std::endl;

cout << "barely start" << endl;


    int nx = image.width();
    int ny = image.height();
    float d = distance(eye,view);

    total_pixels = nx*ny;


    vec3 translate_vec = vec3(-nx/2,-ny/2,d);
    mat4 T1 = translate(mat4(),translate_vec);

    float h = 2*d*tan(radians(fovy)/2);
    float w = nx*h/ny;
    vec3 scale_vec = vec3(-w/nx,-h/ny,1);
    mat4 S2 = scale(mat4(),scale_vec);

    vec3 w_view = normalize(view-eye);
    vec3 u = normalize(cross(up,w_view));
    vec3 v = cross(w_view,u);
    mat4 R3 = mat4(vec4(u,0.0f),vec4(v,0.0f),vec4(w_view,0.0f),vec4(0.0f,0.0f,0.0f,1.0f));


    mat4 T4 = mat4(vec4(1,0,0,0),
                   vec4(0,1,0,0),
                   vec4(0,0,1,0),
                   vec4(eye,1));



    //concurrency
    thread my_threads[nx];
    vector<vector<vec3> > my_colors;
    for(int i = 0;i<nx;++i){
        vec3 one_color(0,0,0);
        vector<vec3> row(ny,one_color);
        my_colors.push_back(row);
    }

    mat4 p_temp_world = T4*R3*S2*T1;

    //print_mat4(p_temp_world);
    int pixels_done = 0;
    int total_pixels = nx*ny;

    find_brother(root);


    for (int x = 0; x < nx; ++x) {
      //if(x==155){
        my_threads[x] = thread(wrap_ray_color,root,eye,p_temp_world,lights,ambient,0,ref(my_colors),x,ny,nx,view);



      //}
    }


    cout << "complete" << endl;



    for(int x = 0;x<nx;++x){
        if(my_threads[x].joinable()){
            my_threads[x].join();
        }
        //pixels_done+=ny;
        //std::cout << "\r" << std::fixed << std::setprecision(2) << ((float(pixels_done) / total_pixels) * 100) << "% done";
    }

    for(int y = 0;y<ny;++y){
        for(int x = 0;x<nx;++x){
            //cout << "col r is " <<my_colors[x][y].r << " col g is " << my_colors[x][y].g << " col b is " << my_colors[x][y].b << endl;
            image(x, y, 0) = my_colors[x][y].r;
            // Green:
            image(x, y, 1) = my_colors[x][y].g;
            // Blue:
            image(x, y, 2) = my_colors[x][y].b;
        }
    }
    //cout << "sum is " << sum << endl;
    /*std::time_t start_time_t = std::chrono::system_clock::to_time_t(start_time);
    std::time_t end_time_t = std::chrono::system_clock::to_time_t(end_time);
    std::cout << "Start time: " << std::ctime(&start_time_t) << std::endl;
    std::cout << "End time: " << std::ctime(&end_time_t) << std::endl;*/
			// Red:

           // }
}
