// Winter 2019

#include "Material.hpp"
using namespace glm;

Material::Material()
{}

Material::~Material()
{}


  vec3 Texture::uv_cal(double u,double v){
    vec3 col(0,0,0);
    //Image cur_image = image;
    //cout << "u is " << u << " v is " << v << endl;
    double di = (image.width() - 1)*u;
    //cout << "di is " << di << endl;
  	double dj = (image.height() - 1)*v;
    //cout << "dj is " << dj << endl;
    int i = int(di);
    double up = di-i;
    int j = int(dj);
    double vp = dj-j;
    //cout << "image i j 0 is " << image(i,j,0) << endl;
    col.r = image(i,j,0)*(1-up)*(1-vp) + image(i,j+1,0)*(1-up)*vp
            + image(i+1,j,0)*up*(1-vp) + image(i+1,j+1,0)*up*vp;
    col.g =  image(i,j,1)*(1-up)*(1-vp) + image(i,j+1,1)*(1-up)*vp
            + image(i+1,j,1)*up*(1-vp) + image(i+1,j+1,1)*up*vp;
    col.b =  image(i,j,2)*(1-up)*(1-vp) + image(i,j+1,2)*(1-up)*vp
            + image(i+1,j,2)*up*(1-vp) + image(i+1,j+1,2)*up*vp;
    return col;
  }


  vec3 Bump::uv_cal(double u,double v){
    vec3 col(0,0,0);
    //Image cur_image = image;
    //cout << "u is " << u << " v is " << v << endl;
    double di = (image.width() - 1)*u;
    //cout << "di is " << di << endl;
  	double dj = (image.height() - 1)*v;
    //cout << "dj is " << dj << endl;
    int i = int(di);
    double up = di-i;
    int j = int(dj);
    double vp = dj-j;
    //cout << "image i j 0 is " << image(i,j,0) << endl;
    col.r = image(i,j,0)*(1-up)*(1-vp) + image(i,j+1,0)*(1-up)*vp
            + image(i+1,j,0)*up*(1-vp) + image(i+1,j+1,0)*up*vp;
    col.g =  image(i,j,1)*(1-up)*(1-vp) + image(i,j+1,1)*(1-up)*vp
            + image(i+1,j,1)*up*(1-vp) + image(i+1,j+1,1)*up*vp;
    col.b =  image(i,j,2)*(1-up)*(1-vp) + image(i,j+1,2)*(1-up)*vp
            + image(i+1,j,2)*up*(1-vp) + image(i+1,j+1,2)*up*vp;
    return col;
  }
