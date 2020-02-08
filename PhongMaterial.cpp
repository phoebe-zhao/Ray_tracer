// Winter 2019

#include "PhongMaterial.hpp"

PhongMaterial::PhongMaterial(
	const glm::vec3& kd, const glm::vec3& ks, double shininess )
	: m_kd(kd)
	, m_ks(ks)
	, m_shininess(shininess)
{
	m_dielectric = false;
}

PhongMaterial::PhongMaterial(
	const glm::vec3& kd, const glm::vec3& ks, double shininess,bool transparent,float ior )
	: m_kd(kd)
	, m_ks(ks)
	, m_shininess(shininess)
{
	m_dielectric = transparent;
	m_ior = ior;
}

PhongMaterial::~PhongMaterial()
{}
