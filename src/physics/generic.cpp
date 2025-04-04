//
// Created by PC on 13.03.2025.
//

#include "generic.hpp"

vec3 StressTensor::traction_vector(vec3 direction) const {
	return (*this) * direction;
}

float StressTensor::stress_component(int i, int j) const {
	return (*this)[i][j];
}

float StressTensor::normal_stress(int i) const {
	return (*this)[i][i];
}

float StressTensor::shear_stress(int i, int j) const {
	if (i==j) throw std::format_error("diagonal elements represent normal stress");
	return stress_component(i, j); }
