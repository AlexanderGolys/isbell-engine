#pragma once
#include "exceptions.hpp"
#include "concepts.hpp"
#include "abstractNonsense.hpp"

inline string removeWhitespace(const string &s){
	string res;
	for (char c : s)
		if (!isspace(c))
			res += c;
	return res;
}




inline string vecToString(vec2 v) {
	return "(" + to_string(v.x) + ", " + to_string(v.y) + ")"; }

inline string vecToString(vec3 v) {
	return "(" + to_string(v.x) + ", " + to_string(v.y) + ", " + to_string(v.z) + ")"; }

inline string vecToString(vec4 v) {
	return  "(" + to_string(v.x) + ", " + to_string(v.y) + ", " + to_string(v.z) + ", " + to_string(v.w) + ")"; }

inline string vecToString(ivec2 v) {
	return "(" + to_string(v.x) + ", " + to_string(v.y) + ")"; }

inline string vecToString(ivec3 v) {
	return "(" + to_string(v.x) + ", " + to_string(v.y) + ", " + to_string(v.z) + ")"; }

inline string vecToString(ivec4 v) {
	return  "(" + to_string(v.x) + ", " + to_string(v.y) + ", " + to_string(v.z) + ", " + to_string(v.w) + ")"; }

template<typename T>
void printVector(vector<T> v, string title="vector")
{
	std::cout << title << " [" << v.size() << "]: ";
	for (int i = 0; i < v.size(); i++)
	{
		std::cout << v[i] << ", ";
	}
	std::cout << std::endl;
}


inline vector<string> split(const string &s, char delim) {
	vector<string> res;
	std::stringstream ss(s);
	string item;
	while (std::getline(ss, item, delim)) {
		res.push_back(item);
	}
	return res;
}

inline vector<string> split(const string &s, const string &delim) {
	vector<string> res;
	size_t start = 0;
	size_t end = s.find(delim);
	while (end != string::npos) {
		res.push_back(s.substr(start, end - start));
		start = end + delim.length();
		end = s.find(delim, start);
	}
	res.push_back(s.substr(start, end));
	return res;
}
