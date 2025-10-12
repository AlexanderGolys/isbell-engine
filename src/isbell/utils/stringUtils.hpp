#pragma once
#include "exceptions.hpp"
#include "concepts.hpp"
#include "abstractNonsense.hpp"

inline string removeWhitespace(const string &s) {
	string res;
	for (char c: s)
		if (!isspace(c)) res += c;
	return res;
}

inline bool isIdentStart(char c) {
	return (c == '_') || (c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z');
}


inline bool isIdentChar(char c) {
	return isIdentStart(c) || (c >= '0' && c <= '9');
}


inline string removeComments(const string &s) {
	string out;
	out.reserve(s.size());

	bool inBlock = false;
	bool inLine = false;
	bool inString = false;

	for (size_t i = 0; i < s.size(); ) {
		char c = s[i];

		if (inLine) {
			if (c == '\n') {
				inLine = false;
				out += c;
				i++;
				continue;
			}
			i++;
			continue;
		}

		if (inBlock) {
			if (c == '*' && i + 1 < s.size() && s[i + 1] == '/') {
				inBlock = false;
				i += 2;
				continue;
			}
			i++;
			continue;
		}

		if (inString) {
			if (c == '\\' && i + 1 < s.size()) {
				out += s[i];
				out += s[i + 1];
				i += 2;
				continue;
			}
			out += c;
			if (c == '"')
				inString = false;
			i++;
			continue;
		}

		if (c == '/' && i + 1 < s.size() && s[i + 1] == '/') {
			inLine = true;
			i += 2;
			continue;
		}
		if (c == '/' && i + 1 < s.size() && s[i + 1] == '*') {
			inBlock = true;
			i += 2;
			continue;
		}
		if (c == '"') {
			inString = true;
			out += c;
			i++;
			continue;
		}

		out += c;
		i++;
	}

	return out;
}


inline string stripCurlyBlocks(const string &s) {
	string out;
	out.reserve(s.size());
	int depth = 0;

	for (size_t i = 0; i < s.size(); i++) {
		char c = s[i];
		if (c == '{') {
			depth++;
			out += ' ';
			i++;
			while (i < s.size() && depth > 0) {
				if (s[i] == '{')
					depth++;
				else if (s[i] == '}')
					depth--;
				i++;
			}
			if (i > 0)
				i--;
			continue;
		}
		if (depth == 0)
			out += c;
	}

	return out;
}


inline string trim(const string &line) {
	size_t a = 0;
	size_t b = line.size();
	while (a < b && (line[a] == ' ' || line[a] == '\t' || line[a] == '\r'))
		a++;
	while (b > a && (line[b - 1] == ' ' || line[b - 1] == '\t' || line[b - 1] == '\r' || line[b - 1] == '\n'))
		b--;
	return line.substr(a, b - a);
}



inline string vecToString(vec2 v) {
	return "(" + to_string(v.x) + ", " + to_string(v.y) + ")";
}

inline string vecToString(vec3 v) {
	return "(" + to_string(v.x) + ", " + to_string(v.y) + ", " + to_string(v.z) + ")";
}

inline string vecToString(vec4 v) {
	return "(" + to_string(v.x) + ", " + to_string(v.y) + ", " + to_string(v.z) + ", " + to_string(v.w) + ")";
}

inline string vecToString(ivec2 v) {
	return "(" + to_string(v.x) + ", " + to_string(v.y) + ")";
}

inline string vecToString(ivec3 v) {
	return "(" + to_string(v.x) + ", " + to_string(v.y) + ", " + to_string(v.z) + ")";
}

inline string vecToString(ivec4 v) {
	return "(" + to_string(v.x) + ", " + to_string(v.y) + ", " + to_string(v.z) + ", " + to_string(v.w) + ")";
}

template<typename T>
void printVector(vector<T> v, string title = "vector") {
	std::cout << title << " [" << v.size() << "]: ";
	for (int i = 0; i < v.size(); i++) {
		std::cout << v[i] << ", ";
	}
	std::cout << std::endl;
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

inline vector<string> split(const string &s, char delim) {
	return split(s, string(1, delim));
}

inline string replace(const string &s, const string &from, const string &to) {
	string res = s;
	size_t start_pos = 0;

	while ((start_pos = res.find(from, start_pos)) != string::npos) {
		res.replace(start_pos, from.length(), to);
		start_pos += to.length();
	}
	return res;
}

inline string replaceRecursive(const string &s, const string &from, const string &to) {
	string res = replace(s, from, to);
	if (res == s) return res;
	return replaceRecursive(res, from, to);
}

inline string remove(const string &s, const string &key) {
	return replace(s, key, "");
}

inline string removeRecursive(const string &s, const string &key) {
	return replaceRecursive(s, key, "");
}

inline string join(const string &delim, const vector<string> &parts) {
	string res;
	for (int i = 0; i < parts.size(); i++) {
		res += parts[i];
		if (i < parts.size() - 1) res += delim;
	}
	return res;
}
