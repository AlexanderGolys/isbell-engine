#pragma once
#include "macros.hpp"


inline string plural (string word) {
	if (word.size() == 0)
		return word;
	bool caps = word.back() >= 'A' && word.back() <= 'Z' && word.size() > 1;
	bool es = word.back() == 's'||
			  word.back() == 'x' ||
			  word.back() == 'z' ||
			  (word.size() > 1 && word[word.size()-2] == 'c' && word.back() == 'h') ||
			  (word.size() > 1 && word[word.size()-2] == 's' && word.back() == 'h') ||
			  word.back() == 'S'||
			  word.back() == 'X' ||
			  word.back() == 'Z' ||
			  (word.size() > 1 && word[word.size()-2] == 'C' && word.back() == 'H') ||
			  (word.size() > 1 && word[word.size()-2] == 'S' && word.back() == 'H');
	if (es && !caps)
		return word + "es";
	if (es)
		return word + "ES";
	if (caps)
		return word + "S";
	return word + "s";
}


class ErrorClassWrapper : public std::exception {
	string msg_;
public:
	explicit ErrorClassWrapper(const string& msg) : msg_(msg) {}
	const char* what() const noexcept override {
		return msg_.c_str();
	}
	string message() const { return msg_; }
};

class NotImplementedError : public ErrorClassWrapper {
public:
    explicit NotImplementedError(const string& notImplementedMethodName, const string& lackingType="Method")
        : ErrorClassWrapper(lackingType + " " + notImplementedMethodName + " is not implemented yet.") {}
};

class RecursionLimitExceeded : public ErrorClassWrapper {
public:
	explicit RecursionLimitExceeded(int limit, const string &where)
		: ErrorClassWrapper("Recursion limit (" + std::to_string(limit) + "levels) exceeded during " + where) {}
};

class IndexOutOfBounds : public ErrorClassWrapper {
public:
	IndexOutOfBounds(int index, int size, const string &indexName="i")
		: ErrorClassWrapper("Index " + indexName + " is out of bounds [" + std::to_string(index) + "/" + std::to_string(size) + "].") {}
	IndexOutOfBounds(const string &index, const string &size, const string &indexName="i")
		: ErrorClassWrapper("Index " + indexName + " is out of bounds [" + index + "/" + size + "].") {}
	IndexOutOfBounds(const ivec2 &index, const ivec2 &size, const string &indexName="i")
		: ErrorClassWrapper("Index " + indexName + " is out of bounds [" + format("({}, {})", index.x, index.y) + "/" + format("({}, {})", size.x, size.y) + "].") {}
	IndexOutOfBounds(const ivec3 &index, const ivec3 &size, const string &indexName="i")
	: ErrorClassWrapper("Index " + indexName + " is out of bounds [" + format("({}, {}, {})", index.x, index.y, index.z) + "/" + format("({}, {}, {})", size.x, size.y, size.z) + "].") {}
	IndexOutOfBounds(const ivec4 &index, const ivec4 &size, const string &indexName="i")
	: ErrorClassWrapper("Index " + indexName + " is out of bounds [" + format("({}, {}, {}, {})", index.x, index.y, index.z, index.w) + "/" + format("({}, {}, {}, {})", size.x, size.y, size.z, size.w) + "].") {}
};

class NotImplementedMethodError : public NotImplementedError {
public:
  explicit NotImplementedMethodError(const string& methodName)
      :NotImplementedError(methodName, "Method") {}
};

class NotImplementedFunctionError : public NotImplementedError {
public:
  explicit NotImplementedFunctionError(const string& name)
      :NotImplementedError(name, "Function") {}
};

class NotImplementedVariantError : public NotImplementedError {
public:
  NotImplementedVariantError(const string& variant, const string& ofWhat)
      :NotImplementedError(variant + " of " + ofWhat, "Variant") {}
};

class UnknownVariantError : public ErrorClassWrapper {
public:
  UnknownVariantError(const string& variant, const string& ofWhat)
      : ErrorClassWrapper("Variant " + variant + " of  " + ofWhat + " is an unknown type in this context.") {}
  explicit UnknownVariantError(const string& msg) : ErrorClassWrapper(msg) {}

};

class IllegalVariantError : public ErrorClassWrapper {
public:
  IllegalVariantError(const string& variant, const string& ofWhat, const string& rejectingMethod)
      : ErrorClassWrapper(plural(ofWhat) + " in variant " + variant + " are considered invalid by method " + rejectingMethod + ".") {}
  explicit IllegalVariantError(const string& msg) : ErrorClassWrapper(msg) {}

};

class IllegalArgumentError : public ErrorClassWrapper {
public:
	explicit IllegalArgumentError(const string& msg) : ErrorClassWrapper(msg) {}
};

class ValueError : public ErrorClassWrapper {
public:
	explicit ValueError(const string& msg) : ErrorClassWrapper(msg) {}
};

class SystemError : public ErrorClassWrapper {
public:
	explicit SystemError(const string& msg) : ErrorClassWrapper(msg) {}
};

class FileNotFoundError : public SystemError {
public:
	explicit FileNotFoundError(const string& filename) : SystemError("File " + filename + " not found.") {}
	FileNotFoundError(const string& filename, const string& dir) : SystemError("File " + filename + " not found in " + dir) {}
};

class InvalidFileError : public SystemError {
public:
	explicit InvalidFileError(const string& filename, const string& reason="")
	: SystemError(format("File {0} is invalid{1}{2}{3}.", filename, reason == "" ? " (" : "", reason, reason == "" ? ")" : "")) {}
};

class ZeroDivisionError : public ValueError {
public:
	explicit ZeroDivisionError(const string& msg="") : ValueError(msg) {}
};
