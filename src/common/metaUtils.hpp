#ifndef METAUTILS_HPP
#define METAUTILS_HPP

#include <exception>
#include <string>

inline std::string plural (std::string word) {
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

class NotImplementedError : public std::exception {
    std::string msg_;
public:
    NotImplementedError(const std::string& notImplementedMethodName, const std::string& lackingType="Method")
        : msg_(lackingType + " " + notImplementedMethodName + " is not implemented yet.") {}

    const char* what() const noexcept override {
        return msg_.c_str();
    }
};

class NotImplementedMethodError : public NotImplementedError {
public:
  NotImplementedMethodError(const std::string& methodName)
      :NotImplementedError(methodName, "Method") {}
};

class NotImplementedFunctionError : public NotImplementedError {
public:
  NotImplementedFunctionError(const std::string& name)
      :NotImplementedError(name, "Function") {}
};

class NotImplementedVariantError : public NotImplementedError {
public:
  NotImplementedVariantError(const std::string& variant, const std::string& ofWhat)
      :NotImplementedError(variant + " of " + ofWhat, "Variant") {}
};

class UnknownVariantError : public std::exception {
  std::string msg_;
public:
  UnknownVariantError(const std::string& variant, const std::string& ofWhat)
      : msg_("Variant " + variant + " of  " + ofWhat + " is an unknown type in this context.") {}
  UnknownVariantError(const std::string& msg)
      : msg_(msg) {}


  const char* what() const noexcept override {
    return msg_.c_str();
  }
};

class IllegalVariantError : public std::exception {
  std::string msg_;
public:
  IllegalVariantError(const std::string& variant, const std::string& ofWhat, const std::string& rejectingMethod)
      : msg_(plural(ofWhat) + " in variant " + variant + " are considered invalid by method " + rejectingMethod + ".") {}
  IllegalVariantError(const std::string& msg)
      : msg_(msg) {}

  const char* what() const noexcept override {
    return msg_.c_str();
  }
};

#endif