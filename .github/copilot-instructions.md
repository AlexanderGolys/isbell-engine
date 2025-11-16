# Basic guidlines
- when I say "length of" I usually mean the number of elements in some sequential data structure (like vector size). By size I usually mean size in bytes


# Code Style
- use macro `HOM(X, Y)` for `std::function<Y(X)>` and `BIHOM(X, Y, Z)` for `std::function<Z(X, Y)>`
- Try to mimick the existing style in the file
- use `and`, `or` and `not` instead of `&&`, `||` and `!` operators
- Always split definitions and declarations unless the source file associated to the header doesnt exist or it is inline function outsdide the body of class/struct
- never use brackets for if/for/while statements when you are not forced by language syntax. Note that body may be nested and not single line technically speaking 
```
    for (int i = 0; i < 10; i++) // no brackets, but body has 2 lines
        if (condition)
            doSomething();
```
- always start header files with `#pragma once `
- for pointers used only as memory address to be dereferenced later use the type `raw_data_ptr` (alias for `const void*`)
- the pointers to memory buffers that preserve the data type information should be used with type `data_ptr<T> `(alias for `const T*`)
- mutable pointers to buffers are also have own alias `data_ptr_mut<T>` (alias for `T*`)
- Always use some alias when pointer is used for pointing to to data of some memory buffer, like when you bufferData in opengl from vector of floats (one from the 3 above)
- when a number is used to represent size in bytes, used type `byte_size` (alias for `size_t`)
- when a number is used to represent the length of some sequential data (usually some vector length) use `array_len` (alias for `size_t`)
- when a number is used to represent number of components of a vector (2 for `vec2`, 3 for `vec3` etc) use `vs_dim` (alias for `unsigned short`)
- Besides pointers to buffers and text, **NEVER** use raw pointers (I usually use `shared_pointer`, and somewtimes `unique_pointer` when ownership is clear)
- **NEVER** use `new` or `delete` keywords (unless in some copy-pasted code for example reading a file)
- for logging use macros from _logging.hpp_, also listed below with usage explained in comments:
```
    LOG("msg");               // default log you'll use most of the time
    LOG_ERROR("msg");        // error level log, used for adding extra info about exceptions and in testing. It should be used together with throw (besides unittests)
    LOG_PURE("i");         // pure log, use it only when you want to log a lot of lines at once, for example when listing files or representing a matrix
    
    // logs that I haven't use so far, but they mayt be suggested if some type of warning will be appropriate to communicate 
    LOG_WARN("warning");
    LOG_ERROR_PURE("error pure");
    LOG_WARN_PURE("warning pure");
    LOG_ERROR("formatted string {}", arg);
```
- never use namespace `std::` prefix unless in following cases:
- - `std::function` (that should be avoided anyway by using `HOM`/`BIHOM`)`
- - `std::map`
- - `std::set`
- - when using this structure for the first time and it is not enlisted in _macros.hpp_ as names from `std::` used in this namespace
- - it was used before in the same file with `std::` prefix (always follow the existing style in the file)

# Handling Exceptions
- **NEVER** use `throw` keyword, use only my macros
- **NEVER** use exceptions that are defined outside _exceptions.h_
- in particular, do not use exceptions from `std::`
- dont put `__FILE__` and `__LINE__` manually as arguments (macros do that for you)
- if the exception is thrown after checking some sondition, i.e. it is effectively assert, use macro `THROW_IF(condition, ExceptionClass, args...)`
``` 
    THROW_IF(not assertion_checker(), SystemError, "Something bad happened");
```
- use macro `THROW(ExceptionClass, args...)` for usual throws
```
    THROW(SystemError, "Something bad happened");
```
- if the error is caused by index out of bounds, use class `IndexOutOfBounds` with arguments `(index, size, container_name)` or `(index, size)` if container name is not known, obvious or it is a vector
- check if index is out of bounds automatically with macro `CHECK_OUT_OF_BOUNDS(index, size)`. It will throw appropriate exception if needed
- if the name of data structure where ther index is checked is important (useful for debugging), use macro `CHECK_OUT_OF_BOUNDS_NAME(index, size, container_name)`

