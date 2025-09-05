#define NOMINMAX
#include <windows.h>
#include "../filesUtils.hpp"



FileDescriptor::~FileDescriptor() {
    if (address)
        UnmapViewOfFile(address);
    if (mappingHandle)
        CloseHandle(mappingHandle);
    if (fileHandle)
        CloseHandle(fileHandle);
}


void FileDescriptor::mapFile() {
    string fullPath = path.to_str() + "\\" + filename;
    fileHandle = CreateFileA(fullPath.c_str(), GENERIC_READ, FILE_SHARE_READ, 0, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, 0);
    if (fileHandle == INVALID_HANDLE_VALUE)
        throw InvalidFileError(fullPath, "Invalid file handle value", __FILE__, __LINE__);
    LARGE_INTEGER size;
    if (!GetFileSizeEx(fileHandle, &size))
        throw InvalidFileError(fullPath, "File size retrieval failed", __FILE__, __LINE__);
    bytesize = size.QuadPart;
    mappingHandle = CreateFileMappingA(fileHandle, nullptr, PAGE_READONLY, 0, 0, nullptr);
    if (!mappingHandle)
        throw InvalidFileError(fullPath, "File mapping creation failed", __FILE__, __LINE__);
    address = MapViewOfFile(mappingHandle, FILE_MAP_READ, 0, 0, 0);
    if (!address)
        throw InvalidFileError(fullPath, "Mapping view creation failed", __FILE__, __LINE__);
}
