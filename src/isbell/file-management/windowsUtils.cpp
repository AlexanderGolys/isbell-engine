#include <windows.h>
#include "filesUtils.hpp"




void FileDescriptor::mapFile() {
    if (isMapped())
        return;
    string fullPath = path.string();
    if (!exists())
        THROW(FileNotFoundError, fullPath);
    fileHandle = CreateFileA(fullPath.c_str(), GENERIC_READ | GENERIC_WRITE, FILE_SHARE_READ, 0, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, 0);
    if (fileHandle == INVALID_HANDLE_VALUE)
        throw InvalidFileError(fullPath, "Invalid file handle value", __FILE__, __LINE__);
    mappingHandle = CreateFileMappingA(fileHandle, nullptr, PAGE_READWRITE, 0, 0, nullptr);
    if (!mappingHandle)
        throw InvalidFileError(fullPath, "File mapping creation failed", __FILE__, __LINE__);
    address = MapViewOfFile(mappingHandle, FILE_MAP_WRITE | FILE_MAP_READ, 0, 0, 0);
    if (!address)
        throw InvalidFileError(fullPath, "Mapping view creation failed", __FILE__, __LINE__);
}

void FileDescriptor::flush() const {
    if (address == nullptr)
        THROW(FileSystemError, "File not mapped, cannot be flushed");
    FlushViewOfFile(address, bytesize);
}



void FileDescriptor::moveTo(const Path &destination) {
    bool wasMapped = isMapped();
    if (wasMapped)
        closeFile();
    DirectoryEntry::moveTo(destination);
    if (wasMapped)
        mapFile();
}

void FileDescriptor::remove() {
    closeFile();
    filesystem::remove(path);
    bytesize = 0;
    type = FILE_NOT_FOUND;
}

size_t FileDescriptor::getSize() const {
    return bytesize;
}



void FileDescriptor::closeFile() {
    if (address == nullptr)
        return;
    if (address)
        UnmapViewOfFile(address);
    address = nullptr;
    if (mappingHandle)
        CloseHandle(mappingHandle);
    mappingHandle = nullptr;
    if (fileHandle)
        CloseHandle(fileHandle);
    fileHandle = nullptr;
}
