#pragma once
#include <expected>

#include "metaUtils.hpp"

typedef std::filesystem::path Path;


enum RecognisedFileType {
	REG_FILE,
	DIR,
	SYM_LINK,
	UNRECOGNISED,
	EMPTY_PATH,
	FILE_NOT_FOUND
};



class PathDereference {
	Path path;
	RecognisedFileType file_status;


public:

	PathDereference(const Path &p);

	Path getPath() const;
	string getFilename() const;
	string getExtension() const;
	bool exists() const;
	void copyTo(const Path &destination) const;
	void moveTo(const Path &destination) const;
	void rename(const string &newName) const;
	void remove() const;
	void resize(size_t newSize) const;
	size_t size() const;
	string to_str() const;
	bool isPrimitive() const;
};


class CodeFileDescriptor {
	PathDereference path;

public:
	CodeFileDescriptor(const string &filename, const string &directory, bool rootRelative = true);
	CodeFileDescriptor(const string &filename, const Path &directory, bool rootRelative = true);
	explicit CodeFileDescriptor(const Path &path, bool rootRelative = true);
	explicit CodeFileDescriptor(const string &path, bool rootRelative = true);
	virtual ~CodeFileDescriptor() = default;

	Path getPath() const;
	string getFilename() const;
	Path getDirectory() const;
	string extension() const;
	virtual string getCode();
	string readCode() const;
	bool exists() const;
	void writeCode(const string &code) const;
	void modifyCode(const string &code) const;
	void saveCopyToNewFile(const string &code) const;
	bool recogniseDirectoryNamingStyle();
	void changeDirectoryNamingStyle(bool unix = true);
	void changeLine(const string &line, int lineNumber);
	string readLine(int lineNumber) const;
};


class FileDescriptor {
	void *address = nullptr;
	void *fileHandle = nullptr;
	void *mappingHandle = nullptr;
	size_t bytesize;
	PathDereference path;


public:
	FileDescriptor(const PathDereference &filePath);
	FileDescriptor(const string &path) : FileDescriptor(PathDereference(path)) {}
	FileDescriptor(const char* path) : FileDescriptor(string(path)) {}

	FileDescriptor(const FileDescriptor &other) : FileDescriptor(other.path) {}
	FileDescriptor(FileDescriptor &&other) noexcept : FileDescriptor(std::move(other.path)) {}

	FileDescriptor &operator=(const FileDescriptor &other);
	FileDescriptor &operator=(FileDescriptor &&other) noexcept;

	~FileDescriptor();

	void mapFile();
	void closeFile();
	size_t getSize();
	void* getAddress();
	void flush() const;
	string getPath() const;

	string getExtension() const;
	string getFilename() const;
	void resize(size_t newSize);
	void addPaddingAtStart(size_t padding);
	void writeData(const void* data, size_t size, size_t offset);
	void rename(const string &newName);
	void removeDataFromStart(size_t deletedChunkSize);
	void readData(void* destination, size_t size, size_t offset);
};


bool isValidFilename(const string &filename);
bool isValidFilenameCharacter(char c);



class DirectoryDescriptor {
	Path path;

public:
	DirectoryDescriptor(const Path &p);
	DirectoryDescriptor(const string &p);
	Path getPath() const;
	bool exists() const;
	vector<PathDereference> listFiles() const;
	vector<DirectoryDescriptor> listDirectories() const;
	vector<PathDereference> listAllFilesRecursively() const;
	size_t directorySize() const;
	bool operator==(const DirectoryDescriptor &other) const;
};
