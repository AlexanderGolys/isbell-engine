#include "filesUtils.hpp"

#include <fstream>
#include <io.h>
#include <iostream>
#include <sstream>

#include "configFiles.hpp"
#include "exceptions.hpp"
#include "metaUtils.hpp"
#include "randomUtils.hpp"

using namespace glm;



CodeFileDescriptor::CodeFileDescriptor(const string &filename, const Path &directory, bool rootRelative)
: path(directory, filename) {
	if (rootRelative)
		path = PathDereference(path.getPath().makeAbsolute());
}

CodeFileDescriptor::CodeFileDescriptor(const Path &path, bool rootRelative)
: CodeFileDescriptor(path.filename(), path.currentDirectory(), rootRelative) {}

CodeFileDescriptor::CodeFileDescriptor(const string &filename, const string &directory, bool rootRelative)
: CodeFileDescriptor(filename, Path(directory), rootRelative) {}

CodeFileDescriptor::CodeFileDescriptor(const string &path, bool rootRelative)
: CodeFileDescriptor(Path(path), rootRelative) {}


Path CodeFileDescriptor::getPath() const {
	return path.getPath();
}
string CodeFileDescriptor::getFilename() const {
	return path.getFilename();
}
Path CodeFileDescriptor::getDirectory() const {
	return path.getPath().currentDirectory();
}






Path::Path(const filesystem::path &path, bool extensionlessFile)
: path(path) {
	makePrefferedStyle();

	if (is_empty(path))
		throw ValueError("Empty paths are not supported", __FILE__, __LINE__);

	if (extensionlessFile) {
		if constexpr (IS_WINDOWS)
			THROW(FileSystemError, "Windows does not support extensionless files like (" + path.string() + ").");
		if (path.has_extension())
			THROW(FileSystemError, "Path is not extensionless: " + path.string());
		isExtensionlessFile = true;
	}
}




void Path::makePrefferedStyle()
{
	path = path.make_preferred();
}

Path Path::operator/(const string &other) const {
	return Path(path / other);
}

Path Path::operator+(const Path &other) const {
	return Path(path / other.path);
}


bool Path::isRoot() const
{
	return is_empty(path.parent_path());
}

Path Path::currentDirectory() const
{
	if (path.has_filename())
		return Path(path.parent_path());

	return *this;
}

string Path::filename() const
{
	if (path.has_filename())
		return path.filename().string();
	throw ValueError("Path does not end with a file", __FILE__, __LINE__);
}

string Path::fileExtension() const
{
	if (path.has_extension())
		return path.extension().string();

	throw ValueError("Path does not end with a file", __FILE__, __LINE__);
}

Path Path::makeRelative(const Path &other) const
{
	return Path(filesystem::relative(path, other.path));
}

Path Path::makeAbsolute() const
{
	return Path(filesystem::absolute(path));
}

Path Path::makeRelative(const string &other) const
{
	return makeRelative(Path(other));
}

Path Path::goUp() const
{
	if (isRoot())
		throw ValueError("Parent path is empty, as the path is already primitive (" + to_str() + ")", __FILE__, __LINE__);

	return Path(path.parent_path());
}

string Path::to_str() const
{
	return path.string();
}

bool Path::operator==(const Path &other) const
{
	return path == other.path;
}

PathDereference::PathDereference(const string &path)
: path(path)
{
	// 	throw ValueError("Path does not end with filename: " + path, __FILE__, __LINE__);
}

PathDereference::PathDereference(const Path &path)
: PathDereference(path.to_str())
{
}

PathDereference::PathDereference(const Path &directory, const string &filename)
: PathDereference(directory / filename)
{
}

PathDereference::PathDereference(const string &directory, const string &filename)
: PathDereference(Path(directory) / filename)
{
}

PathDereference::PathDereference(const filesystem::path &path)
: PathDereference(path.string())
{
}

PathDereference::PathDereference(const Path &p): path(p) {
	if (path == Path())
		file_status = EMPTY_PATH;
	else if (!filesystem::exists(path))
		file_status = FILE_NOT_FOUND;
	else {
		auto st = filesystem::symlink_status(path);
		if (st.type() == filesystem::file_type::symlink)
			file_status = REG_FILE;
		else if (st.type() == filesystem::file_type::directory)
			file_status = DIR;
		else if (st.type() == filesystem::file_type::regular)
			file_status = REG_FILE;
		else
			file_status = UNRECOGNISED;
	}
}

Path PathDereference::getPath() const
{
	return path;
}

string PathDereference::getFilename() const
{
	return path.filename();
}

string PathDereference::getExtension() const
{
	return path.fileExtension();
}

bool PathDereference::exists() const
{
	return filesystem::exists(path.path);
}

void PathDereference::copyTo(const Path &destination) const
{
	filesystem::copy(path.path, destination.path, filesystem::copy_options::overwrite_existing);
}

void PathDereference::moveTo(const Path &destination) const
{
	filesystem::rename(path.path, destination.path);
}

void PathDereference::rename(const string &newName) const
{
	filesystem::rename(path.path, path.currentDirectory().path / newName);
}

void PathDereference::remove() const
{
	filesystem::remove(path.path);
}

void PathDereference::resize(size_t newSize) const
{
	filesystem::resize_file(path.path, newSize);
}

size_t PathDereference::size() const
{
	return filesystem::file_size(path.path);
}

string PathDereference::to_str() const
{
	return path.to_str();
}

bool PathDereference::isPrimitive() const
{
	return path.isRoot();
}

// Path::operator string() const { return path; }
string FileDescriptor::getPath() const
{
	return path.to_str();
}

string FileDescriptor::getExtension() const
{
	return path.getExtension();
}

string FileDescriptor::getFilename() const
{
	return path.getFilename();
}

void FileDescriptor::resize(size_t newSize)
{
	bool wasMapped = address != nullptr;
	if (wasMapped)
		closeFile();
	path.resize(newSize);
	bytesize = newSize;
	if (wasMapped)
		mapFile();
}

void FileDescriptor::addPaddingAtStart(size_t padding)
{
	bool wasMapped = address != nullptr;
	size_t old_size = getSize();
	size_t new_size = old_size + padding;
	resize(new_size);
	if (!wasMapped)
		mapFile();
	memmove((char*)address + padding, address, old_size);
	memset(address, 0, padding);
	flush();
	if (!wasMapped)
		closeFile();
}

void FileDescriptor::writeData(const void *data, size_t size, size_t offset)
{
	bool wasMapped = address != nullptr;
	if (!wasMapped)
		mapFile();
	if (offset + size > bytesize)
		throw ValueError("Writing data out of bounds: " + std::to_string(offset + size) + " > " + std::to_string(bytesize), __FILE__, __LINE__);
	memcpy((char*)address + offset, data, size);
	flush();
	if (!wasMapped)
		closeFile();
}

void FileDescriptor::rename(const string &newName)
{
	bool wasMapped = address != nullptr;
	if (wasMapped)
		closeFile();
	path.rename(newName);
	if (wasMapped)
		mapFile();
}

void FileDescriptor::removeDataFromStart(size_t deletedChunkSize)
{
	bool wasMapped = address != nullptr;
	size_t old_size = getSize();
	size_t new_size = old_size - deletedChunkSize;
	if (!wasMapped)
		mapFile();
	memmove(address, (char*)address + deletedChunkSize, new_size);
	flush();
	if (!wasMapped)
		closeFile();
	resize(new_size);
}

void FileDescriptor::readData(void *destination, size_t size, size_t offset)
{
	bool wasMapped = address != nullptr;
	if (!wasMapped)
		mapFile();
	if (offset + size > bytesize)
		throw ValueError("Reading data out of bounds: " + std::to_string(offset + size) + " > " + std::to_string(bytesize), __FILE__, __LINE__);
	memcpy(destination, (char*)address + offset, size);
	if (!wasMapped)
		closeFile();
}

bool isValidFilenameCharacter(char c)
{
	const string valid_nonalphanumeric = "._- ";
	if (isalnum(c))
		return true;
	for (char valid : valid_nonalphanumeric)
		if (c == valid)
			return true;
	return false;
}

DirectoryDescriptor::DirectoryDescriptor(const Path &p)
: path(p)
{
	if (is_directory(path.path))
		THROW(FileSystemError, "Path " + path.to_str() + " points to a file, not a directory");
}

DirectoryDescriptor::DirectoryDescriptor(const string &p)
: DirectoryDescriptor(Path(p))
{
}

DirectoryDescriptor::DirectoryDescriptor(const filesystem::path &p)
: DirectoryDescriptor(Path(p))
{
}

Path DirectoryDescriptor::getPath() const
{
	return path;
}

	if (filesystem::is_directory(path.path))
{
	return filesystem::exists(path.path) && filesystem::is_directory(path.path);
}
bool DirectoryDescriptor::exists() const
vector<PathDereference> DirectoryDescriptor::listFiles() const
{
	vector<PathDereference> files;
	if (!exists())
		throw FileSystemError("Directory " + path.to_str() + " not found at this device.", __FILE__, __LINE__);

	for (const auto & entry : filesystem::directory_iterator(path.path))
		if (filesystem::is_regular_file(entry.status()) && !filesystem::is_empty(entry.path()))
			files.emplace_back(entry.path());

	return files;
}

vector<DirectoryDescriptor> DirectoryDescriptor::listDirectories() const
{
	vector<DirectoryDescriptor> dirs;
	if (path.isRoot())
		return dirs;
	if (!exists())
		THROW(FileSystemError, "Directory " + path.to_str() + " not found at this device.");
	for (const auto & entry : filesystem::directory_iterator(path.path))
		if (filesystem::is_directory(entry.status()) && !filesystem::is_empty(entry.path()))
			dirs.emplace_back(DirectoryDescriptor(Path(entry.path())));

	return dirs;
}

vector<PathDereference> DirectoryDescriptor::listAllFilesRecursively() const
{
	vector<PathDereference> files = listFiles();
	vector<DirectoryDescriptor> dirs = listDirectories();
	for (const DirectoryDescriptor &dir : dirs)
		addAll(files, dir.listAllFilesRecursively());

	return files;
}

size_t DirectoryDescriptor::directorySize() const
{
	vector<PathDereference> files = listAllFilesRecursively();
	return sum<size_t, vector<size_t>>(map<PathDereference, size_t>(files, [](const PathDereference &f) { return f.size(); }));
}

bool DirectoryDescriptor::operator==(const DirectoryDescriptor &other) const
{
	return path == other.path;
}

bool isValidFilename(const string &filename)
{
	if (filename.empty())
		return false;
	if (filename.size() > 255)
		return false;

	bool contains_dot = false;
	bool contains_sth_after_dot = false;

	for (char c : filename)
	{
		if (!isValidFilenameCharacter(c))
			return false;

		if (contains_dot)
			contains_sth_after_dot = true;

		if (c == '.')
			contains_dot = true;
	}

	return contains_sth_after_dot;
}

void CodeFileDescriptor::changeLine(const string &line, int lineNumber)
{
	std::ifstream file(getPath().to_str());
	if (!exists())
		throw FileNotFoundError(getFilename(), __FILE__, __LINE__);
	std::string code;
	int i = 0;
	while (std::getline(file, code))
	{
		if (i == lineNumber)
		{
			code = line;
			break;
		}
		i++;
	}
	file.close();
	writeCode(code);
}

string CodeFileDescriptor::readLine(int lineNumber) const
{
	std::ifstream file(getPath().to_str());
	if (!exists())
		throw FileNotFoundError(getFilename(), __FILE__, __LINE__);
	std::string code;
	int i = 0;
	while (std::getline(file, code))
	{
		if (i == lineNumber)
			break;
		i++;
	}
	file.close();
	return code;
}

FileDescriptor::FileDescriptor(const PathDereference &filePath)
: bytesize(0), path(filePath)
{
	if (!path.exists())
		throw FileNotFoundError(filePath.to_str(), __FILE__, __LINE__);
	bytesize = path.size();
}


string CodeFileDescriptor::extension() const
{
	return path.getExtension();
}

string CodeFileDescriptor::getCode()
{
	std::ifstream shaderStream(getPath().to_str(), std::ios::in);
	if (shaderStream.is_open())
	{
		string code;
		std::stringstream sstr;
		sstr << shaderStream.rdbuf();
		code = sstr.str();
		shaderStream.close();
		return code;
	}

	throw FileNotFoundError(getPath().to_str(), __FILE__, __LINE__);
}

string CodeFileDescriptor::readCode() const
{
	std::ifstream shaderStream(getPath().to_str(), std::ios::in);
	if (shaderStream.is_open())
	{
		string code;
		std::stringstream sstr;
		sstr << shaderStream.rdbuf();
		code = sstr.str();
		shaderStream.close();
		return code;
	}

	throw FileNotFoundError(getPath().to_str(), __FILE__, __LINE__);
}

bool CodeFileDescriptor::exists() const
{
	std::ifstream shaderStream(getPath().to_str(), std::ios::in);
	return shaderStream.is_open();
}

void CodeFileDescriptor::writeCode(const string &code) const
{
	std::ofstream shaderStream(getPath().to_str(), std::ios::out);
	if (shaderStream.is_open())
	{
		shaderStream << code;
		shaderStream.close();
		return;
	}

	throw FileNotFoundError(getFilename(), __FILE__, __LINE__);
}

void CodeFileDescriptor::modifyCode(const string &code) const
{
	if (!exists())
		throw FileNotFoundError(getFilename(), __FILE__, __LINE__);
	writeCode(code);
}

void CodeFileDescriptor::saveCopyToNewFile(const string &code) const
{
	if (exists())
		throw SystemError("File " + getFilename() + " already exists in " + getDirectory().to_str() + ".", __FILE__, __LINE__);
	writeCode(code);
}

bool CodeFileDescriptor::recogniseDirectoryNamingStyle()
{
	return getDirectory().to_str().find('/') != string::npos;
}

size_t FileDescriptor::getSize()
{
	return bytesize;
}

void* FileDescriptor::getAddress()
{
	if (address == nullptr)
		mapFile();
	return address;
}


FileDescriptor & FileDescriptor::operator=(const FileDescriptor &other) {
	if (this != &other) {
		closeFile();
		path = other.path;
		bytesize = other.bytesize;
		address = nullptr;
		fileHandle = nullptr;
		mappingHandle = nullptr;
	}
	return *this;
}

FileDescriptor & FileDescriptor::operator=(FileDescriptor &&other) noexcept {
	if (this != &other) {
		closeFile();
		path = std::move(other.path);
		bytesize = other.bytesize;
		address = nullptr;
		fileHandle = nullptr;
		mappingHandle = nullptr;
	}
	return *this;
}

FileDescriptor::~FileDescriptor() {
	closeFile();
}
