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




DirectoryEntry::DirectoryEntry(const string &p)
: path(filesystem::absolute(Path(p))) {
	if (p == "")
		type = EMPTY_PATH;
	else if (!filesystem::exists(p))
		type = FILE_NOT_FOUND;
	else {
		auto st = filesystem::symlink_status(p);
		if (st.type() == filesystem::file_type::symlink)
			type = SYM_LINK;
		else if (st.type() == filesystem::file_type::directory)
			type = DIR;
		else if (st.type() == filesystem::file_type::regular)
			type = REG_FILE;
		else type = UNRECOGNISED;
	}
}

DirectoryEntry::DirectoryEntry(const DirectoryEntry &other)
: path(other.path), type(other.type) {}

DirectoryEntry::DirectoryEntry(DirectoryEntry &&other) noexcept
: path(std::move(other.path)), type(other.type) {}

DirectoryEntry &DirectoryEntry::operator=(const DirectoryEntry &other) {
	if (this == &other) return *this;
	path = other.path;
	type = other.type;
	return *this;
}

DirectoryEntry &DirectoryEntry::operator=(DirectoryEntry &&other) noexcept {
	if (this == &other) return *this;
	path = std::move(other.path);
	type = other.type;
	return *this;
}

string DirectoryEntry::getName() const {
	return path.filename().string();
}

DirectoryEntry::DirectoryEntry(const Path &p)
: DirectoryEntry(p.string()) {}

Path DirectoryEntry::getPath() const {
	return path;
}

DirectoryEntryType DirectoryEntry::getType() const {
	return type;
}

bool DirectoryEntry::exists() {
	return filesystem::exists(path);
}

DirectoryEntry DirectoryEntry::copyTo(const Path &destination) {
	filesystem::copy(path, destination, filesystem::copy_options::overwrite_existing);
	return DirectoryEntry(destination);
}

void DirectoryEntry::moveTo(const Path &destination) {
	filesystem::rename(path, destination);
	path = destination;
}

void DirectoryEntry::rename(const string &newName) {
	filesystem::rename(path, path.parent_path() / newName);
	path = path.parent_path() / newName;
}

void DirectoryEntry::remove() {
	filesystem::remove(path);
	type = FILE_NOT_FOUND;
}

void DirectoryEntry::refresh() {
	if (path == "")
		type = EMPTY_PATH;
	else if (!filesystem::exists(path))
		type = FILE_NOT_FOUND;
	else {
		auto st = filesystem::symlink_status(path);
		if (st.type() == filesystem::file_type::symlink) type = SYM_LINK;
		else if (st.type() == filesystem::file_type::directory) type = DIR;
		else if (st.type() == filesystem::file_type::regular) type = REG_FILE;
		else type = UNRECOGNISED;
	}
}

size_t DirectoryEntry::getSize() const {
	THROW(NotImplementedError, "getSize is abstract.");
}

bool DirectoryEntry::operator==(const DirectoryEntry &other) const {
	return path == other.path;
}

string DirectoryEntry::to_str() const {
	return path.string();
}

DirectoryDescriptor DirectoryEntry::getParent() const {
	if (!hasParent())
		THROW(FileSystemError, "This path has no parent.");
	return DirectoryDescriptor(path.parent_path());
}

bool DirectoryEntry::hasParent() const {
	return path.has_parent_path();
}

bool DirectoryEntry::isFile() const {
	return type == REG_FILE or type == SYM_LINK;
}

bool DirectoryEntry::isDirectory() const {
	return type == DIR;
}

filesystem::directory_entry DirectoryEntry::getStdDirectoryEntry() const {
	return filesystem::directory_entry(path);
}


filesystem::file_status DirectoryEntry::getStdFileStatus() const {
	return filesystem::status(path);
}

filesystem::perms DirectoryEntry::getPermissions() const {
	return filesystem::status(path).permissions();
}






bool FileDescriptor::exists() {
	if (isMapped())
		return true;
	return DirectoryEntry::exists();
}

bool FileDescriptor::isMapped() const {
	return address != nullptr;
}

string FileDescriptor::getExtension() const {
	return getPath().extension().string();
}

void FileDescriptor::resize(size_t newSize) {
	bool wasMapped = address != nullptr;
	if (wasMapped)
		closeFile();
	filesystem::resize_file(getPath(), newSize);
	bytesize = newSize;
	if (wasMapped)
		mapFile();
}

void FileDescriptor::refresh() {
	bool mapped = isMapped();
	if (mapped)
		closeFile();
	DirectoryEntry::refresh();
	bytesize = exists() ? getSize() : 0;
	if (mapped)
		mapFile();
}

void FileDescriptor::addPaddingAtStart(size_t padding) {
	bool wasMapped = address != nullptr;
	size_t old_size = getSize();
	size_t new_size = old_size + padding;
	resize(new_size);
	if (!wasMapped) mapFile();
	memmove((char *) address + padding, address, old_size);
	memset(address, 0, padding);
	flush();
	if (!wasMapped) closeFile();
}

void FileDescriptor::writeData(const void *data, size_t size, size_t offset) {
	bool wasMapped = address != nullptr;
	if (!wasMapped) mapFile();
	if (offset + size > bytesize)
		throw ValueError("Writing data out of bounds: " + std::to_string(offset + size) + " > " + std::to_string(bytesize), __FILE__, __LINE__);
	memcpy((char *) address + offset, data, size);
	flush();
	if (!wasMapped) closeFile();
}

void FileDescriptor::rename(const string &newName) {
	bool wasMapped = address != nullptr;
	if (wasMapped) closeFile();
	DirectoryEntry::rename(newName);
	if (wasMapped) mapFile();
}

void FileDescriptor::removeDataFromStart(size_t deletedChunkSize) {
	bool wasMapped = address != nullptr;
	size_t old_size = getSize();
	size_t new_size = old_size - deletedChunkSize;
	if (!wasMapped) mapFile();
	memmove(address, (char *) address + deletedChunkSize, new_size);
	flush();
	if (!wasMapped) closeFile();
	resize(new_size);
}

void FileDescriptor::readData(void *destination, size_t size, size_t offset) {
	bool wasMapped = address != nullptr;
	if (!wasMapped) mapFile();
	if (offset + size > bytesize) throw ValueError("Reading data out of bounds: " + std::to_string(offset + size) + " > " + std::to_string(bytesize), __FILE__, __LINE__);
	memcpy(destination, (char *) address + offset, size);
	if (!wasMapped) closeFile();
}

bool isValidFilenameCharacter(char c) {
	const string valid_nonalphanumeric = "._- ";
	if (isalnum(c)) return true;
	for (char valid: valid_nonalphanumeric)
		if (c == valid) return true;
	return false;
}



bool DirectoryDescriptor::FileIterator::_points_to_file() const {
	return it->is_regular_file() or it->is_symlink();
}

void DirectoryDescriptor::FileIterator::_make_end_independent() {
	it = filesystem::end(filesystem::directory_iterator());
}

bool DirectoryDescriptor::FileIterator::isEnd() const {
	return it == filesystem::end(it);
}

DirectoryDescriptor::FileIterator::FileIterator(const DirectoryDescriptor &dir)
: it(dir.getPath()) {
	while (it != filesystem::end(it) and !_points_to_file())
		++it;
}

DirectoryDescriptor::FileIterator::FileIterator()
: it(std::filesystem::end(std::filesystem::directory_iterator())) {}

DirectoryDescriptor::FileIterator & DirectoryDescriptor::FileIterator::operator++() {
	if (isEnd())
		return *this;
	while (++it != filesystem::end(it) and !_points_to_file()) {}
	if (isEnd())
		_make_end_independent();
	return *this;
}

FileDescriptor DirectoryDescriptor::FileIterator::operator*() const {
	THROW_IF(isEnd(), IteratorEndReferenceError, "Dereferencing end iterator");
	return FileDescriptor(*it);
}

bool DirectoryDescriptor::FileIterator::operator!=(const FileIterator &other) const {
	return it != other.it;
}

DirectoryDescriptor::FileIterator DirectoryDescriptor::beginFiles() const {
	return FileIterator(*this);
}

DirectoryDescriptor::FileIterator DirectoryDescriptor::EndFiles() const {
	return FileIterator();
}

bool DirectoryDescriptor::SubdirectoryIterator::_points_to_dir() const {
	return it->is_directory();
}

void DirectoryDescriptor::SubdirectoryIterator::_make_end_independent() {
	it = filesystem::end(filesystem::directory_iterator());
}

bool DirectoryDescriptor::SubdirectoryIterator::isEnd() const {
	return it == filesystem::end(it);
}

DirectoryDescriptor::SubdirectoryIterator::SubdirectoryIterator(const DirectoryDescriptor &dir)
: it(dir.getPath()) {
	while (it != filesystem::end(it) and !_points_to_dir())
		++it;
	if (isEnd())
		_make_end_independent();
}

DirectoryDescriptor::SubdirectoryIterator::SubdirectoryIterator()
: it(std::filesystem::end(std::filesystem::directory_iterator())) {}

DirectoryDescriptor::SubdirectoryIterator & DirectoryDescriptor::SubdirectoryIterator::operator++() {
	while (++it != filesystem::end(it) and !_points_to_dir()) {}
	if (it == filesystem::end(it))
		_make_end_independent();
	return *this;
}

DirectoryDescriptor DirectoryDescriptor::SubdirectoryIterator::operator*() const {
	return DirectoryDescriptor(*it);
}

DirectoryDescriptor::SubdirectoryIterator & DirectoryDescriptor::SubdirectoryIterator::reset() {
	it = filesystem::begin(it);
	return *this;
}

void DirectoryDescriptor::SubdirectoryIterator::inc() { ++it; }


bool DirectoryDescriptor::SubdirectoryIterator::operator!=(const SubdirectoryIterator &other) const {
	return it != other.it;
}

DirectoryDescriptor::SubdirectoryIterator DirectoryDescriptor::beginDir() const {
	return SubdirectoryIterator(*this);
}

DirectoryDescriptor::SubdirectoryIterator DirectoryDescriptor::endDir() const {
	return SubdirectoryIterator();
}

void DirectoryDescriptor::RecursiveSubdirectoryIterator::_make_new_child() {
	child = make_shared<RecursiveSubdirectoryIterator>(*it);
}

void DirectoryDescriptor::RecursiveSubdirectoryIterator::_inc_child() {
	THROW_IF(child == nullptr, IteratorError, "Child not set, nullptr cannot be increased");
	child->inc();
}

void DirectoryDescriptor::RecursiveSubdirectoryIterator::_inc_it() {
	THROW_IF(it.isEnd(), IteratorEndIncrementError, "Incrementing terminal iterator is not allowed");
	++it;
}


DirectoryDescriptor::RecursiveSubdirectoryIterator::RecursiveSubdirectoryIterator(const DirectoryDescriptor &dir): it(dir), child(nullptr) {}

DirectoryDescriptor::RecursiveSubdirectoryIterator & DirectoryDescriptor::RecursiveSubdirectoryIterator::operator++() {
	THROW_IF(isEnd(), IteratorEndIncrementError, "Incrementing terminal iterator is not allowed");

	if (pointingToSelf())
		_make_new_child();

	else if (not childEnded())
		_inc_child();

	if (childEnded() and not isEnd()) {
		_inc_it();
		_make_new_child();
	}

	return *this;
}

void DirectoryDescriptor::RecursiveSubdirectoryIterator::inc() { operator++(); }

DirectoryDescriptor DirectoryDescriptor::RecursiveSubdirectoryIterator::operator*() const {
	THROW_IF(isEnd(), IteratorEndReferenceError);
	if (pointingToSelf())
		return *it;
	return **child;
}

bool DirectoryDescriptor::RecursiveSubdirectoryIterator::operator!=(const RecursiveSubdirectoryIterator &other) const {
	if (pointingToSelf() != other.pointingToSelf() or it != other.it)
		return true;
	if (pointingToSelf())
		return false;
	return *child != *(other.child);
}

bool DirectoryDescriptor::RecursiveSubdirectoryIterator::isEnd() const {
	return it.isEnd() and childEnded();
}

bool DirectoryDescriptor::RecursiveSubdirectoryIterator::pointingToSelf() const {
	return child == nullptr;
}

bool DirectoryDescriptor::RecursiveSubdirectoryIterator::childEnded() const {
	return not pointingToSelf() and child->isEnd();
}

DirectoryDescriptor::RecursiveSubdirectoryIterator DirectoryDescriptor::beginRecursiveDir() const { return RecursiveSubdirectoryIterator(*this); }

DirectoryDescriptor::RecursiveSubdirectoryIterator DirectoryDescriptor::endRecursiveDir() const { return RecursiveSubdirectoryIterator(); }

DirectoryDescriptor::RecursiveFileIterator::RecursiveFileIterator(const DirectoryDescriptor &dir): dir_it(dir), file_it(dir) {
	while (not dir_it.isEnd() and file_it.isEnd())
		if (not (++dir_it).isEnd())
			file_it = FileIterator(*dir_it);
}

DirectoryDescriptor::RecursiveFileIterator & DirectoryDescriptor::RecursiveFileIterator::operator++() {
	THROW_IF(isEnd(), IteratorEndIncrementError);
	++file_it;
	while (not dir_it.isEnd() and file_it.isEnd())
		if (not (++dir_it).isEnd())
			file_it = FileIterator(*dir_it);
	return *this;
}

FileDescriptor DirectoryDescriptor::RecursiveFileIterator::operator*() const {
	THROW_IF(isEnd(), IteratorEndReferenceError);
	return *file_it;
}

bool DirectoryDescriptor::RecursiveFileIterator::operator!=(const RecursiveFileIterator &other) const {
	return dir_it != other.dir_it or file_it != other.file_it;
}

bool DirectoryDescriptor::RecursiveFileIterator::isEnd() const {
	return file_it.isEnd();
}

DirectoryDescriptor::RecursiveFileIterator DirectoryDescriptor::beginRecursiveFiles() const {
	return RecursiveFileIterator(*this);
}

DirectoryDescriptor::RecursiveFileIterator DirectoryDescriptor::endRecursiveFiles() const {
	return RecursiveFileIterator();
}

DirectoryDescriptor::RecursiveFileIterator DirectoryDescriptor::begin() const {
	return beginRecursiveFiles();
}

DirectoryDescriptor::RecursiveFileIterator DirectoryDescriptor::end() const {
	return endRecursiveFiles();
}

size_t DirectoryDescriptor::getSize() const {
	return std::accumulate(begin(), end(), size_t(0), [](size_t sum, const FileDescriptor &f) { return sum + f.getSize(); });
}




bool isValidFilename(const string &filename) {
	if (filename.empty()) return false;
	if (filename.size() > 255) return false;

	bool contains_dot = false;
	bool contains_sth_after_dot = false;

	for (char c: filename) {
		if (!isValidFilenameCharacter(c)) return false;

		if (contains_dot) contains_sth_after_dot = true;

		if (c == '.') contains_dot = true;
	}

	return contains_sth_after_dot;
}

FileDescriptor::FileDescriptor(const DirectoryEntry &filePath)
: DirectoryEntry(filePath), bytesize(0) {
	if (type == REG_FILE or type == SYM_LINK)
		bytesize = filesystem::file_size(path);
	else if (type != FILE_NOT_FOUND)
		THROW(FileSystemError, "Invalid file type");
}

DirectoryEntry FileDescriptor::copyTo(const Path &destination) {
	bool wasMapped = isMapped();
	if (wasMapped)
		closeFile();
	DirectoryEntry copied = DirectoryEntry::copyTo(destination);
	if (wasMapped)
		mapFile();
	return copied;
}

FileDescriptor::FileDescriptor(const string &path) : FileDescriptor(DirectoryEntry(path)) {}

FileDescriptor::FileDescriptor(const char *path) : FileDescriptor(string(path)) {}

void *FileDescriptor::getAddress() {
	if (address == nullptr)
		mapFile();
	return address;
}


FileDescriptor &FileDescriptor::operator=(const FileDescriptor &other) {
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

FileDescriptor &FileDescriptor::operator=(FileDescriptor &&other) noexcept {
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








string CodeFileDescriptor::extension() const {
	return file.getExtension();
}

string CodeFileDescriptor::getCode() {
	std::ifstream shaderStream(getPath().string(), std::ios::in);
	if (shaderStream.is_open()) {
		string code;
		std::stringstream sstr;
		sstr << shaderStream.rdbuf();
		code = sstr.str();
		shaderStream.close();
		return code;
	}

	THROW(FileNotFoundError, getPath().string());
}

string CodeFileDescriptor::readCode() const {
	std::ifstream shaderStream(getPath().string(), std::ios::in);
	if (shaderStream.is_open()) {
		string code;
		std::stringstream sstr;
		sstr << shaderStream.rdbuf();
		code = sstr.str();
		shaderStream.close();
		return code;
	}

	THROW(FileNotFoundError, getPath().string());
}

bool CodeFileDescriptor::exists() const {
	return getPath().empty() ? false : filesystem::exists(getPath());
}

void CodeFileDescriptor::writeCode(const string &code) const {
	std::ofstream shaderStream(getPath().string(), std::ios::out);
	if (shaderStream.is_open()) {
		shaderStream << code;
		shaderStream.close();
		return;
	}

	THROW(FileNotFoundError, getFilename());
}

void CodeFileDescriptor::modifyCode(const string &code) const {
	if (!exists())
		THROW(FileNotFoundError, getFilename());
	writeCode(code);
}

void CodeFileDescriptor::saveCopyToNewFile(const string &code) const {
	if (exists())
		THROW(SystemError, "File " + getFilename() + " already exists in " + getDirectory().string() + ".");
	writeCode(code);
}

bool CodeFileDescriptor::recogniseDirectoryNamingStyle() {
	return getDirectory().string().find('/') != string::npos;
}

CodeFileDescriptor::CodeFileDescriptor(const string &filename, const Path &directory, bool rootRelative)
: file(directory / filename) {
	if (rootRelative)
		file = FileDescriptor(filesystem::absolute(directory / filename));
}

CodeFileDescriptor::CodeFileDescriptor(const Path &path, bool rootRelative)
: CodeFileDescriptor(path.filename().string(), path.parent_path(), rootRelative) {}

CodeFileDescriptor::CodeFileDescriptor(const string &filename, const string &directory, bool rootRelative)
: CodeFileDescriptor(filename, Path(directory), rootRelative) {}

CodeFileDescriptor::CodeFileDescriptor(const string &path, bool rootRelative)
: CodeFileDescriptor(Path(path), rootRelative) {}

void CodeFileDescriptor::changeLine(const string &line, int lineNumber) {
	std::ifstream file_stream(getPath().string());
	if (!exists())
		THROW(FileNotFoundError, getFilename());

	vector<string> lines;
	string current_line;
	while (std::getline(file_stream, current_line)) {
		lines.push_back(current_line);
	}
	file_stream.close();

	if (lineNumber < 0 || lineNumber >= lines.size())
		THROW(ValueError, "Line number out of range");

	lines[lineNumber] = line;

	string full_code;
	for (const auto& l : lines) {
		full_code += l + "\n";
	}

	writeCode(full_code);
}

string CodeFileDescriptor::readLine(int lineNumber) const {
	std::ifstream file_stream(getPath().string());
	if (!exists())
		THROW(FileNotFoundError, getFilename());

	string current_line;
	int i = 0;
	while (std::getline(file_stream, current_line)) {
		if (i == lineNumber) {
			file_stream.close();
			return current_line;
		}
		i++;
	}
	file_stream.close();
	THROW(ValueError, "Line number out of range");
}

Path CodeFileDescriptor::getPath() const {
	return file.getPath();
}

string CodeFileDescriptor::getFilename() const {
	return file.getName();
}

Path CodeFileDescriptor::getDirectory() const {
	return file.getPath().parent_path();
}
