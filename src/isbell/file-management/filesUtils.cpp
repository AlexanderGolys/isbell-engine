#include "filesUtils.hpp"

#include <fstream>
#include <io.h>
#include <iostream>
#include <sstream>

#include "configFiles.hpp"
#include "exceptions.hpp"
#include "logging.hpp"
#include "metaUtils.hpp"
#include "randomUtils.hpp"

using namespace glm;




DirectoryEntry::DirectoryEntry(const Path &p)
: path(filesystem::absolute(p)) {
	path = path.make_preferred();
	DirectoryEntry::refresh();
}

DirectoryEntry::DirectoryEntry(const DirectoryEntry &other)
: path(other.path), type(other.type) {}

DirectoryEntry::DirectoryEntry(DirectoryEntry &&other) noexcept
: path(std::move(other.path)), type(other.type) {}

DirectoryEntry &DirectoryEntry::operator=(const DirectoryEntry &other) {
	if (this == &other) return *this;
	path = other.path;
	type = other.type;
	refresh();
	return *this;
}

DirectoryEntry &DirectoryEntry::operator=(DirectoryEntry &&other) noexcept {
	if (this == &other) return *this;
	path = std::move(other.path);
	type = other.type;
	refresh();
	return *this;
}

string DirectoryEntry::getName() const {
	if (type == EMPTY_PATH)
		return "";

	THROW_IF(not path.has_filename(), FileSystemError, "Path has no filename: " + path.string());

	return path.filename().string();
}

DirectoryEntry::DirectoryEntry(const string &p)
: DirectoryEntry(Path(p)) {}

DirectoryEntry::DirectoryEntry(): DirectoryEntry(filesystem::current_path()) {
	try {
		path = path.parent_path().parent_path().parent_path().parent_path();
		DirectoryEntry::refresh();
	} catch (...) {}
}

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
	filesystem::rename(path, destination / getName());
	setPath(destination / getName());
}

void DirectoryEntry::rename(const string &newName) {
	filesystem::rename(path, path.parent_path() / newName);
	setPath(path.parent_path() / newName);
}

void DirectoryEntry::remove() {
	filesystem::remove(path);
	type = FILE_NOT_FOUND;
}

void DirectoryEntry::refresh() {
	if (path == Path())
		type = EMPTY_PATH;
	else if (!filesystem::exists(path))
		type = FILE_NOT_FOUND;
	else {
		auto st = filesystem::symlink_status(path);
		if (st.type() == filesystem::file_type::symlink)
			type = SYM_LINK;
		else if (st.type() == filesystem::file_type::directory)
			type = DIR;
		else if (st.type() == filesystem::file_type::regular)
			type = REG_FILE;
		else
			type = UNRECOGNISED;
	}
}

size_t DirectoryEntry::getSize() const {
	THROW(NotImplementedError, "getSize is abstract.");
}

bool DirectoryEntry::operator==(const DirectoryEntry &other) const {
	return path == other.path;
}

bool DirectoryEntry::operator<=(const DirectoryEntry &other) const {
	return getPath().string().contains(other.getPath().string());
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

void DirectoryEntry::setPath(const Path &newPath) {
	path = newPath;
	path = filesystem::absolute(path).make_preferred();
	refresh();
}

bool DirectoryEntry::isFile() const {
	return type == REG_FILE or type == SYM_LINK;
}

bool DirectoryEntry::isDirectory() const {
	return type == DIR;
}

DirectoryEntry DirectoryEntry::operator/(const string &subpath) const {
	return DirectoryEntry(path / subpath);
}

DirectoryEntry DirectoryEntry::operator/(const char* subpath) const {
	return DirectoryEntry(path / string(subpath));
}

DirectoryEntry& DirectoryEntry::operator/=(const string &subpath) {
	return operator/=(Path(subpath));

}

DirectoryEntry& DirectoryEntry::operator/=(const Path &subpath) {
	setPath(path / subpath);
	return *this;}

DirectoryEntry& DirectoryEntry::operator/=(const char *subpath) {
	return operator/=(Path(subpath));
}

DirectoryEntry DirectoryEntry::operator/(const Path &subpath) const {
	return DirectoryEntry(path / subpath);
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

DirectoryDescriptor DirectoryDescriptor::operator/(const string &subpath) const {
	return DirectoryDescriptor(path / subpath);
}

DirectoryDescriptor DirectoryDescriptor::operator/(const Path &subpath) const {
	return DirectoryDescriptor(path / subpath);
}

DirectoryDescriptor DirectoryDescriptor::operator/(const char *subpath) const {
	return *this / string(subpath);
}

DirectoryDescriptor & DirectoryDescriptor::operator/=(const Path &subpath) {
	setPath(path / subpath);
	return *this;
}

DirectoryDescriptor & DirectoryDescriptor::operator/=(const string &subpath) {
	return *this /= Path(subpath);
}

DirectoryDescriptor & DirectoryDescriptor::operator/=(const char *subpath) {
	return *this /= string(subpath);
}

FileDescriptor DirectoryDescriptor::operator+(const Path &subpath) const {
	return FileDescriptor(*this / subpath);
}

FileDescriptor DirectoryDescriptor::operator+(const string &subpath) const {
	return *this + Path(subpath);
}

FileDescriptor DirectoryDescriptor::operator+(const char *subpath) const {
	return *this + string(subpath);
}



bool DirectoryDescriptor::FileIterator::_points_wrong_thing() const {
	return not isEnd() and not (it->is_regular_file() or it->is_symlink());
}

void DirectoryDescriptor::FileIterator::_make_end_independent() {
	it = filesystem::end(filesystem::directory_iterator());
}

void DirectoryDescriptor::FileIterator::_end_check() {
	if (isEnd())
		_make_end_independent();
}

bool DirectoryDescriptor::FileIterator::isEnd() const {
	return it == filesystem::end(it);
}

DirectoryDescriptor::FileIterator::FileIterator(const DirectoryDescriptor &dir)
: it(dir.getPath()) {
	while (_points_wrong_thing())
		++it;
	_end_check();
}

DirectoryDescriptor::FileIterator::FileIterator()
: it(filesystem::end(filesystem::directory_iterator())) {}

DirectoryDescriptor::FileIterator & DirectoryDescriptor::FileIterator::operator++() {
	THROW_IF(isEnd(), IteratorEndIncrementError);
	do ++it; while (_points_wrong_thing());
	_end_check();
	return *this;
}

void DirectoryDescriptor::FileIterator::inc() { operator++(); }

FileDescriptor DirectoryDescriptor::FileIterator::operator*() const {
	THROW_IF(isEnd(), IteratorEndReferenceError, "Dereferencing end iterator");
	return FileDescriptor(*it);
}

bool DirectoryDescriptor::FileIterator::operator!=(const FileIterator &other) const {
	if (isEnd() and other.isEnd())
		return false;
	return it != other.it;
}

DirectoryDescriptor::FileIterator DirectoryDescriptor::beginFiles() const {
	return FileIterator(*this);
}

DirectoryDescriptor::FileIterator DirectoryDescriptor::endFiles() const {
	return FileIterator();
}

bool DirectoryDescriptor::SubdirectoryIterator::_points_wrong_thing() const {
	return not isEnd() and not it->is_directory();
}

void DirectoryDescriptor::SubdirectoryIterator::_make_end_independent() {
	it = filesystem::end(filesystem::directory_iterator());
}

void DirectoryDescriptor::SubdirectoryIterator::_end_check() {
	if (isEnd())
		_make_end_independent();
}

bool DirectoryDescriptor::SubdirectoryIterator::isEnd() const {
	return it == filesystem::end(it);
}

DirectoryDescriptor::SubdirectoryIterator::SubdirectoryIterator(const DirectoryDescriptor &dir)
: it(dir.getPath())
{
	while (_points_wrong_thing())
		++it;
	_end_check();
}

DirectoryDescriptor::SubdirectoryIterator::SubdirectoryIterator()
: it(filesystem::end(filesystem::directory_iterator())) {}

DirectoryDescriptor::SubdirectoryIterator& DirectoryDescriptor::SubdirectoryIterator::operator++() {
	THROW_IF(isEnd(), IteratorEndIncrementError);
	++it;
	while (_points_wrong_thing())
		++it;
	_end_check();
	return *this;
}

DirectoryDescriptor DirectoryDescriptor::SubdirectoryIterator::operator*() const {
	THROW_IF(isEnd(), IteratorEndReferenceError);
	return DirectoryDescriptor(*it);
}

void DirectoryDescriptor::SubdirectoryIterator::inc() {
	operator++();
}


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
	THROW_IF(isEnd(), IteratorEndReferenceError, "Creating new child from end iterator");
	child = make_shared<RecursiveSubdirectoryIterator>(*it);
}


DirectoryDescriptor::RecursiveSubdirectoryIterator::RecursiveSubdirectoryIterator(const DirectoryDescriptor &dir)
: it(dir), child(nullptr) {}

DirectoryDescriptor::RecursiveSubdirectoryIterator & DirectoryDescriptor::RecursiveSubdirectoryIterator::operator++() {
	THROW_IF(isEnd(), IteratorEndIncrementError, "Incrementing terminal iterator is not allowed");

	if (child == nullptr)
		_make_new_child();

	else if (not child->isEnd())
		child->inc();

	if (child->isEnd()) {
		++it;
		child = nullptr;
	}
	return *this;
}

void DirectoryDescriptor::RecursiveSubdirectoryIterator::inc() {
	operator++();
}

DirectoryDescriptor DirectoryDescriptor::RecursiveSubdirectoryIterator::operator*() const {
	THROW_IF(isEnd(), IteratorEndReferenceError);
	if (child == nullptr)
		return *it;
	THROW_IF(child->isEnd(), IteratorEndReferenceError, "Child iterator is at end");
	return **child;
}

bool DirectoryDescriptor::RecursiveSubdirectoryIterator::operator!=(const RecursiveSubdirectoryIterator &other) const {
	if (child == nullptr and other.child == nullptr)
		return it != other.it;
	if (child == nullptr xor other.child == nullptr)
		return true;
	return *child != *other.child or it != other.it;
}

bool DirectoryDescriptor::RecursiveSubdirectoryIterator::isEnd() const {
	return it.isEnd();
}





DirectoryDescriptor::RecursiveSubdirectoryIterator DirectoryDescriptor::beginRecursiveDir() const {
	return RecursiveSubdirectoryIterator(*this);
}

DirectoryDescriptor::RecursiveSubdirectoryIterator DirectoryDescriptor::endRecursiveDir() const {
	return RecursiveSubdirectoryIterator();
}

void DirectoryDescriptor::RecursiveFileIterator::_next_dir() {
	THROW_IF(dir_it.isEnd(), IteratorEndIncrementError, "No more directories to iterate");

	if (rootDirIter)
		rootDirIter = false;
	else
		++dir_it;

	if (dir_it.isEnd()) {
		file_it = FileIterator();
		return;
	}

	file_it = FileIterator(*dir_it);

	if (file_it.isEnd())
		_next_dir();

}

DirectoryDescriptor::RecursiveFileIterator::RecursiveFileIterator(const DirectoryDescriptor &dir)
: dir_it(dir), file_it(dir), rootDirIter(true)
{
	if (file_it.isEnd() and not dir_it.isEnd())
		_next_dir();
}

DirectoryDescriptor::RecursiveFileIterator& DirectoryDescriptor::RecursiveFileIterator::operator++() {
	THROW_IF(isEnd(), IteratorEndIncrementError);
	++file_it;
	if (file_it.isEnd() and not dir_it.isEnd())
		_next_dir();
	return *this;
}

FileDescriptor DirectoryDescriptor::RecursiveFileIterator::operator*() const {
	THROW_IF(isEnd(), IteratorEndReferenceError);
	return *file_it;
}

bool DirectoryDescriptor::RecursiveFileIterator::operator!=(const RecursiveFileIterator &other) const {
	return file_it != other.file_it or dir_it != other.dir_it;
}

bool DirectoryDescriptor::RecursiveFileIterator::isEnd() const {
	return dir_it.isEnd() and file_it.isEnd();
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
	return std::accumulate(begin(), end(), size_t(0), [](size_t sum, const FileDescriptor &f){
		return sum + f.getSize();
	});
}




bool isValidFilename(const string &filename) {
	if (filename.empty()) return false;
	if (filename.size() > 255) return false;

	bool contains_dot = false;
	bool contains_sth_after_dot = false;

	for (char c: filename) {
		if (!isValidFilenameCharacter(c))
			return false;
		if (contains_dot)
			contains_sth_after_dot = true;
		if (c == '.')
			contains_dot = true;
	}

	return contains_sth_after_dot;
}
FileDescriptor::FileDescriptor(const Path &filePath)
: FileDescriptor(DirectoryEntry(filePath)) {}

FileDescriptor::FileDescriptor(const DirectoryEntry &filePath)
: DirectoryEntry(filePath), bytesize(0) {
	if (type == REG_FILE or type == SYM_LINK)
		bytesize = filesystem::file_size(path);
	else if (type != FILE_NOT_FOUND)
		THROW(FileSystemError, "Invalid file type");
	// if (not DirectoryEntry::exists())
	// 	THROW(FileNotFoundError, getPath().string());

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

FileDescriptor::FileDescriptor(const string &path)
: FileDescriptor(DirectoryEntry(path)) {}

FileDescriptor::FileDescriptor(const char *path)
: FileDescriptor(string(path)) {}

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
		refresh();
		if (other.isMapped())
			mapFile();
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
	return readCode();
}

string CodeFileDescriptor::readCode() const {
	std::ifstream stream(getPath().string(), std::ios::in);
	if (stream.is_open()) {
		string code;
		std::stringstream sstr;
		sstr << stream.rdbuf();
		code = sstr.str();
		stream.close();
		return code;
	}
	THROW(FileNotFoundError, getPath().string());
}

bool CodeFileDescriptor::exists() const {
	return filesystem::exists(getPath());
}

void CodeFileDescriptor::writeCode(const string &code) const {
	std::ofstream stream(getPath().string(), std::ios::out);
	if (stream.is_open()) {
		stream << code;
		stream.close();
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

CodeFileDescriptor::CodeFileDescriptor(const CodeFileDescriptor &other)
: file(other.file.getPath()) {}

CodeFileDescriptor::CodeFileDescriptor(CodeFileDescriptor &&other) noexcept
: file(std::move(other.file.getPath())) {}

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

CodeFileDescriptor & CodeFileDescriptor::operator=(const CodeFileDescriptor &other) {
	if (this == &other) return *this;
	file = other.file;
	return *this;
}

CodeFileDescriptor & CodeFileDescriptor::operator=(CodeFileDescriptor &&other) noexcept {
	if (this == &other) return *this;
	file = std::move(other.file);
	return *this;
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
