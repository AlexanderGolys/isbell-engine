#pragma once

#include <numeric>

#include "metaUtils.hpp"

typedef std::filesystem::path Path;


enum DirectoryEntryType {
	REG_FILE,
	DIR,
	SYM_LINK,
	UNRECOGNISED,
	EMPTY_PATH,
	FILE_NOT_FOUND
};


class DirectoryDescriptor;


class DirectoryEntry {
protected:
	Path path;
	DirectoryEntryType type;

public:
	DirectoryEntry(const Path &p);
	DirectoryEntry(const string &p);
	virtual ~DirectoryEntry() = default;
	DirectoryEntry(const DirectoryEntry &other);
	DirectoryEntry(DirectoryEntry &&other) noexcept;
	DirectoryEntry & operator=(const DirectoryEntry &other);
	DirectoryEntry & operator=(DirectoryEntry &&other) noexcept;

	virtual string getName() const;
	virtual bool exists();
	virtual DirectoryEntry copyTo(const Path &destination);
	virtual void moveTo(const Path &destination);
	virtual void rename(const string &newName);
	virtual void remove();

	virtual void refresh();
	virtual size_t getSize() const;

	bool operator==(const DirectoryEntry &other) const;
	string to_str() const;
	DirectoryDescriptor getParent() const;
	bool hasParent() const;

	bool isFile() const;
	bool isDirectory() const;

	filesystem::directory_entry getStdDirectoryEntry() const;
	filesystem::file_status getStdFileStatus() const;
	filesystem::perms getPermissions() const;

	Path getPath() const;
	DirectoryEntryType getType() const;
};


class FileDescriptor : public DirectoryEntry {
	void* address = nullptr;
	void* fileHandle = nullptr;
	void* mappingHandle = nullptr;
	size_t bytesize;

public:
	explicit FileDescriptor(const DirectoryEntry &filePath);
	explicit FileDescriptor(const Path &path) : FileDescriptor(DirectoryEntry(path)) {}
	explicit FileDescriptor(const string &path);
	explicit FileDescriptor(const char* path);
	~FileDescriptor() override;
	FileDescriptor &operator=(const FileDescriptor &other);
	FileDescriptor &operator=(FileDescriptor &&other) noexcept;

	void mapFile();
	void closeFile();
	void* getAddress();
	void flush() const;
	bool isMapped() const;


	bool exists() override;
	DirectoryEntry copyTo(const Path &destination) override;
	void moveTo(const Path &destination) override;
	void remove() override;
	size_t getSize() const override;

	string getExtension() const;
	void resize(size_t newSize);
	void refresh() override;

	void rename(const string &newName) override;

	void addPaddingAtStart(size_t padding);
	void writeData(const void* data, size_t size, size_t offset);
	void removeDataFromStart(size_t deletedChunkSize);
	void readData(void* destination, size_t size, size_t offset);
};


bool isValidFilename(const string &filename);
bool isValidFilenameCharacter(char c);


class DirectoryDescriptor : public DirectoryEntry {
public:
	using DirectoryEntry::DirectoryEntry;

	class FileIterator {
		filesystem::directory_iterator it;
		bool _points_to_file() const;
		void _make_end_independent();

	public:
		explicit FileIterator(); // end
		explicit FileIterator(const DirectoryDescriptor &dir);
		FileIterator& operator++();
		FileDescriptor operator*() const;
		bool operator!=(const FileIterator &other) const;
		bool isEnd() const;
	};
	FileIterator beginFiles() const;
	FileIterator EndFiles() const;

	class SubdirectoryIterator {
		filesystem::directory_iterator it;
		bool _points_to_dir() const;
		void _make_end_independent();

	public:
		explicit SubdirectoryIterator(); // end
		explicit SubdirectoryIterator(const DirectoryDescriptor &dir);
		SubdirectoryIterator &operator++();
		DirectoryDescriptor operator*() const;
		SubdirectoryIterator& reset();
		void inc();

		bool operator!=(const SubdirectoryIterator &other) const;
		bool isEnd() const;
	};
	SubdirectoryIterator beginDir() const;
	SubdirectoryIterator endDir() const;

	class RecursiveSubdirectoryIterator {
		shared_ptr<RecursiveSubdirectoryIterator> child;
		SubdirectoryIterator it;

		void _make_new_child();
		void _inc_child();
		void _inc_it();

	public:
		RecursiveSubdirectoryIterator() = default; // end
		explicit RecursiveSubdirectoryIterator(const DirectoryDescriptor &dir);

		RecursiveSubdirectoryIterator& operator++();
		void inc();
		DirectoryDescriptor operator*() const;
		bool operator!=(const RecursiveSubdirectoryIterator &other) const;

		bool isEnd() const;
		bool pointingToSelf() const;
		bool childEnded() const;
	};
	RecursiveSubdirectoryIterator beginRecursiveDir() const;
	RecursiveSubdirectoryIterator endRecursiveDir() const;

	class RecursiveFileIterator {
		RecursiveSubdirectoryIterator dir_it;
		FileIterator file_it;

	public:
		RecursiveFileIterator() = default;
		explicit RecursiveFileIterator(const DirectoryDescriptor &dir);

		RecursiveFileIterator &operator++();
		FileDescriptor operator*() const;
		bool operator!=(const RecursiveFileIterator &other) const;
		bool isEnd() const;
	};
	RecursiveFileIterator beginRecursiveFiles() const;
	RecursiveFileIterator endRecursiveFiles() const;

	RecursiveFileIterator begin() const;
	RecursiveFileIterator end() const;

	size_t getSize() const override;
};


class CodeFileDescriptor {
	FileDescriptor file;

public:
	CodeFileDescriptor(const string &filename, const string &directory, bool rootRelative = true);
	CodeFileDescriptor(const string &filename, const Path &directory, bool rootRelative = true);
	explicit CodeFileDescriptor(const Path &path, bool rootRelative = true);
	explicit CodeFileDescriptor(const string &path, bool rootRelative = true);
	virtual ~CodeFileDescriptor() = default;

	CodeFileDescriptor(const CodeFileDescriptor &other);
	CodeFileDescriptor(CodeFileDescriptor &&other) noexcept;
	CodeFileDescriptor & operator=(const CodeFileDescriptor &other);
	CodeFileDescriptor & operator=(CodeFileDescriptor &&other) noexcept;

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
	void changeLine(const string &line, int lineNumber);
	string readLine(int lineNumber) const;
};
