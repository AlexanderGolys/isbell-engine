#pragma once
#include "unittests.hpp"
#include "isbell.hpp"


inline bool descriptorsTest() {
	bool passed = true;
	DirectoryDescriptor d = DirectoryDescriptor();
	passed &= assertEqual_UT(d.getName(), "isbell");
	d /= "src/tests/test_dir/d1";
	FileDescriptor f2 = d + "f2.txt";
	return passed;
}

inline bool filesystemIteratorsTest()
{
	bool passed = true;
	DirectoryDescriptor d = DirectoryDescriptor() / "src" / "tests" / "test_dir";

	int files = 0;
	int dirs = 0;
	int rec_dirs = 0;
	int rec_files = 0;

	for (auto it = d.beginFiles(); it != d.endFiles(); ++it)
		++files;
	for (auto it = d.beginDir(); it != d.endDir(); ++it)
		++dirs;
	for (auto it = d.beginRecursiveDir(); it != d.endRecursiveDir(); ++it)
		++rec_dirs;
	for (auto _ : d)
		++rec_files;


	passed &= assertEqual_UT(files, 1);
	passed &= assertEqual_UT(dirs, 2);
	passed &= assertEqual_UT(rec_dirs, 7);
	passed &= assertEqual_UT(rec_files, 10);

	return passed;
}

inline UnitTestResult filesystemTests__all()
{

	UnitTestResult result;

	result.runTest(descriptorsTest);
	result.runTest(filesystemIteratorsTest);

	return result;
}
