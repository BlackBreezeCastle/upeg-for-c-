#pragma once
#ifndef PY_BUILD
#include "json.h"
#else
#include"json.cpp"
#endif
class Config
{
private:
	JsonObject root;
private:
	bool ReadFile(std::string FileName, std::string &out)
	{
		out = "";
		FILE *fp = NULL;
		fopen_s(&fp, FileName.c_str(), "rb");
		if (NULL == fp)
		{
			return false;
		}
		fseek(fp, 0L, SEEK_END);
		long length = ftell(fp);
		char *pContent = new char[length];
		fseek(fp, 0L, SEEK_SET);
		fread(pContent, length, 1, fp);
		fclose(fp);
		std::string tmp(pContent, length);
		out += tmp;
		delete pContent;
		return true;
	}
public:
	Config(std::string config_file_name)
	{
		std::string configStr;
		ReadFile(config_file_name, configStr);
		root.parse(configStr.c_str());
	}
	inline JsonObject operator[](std::string str)
	{
		return root[str];
	}

	inline JsonObject operator[](int i)
	{
		return root[i];
	}
};