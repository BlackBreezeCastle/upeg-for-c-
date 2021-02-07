#include<map>
#include<vector>
#include<string>

class JsonObject
{
public:
	JsonObject();

	JsonObject(std::string str);

	bool parse(std::string str);

	JsonObject& operator[](int i);

	JsonObject& operator[](std::string id);

	double asDouble();

	std::string asString();

	int size();
private:
	void free(JsonObject*p);

	void parseObject(const std::string &str, int &i, JsonObject*object);

	void parseArray(const std::string &str, int &i, JsonObject*object);
private:
	std::map<std::string, JsonObject*> m_children;

	std::vector<JsonObject*> m_array;

	std::string m_content;

};
