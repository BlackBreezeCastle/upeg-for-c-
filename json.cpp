#include "json.h"

void parseSkip(const std::string &str, int &i)
{
	while (str[i] == ' ' ||str[i]== '\n' ||str[i]=='\r')
	{
		i++;
	}
}

std::string parseString(const std::string &str, int &i)
{
	i++;
	std::string ret="";
	while (str[i] != '"')
	{
		ret += str[i];
		i++;
	}
	i++;
	return ret;
}

std::string parseValue(const std::string &str, int &i)
{
	if (str[i] == '"')
	{
		return parseString(str, i);
	}
	else
	{
		std::string ret = "";
		while (str[i] != ',' && str[i] != ' ' && str[i] != '\n'&&str[i] != '}'&&str[i] != ']')
		{
			ret += str[i];
			i++;
		}
		return ret;
	}
}

void JsonObject::parseObject(const std::string &str, int &i, JsonObject*object)
{
	i++;
	while (str[i]!='}')
	{
		parseSkip(str, i);
		std::string key;
		if (str[i] == '"')
		{
			key = parseString(str, i);
		}
		else
		{
			//error
		}

		parseSkip(str, i);
		if (str[i] == ':')
		{
			i++;
		}
		else
		{
			//error
		}
		parseSkip(str, i);

		if (str[i] == '{')
		{
			JsonObject *child = new JsonObject();
			child->parseObject(str, i, child);
			object->m_children.insert(std::map<std::string, JsonObject *>::value_type(key, child));
		}
		else if (str[i] == '[')
		{
			JsonObject *child = new JsonObject();
			object->m_children.insert(std::map<std::string, JsonObject *>::value_type(key, child));
			child->parseArray(str, i, child);
		}
		else
		{
			JsonObject *child = new JsonObject();
			child->m_content = parseValue(str, i);
			object->m_children.insert(std::map<std::string, JsonObject *>::value_type(key, child));
		}

		parseSkip(str, i);
		if (str[i] == ',')
		{
			i++;
		}
	}
	i++;
}


void JsonObject::parseArray(const std::string &str, int &i, JsonObject*object)
{
	i++;
	parseSkip(str, i);
	if (str[i] == '{')
	{
		while (str[i] != ']')
		{
			parseSkip(str, i);
			if (str[i] != '{')
			{
				//error
			}
			JsonObject *child = new JsonObject();
			m_array.push_back(child);
			parseObject(str, i, child);
			parseSkip(str, i);
			if (str[i] == ',')
			{
				i++;
			}
		}
	}
	else
	{
		while (str[i] != ']')
		{
			parseSkip(str, i);
			JsonObject *child = new JsonObject();
			m_array.push_back(child);
			child->m_content=parseValue(str, i);
			parseSkip(str, i);
			if (str[i] == ',')
			{
				i++;
			}
		}
	}
	i++;
}

JsonObject::JsonObject()
{

};

JsonObject::JsonObject(std::string str)
{
	parse(str);
}

bool JsonObject::parse(std::string str)
{
	int i = 0;
	while (i < str.size())
	{
		parseSkip(str, i);
		if (str[i] == '{')
		{
			parseObject(str, i,this);
		}
	}
	return true;
}

JsonObject& JsonObject::operator[](int i)
{
	return *m_array[i];
}

JsonObject& JsonObject::operator[](std::string id)
{
	return *m_children[id];
}

double JsonObject::asDouble()
{
	return atof(m_content.c_str());
}

std::string JsonObject::asString()
{
	return m_content;
}

int JsonObject::size()
{
	return m_array.size();
}

void JsonObject::free(JsonObject*object)
{
	if (object != NULL)
	{
		for (auto mapIt = object->m_children.begin();
			mapIt != object->m_children.end();
			mapIt++)
		{
			free(mapIt->second);
		}

		for (int i = 0; i < object->m_array.size(); i++)
		{
			free(m_array[i]);
		}
		delete object;
		object = NULL;
	}
}