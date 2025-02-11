#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <exception>
#include "ini_reader.hpp"

using std::string;
using namespace ufd;

IniReader::IniReader(const string& filename, const string& comment) {
    error_ = parseIni(filename, comment);
    if (error_ != 0){
        string msg = "Error parsing file '" + filename + "' in line " + std::to_string(error_);
        throw std::runtime_error(msg);
    }
}   

IniReader::~IniReader(){
}

int IniReader::parseError() const {
    return error_;
}

string IniReader::getRawString(const string& section, const string& name, const string& default_value) const {
    string key = makeKey(section, name);
    return values_.count(key) ? values_.find(key)->second : default_value;
}

string IniReader::getString(const string& section, const string& name, const string& default_value) const {
    const string str = getRawString(section, name, "");
    return str.empty() ? default_value : str;
}

string IniReader::getLowerString(const string& section, const string& name, const string& default_value) const {
    string str = getRawString(section, name, "");
    std::transform(str.begin(), str.end(), str.begin(), ::tolower);
    return str.empty() ? default_value : str;
}

int IniReader::getInteger(const string& section, const string& name, const int default_value) const {
    string valstr =  getRawString(section, name, "");
    int value;
    try{
        value = std::stoi(valstr);
    }
    catch (...) {
        value = default_value;
    }
    return value;
}

double IniReader::getDouble(const string& section, const string& name, const double default_value) const {
    string valstr =  getRawString(section, name, "");
    double value;
    try{
        value = std::stod(valstr);
    }
    catch (...) {
        value = default_value;
    }
    return value;
}

bool IniReader::getBoolean(const string& section, const string& name, const bool default_value) const {
    string valstr = getRawString(section, name, "");
    // convert to lower case to make comparison case insensitive
    std::transform(valstr.begin(), valstr.end(), valstr.begin(), ::tolower);
    if (valstr == "true" || valstr == "yes" || valstr == "on" || valstr == "1")
        return true;
    else if (valstr == "false" || valstr == "no" || valstr == "off" || valstr == "0") 
        return false;
    else 
        return default_value;
}

bool IniReader::valueExists(const string& section, const string& name) const {
    string key = makeKey(section, name);
    return values_.count(key);
}

string IniReader::makeKey(const string& section, const string& name){
    string key = section + "::" + name;
    // Convert to lower case to make section/name lookup case-insensitive
    std::transform(key.begin(), key.end(), key.begin(), ::tolower);
    return key;
}

int IniReader::parseIni(const string& filename, const string& comment){
    int error = 0;

    std::ifstream file(filename);
    string line, name = "", value, section = "", key;
    size_t pos, pos2;
    uint line_number = 0;

    if (file.is_open()){
        try{
            while (std::getline(file, line)){ 
                line_number++;
                pos = line.find(comment); // find position of comment
                if (pos != string::npos)
                    line.erase(pos, string::npos); // erase comment
               
                if (!line.empty() && line.find_first_not_of(" \t") != string::npos){
                    // search for sections
                    pos = line.find("[");
                    pos2 = line.find("]", pos);
                    if (pos != string::npos && pos2 != string::npos){
                        pos += 1; // skip '['
                        pos2 -= 1; // skip ']'
                        section = line.substr(pos, pos2 - pos + 1);
                        section = trim(section);
                    }
                    // if a section was found there won't be searched for name value pairs in the same line
                    else {
                        if (line.find("=") == string::npos){
                            // if no '=' is found this is interpreted as multiline value
                            values_[key] += "\n" + trim(line);
                        }
                        else {
                            std::istringstream ss(line);
                            std::getline(ss, name, '=');
                            name = trim(name);
                            std::getline(ss, value);
                            if (value.find("=") != string::npos){
                                // there is a second name-value pair in the same line or the value contains a '=' which is not allowed
                                // so the line number where this occured is returned as error
                                return line_number;
                            }
                            value = trim(value);

                            key = makeKey(section, name);
                            values_[key] = value; // potential previous values are discarded
                        }
                    }
                }
            }
        }
        catch (...){
            // if any exception is thrown the line number of the line where the error occured is returned
            return line_number;
        }
        file.close();
    }
    else {
        error = -1;
    }

    return error;
}

std::string IniReader::trim(const string& str, const string& whitespace){
    const auto begin = str.find_first_not_of(whitespace);
    if (begin == string::npos)
        return ""; // no content

    const auto end = str.find_last_not_of(whitespace);
    return str.substr(begin, end - begin + 1);
}



