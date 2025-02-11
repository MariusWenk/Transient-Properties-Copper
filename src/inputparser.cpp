#include "inputparser.hpp"
using namespace ufd;

Inputparser::Inputparser() {
    //__NOT_DEFINED_YET__;
}

Inputparser::Inputparser(std::string inifile, std::string default_inifile="") {
    reader_ = new IniReader(inifile);
    if(default_inifile==""){
        default_reader_=NULL;
    }
    else{
        default_reader_ = new IniReader(default_inifile);
    }
}

Inputparser::~Inputparser() {
    delete reader_;
    delete default_reader_;
}

double Inputparser::getDouble(const std::string& section, const std::string& name) const{
    //reads double variables from the .INI file and returns the entry in the given default file, if no entry is found in the starting.INI file
    double default_value = NAN;
    if (default_reader_ != NULL) {
        default_value = default_reader_ ->getDouble(section, name, default_value);
    }
    return reader_->getDouble(section, name, default_value);
}

std::string Inputparser::getString(const std::string& section, const std::string& name) const{
    //reads string variables from the .INI file and returns the entry in the given default file, if no entry is found in the starting.INI file
    std::string default_value = "";
    if (default_reader_ != NULL) {
        default_value = default_reader_ -> getRawString(section, name, default_value);
    }
    std::string value = reader_ -> getRawString(section, name, default_value);
    //converting respective string to lowercase
    std::transform(value.begin(), value.end(), value.begin(), ::tolower);
    return value;
}

std::string Inputparser::getOrigString(const std::string& section, const std::string& name) const{
    //reads string variables from the .INI file and returns the entry in the given default file, if no entry is found in the starting.INI file
    std::string default_value = "";
    if (default_reader_ != NULL) {
        default_value = default_reader_ -> getRawString(section, name, default_value);
    }
    std::string value = reader_ -> getRawString(section, name, default_value);
    return value;
}

int Inputparser::getInteger(const std::string& section, const std::string& name) const{
    //reads integer variables from the .INI file and returns the entry in the given default file, if no entry is found in the starting.INI file
    int default_value = -1;
    if (default_reader_ != NULL) {
        default_value = default_reader_ ->getInteger(section, name, default_value);
    }
    return reader_->getInteger(section, name, default_value);
}

bool Inputparser::getBoolean(const std::string& section, const std::string& name) const{
    //reads boolean variables from the .INI file and returns the entry in the given default file, if no entry is found in the starting.INI file
    bool default_value = false;
    if (default_reader_ != NULL) {
        default_value = default_reader_->getBoolean(section, name, default_value);
    }
    return reader_->getBoolean(section, name, default_value);
}
