#include "basic_inputparser.hpp"
//using namespace ufd::BasicInputparser;
namespace ufd::BasicInputparser{
std::string deleteComment(std::string str){
    uint found = str.find_first_of("#");\
                        return str.substr(0,found);
}

std::string deleteWhitespaces(std::string str){
    str.erase( std::remove_if( str.begin(), str.end(), ::isspace ), str.end() );
    return str;
}

std::vector<double> readColumn(std::string inputname, uint column_number){
    std::vector<double> output;
    std::string line;
    std::ifstream inputfile(inputname.c_str());
    if (inputfile.is_open()){
        while ( getline(inputfile,line))
        {
            std::istringstream ss(line);
            double value;
            uint counter = 0;
            while ((ss >> value) && (counter < column_number)){
                counter++;
            }
            output.push_back(value);
        }
    }

    return output;
}

bool stringToBool(std::string s){
    if (s.compare("false")==0 || s.compare("off")==0 || s.compare("0")==0){
        return false;
    }
    if (s.compare("true")==0 || s.compare("on")==0 || s.compare("1")==0){
        return true;
    }

    std::cerr << "Error in stringToBool: " << s << "has non-valid type" << std::endl;
    exit(999);
}

double stringToDouble(std::string s){
    std::stringstream ss(s);
    double number;
    ss >> number;
    return number;
}

int stringToInt(std::string s){
    std::stringstream ss(s);
    int number;
    ss >> number;
    return number;
}

std::vector<std::string> split(std::string s, char delim){
    std::vector<std::string> elems;
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}
}
