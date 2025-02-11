#include "debug.hpp"

using namespace ufd;

MessageHandler::MessageHandler(){
}

MessageHandler::~MessageHandler(){
}

Debug MessageHandler::fine() const {
    return Debug(Debug::MsgType::fine); 
}

void MessageHandler::fine(const char* msg) const {
    fine() << msg; 
}

Debug MessageHandler::debug() const {
    return Debug(Debug::MsgType::debug); 
}

void MessageHandler::debug(const char* msg) const {
    debug() << msg; 
}

Debug MessageHandler::info() const {
    return Debug(Debug::MsgType::info); 
}

void MessageHandler::info(const char* msg) const {
    info() << msg; 
}
Debug MessageHandler::warning() const {
    return Debug(Debug::MsgType::warning); 
}

void MessageHandler::warning(const char* msg) const {
    warning() << msg; 
}
Debug MessageHandler::error() const {
    return Debug(Debug::MsgType::error); 
}

void MessageHandler::error(const char* msg) const {
    error() << msg; 
}
Debug MessageHandler::fatal() const {
    return Debug(Debug::MsgType::fatal); 
}

void MessageHandler::fatal(const char* msg) const {
    fatal() << msg; 
}
NoDebug MessageHandler::noOutput() const {
    return NoDebug();
}

Debug::Debug(MsgType type){
    type_ = type;
    stream = new std::stringstream();
}

Debug::~Debug(){
    std::string type_string = "";
    std::string color_string = "";
    switch (type_){
        case fine: type_string = "[ FINE ] "; color_string = str_white; break;
        case debug: type_string = "[ DEBUG ] "; color_string = str_cyan; break;
        case info: type_string = "[ INFO  ] "; color_string = str_green; break;
        case warning: type_string = "[ WARNING ] "; color_string = str_yellow; break;
        case error: type_string = "[ ERROR ] "; color_string = str_red; break;
        case fatal: type_string = "[ FATAL ] "; color_string = str_red; break;
    }
    std::string msg = stream->str();
    std::cout << color_string << type_string << str_clear << msg << std::endl;
    delete(stream);
    #if defined (DEBUG_TO_FILE)
        // save debug messages to log file
        std::fstream file("var/log/debug.log", std::ios::out | std::ios::app);
        file << type_string << msg << std::endl;
        file.close();
    #endif
}

