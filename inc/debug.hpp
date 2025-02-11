#ifndef DEBUG_HPP_
#define DEBUG_HPP_

#include <iostream>
#include <sstream>
#include <string>
#include <fstream>

/*!
 * \file
 * \brief Provides debug functions with different levels.
 * \details To use the debuging functionality of the librethfeld the methods
 * \a arDebug, \a arInfo, \a arWarning, \a arError and \a arFatal are available. 
 * Any of them can be used to stream debug messages into or pass a message as parameter like in 
 * a normal function. 
 * \warning Note that \a std::endl is not accepted in streams! All messages will automatically be terminated with std::endl.
 * \details To enable or disable a certain debug level pass one of the following arguments to the compiler:
 * FINE, DEBUG, INFO, WARNING, ERROR, FATAL, NO_DEBUG_OUTPUT, NO_INFO_OUTPUT, NO_WARNING_OUTPUT, NO_E
 * , NO_FATAL_OUTPUT.
 * These levels are hierachically ordered wich leads to the suppression of output for all levels
 * below the specified one. E.g. if \a -DWARNING was passed to the compiler, output of WARNING, ERROR and FATAL will be printed but not of INFO, DEBUG and FINE. 
 * The NO_<X>_OUTPUT commands work similarly.
 * If no level was passed to the compiler always the lowest (FINE) is assumed. 
 * If additionally the DEBUG_TO_FILE flag is passed to the compiler, all messages 
 * will be written into a file called \a debug.hpp located at var/log . This file is never 
 * overwritten but all content will just be appended. 
 *
 * This file also provides additional macros for debuging. 
 * These are __INTERFACE_CLASS__, __WRONG_FUNCTION__, __NOT_DEFINED_YET__, EXIT(X), ISSUE(X) and WHERE_AM_I. 
 * The WHERE_AM_I macro prints the name and signature of a function, but only when the FUNC_INFO flag was passed to the compiler.
 *
 * This is basically a simple version of the Qt debugger so some ideas were taken from qlooging.h/.cpp and qdebug.h/.cpp.
 */

namespace ufd {

    class Debug  {
        public:
            enum MsgType{fatal, error, warning, info, debug, fine};
            Debug(MsgType type);
            ~Debug();

            template<typename T>
            inline Debug& operator<<(T t) {*stream << t; return *this;}; 

        private:
            MsgType type_;
            std::stringstream* stream;

    };

    class NoDebug {
        public:
            NoDebug(){}
            ~NoDebug(){}

            template<typename T>
            inline NoDebug& operator<<(T t){return *this;};
    };

    class MessageHandler {
        public:
            //! Constructor
            MessageHandler();
                
            //! Destructor
            ~MessageHandler();

            /*!
             * \brief creates a Debug object to stream messages into
             * If the level FINE is enabled the stream will be printed
             */
            Debug fine() const;

            //! prints \a msg if the level FINE is enabled
            void fine(const char* msg) const;

            /*!
             * \brief creates a Debug object to stream messages into
             * If the level DEBUG is enabled the stream will be printed
             */
            Debug debug() const;

            //! prints \a msg if the level DEBUG is enabled
            void debug(const char* msg) const;

            /*!
             * \brief creates a Debug object to stream messages into
             * If the level INFO is enabled the stream will be printed
             */
            Debug info() const;
            
            //! prints \msg if the level INFO is enabled
            void info(const char* msg) const;

            /*!
             * \brief creates a Debug object to stream messages into
             * If the level WARNING is enabled the stream will be printed
             */
            Debug warning() const;
            
            //! prints \msg if the level WARNING is enabled
            void warning(const char* msg) const;

            /*!
             * \brief creates a Debug object to stream messages into
             * If the level ERROR is enabled the stream will be printed
             */
            Debug error() const;
            
            //! prints \a msg if the level ERROR is enabled
            void error(const char* msg) const;


            /*!
             * \brief creates a Debug object to stream messages into
             * If the level FATAL is enabled the stream will be printed
             */
            Debug fatal() const;
            
            //! prints \msg if the level FATAL is enabled
            void fatal(const char* msg) const;

            /*!
             * \brief returns a noOutput object to stream messages into
             * Nothing will be printed so this is used for all debug functions below the current debug level
             */
            NoDebug noOutput() const;

            //! the message \a msg is discarded
            void noOutput(const char* msg) const {}


    };



}

#define arFine ufd::MessageHandler().fine
#define arDebug ufd::MessageHandler().debug
#define arInfo ufd::MessageHandler().info
#define arWarning ufd::MessageHandler().warning
#define arError ufd::MessageHandler().error
#define arFatal ufd::MessageHandler().fatal

#define NO_OUTPUT_MACRO while (false) ufd::MessageHandler().noOutput

#if defined (FATAL)
    #define NO_ERROR_OUTPUT 
#endif

#if defined (ERROR)
    #define NO_WARNING_OUTPUT 
#endif

#if defined (WARNING)
    #define NO_INFO_OUTPUT 
#endif

#if defined (INFO)
    #define NO_DEBUG_OUTPUT 
#endif

#if defined (DEBUG)
    #define NO_FINE_OUTPUT
#endif 

#if defined (NO_FATAL_OUTPUT)
    #define NO_ERROR_OUTPUT
    #undef arFatal
    #define arFatal NO_OUTPUT_MACRO
#endif

#if defined (NO_ERROR_OUTPUT)
    #define NO_WARNING_OUTPUT
    #undef arError
    #define arError NO_OUTPUT_MACRO
#endif

#if defined (NO_WARNING_OUTPUT)
    #define NO_INFO_OUTPUT
    #undef arWarning
    #define arWarning NO_OUTPUT_MACRO
#endif

#if defined (NO_INFO_OUTPUT)
    #define NO_DEBUG_OUTPUT
    #undef arInfo
    #define arInfo NO_OUTPUT_MACRO
#endif

#if defined (NO_DEBUG_OUTPUT)
    #define NO_FINE_OUTPUT
    #undef arDebug
    #define arDebug NO_OUTPUT_MACRO
#endif

#if defined (NO_FINE_OUTPUT)
    #undef arFine
    #define arFine NO_OUTPUT_MACRO
#endif
 
// BASH PS1 colors 
// It changes/resets the color of all following output
// See https://wiki.archlinux.org/index.php/Bash/Prompt_customization
// Better explanation (but german) https://wiki.archlinux.de/title/Bash-Prompt_anpassen

const char str_clear  []="\033[0m";
const char str_black  []="\033[30m";
const char str_red    []="\033[31m";
const char str_green  []="\033[32m";
const char str_yellow []="\033[33m";
const char str_blue   []="\033[34m";
const char str_purple []="\033[35m";
const char str_cyan   []="\033[36m";
const char str_white  []="\033[37m";

// Predefined error messages with exit code
// Should be used when class is an interface and shouldn't be called directly
// __PRETTY_FUNCTION__ shows the name and signature of the function ("Where to find it")
#define __INTERFACE_CLASS__\
    (\
        (arError() << str_red << "The code shouldn't call function "  << __PRETTY_FUNCTION__  << ".\n It is an interface class. Call one of the derived classes." << str_clear),\
        (exit(1))\
    )
// Predefined error messages with exit code
// Error message, that you have called a function, that should never be called.
// __PRETTY_FUNCTION__ shows the name and signature of the function ("Where to find it")
#define __WRONG_FUNCTION__\
    (\
        (arError() << str_red << "The code shouldn't call function "  << __PRETTY_FUNCTION__  << str_clear),\
        (exit(1))\
    )
// Predefined error messages with exit code
// Error message for functions, that have to exist, but are not implemented yet.
// __PRETTY_FUNCTION__ shows the name and signature of the function ("Where to find it")
#define __NOT_DEFINED_YET__\
    (\
        (arError() << str_red << "The function "  << __PRETTY_FUNCTION__  << " is not defined yet." << str_clear),\
        (exit(2))\
    )
// Kills code with an error code and a message, where it died
// __PRETTY_FUNCTION__ shows the name and signature of the function ("Where to find it")
#define EXIT(X)\
    (\
        (arFatal() << str_red << X << " \n in function "  << __PRETTY_FUNCTION__ << str_clear),\
        (exit(3))\
    )

// the gitlab issue is linked. 
// __PRETTY_FUNCTION__ shows the name and signature of the function ("Where to find it") 
#define ISSUE(X)\
	(\
         (arInfo() << str_green << "Further reading in https://gitlab.rlp.net/ag-rethfeld/programs/monstr/monstr/issues/" << X << str_clear)\
     )    
//These function are only defined, if the FUNC_INFO flag was set
#ifdef FUNC_INFO 
    // Writes the name of the function that was called.
    // __PRETTY_FUNCTION__ shows the name and signature of the function ("Where to find it")
    #define WHERE_AM_I\
        (\
            (arDebug() << str_purple << "I'm here: " << __PRETTY_FUNCTION__ << str_clear)\
        )
#else
    //If the FUNC_INFO flag is not set, the function should be quiet
    #define WHERE_AM_I
#endif
 
 
 
 
 
 
 #endif // DEBUG_HPP_ 
 
 
 
 
 
 
 
 
