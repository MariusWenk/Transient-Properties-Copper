#ifndef INPUTPARSER_H
#define INPUTPARSER_H

#include <iostream>
#include <string>
#include <algorithm>
#include <cmath>
#include "ini_reader.hpp"

namespace ufd{
    /*! \~ \brief class to read .ini files for further use in the simulation

        uses a starting and a default .ini-File
        uses INIReader.h
        values are saved in sections with names
        does not differentiate between lower or upper case letters for section and names
    */
    class Inputparser{
        public:

            //! \~ constructor (does nothing)
            Inputparser();
 
            /*! \~ \brief constructor

                if no default.ini exists use "" and hardcodet default_values are used
                @param inifile .ini file with desired values
                @param default_inifile .ini file with default values
                \sa reader_
                \sa default_reader_
             */
            Inputparser(std::string inifile, std::string default_inifile);

            //! \~ destructor (does nothing)
            ~Inputparser();

            /*! \~ \brief reads a double 

                reads an entry from the \a inifile with the given \a name in the
                given \a section and returns it as a double variable
                if no entry was found, uses the entry in the \a default_inifile file in the 
                \a default_reader_
                if no \a default_reader_ was found, uses NAN 
                @param section
                @param name
                @return double from .ini file
             */
            double getDouble(const std::string& section, const std::string& name) const;

             /*! \~ \brief reads a string 

                reads an entry from the \a inifile with the given \a name in the
                given \a section and returns it as a string variable
                if no entry was found, uses the entry in the \a default_inifile file in the 
                \a default_reader_
                if no \a default_reader_ was found, uses ""
                converts upper to lower case
                @param section
                @param name
                @return string from .ini file
                \sa getOrigString
             */
            std::string getString(const std::string& section, const std::string& name) const;

             /*! \~ \brief reads a string 

                reads an entry from the \a inifile with the given \a name in the
                given \a section and returns it as a string variable
                if no entry was found, uses the entry in the \a default_inifile file in the 
                \a default_reader_
                if no \a default_reader_ was found, uses ""
                does not convert upper to lower case
                @param section
                @param name
                @return string from .ini file
                \sa getString
             */
            std::string getOrigString(const std::string& section, const std::string& name) const;

             /*! \~ \brief reads an Integer 

                reads an entry from the \a inifile with the given \a name in the
                given \a section and returns it as an integer variable
                if no entry was found, uses the entry in the \a default_inifile file in the 
                \a default_reader_
                if no \a default_reader_ was found, uses -1
                @param section
                @param name
                @return Integer from .ini file
             */
            int getInteger(const std::string& section, const std::string& name) const;

             /*! \~ \brief reads a Boolean 

                reads an entry from the \a inifile with the given \a name in the
                given \a section and returns it as a string variable
                if no entry was found, uses the entry in the \a default_inifile file in the 
                \a default_reader_
                if no \a default_reader_ was found, uses \a false
                @param section
                @param name
                @return boolean from .ini file
             */
            bool getBoolean(const std::string& section, const std::string& name) const;


        private:
            /*! \~ \brief reader which reads the \a inifile

                reads the inifile that is specified in the constructor
                \sa Inputparser(std::string inifile, std::string default_inifile)
            */
            IniReader* reader_;
            
            /*! \~ \brief reader which reads the \a default_inifile

                reads the default_inifile that is specified in the constructor
                \sa Inputparser(std::string inifile, std::string default_inifile)
            */
            IniReader* default_reader_;

    };
}

#endif /* INPUTPARSER_H */
