#ifndef INI_READER_HPP_
#define INI_READER_HPP_

#include <string>
#include <map>


namespace ufd{
    /*!
     * \brief Class that reads an ini file and provides basic methods to extract values
     *
     * All sections and names are not case sensitive. But the values are case sensitive except of the ones used for the \a getBoolean method. If the same name occurs multiple times in a section, the value of the last occurance is stored
     * and all others are discarded. If a line does not contain a '=' character this line is interpreted as part of a multiline value and therefore
     * appended to the value of the last name separated by a '\n'.
     * Only one name-value pair per line is allowed.
     * Do not use '=' as part of a value because this will be recognized as a parsing error and therefore an exception will be thrown. Currently there is no way to escape this. Of course the comment symbol neither can be part of a value.
     * If a value should be interpreted as number type, as many characters as possible are taken to form a valid number of that type.
     *
     * This class is at least partially based on the project on https://github.com/benhoyt/inih
     */
    class IniReader {
        public:
            /*!
             * \brief Constructor
             * @param filename the file to read from
             * @param comment symbol indicating comments in the file. Default ';'
             * @throws std::runtime_error Will be thrown if any error occurred while parsing the file
             */
            explicit IniReader(const std::string& filename, const std::string& comment=";");

            //! Destructor
            ~IniReader();

            /*
             * \brief returns error codes for occured errors
             * @return 0 on success, -1 on file open error or line number of first error on parse error
             */
            int parseError() const;

            /*!
             * \brief get a string value from an ini file exactly as it was found. 
             * @param section section of the name
             * @param name name whose value to get 
             * @param default_value value to return if name is not found
             * @return value of \a name in ini file or \default_value if it is not found
             */
            std::string getRawString(const std::string& section, const std::string& name, const std::string& default_value) const;

            /*!
             * \brief get a string value from an ini file. 
             * @param section section of the name
             * @param name name whose value to get 
             * @param default_value value to return if name is not found
             * @return value of \a name in ini file or \default_value if it is not found or is empty or contains only whitespaces
             */
            std::string getString(const std::string& section, const std::string& name, const std::string& default_value) const;

            /*!
             * \brief get a string value from an ini file in lower case. 
             * @param section section of the name
             * @param name name whose value to get 
             * @param default_value value to return if name is not found
             * @return value of \a name in ini file in lower case or \default_value if it is not found or is empty or contains only whitespaces
             */
            std::string getLowerString(const std::string& section, const std::string& name, const std::string& default_value) const;

            /*!
             * \brief get a integer value from an ini file. 
             * @param section section of the name
             * @param name name whose value to get 
             * @param default_value value to return if name is not found
             * @return value of \a name in ini file or \default_value if it is not found or not a valid integer
             */
            int getInteger(const std::string& section, const std::string& name, const int default_value) const;

            /*!
             * \brief get a double value from an ini file. 
             * @param section section of the name
             * @param name name whose value to get 
             * @param default_value value to return if name is not found
             * @return value of \a name in ini file or \default_value if it is not found or not a valid double
             */
            double getDouble(const std::string& section, const std::string& name, const double default_value) const;

            /*!
             * \brief get a boolean value from an ini file. 
             * @param section section of the name
             * @param name name whose value to get 
             * @param default_value value to return if name is not found
             * @return value of \a name in ini file or \default_value if it is not found or not a valid true/false value. 
             * Valid true values are "true", "on", "yes", "1" and valid false values are "false", "off", "no", "0"
             */
            bool getBoolean(const std::string& section, const std::string& name, const bool default_value) const;

            /*! 
             * \brief Checks whether the given name exists in the ini file
             * @param section section of the name
             * @param name name whose value to get 
             * @return true if the name exists 
             */
            bool valueExists(const std::string& section, const std::string& name) const;

        private: 
            /*!
             * \brief makes a key to store the data from the section and the name
             * @param section section of the name
             * @param name name  
             * @return "<section>::<name>"
             */
            static std::string makeKey(const std::string& section, const std::string& name);

            /*!
             * \brief parses the content of the ini file 
             * @param filename name of the ini file
             * @return error code (see parseError for more information)
             */
            int parseIni(const std::string& filename, const std::string& comment);

            /*!
             * \brief removes leading and trailing whitespaces from the str
             */
            std::string trim(const std::string& str, const std::string& whitespace=" \t");

            //! Stores the error code returned by the \a parseIni method
            int error_;

            //! Stores the name-value pairs read from the file in the form key-value with the key as returned from \a makeKey
            std::map<std::string, std::string> values_;
    };
}

#endif //INI_READER_HPP_
