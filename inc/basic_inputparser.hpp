#ifndef BASIC_INPUTPARSER_H
#define BASIC_INPUTPARSER_H

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>

namespace ufd{
    namespace BasicInputparser{

    /*! \file basic_inputparser.hpp
        \~english \brief fundamental interface between text files and c++

        This class contains a fundamental interface between text files
        and c++ (delete comments, delete whitespaces , read in columns).
        Additionally this class can parse and split strings

        \~german \brief grundlegende Schnittstelle zwischen Textdateien und c++

        Diese Klasse beinhaltet eine grundlegende Schnittstelle zwischen Textdateien
        und c++ (kommentare löschen, Whitespaces löschen, Spalten einlesen).
        Zudem kann die Klasse strings parsen und splitten
       */
    
            /*! \~english \brief erases comments from a string

                deletes everything after '#' including '#' from \a str
                @param str input string
                @return manipulated string, from which all comments were deleted

                \~german \brief löscht Kommentare aus String

                löscht alles was nach '#' kommt inclusive '#' aus \a str
                @param str Eingabe-String
                @return manipulierten String, bei dem die Kommentare gelöscht wurden
            */
            std::string deleteComment(std::string str);

            /*! \~english \brief erases whitespaces from a string

                erases all whitespaces from \a str
                @param str input string
                @return manipulated string, from which all whitespaces were deleted

                \~german \brief löscht Whitespaces aus String

                löscht alle Whitespaces aus \a str
                @param str Eingabe-String
                @return manipulierten String, bei dem die Kommentare gelöscht wurden
            */
            std::string deleteWhitespaces(std::string str);

            /*! \~english \brief reads in a column from a text file

                reads in column \a column_number from the text file \a inputname
                columns are separated by tabs or whitespaces
                @param inputname input file
                @param column_number column number of the data
                @return vector of the column length filled with values from the file

                \~german \brief liest Spalte aus Textdatei

                liest Spalte \a column_number aus der Textdatei \a inputname
                Spalten sind durch tabs oder Leerzeichen getrennt
                @param inputname Eingabedatei
                @param column_number Spaltennummer der Daten
                @return Vektor der Länge der Spalte mit Werten aus Datei
            */
            std::vector<double> readColumn(std::string inputname, uint column_number);

            /*! \~english \brief converts a string into a bool

                accepts "false", "off", "0", "true", "on", "1"
                exit code 999, if no valid value was used
                @param s input string
                @return parsed boolean

                \~german \brief convertiert einen String zu einem Bool

                akzeptiert "flase", "off", "0", "true", "on", "1"
                exit code 999, wenn kein gültiger Wert verwendet wurde
                @param s Eingabe-String
                @return geparsten boolean
            */
            bool stringToBool(std::string s);

            /*! \~english \brief converts a string into a double

                @param s input string
                @return parsed double

                \~german \brief convertiert einen String zu einem Double

                @param s Eingabe-String
                @return geparsten double
            */
            double stringToDouble(std::string s);

            /*! \~english \brief converts a string into an int

                @param s input string
                @return parsed Integer

                \~german \brief convertiert einen String zu einem Int

                @param s Eingabe-String
                @return geparsten Integer
            */
            int stringToInt(std::string s);

            /*! \~english \brief divides a string into a vector

                divides the components of a string, separated by \a delim and saves them in a vector
                @param s input string
                @param delim sign after which string is splitted
                @return vector filled with components of the string

                \~german \brief teilt einen String in einen Vektor

                teilt die Bestandteile eines Strings, die durch \a delim getrennt sind und speichert sie in einem Vektor
                @param s Eingabe-String
                @param delim Zeichen, nach dem gesplittet wird
                @return Vektor mit Teilstrings
            */
            std::vector<std::string> split(std::string s, char delim);
        }
}

#endif /* BASIC_INPUTPARSER_H */
