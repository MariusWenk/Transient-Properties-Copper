/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   read_file.cpp
 * Author: ndione
 * 
 * Created on 13. Juni 2018, 20:34
 */

#include "read_file.h"

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>

using namespace std;

struct param_line{
	string param;
	double value;
};

class File_reader{
public:
	File_reader(string filename); //Constructor, reads the file "filename" into "data"
	~File_reader(); //Destructor
	double read(string param); //returns the value of the parameter "param", if the parameter is not found -1 is returned
private:
	vector<param_line> data; //stores all the data from the file
	string filename;

};



