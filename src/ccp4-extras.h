/*! \file ccp4-extras.h CCP4 C++ program extras */
/* Copyright 2003-2004 Kevin Cowtan & University of York all rights reserved */

//L   This code is distributed under the terms and conditions of the
//L   CCP4 Program Suite Licence Agreement as a CCP4 Application.
//L   A copy of the CCP4 licence can be obtained by writing to the
//L   CCP4 Secretary, Daresbury Laboratory, Warrington WA4 4AD, UK.

#include <vector>
#include <iostream>
#include <string>


//! Mini-parser for command line input
/*! This class removes any escape characters from the command-line,
  and if a -stdin option is present reads standard input and adds the
  result to the argument list. */
class CommandInput : public std::vector<std::string>
{
 public:
  CommandInput( int argc, char** argv, bool echo = false );
};


//! class for program start and end
class CCP4program
{
 public:
  CCP4program( const char* name, const char* vers, const char* rcsdate );
  ~CCP4program();
};
