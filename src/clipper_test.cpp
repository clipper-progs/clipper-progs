// Clipper app to run self-tests
/* Copyright 2003-2004 Kevin Cowtan & University of York all rights reserved */

#include <clipper/clipper.h>
#include <clipper/core/test_core.h>
#include <clipper/contrib/test_contrib.h>

#include <fstream>
extern "C" {
#include <stdlib.h>
}


int main( int argc, char** argv )
{
  // do self tests
  std::cout << "Test core:\n";
  clipper::Test_core test_core;
  bool result_core = test_core();
  if ( result_core ) std::cout << "OK\n";
  std::cout << "Test contrib:\n";
  clipper::Test_contrib test_contrib;
  bool result_contrib = test_contrib();
  if ( result_contrib ) std::cout << "OK\n";
  // done self tests

  // error exit code
  if ( !( result_core && result_contrib ) ) exit(1);
}
