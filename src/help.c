#include <iostream>
#include <iomanip>
#include <string>
#include <stdlib.h>
#include <fstream>
#include <cstring>
#include <new>

using namespace std;

//file help.c
/*
+ help - help file interface
    Description :
     This  uses  the  command  string   in  topic  to  print  out  the
     appropriate  part of the  muvi.help file.   If topic  = "TOPICS",
     then all commands  are listed.  Otherwise it prints  help for the
     command given by key.

*/
void  help(char topic[])
{
  int ic;
  bool found;
  char delim, temp, path[50], name[60];
  string keyw, line;
  std::ifstream file;
  strcpy(path, "mg");
  if (getenv(path) == NULL){
      cout << path << " not defined" << '\n';
      return;
  }
  strcpy(name, getenv(path));
  strcat(name, "/src/help/muvi.help");
  file.open(name);
  delim = '?';
  keyw = topic;
  ic = 0;
  if (keyw == "TOPICS"){
    while (!file.eof()){
      temp = file.get();
      if (temp == delim){
	file >> line;
	cout << right << setw(15) << line << "    ";
	ic++;
	if (ic > 3){
	  std::cout << '\n';
	  ic = 0;
	}
	getline(file, line, delim);
	file.get();
      }
    }
  cout << '\n';
  }
  else{
    found = false;
    while (!file.eof()){
      temp = file.get();
      if (temp == delim){
	file >> line;
	if (line == keyw){
	  getline(file, line, delim);
	  cout << line;
	  found = true;
	}
      }
    }
    if (!found)
      cout << "Unknown topic " << keyw << '\n';
  }
  file.close();
  return;
}
//end help.c

 
