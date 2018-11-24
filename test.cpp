#include<iostream>
#include<cmath>
#include<vector>
#include<string>
#include"cholesky.h"



using namespace std;

#define ANSI_COLOR_RED  "\x1b[1;31m"
#define ANSI_COLOR_GREEN   "\x1b[1;32m"
#define ANSI_COLOR_RESET   "\x1b[0m"

static int total = 0;
static int passed = 0;


void test__count(bool res, const char *name) {
  printf("%s : ", name);
  total++;
  passed += res;
  if (res){
    printf(ANSI_COLOR_GREEN "PASSED \n" ANSI_COLOR_RESET);
        }
  else {
    printf(ANSI_COLOR_RED "FAILED \n" ANSI_COLOR_RESET);
  }
}


int main(){
  test__count(true,"JUST A TEST");
  test__count(false,"JUST A TEST");
	return 0;
}