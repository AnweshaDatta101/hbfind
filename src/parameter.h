/*
 * =====================================================================================
 *
 *       Filename:  parameter.h
 *
 *    Description:  Parameters for the program
 *
 *        Version:  1.0
 *        Created:  Friday 31 December 2021 04:34:47  IST
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  PARTHAJIT ROY (PR), roy.parthajit@gmail.com
 *   Organization:  The University of Burdwan
 *
 * =====================================================================================
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

struct parameter{
      struct{
	    char basename[128];
	    char ext[10];
	    char path[512];
	    char type;
	    char full_name[512];
      }file;
      struct{
	    double hbdist;
	    char occu;
	    double HH_dist;
      }bio;
      struct{
	    char os[10];
      }sys;
};


void param_init(struct parameter* args);


void process_argv(int argc, char* argv[], struct parameter* args, int file_index[], int* file_count);



