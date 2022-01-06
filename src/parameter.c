/*
 * =====================================================================================
 *
 *       Filename:  parameter.c
 *
 *    Description:  Parameters descriptions.
 *
 *        Version:  1.0
 *        Created:  Friday 31 December 2021 04:36:14  IST
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  PARTHAJIT ROY (PR), roy.parthajit@gmail.com
 *   Organization:  The University of Burdwan
 *
 * =====================================================================================
 */


#include "parameter.h"

void param_init(struct parameter* args)
{
      args->bio.occu = 'B';
      args->bio.HH_dist = 2.00;
      strcpy(args->sys.os, "linux");
}

void process_argv(int argc, char* argv[], struct parameter* args, int file_index[], int* file_count)
{
      for(int i=1; i<argc; ++i){
	    if(strcmp(argv[i], "-occ") == 0){
		  if(strcmp(argv[i+1], "s") == 0){
			args->bio.occu = 'S';
			++i;
		  }
	    }else{
		  file_index[*file_count] = i;
		  *file_count += 1;
	    }
      }
}

