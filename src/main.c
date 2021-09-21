/*
 * =====================================================================================
 *
 *       Filename:  main.c
 *
 *    Description:  HBFind 
 *
 *        Version:  1.0
 *        Created:  13/09/21 09:43:14 PM IST
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
#include "geom3d.h"
#include "biodefs.h"
#include "bioio.h"
#include "polymer.h"
#include "hbfind.h"



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
      }bio;
      struct{
	    char os[10];
      }sys;
};


void param_init(struct parameter* args)
{
      args->bio.occu = 'B';
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

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  main
 *  Description:  
 * =====================================================================================
 */

int all_residues(char* res){
      if(res != NULL) return 1;
      return 0;
}
int main ( int argc, char *argv[] )
{
      printf("HBFind Starts\n");
      struct parameter args;
      int file_index[1000];
      int file_count = 0;

      param_init(&args);

      
      process_argv(argc, argv, &args, file_index, &file_count);
      
      
      struct atom* atoms;
      int numatoms;
      for(int i=0; i<file_count; ++i){
	    strcpy(args.file.full_name, argv[file_index[i]]);

	    
//	    fname_split(args.file.path, args.file.basename, args.file.ext, args.file.full_name);
	    strcpy(args.file.ext, ".pdb");
	    if(strcmp(args.file.ext, ".cif") == 0){
		  scancif(args.file.full_name, all_residues, NULL, NULL, &atoms, &numatoms, ALL_TYPE, "label", args.bio.occu);
	    }else if(strcmp(args.file.ext, ".pdb") ==0 ){
		  scanpdb(args.file.full_name, all_residues, NULL, NULL, &atoms, &numatoms, ALL_TYPE, args.bio.occu);
	    }else{
		  fprintf(stderr, "Error in function %s()... Unrecognized file type supplied.\n", __func__);
		  exit(EXIT_FAILURE);
	    }
	    exec_hbfind(atoms, numatoms);
//	    struct polymer polymer;
//	    polymer_create(&polymer, atoms, numatoms);
//	    struct residue res = residue_at(&polymer, 0);
//	    cys_addh(&res);
//	    printpdb("abc.pdb", atoms, numatoms);
      }
      return EXIT_SUCCESS;
}				/* ----------  end of function main  ---------- */
