/*
 * =====================================================================================
 *
 *       Filename:  controller.c
 *
 *    Description:   
 *
 *        Version:  1.0
 *        Created:  20/09/21 04:32:08 PM IST
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  PARTHAJIT ROY (PR), roy.parthajit@gmail.com
 *   Organization:  The University of Burdwan
 *
 * =====================================================================================
 */

#include "controller.h"


void site_init(struct site* site)
{
      site->src = NULL;
      site->maxnn = MAX_NN_RES;
}
void site_setsrc(struct site* site, struct residue* src)
{
      site->src = src;
      site->numnn = 0;
}

void site_fill_neighbor(struct site* site, struct kdtree* tree, double incldist)
{
      kdtree_neighbors(site->nnres, &site->numnn, site->maxnn, tree, site->src, incldist *incldist);
}

void exec_hbfind(struct atom* atoms, int numatom)
{
      struct polymer poly;
      
      
      polymer_create(&poly, atoms, numatom);

      
      struct kdtree kdtree;
      printf("numres = %d\n", poly.numres);
      kdtree_init(&kdtree, poly.residues, poly.numres);

      
      kdtree_build(&kdtree);

      struct site site;
      site_init(&site);
      printf("Trace %d\n", poly.numres);
      for(int i=0; i<poly.numres; ++i){
	    site_setsrc(&site, poly.residues+i);
	    printf("------src %d -- of-- %d\n", i, poly.numres);
	    printf("src = \n");
	    print_pdb_line(stdout, site.src->atoms);
	    site_fill_neighbor(&site, &kdtree, 1236.5);

	    
	    printf("nbh = \n");
	    for(int j=0; j<site.numnn; ++j){
		  struct residue* res = site.nnres[j];
		  print_pdb_line(stdout,res->atoms);
	    }
	    printf("end ---------\n");
      }

}

