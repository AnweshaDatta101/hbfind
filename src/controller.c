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

void exec_hbfind(struct polymer* poly)
{
     
           
      for(int i=0; i<poly->numres; ++i){
          if(strcmp(poly->residues[i].name, "CYS") == 0){
              cys_addh(poly->residues + i);
          }
          else if(strcmp(poly->residues[i].name, "ASP") == 0){
              asp_addh(poly->residues + i);
          }
          else if(strcmp(poly->residues[i].name, "GLU") == 0){
              glu_addh(poly->residues + i);
          }
          /*else if(strcmp(poly->residues[i].name, "PHE") == 0){
              phe_addh(poly->residues + i);
          }*/ 
          else if(strcmp(poly->residues[i].name, "GLY") == 0){
              gly_addh(poly->residues + i);
          }
         
           /*else if(strcmp(poly->residues[i].name, "HIS") == 0){
              his_addh(poly->residues + i);
          }*/
         /* else if(strcmp(poly->residues[i].name, "ILE") == 0){
              ile_addh(poly->residues + i);
          }*/
         /* else if(strcmp(poly->residues[i].name, "ASN") == 0){
              asn_addh(poly->residues + i);
          }*/



      }

 /*     
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
      */

}

