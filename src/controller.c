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

void sites_init(struct site** sites, int size)
{
      *sites	= (struct site*) malloc ( size * sizeof(struct site) );
      if ( *sites==NULL ) {
	    fprintf ( stderr, "\ndynamic memory allocation failed in function %s()\n" , __func__);
	    exit (EXIT_FAILURE);
      }
}

void sites_free(struct site** sites)
{
      free ( *sites );
      *sites	= NULL;
}

void sites_build(struct site* sites, struct polymer* poly ){
      struct kdtree kdtree;
      kdtree_init(&kdtree, poly->residues, poly->numres);

      
      kdtree_build(&kdtree);

//      printf("Trace %d\n", poly->numres);
      for(int i=0; i<poly->numres; ++i){
	    site_init(sites + i);
	    site_setsrc(sites + i, poly->residues+i);
//	    printf("------src %d -- of-- %d\n", i, poly->numres);
//	    printf("src = \n");
//	    print_pdb_line(stdout, site.src->atoms);
	    site_fill_neighbor(sites + i, &kdtree, 10.0);

	    
//	    printf("nbh = \n");
//	    for(int j=0; j<site.numnn; ++j){
//		  struct residue* res = site.nnres[j];
//		  print_pdb_line(stdout,res->atoms);
//	    }
//	    printf("end ---------\n");
      }
      kdtree_free(&kdtree);

}

void site_hydro_clash(struct site* site)
{
      struct residue* res1 = site->src;
      struct residue* res2;
      for(int i=0; i<res1->numh; ++i){
	    for(int j=0; j<site->numnn; ++j){

		  res2 = site->nnres[j];
		  if(res1 == res2) continue;
		  for(int k=0; k<res2->numh; ++k){
			double d2 = distsqr(res1->H[i].center, res2->H[k].center) ;
			if( d2 <= 6.00 ){
			      fprintf(stdout, "Clash found %lf\n", d2);
			      print_pdb_line(stdout, res1->H + i);
			      print_pdb_line(stdout, res2->H + k);
			}
		  }
	    }
      }
}
void find_hydro_clash(struct site* sites, int size)
{
      struct residue* res;
      for(int i=0; i<size; ++i){
	    site_hydro_clash(sites + i);
      }
}

void exec_hbfind(struct polymer* poly)
{
     
           
      for(int i=0; i<poly->numres; ++i){
          if(strcmp(poly->residues[i].name, "CYS") == 0){
              cys_addh(poly->residues + i);
          }
          else if(strcmp(poly->residues[i].name, "ALA") == 0){
              ala_addh(poly->residues + i);
          }
          else if(strcmp(poly->residues[i].name, "ASP") == 0){
              asp_addh(poly->residues + i);
          }
          else if(strcmp(poly->residues[i].name, "GLU") == 0){
              glu_addh(poly->residues + i);
          }
          else if(strcmp(poly->residues[i].name, "PHE") == 0){
              phe_addh(poly->residues + i);
          }
          else if(strcmp(poly->residues[i].name, "GLY") == 0){
              gly_addh(poly->residues + i);
          }
         
           else if(strcmp(poly->residues[i].name, "HIS") == 0){
              his_addh(poly->residues + i);
          }
           else if(strcmp(poly->residues[i].name, "ILE") == 0){
              ile_addh(poly->residues + i);
          }
          else if(strcmp(poly->residues[i].name, "ASN") == 0){
              asn_addh(poly->residues + i);
          }
          else if(strcmp(poly->residues[i].name, "LYS") == 0){
              lys_addh(poly->residues + i);
          }
          else if(strcmp(poly->residues[i].name, "LEU") == 0){
              leu_addh(poly->residues + i);
          }

        else if(strcmp(poly->residues[i].name, "MET") == 0){
              met_addh(poly->residues + i);
          }
        else if(strcmp(poly->residues[i].name, "PRO") == 0){
              pro_addh(poly->residues + i);
          }

         else if(strcmp(poly->residues[i].name, "GLN") == 0){
              gln_addh(poly->residues + i);
          }
          else if(strcmp(poly->residues[i].name, "ARG") == 0){
              arg_addh(poly->residues + i);
          }
         else if(strcmp(poly->residues[i].name, "SER") == 0){
              ser_addh(poly->residues + i);
          }

        else if(strcmp(poly->residues[i].name, "THR") == 0){
              thr_addh(poly->residues + i);
          }
          else if(strcmp(poly->residues[i].name, "VAL") == 0){
              val_addh(poly->residues + i);
          }
        else if(strcmp(poly->residues[i].name, "TRP") == 0){
              trp_addh(poly->residues + i);
          }
        else if(strcmp(poly->residues[i].name, "TYR") == 0){
              tyr_addh(poly->residues + i);
          }
        else if(strcmp(poly->residues[i].name, "A") == 0){
              ade_addh(poly->residues + i);
          }
        /*else if(strcmp(poly->residues[i].name, "G") == 0){
              dg_addh(poly->residues + i);
          }*/
          else if(strcmp(poly->residues[i].name, "G") == 0){
              gua_addh(poly->residues + i);
          }
          else if(strcmp(poly->residues[i].name, "C") == 0){
              cyt_addh(poly->residues + i);
          }
          else if(strcmp(poly->residues[i].name, "T") == 0){
              thy_addh(poly->residues + i);
          }
           else if(strcmp(poly->residues[i].name, "U") == 0){
              ura_addh(poly->residues + i);
          }  

 }

      
      struct kdtree kdtree;
      printf("numres = %d\n", poly->numres);
      kdtree_init(&kdtree, poly->residues, poly->numres);

      
      kdtree_build(&kdtree);

      struct site site;
      site_init(&site);
//      printf("Trace %d\n", poly->numres);
      for(int i=0; i<poly->numres; ++i){
	    site_setsrc(&site, poly->residues+i);
//	    printf("------src %d -- of-- %d\n", i, poly->numres);
//	    printf("src = \n");
//	    print_pdb_line(stdout, site.src->atoms);
	    site_fill_neighbor(&site, &kdtree, 10.0);

	    
//	    printf("nbh = \n");
//	    for(int j=0; j<site.numnn; ++j){
//		  struct residue* res = site.nnres[j];
//		  print_pdb_line(stdout,res->atoms);
//	    }
//	    printf("end ---------\n");
      }
      

}

