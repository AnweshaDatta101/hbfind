/*
 * =====================================================================================
 *
 *       Filename:  hbutil.c
 *
 *    Description:  Utility functions 
 *
 *        Version:  1.0
 *        Created:  19/09/21 08:17:47 PM IST
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  PARTHAJIT ROY (PR), roy.parthajit@gmail.com
 *   Organization:  The University of Burdwan
 *
 * =====================================================================================
 */


#include "biodefs.h"
#include "polymer.h"
#include "bioio.h"
double res_distsqr(struct residue* residue1, struct residue* residue2)
{
      struct residue* res1 = (struct residue*) residue1;
      struct residue* res2 = (struct residue*) residue2;
      double d2 = distsqr(res1->atoms[0].center, res2->atoms[0].center);
      return d2;

}

int res_comp(struct residue* residue1, struct residue* residue2, int index)
{
      
      
      struct residue* res1 = (struct residue*) residue1;
      struct residue* res2 = (struct residue*) residue2;

      if(index == 0){
	    if(res1->atoms[0].center.x < res2->atoms[0].center.x) return -1;
	    if(res1->atoms[0].center.x > res2->atoms[0].center.x) return +1;
	    return 0;
      }
      if(index == 1){
	    if(res1->atoms[0].center.y < res2->atoms[0].center.y) return -1;
	    if(res1->atoms[0].center.y > res2->atoms[0].center.y) return +1;
	    return 0;
      }
      if(index == 2){
	    if(res1->atoms[0].center.z < res2->atoms[0].center.z) return -1;
	    if(res1->atoms[0].center.z > res2->atoms[0].center.z) return +1;
	    return 0;
      }
      fprintf(stderr, "Error in function %s()... Invalid dimension.\n", __func__);
      exit(EXIT_FAILURE);
      return -999;
}
double res_value_at(struct residue* residue, int index)
{

      struct residue* res = (struct residue*) residue;
      if(index == 0) return res->atoms[0].center.x;
      if(index == 1) return res->atoms[0].center.y;
      if(index == 2) return res->atoms[0].center.z;
      fprintf(stderr, "Error in function %s()... Invalid dimension.\n", __func__);
      exit(EXIT_FAILURE);
      return -999;

}
