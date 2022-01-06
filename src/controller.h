		/*
 * =====================================================================================
 *
 *       Filename:  controller.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  20/09/21 09:45:30 AM IST
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  PARTHAJIT ROY (PR), roy.parthajit@gmail.com
 *   Organization:  The University of Burdwan
 *
 * =====================================================================================
 */

#ifndef  __controller_H__
#define  __controller_H__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "biodefs.h"
#include "bioio.h"
#include "polymer.h"
#include "kdtree.h"
#include "hbutil.h"
#include "hbfind.h"

#define MAX_NN_RES (40)
struct site{
      struct residue* src;
      char name[5];
      char type;
      struct residue* nnres[MAX_NN_RES];
      int numnn;
      int maxnn;
};

//void site_init(struct site* site);
//void site_setsrc(struct site* site, struct residue* src);

//void site_fill_neighbor(struct site* site, struct kdtree* tree, double incldist);

void sites_init(struct site** sites, int size);
void sites_build(struct site* sites, struct polymer* poly);
void sites_free(struct site** sites);




void exec_hbfind(struct polymer* poly);

void find_hydro_clash(struct site* sites, int size);


#endif   /* ----- #ifndef __controller_H__  ----- */
