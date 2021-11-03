/*
 * =====================================================================================
 *
 *       Filename:  kdtree.h
 *
 *    Description:  K-D Tree 
 *
 *        Version:  1.0
 *        Created:  19/09/21 09:56:49 AM IST
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  PARTHAJIT ROY (PR), roy.parthajit@gmail.com
 *   Organization:  The University of Burdwan
 *
 * =====================================================================================
 */



#ifndef  __kdtree_H__
#define  __kdtree_H__


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "bioio.h"
#include "polymer.h"
#include "hbutil.h"

struct kdt_node;

struct kdtree{
      struct kdt_node* nodes;
      struct kdt_node* root;
//      int node_size;
      int size;
//      double (*value_at)(void*, int dim);
//      int (*objcmp)(void*, void*, int);
//      double (*pfdist2)(void*, void*);
};

void kdtree_init(struct kdtree* tree, struct residue* objs, int size);
void kdtree_free(struct kdtree* tree);
void kdtree_build(struct kdtree* tree);
void kdtree_neighbors(struct residue* neighbors[], int* count, int maxnn, struct kdtree* tree, struct residue* obj, double range_sqr);
#endif   /* ----- #ifndef __kdtree_H__  ----- */
