/*
 * =====================================================================================
 *
 *       Filename:  hbutil.h
 *
 *    Description:  hb unil 
 *
 *        Version:  1.0
 *        Created:  19/09/21 08:37:15 PM IST
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  PARTHAJIT ROY (PR), roy.parthajit@gmail.com
 *   Organization:  The University of Burdwan
 *
 * =====================================================================================
 */

#ifndef  __hbutil_H__
#define  __hbutil_H__
#include "biodefs.h"
#include "polymer.h"
/* These functions are for kd-tree*/
double res_distsqr(struct residue* residue1, struct residue* residue2);
int res_comp(struct residue* residue1, struct residue* residue2, int index);
double res_value_at(struct residue* residue, int index);


#endif   /* ----- #ifndef __hbutil_H__  ----- */
