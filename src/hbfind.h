/*
 * =====================================================================================
 *
 *       Filename:  hbfind.h
 *
 *    Description: hydrogen addition. 
 *
 *        Version:  1.0
 *        Created:  13/09/21 09:46:15 PM IST
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  PARTHAJIT ROY (PR), roy.parthajit@gmail.com
 *   Organization:  The University of Burdwan
 *
 * =====================================================================================
 */

 
#ifndef  __hbfind_H__
#define  __hbfind_H__



#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "biodefs.h"
#include "polymer.h"

 
void cys_addh(struct residue* res);

#endif   /* ----- #ifndef __hbfind_H__  ----- */
