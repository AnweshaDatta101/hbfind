/*
 * =====================================================================================
 *
 *       Filename:  kdtree.c
 *
 *    Description:  KD-Tree 
 *
 *        Version:  1.0
 *        Created:  19/09/21 07:51:31 PM IST
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  PARTHAJIT ROY (PR), roy.parthajit@gmail.com
 *   Organization:  The University of Burdwan
 *
 * =====================================================================================
 */


#include "kdtree.h"
#define KDTREE_DIM (3)

struct kdt_node{
      struct residue*  key;
      struct kdt_node* left;
      struct kdt_node* right;
};



void kdtree_init(struct kdtree* tree, struct residue* objs, int size)
{
      tree->nodes = (struct kdt_node*) malloc ( size * sizeof(struct kdt_node) );
      if ( tree->nodes==NULL ) {
	    fprintf ( stderr, "\ndynamic memory allocation failed in function %s()\n" , __func__);
	    exit (EXIT_FAILURE);
      }

      tree->size = size;
      for(int i=0; i<size; ++i){
	    tree->nodes[i].key = objs + i;
	    tree->nodes[i].left = NULL;
	    tree->nodes[i].right= NULL;
      }
}


void kdtree_free(struct kdtree* tree){
      free(tree->nodes);
      tree->nodes = NULL;
}

static void swap(struct kdt_node* n1, struct kdt_node* n2)
{
      struct residue* tmp = NULL;// = malloc(size * sizeof(char));

      tmp = n1->key;
      n1->key = n2->key;
      n2->key = tmp;
}


static int partition(struct kdt_node* nodes, int beg, int end, int idx)
{
      if(end < beg) return -1;
      if( end == beg ) return beg;
      int mid = (end - beg) / 2;
      swap(nodes+mid, nodes+end);
      struct kdt_node* pivot = nodes+end;
      int pos = beg - 1;
      for(int i=beg; i<end; ++i){
	    if(res_comp(nodes[i].key, pivot->key, idx) <= 0){
		  pos ++;
		  swap(nodes+pos, nodes+i);
	    }
      }
      swap(nodes+pos+1, nodes+end);
      return pos + 1;
}

static struct kdt_node* median(struct kdt_node* nodes, int beg, int end, int i, int idx)
{
      if(beg > end) return NULL;
      if(beg == end) return nodes+beg;
      int pos = partition(nodes, beg, end, idx);
      int k = pos - beg + 1;
      if(i == k){
	    return nodes + pos;
      }else if(i < k){
	    return median(nodes, beg, pos-1, i, idx);
      }else{
	    return median(nodes, pos+1, end, i - k, idx);
      }
}

static struct kdt_node* kdt_make(struct kdt_node* base, int len, int i)
{
      int middle = (0 + len) / 2;
      struct kdt_node* med = median(base, 0, len-1, middle, i);
      if(med != NULL){
	    i = (i + 1) % KDTREE_DIM; 
	    med->left  = kdt_make(base, med - base, i);
	    med->right = kdt_make(med + 1, (base + len) - (med + 1), i);
      }
      return med;
}


void kdtree_build(struct kdtree* tree)
{
      tree->root = kdt_make(tree->nodes, tree->size, 0);
}
static void kdt_node_neighbor(struct kdt_node* node, struct residue* key, int i, double range2, struct residue* neighbors[], int* count, int maxnn)
{
      if(node == NULL) return;
      
      double dist2 = res_distsqr(node->key, key);
      double delta = res_value_at(node->key, i) - res_value_at(key, i);
      double delta2 = delta * delta;

      if(dist2 <= range2){
	    if( *count >= maxnn ){    /* Exception Handling */ 
		  fprintf(stderr, "Error in function %s()... too many neighbors found.\n", __func__);
		  exit(EXIT_FAILURE);
	    }
	    neighbors[*count] = node->key;
	    *count = (*count) + 1;
      }
      i = (i + 1) % KDTREE_DIM;
      kdt_node_neighbor(delta > 0.0 ? node->left : node->right, key, i, range2,  neighbors, count, maxnn);
      if(delta2 >= range2) return;
      kdt_node_neighbor(delta > 0.0 ? node->right : node->left, key, i, range2,  neighbors, count, maxnn);
}
void kdtree_neighbors(struct residue* neighbors[], int* count, int maxnn, struct kdtree* tree, struct residue* obj, double range2)
{
      *count = 0;
      kdt_node_neighbor(tree->root, obj, 0, range2, neighbors, count, maxnn);
}
