/*
 * =====================================================================================
 *
 *       Filename:  hbfind.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  16/09/21 08:05:44 PM IST
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  PARTHAJIT ROY (PR), roy.parthajit@gmail.com
 *   Organization:  The University of Burdwan
 *
 * =====================================================================================
 */

#include "hbfind.h"

#include "bioio.h"



static Point3d find_location(Vector3d axis, Vector3d ref_normal, Point3d reference_point, double bond_len, double bond_angle, double torsion_angle){
      Vector3d bang_add = vec3d_polar_rotation(ref_normal, axis, PI - bond_angle);
      Vector3d tors_add = vec3d_polar_rotation(axis, bang_add, torsion_angle);
      Vector3d blen_add = vec3d_scal_mult(tors_add, bond_len);
      Point3d  point_add= (Point3d) vec3d_add(reference_point, blen_add);
      return point_add;
}

static Point3d fix_atom(Point3d pcur1, Point3d pcur2, Point3d pcur3, double bond_len, double bond_angle, double torsion_angle){
      Plane pcur_plane = plane_create(pcur1, pcur2, pcur3);
      Vector3d pcur_unit_vector = vec3d_unit(vec3d_sub(pcur3, pcur2));
      Point3d hloc = find_location(pcur_unit_vector, pcur_plane.unit_normal, pcur3, bond_len, bond_angle, torsion_angle);
      return hloc;
}

static void update_atom(struct atom* atom, int atomid, char* loc_name, char* symb, Point3d* point){
      atom->id = atomid;
      strcpy(atom->loc, loc_name);
      strcpy(atom->symbol, symb);
      atom->center = *point;
}


void cys_addh(struct residue* res)
{
      Point3d N ;
      Point3d CA;
      Point3d CB;
      Point3d SG;
      int count = 0;
      for(int i=0; i<res->size; ++i){
	    print_pdb_line(stdout, res->atoms+i);
	    if(strcmp(res->atoms[i].loc, "N") == 0){
		  N = res->atoms[i].center;
		  count++;
	    }else if(strcmp(res->atoms[i].loc, "CA") == 0){
		  CA = res->atoms[i].center;
		  count++;
	    }else if(strcmp(res->atoms[i].loc, "CB") == 0){
		  CB = res->atoms[i].center;
		  count++;
	    }else if(strcmp(res->atoms[i].loc, "SG") == 0){
		  SG = res->atoms[i].center;
		  count++;
	    }
      }
      if( count != 4 ){    /* Exception Handling */ 
	    fprintf(stderr, "Error in function %s()... All required atoms not found %d\n", __func__, count);
	    exit(EXIT_FAILURE);
      }
      double tors_ang = torsion_angle(N, CA, CB, SG);

      double adjust = torad(120.0);
      double hbangle = torad(106.97);
      double hbdist  = 0.97;
      Point3d HB1 = fix_atom(N, CA, CB, hbdist, hbangle, tors_ang + adjust );
      Point3d HB2 = fix_atom(N, CA, CB, hbdist, hbangle, tors_ang - adjust );
      Point3d HG  = fix_atom(CA, CB, SG, 1.20, torad(109.0), torad(180.0));
      Point3d HA  = fix_atom(SG, CB, CA, 0.97, torad(105.0), torad(174.0));
      struct atom newatom = res->atoms[0];
      update_atom(&newatom, 0, "HB1", "H", &HB1);
      printf("--from here--\n");
      print_pdb_line(stdout, &newatom);
      update_atom(&newatom, 0, "HB2", "H", &HB2);
      print_pdb_line(stdout, &newatom);
      update_atom(&newatom, 0, "HG", "H", &HG);
      print_pdb_line(stdout, &newatom);
      update_atom(&newatom, 0, "HA", "H", &HA);
      print_pdb_line(stdout, &newatom);
}

