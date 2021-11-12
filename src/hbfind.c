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

//static void update_atom(struct atom* atom, int atomid, char* loc_name, char* symb, Point3d* point){
//      atom->id = atomid;
//      strcpy(atom->loc, loc_name);
//      strcpy(atom->symbol, symb);
//      atom->center = *point;
//}



void cys_addh(struct residue* res)
{
      printf("Trace... cys exec\n");
      Point3d N ;
      Point3d CA;
      Point3d CB;
      Point3d SG;
      int count = 0;
      for(int i=0; i<res->size; ++i){
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
      residue_addh(res, HB1, "HB1");
      Point3d HB2 = fix_atom(N, CA, CB, hbdist, hbangle, tors_ang - adjust );
      residue_addh(res, HB2, "HB2");
      Point3d HG  = fix_atom(CA, CB, SG, 1.20, torad(90.0), torad(180.0));

      residue_addh(res, HG, "HG");
      Point3d HA  = fix_atom(SG, CB, CA, 0.97, torad(90.0), tors_ang+adjust);
      residue_addh(res, HA, "HA");
}
void asp_addh(struct residue* res)
{
      printf("Trace... asp exec\n");
      Point3d N ;
      Point3d CA;
      Point3d CB;
      Point3d CG;
      int count = 0;
      for(int i=0; i<res->size; ++i){
	    if(strcmp(res->atoms[i].loc, "N") == 0){
		  N = res->atoms[i].center;
		  count++;
	    }else if(strcmp(res->atoms[i].loc, "CA") == 0){
		  CA = res->atoms[i].center;
		  count++;
	    }else if(strcmp(res->atoms[i].loc, "CB") == 0){
		  CB = res->atoms[i].center;
		  count++;
	    }else if(strcmp(res->atoms[i].loc, "CG") == 0){
		  CG = res->atoms[i].center;
		  count++;
	    }
      }
      if( count != 4 ){    /* Exception Handling */ 
	    fprintf(stderr, "Error in function %s()... All required atoms not found %d\n", __func__, count);
	    exit(EXIT_FAILURE);
      }
      double tors_ang = torsion_angle(N, CA, CB, CG);

      double adjust = torad(120.0);
      double hbangle = torad(106.97);
      double hbdist  = 0.97;
      Point3d HB1 = fix_atom(N, CA, CB, hbdist, hbangle, tors_ang + adjust );
      residue_addh(res, HB1, "HB1");
      Point3d HB2 = fix_atom(N, CA, CB, hbdist, hbangle, tors_ang - adjust );
      residue_addh(res, HB2, "HB2");

}

void glu_addh(struct residue* res)
{
      printf("Trace... asp exec\n");
      Point3d CA;
      Point3d CB;
      Point3d CG;
      Point3d CD;
      Point3d  N;
      Point3d  O;
      Point3d  C;
      int count = 0;
      for(int i=0; i<res->size; ++i){
	    if(strcmp(res->atoms[i].loc, "CA") == 0){
		  CA= res->atoms[i].center;
		  count++;
	    }else if(strcmp(res->atoms[i].loc, "CB") == 0){
		  CB = res->atoms[i].center;
		  count++;
	    }else if(strcmp(res->atoms[i].loc, "CG") == 0){
		  CG = res->atoms[i].center;
		  count++;
	    }else if(strcmp(res->atoms[i].loc, "CD") == 0){
		  CD= res->atoms[i].center;
		  count++;
	    }
      }
      if( count != 4 ){    /* Exception Handling */ 
	    fprintf(stderr, "Error in function %s()... All required atoms not found %d\n", __func__, count);
	    exit(EXIT_FAILURE);
      }
      double tors_ang = torsion_angle(N, CA, CB, CG);

      double adjust = torad(120.0);
      double hbangle = torad(106.97);
      double hbdist  = 0.97;
      Point3d HB1 = fix_atom(CD,CG,CB ,hbdist, hbangle, tors_ang + adjust );
      residue_addh(res, HB1, "HB1");
      Point3d HB2 = fix_atom(CD, CG, CB, hbdist, hbangle, tors_ang - adjust );
      residue_addh(res, HB2, "HB2");
      double tors_ang1 = torsion_angle( CA, CB, CG,CD);
      Point3d HE1 = fix_atom(CA,CB,CG ,hbdist, hbangle, tors_ang1 + adjust );
      residue_addh(res, HE1, "HE1");
      Point3d HE2 = fix_atom(CA, CB, CG, hbdist, hbangle, tors_ang1 - adjust );
      residue_addh(res, HE2, "HE2");
      double tors_ang2= torsion_angle(O,C,CA,CB);
      Point3d HA = fix_atom(O,C,CA,hbdist, hbangle, tors_ang2);
      residue_addh(res, HA, "HA");
}


void phe_addh(struct residue* res)
{
      printf("Trace... asp exec\n");
      Point3d C;
      Point3d O;
      Point3d N;
      Point3d CA;
      Point3d CB;
      Point3d CG;
      Point3d CD1;
      Point3d CD2;
      Point3d CE1;
      Point3d CE2;
      Point3d CZ;
      int count = 0;
      for(int i=0; i<res->size; ++i){
          if(strcmp(res->atoms[i].loc, "C") == 0){
		  C= res->atoms[i].center;
		  count++;
	    } 
          else if(strcmp(res->atoms[i].loc, "O") == 0){
		  O= res->atoms[i].center;
		  count++;
	    } 
          else if(strcmp(res->atoms[i].loc, "N") == 0){
		  N= res->atoms[i].center;
		  count++;
	    }
	    else if(strcmp(res->atoms[i].loc, "CA") == 0){
		  CA= res->atoms[i].center;
		  count++;
	    }else if(strcmp(res->atoms[i].loc, "CB") == 0){
		  CB = res->atoms[i].center;
		  count++;
	    }else if(strcmp(res->atoms[i].loc, "CG") == 0){
		  CG = res->atoms[i].center;
		  count++;
	    }else if(strcmp(res->atoms[i].loc, "CD1") == 0){
		  CD1= res->atoms[i].center;
		  count++;
	    }
          else if(strcmp(res->atoms[i].loc, "CD1") == 0){
		  CD2= res->atoms[i].center;
		  count++;
	    }
          else if(strcmp(res->atoms[i].loc, "CE1") == 0){
		  CE1= res->atoms[i].center;
		  count++;
	    }
          else if(strcmp(res->atoms[i].loc, "CE2") == 0){
		  CE2= res->atoms[i].center;
		  count++;
	    }
          else if(strcmp(res->atoms[i].loc, "CZ") == 0){
		  CZ= res->atoms[i].center;
		  count++;
	    }
      }
      if( count != 4 ){    /* Exception Handling */ 
	    fprintf(stderr, "Error in function %s()... All required atoms not found %d\n", __func__, count);
	    exit(EXIT_FAILURE);
      }
      double tors_ang = torsion_angle(N, CA, CB, CG);

      double adjust = torad(120.0);
      double hbangle = torad(106.97);
      double hbdist  = 0.97;
      Point3d HB1 = fix_atom(N,CA,CB,hbdist, hbangle, tors_ang + adjust );
      residue_addh(res, HB1, "HB1");
      Point3d HB2 = fix_atom(N, CA, CB, hbdist, hbangle, tors_ang - adjust );
      residue_addh(res, HB2, "HB2");

      double tors_ang1 = torsion_angle(CB, CG, CD1,CE1);
      Point3d HD1 = fix_atom(CB,CG,CD1 ,hbdist, hbangle, tors_ang1 + adjust );
      residue_addh(res, HD1, "HD1");

      double tors_ang2 = torsion_angle( CG, CD1,CE1,CZ);
      Point3d HE1 = fix_atom(CG, CD1, CE1, hbdist, hbangle, tors_ang2 - adjust );
      residue_addh(res, HE1, "HE1");

      double tors_ang3= torsion_angle( CD1, CE1, CZ,CE2);
      Point3d HZ = fix_atom(CD1,CE1,CZ ,hbdist, hbangle, tors_ang3 + adjust );
      residue_addh(res, HZ, "HZ");

      double tors_ang4= torsion_angle(  CE1, CZ,CE2,CD2);
      Point3d HE2= fix_atom(CE1,CZ,CE2 ,hbdist, hbangle, tors_ang4 + adjust );
      residue_addh(res, HE2, "HE2");

      double tors_ang5= torsion_angle(O,C,CA,CB);
      Point3d HA = fix_atom(O,C,CA ,hbdist, hbangle, tors_ang5 + adjust );
      residue_addh(res, HA, "HA");
}



