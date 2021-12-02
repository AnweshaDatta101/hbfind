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

static Point3d find_location(Vector3d axis, Vector3d ref_normal, Point3d reference_point, double bond_len, double bond_angle, double torsion_angle)
{
      Vector3d bang_add = vec3d_polar_rotation(ref_normal, axis, PI - bond_angle);
      Vector3d tors_add = vec3d_polar_rotation(axis, bang_add, torsion_angle);
      Vector3d blen_add = vec3d_scal_mult(tors_add, bond_len);
      Point3d point_add = (Point3d)vec3d_add(reference_point, blen_add);
      return point_add;
}

static Point3d fix_atom(Point3d pcur1, Point3d pcur2, Point3d pcur3, double bond_len, double bond_angle, double torsion_angle)
{
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

void cys_addh(struct residue *res)
{
      printf("Trace... cys exec\n");
      Point3d N;
      Point3d CA;
      Point3d CB;
      Point3d SG;
      int count = 0;
      for (int i = 0; i < res->size; ++i)
      {
            if (strcmp(res->atoms[i].loc, "N") == 0)
            {
                  N = res->atoms[i].center;
                  count++;
            }
            else if (strcmp(res->atoms[i].loc, "CA") == 0)
            {
                  CA = res->atoms[i].center;
                  count++;
            }
            else if (strcmp(res->atoms[i].loc, "CB") == 0)
            {
                  CB = res->atoms[i].center;
                  count++;
            }
            else if (strcmp(res->atoms[i].loc, "SG") == 0)
            {
                  SG = res->atoms[i].center;
                  count++;
            }
      }
      if (count != 4)
      { /* Exception Handling */
           fprintf(stderr, "Error in function %s()... All required atoms not found %d\n", __func__, count);
            exit(EXIT_FAILURE);
      }
      double tors_ang = torsion_angle(N, CA, CB, SG);

      double adjust = torad(120.0);
      double hbangle = torad(106.97);
      double hbdist = 0.97;
      Point3d HB1 = fix_atom(N, CA, CB, hbdist, hbangle, tors_ang + adjust);
      residue_addh(res, HB1, "HB1");
      Point3d HB2 = fix_atom(N, CA, CB, hbdist, hbangle, tors_ang - adjust);
      residue_addh(res, HB2, "HB2");
      Point3d HG = fix_atom(CA, CB, SG, 1.20, torad(90.0), torad(180.0));

      residue_addh(res, HG, "HG");
      Point3d HA = fix_atom(SG, CB, CA, 0.97, torad(90.0), tors_ang + adjust);
      residue_addh(res, HA, "HA");
}
void asp_addh(struct residue *res)
{
      printf("Trace... asp exec\n");
      Point3d N;
      Point3d CA;
      Point3d CB;
      Point3d CG;
      int count = 0;
      for (int i = 0; i < res->size; ++i)
      {
            if (strcmp(res->atoms[i].loc, "N") == 0)
            {
                  N = res->atoms[i].center;
                  count++;
            }
            else if (strcmp(res->atoms[i].loc, "CA") == 0)
            {
                  CA = res->atoms[i].center;
                  count++;
            }
            else if (strcmp(res->atoms[i].loc, "CB") == 0)
            {
                  CB = res->atoms[i].center;
                  count++;
            }
            else if (strcmp(res->atoms[i].loc, "CG") == 0)
            {
                  CG = res->atoms[i].center;
                  count++;
            }
      }
      if (count != 4)
      { /* Exception Handling */
            fprintf(stderr, "Error in function %s()... All required atoms not found %d\n", __func__, count);
            exit(EXIT_FAILURE);
      }
      double tors_ang = torsion_angle(N, CA, CB, CG);

      double adjust = torad(120.0);
      double hbangle = torad(106.97);
      double hbdist = 0.97;
      Point3d HB1 = fix_atom(N, CA, CB, hbdist, hbangle, tors_ang + adjust);
      residue_addh(res, HB1, "HB1");
      Point3d HB2 = fix_atom(N, CA, CB, hbdist, hbangle, tors_ang - adjust);
      residue_addh(res, HB2, "HB2");

      double tors_ang1 = torsion_angle(N, CA, CB, CG);
      Point3d HA = fix_atom(CG, CB, CA, hbdist, hbangle, tors_ang1 + adjust);
      residue_addh(res, HA, "HA");
      Point3d HN = fix_atom(CB, CA, N, 1.20, torad(111.37), torad(180.0));

      residue_addh(res, HN, "HN");
}

void glu_addh(struct residue *res)
{
      printf("Trace... asp exec\n");
      Point3d CA;
      Point3d CB;
      Point3d CG;
      Point3d CD;
      Point3d N;
     // Point3d O;
     Point3d C;
      int count = 0;
      for (int i = 0; i < res->size; ++i)
      {
            if (strcmp(res->atoms[i].loc, "CA") == 0)
            {
                  CA = res->atoms[i].center;
                  count++;
            }
            else if (strcmp(res->atoms[i].loc, "CB") == 0)
            {
                  CB = res->atoms[i].center;
                  count++;
            }
            else if (strcmp(res->atoms[i].loc, "CG") == 0)
            {
                  CG = res->atoms[i].center;
                  count++;
            }
            else if (strcmp(res->atoms[i].loc, "CD") == 0)
            {
                  CD = res->atoms[i].center;
                  count++;
            }
             else if (strcmp(res->atoms[i].loc, "N") == 0)
            {
                  N= res->atoms[i].center;
                  count++;
            }
             else if (strcmp(res->atoms[i].loc, "C") == 0)
            {
                  C= res->atoms[i].center;
                  count++;
            }
      }

      if (count != 6)
      { /* Exception Handling */
            fprintf(stderr, "Error in function %s()... All required atoms not found %d\n", __func__, count);
            exit(EXIT_FAILURE);
      }
      double tors_ang = torsion_angle(CA, CB, CG, CD);

      double adjust = torad(120.0);
      double hbangle = torad(106.97);
      double hbdist = 0.97;
      Point3d HB1 = fix_atom(CD, CG, CB, hbdist, hbangle, tors_ang + adjust);
      residue_addh(res, HB1, "HB1");
      Point3d HB2 = fix_atom(CD, CG, CB, hbdist, hbangle, tors_ang - adjust);
      residue_addh(res, HB2, "HB2");
      // double tors_ang1 = torsion_angle( CA, CB, CG,CD);
      Point3d HE1 = fix_atom(CA, CB, CG, hbdist, hbangle, tors_ang + adjust);
      residue_addh(res, HE1, "HE1");
      Point3d HE2 = fix_atom(CA, CB, CG, hbdist, hbangle, tors_ang - adjust);
      residue_addh(res, HE2, "HE2");
       tors_ang = torsion_angle(C, CA, CB, CG);
      Point3d HA = fix_atom(CG, CB, CA, hbdist, torad(280), tors_ang + adjust);
      residue_addh(res, HA, "HA");
      Point3d H= fix_atom(CB, CA, N, .86, torad(183.33), torad(90));

      residue_addh(res, H, "H");
}

void phe_addh(struct residue* res)
{
      printf("Trace... asp exec\n");
      Point3d C;
    //  Point3d O;
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
        /*  else if(strcmp(res->atoms[i].loc, "O") == 0){
		  O= res->atoms[i].center;
		  count++;
	    } */
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
          else if(strcmp(res->atoms[i].loc, "CD2") == 0){
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
if( count != 10){    /* Exception Handling */
          fprintf(stderr, "Error in function %s()... All required atoms not found %d\n", __func__, count);
	    exit(EXIT_FAILURE);
      }
      double tors_ang = torsion_angle(C ,CA, CB, CG);

      double adjust = torad(120.0);
      double hbangle = torad(106.97);
      double hbdist  = 0.97;
      Point3d HB2 = fix_atom(C,CA,CB,hbdist, hbangle, tors_ang + adjust );
      residue_addh(res, HB2, "HB2");
      Point3d HB3 = fix_atom(C, CA, CB, hbdist, hbangle, tors_ang - adjust );
      residue_addh(res, HB3, "HB3");
      Point3d HA = fix_atom(CG,CB,CA,hbdist,torad(150.0),tors_ang -adjust);
      residue_addh(res, HA, "HA");

      tors_ang = torsion_angle(CG,CD2,CE2,CZ);
      Point3d HD2 = fix_atom(CZ,CE2,CD2,hbdist, torad(130.0), tors_ang+torad(145.0) );
      residue_addh(res, HD2, "HD2");
      Point3d HE2 = fix_atom(CG,CD2,CE2, hbdist, hbangle, tors_ang + torad(210.0) );
      residue_addh(res, HE2, "HE2");

      tors_ang= torsion_angle( CE2, CZ, CE1,CD1);
      Point3d HZ = fix_atom(CD1,CE1,CZ ,hbdist, torad(120.0), tors_ang + torad(130.0) );
      residue_addh(res, HZ, "HZ");
      Point3d HE1= fix_atom(CE2,CZ,CE1 ,hbdist, hbangle, tors_ang-torad(150.0));
      residue_addh(res, HE1, "HE1");

      tors_ang= torsion_angle(CZ,CE1,CD1,CG);
      Point3d HD1= fix_atom(CZ,CE1 ,CD1,hbdist, torad(235.0), tors_ang  );
      residue_addh(res, HD1, "HD1");

       Point3d H= fix_atom(CB, CA, N, .86, torad(180.00), torad(90));
       residue_addh(res, H, "H");
      
}
     
void gly_addh(struct residue *res)
{
      printf("Trace... gly exec\n");
      Point3d O;
      Point3d C;
      Point3d CA;
      Point3d N;
      int count = 0;
      for (int i = 0; i < res->size; ++i)
      {
            if (strcmp(res->atoms[i].loc, "O") == 0)
            {
                  O = res->atoms[i].center;
                  count++;
            }
            else if (strcmp(res->atoms[i].loc, "C") == 0)
            {
                  C = res->atoms[i].center;
                  count++;
            }
            else if (strcmp(res->atoms[i].loc, "CA") == 0)
            {
                  CA = res->atoms[i].center;
                  count++;
            }
            else if (strcmp(res->atoms[i].loc, "N") == 0)
            {
                  N = res->atoms[i].center;
                  count++;
            }
      }
      if (count != 4)
      { /* Exception Handling */
            fprintf(stderr, "Error in function %s()... All required atoms not found %d\n", __func__, count);
            exit(EXIT_FAILURE);
      }
      double tors_ang = torsion_angle(O, C, CA, N);

      double adjust = torad(120.0);
      double hbangle = torad(106.97);
      double hbdist = 0.97;
      Point3d HB1 = fix_atom(O, C, CA, hbdist, hbangle, tors_ang + adjust);
      residue_addh(res, HB1, "HB1");
      Point3d HB2 = fix_atom(O, C, CA, hbdist, hbangle, tors_ang - adjust);
      residue_addh(res, HB2, "HB2");
      Point3d HG = fix_atom(C, CA, N, 1.20, torad(109.07), torad(180.0));

      residue_addh(res, HG, "HG");
}
void his_addh(struct residue* res)
{
      printf("Trace... his exec\n");
      Point3d O;
      Point3d C;
      Point3d CA;
      Point3d N;
      Point3d CB;
      Point3d CG;
      Point3d ND1;
      Point3d CE1;
      Point3d NE2;
      Point3d CD2;

      int count = 0;
      for(int i=0; i<res->size; ++i){
	    if(strcmp(res->atoms[i].loc, "O") == 0){
		  O = res->atoms[i].center;
		  count++;
	    }else if(strcmp(res->atoms[i].loc, "C") == 0){
		  C = res->atoms[i].center;
		  count++;
	    }else if(strcmp(res->atoms[i].loc, "CA") == 0){
		  CA = res->atoms[i].center;
		  count++;
	    }else if(strcmp(res->atoms[i].loc, "N") == 0){
		  N = res->atoms[i].center;
		  count++;
	    }
          else if(strcmp(res->atoms[i].loc, "CB") == 0){
		  CB = res->atoms[i].center;
		  count++;
	    }
          else if(strcmp(res->atoms[i].loc, "CG") == 0){
		  CG= res->atoms[i].center;
		  count++;
	    }
          else if(strcmp(res->atoms[i].loc, "ND1") == 0){
		  ND1 = res->atoms[i].center;
		  count++;
	    }
          else if(strcmp(res->atoms[i].loc, "CE1") == 0){
		  CE1= res->atoms[i].center;
		  count++;
	    }
          else if(strcmp(res->atoms[i].loc, "NE2") == 0){
		  NE2 = res->atoms[i].center;
		  count++;
	    }
          else if(strcmp(res->atoms[i].loc, "CD2") == 0){
		  CD2 = res->atoms[i].center;
		  count++;
	    }
          
      }
      if( count != 10 ){ /* Exception Handling */
          fprintf(stderr, "Error in function %s()... All required atoms not found %d\n", __func__, count);
	    exit(EXIT_FAILURE);
      }
      double tors_ang = torsion_angle(O,C,CA,N);

      double adjust = torad(120.0);
      double hbangle = torad(106.97);
      double hbdist  = 0.97;
      Point3d HA= fix_atom(O,C ,CA, hbdist, hbangle, tors_ang-adjust);
      residue_addh(res, HA, "HA");
      
      Point3d H  = fix_atom(C,CA,N ,1.20, torad(115.14), torad(180.0)+adjust);

      residue_addh(res, H, "H");

      double tors_ang1 = torsion_angle(C,CA,CB,CG);
      Point3d HB1 = fix_atom(C ,CA,CB, hbdist, hbangle, tors_ang1 + adjust );
      residue_addh(res, HB1, "HB1");
      Point3d HB2 = fix_atom(C, CA,CB, hbdist, hbangle, tors_ang1 - adjust );
      residue_addh(res, HB2, "HB2");

      double tors_ang2 = torsion_angle(CG,ND1,CE1,NE2);
      Point3d HE1= fix_atom(CG,ND1,CE1, hbdist, hbangle, tors_ang2+torad(180.0) );
      residue_addh(res, HE1, "HE1");
      
     double tors_ang3 = torsion_angle(CB,CG,CD2,NE2);
      Point3d HD2= fix_atom(CB,CG,CD2, hbdist, hbangle, tors_ang3+torad(180.0));
      residue_addh(res, HD2, "HD2");
}
/*void ile_addh(struct residue *res)
{
      printf("Trace... cys exec\n");
      Point3d N;
      Point3d CA;
      Point3d CB;
      Point3d CG1;
      Point3d CG2;
      Point3d CD1;

      int count = 0;
      for (int i = 0; i < res->size; ++i)
      {
            if (strcmp(res->atoms[i].loc, "N") == 0)
            {
                  N = res->atoms[i].center;
                  count++;
            }
            else if (strcmp(res->atoms[i].loc, "CA") == 0)
            {
                  CA = res->atoms[i].center;
                  count++;
            }
            else if (strcmp(res->atoms[i].loc, "CB") == 0)
            {
                  CB = res->atoms[i].center;
                  count++;
            }
            else if (strcmp(res->atoms[i].loc, "CG1") == 0)
            {
                  CG1 = res->atoms[i].center;
                  count++;
            }
            else if (strcmp(res->atoms[i].loc, "CG2") == 0)
            {
                  CG2 = res->atoms[i].center;
                  count++;
            }
            else if (strcmp(res->atoms[i].loc, "CD1") == 0)
            {
                  CD1 = res->atoms[i].center;
                  count++;
            }
      }
      if (count != 6)
      { *//* Exception Handling */
        /*   fprintf(stderr, "Error in function %s()... All required atoms not found %d\n", __func__, count);
            exit(EXIT_FAILURE);
      }
      double tors_ang = torsion_angle(CG2,CB,CA,N);

      double adjust = torad(120.0);
      double hbangle = torad(106.97);
      double hbdist = 0.97;
      Point3d HA = fix_atom(CG2, CB, CA, 0.97, torad(90.0), tors_ang + adjust);
       residue_addh(res, HA, "HA");
      Point3d HB = fix_atom(N, CA, CB, hbdist, hbangle, tors_ang + adjust);
       residue_addh(res, HB, "HB");
      Point3d HG22 = fix_atom(CA, CB, CG2, 1.20, torad(90.0), torad(180.0)+adjust);
       residue_addh(res, HG22, "HG22");
      Point3d HG21 = fix_atom(CA, CB, CG2, 1.20, torad(90.0), torad(180.0));
       residue_addh(res, HG21, "HG21");
      Point3d HG23= fix_atom(CA, CB, CG2, 1.20, torad(90.0), torad(180.0)-adjust);
       residue_addh(res, HG23, "HG23");
       Point3d H = fix_atom(CB, CA, N, 1.20, torad(90.0), torad(180.0));
       residue_addh(res, H, "H");
       tors_ang = torsion_angle(CA,CB,CG2,CD1);
      Point3d H12 = fix_atom(CA,CB,CG1, hbdist, hbangle, tors_ang + adjust);
      residue_addh(res, H12, "H12");
      Point3d H13 = fix_atom(CA,CB,CG1, hbdist, hbangle, tors_ang - adjust);
      residue_addh(res, H13, "H13");
       Point3d HD1 = fix_atom(CB, CG1, CD1, 1.20, torad(90.0), torad(162.95)+adjust);
       residue_addh(res, HD1, "HD1");
      Point3d HD2 = fix_atom( CB, CG1,CD1, 1.20, torad(90.0), torad(162.95));
       residue_addh(res, HD2, "HD2");
      Point3d HD3= fix_atom( CB, CG1,CD1, 1.20, torad(90.0), torad(162.95)-adjust);
       residue_addh(res, HD3, "HD3");
    
      
}*/
/*void asn_addh(struct residue *res)
{
      printf("Trace... asp exec\n");
      Point3d N;
      Point3d CA;
      Point3d CB;
      Point3d CG;
      Point3d C;
      Point3d ND2;

      int count = 0;
      for (int i = 0; i < res->size; ++i)
      {
            if (strcmp(res->atoms[i].loc, "N") == 0)
            {
                  N = res->atoms[i].center;
                  count++;
            }
            else if (strcmp(res->atoms[i].loc, "CA") == 0)
            {
                  CA = res->atoms[i].center;
                  count++;
            }
            else if (strcmp(res->atoms[i].loc, "CB") == 0)
            {
                  CB = res->atoms[i].center;
                  count++;
            }
            else if (strcmp(res->atoms[i].loc, "CG") == 0)
            {
                  CG = res->atoms[i].center;
                  count++;
            }
            else if (strcmp(res->atoms[i].loc, "C") == 0)
            {
                  C = res->atoms[i].center;
                  count++;
            }
            else if (strcmp(res->atoms[i].loc, "ND2") == 0)
            {
                  ND2= res->atoms[i].center;
                  count++;
            }
      }
      if (count != 6)
      {*/ /* Exception Handling */
           /* fprintf(stderr, "Error in function %s()... All required atoms not found %d\n", __func__, count);
            exit(EXIT_FAILURE);
      }
      double tors_ang = torsion_angle(C, CA, CB, CG);

      double adjust = torad(120.0);
      double hbangle = torad(106.97);
      double hbdist = 0.97;
      Point3d HB1 = fix_atom(C, CA, CB, hbdist, hbangle, tors_ang + adjust);
      residue_addh(res, HB1, "HB1");
      Point3d HB2 = fix_atom(C, CA, CB, hbdist, hbangle, tors_ang - adjust);
      residue_addh(res, HB2, "HB2");
      Point3d HA = fix_atom(CG, CB, CA, hbdist, hbangle, tors_ang + adjust);
      residue_addh(res, HA, "HA");
      Point3d H = fix_atom(CB, CA, N, 1.20, torad(111.37), torad(180.0)+adjust);
      residue_addh(res, H, "H");
      Point3d HN1 = fix_atom(CB, CG, ND2, 1.20, torad(111.37), torad(180.0)+adjust);

      residue_addh(res, HN1, "HN1");
      Point3d HN2 = fix_atom(CB, CG, ND2, 1.20, torad(111.37), torad(180.0)-adjust);

      residue_addh(res, HN2, "HN2");

}*/
/*void lys_addh(struct residue *res)
{
      printf("Trace... cys exec\n");
      Point3d CA;
      Point3d CB;
      Point3d CG;
      Point3d CD;
      Point3d CE;
      Point3d NZ;
      int count = 0;
      for (int i = 0; i < res->size; ++i)
      {
            if (strcmp(res->atoms[i].loc, "CA") == 0)
            {
                  CA = res->atoms[i].center;
                  count++;
            }
            else if (strcmp(res->atoms[i].loc, "CB") == 0)
            {
                  CB = res->atoms[i].center;
                  count++;
            }
            else if (strcmp(res->atoms[i].loc, "CG") == 0)
            {
                  CG = res->atoms[i].center;
                  count++;
            }
            else if (strcmp(res->atoms[i].loc, "CD") == 0)
            {
                  CD = res->atoms[i].center;
                  count++;
            }
            else if (strcmp(res->atoms[i].loc, "CE") == 0)
            {
                  CE = res->atoms[i].center;
                  count++;
            }
            else if (strcmp(res->atoms[i].loc, "NZ") == 0)
            {
                  NZ = res->atoms[i].center;
                  count++;
            }
      }
      if (count != 6)
      {*/ /* Exception Handling */
         /*  fprintf(stderr, "Error in function %s()... All required atoms not found %d\n", __func__, count);
            exit(EXIT_FAILURE);
      }
      double tors_ang = torsion_angle(CA, CB, CG,CD);

      double adjust = torad(120.0);
      double hbangle = torad(106.97);
      double hbdist = 0.97;
      Point3d HG1 = fix_atom(CA, CB,CG, hbdist, hbangle, tors_ang + adjust);
      residue_addh(res, HG1, "HG1");
      Point3d HG2 = fix_atom(CA, CB,CG, hbdist, hbangle, tors_ang - adjust);
      residue_addh(res, HG2, "HG2");
      Point3d HB1 = fix_atom(CD, CG,CB, hbdist, hbangle, tors_ang + adjust);
      residue_addh(res, HB1, "HB1");
      Point3d HB2 = fix_atom(CD, CG,CB, hbdist, hbangle, tors_ang - adjust);
      residue_addh(res, HB2, "HB2");

      tors_ang = torsion_angle(CG,CD,CE,NZ);
      Point3d HD1 = fix_atom(NZ,CE,CD, hbdist, hbangle, tors_ang + adjust);
      residue_addh(res, HD1, "HD1");
      Point3d HD2 = fix_atom(NZ,CE,CD, hbdist, hbangle, tors_ang - adjust);
      residue_addh(res, HD2, "HD2");
      Point3d HE1 = fix_atom(CG,CD,CE, hbdist, hbangle, tors_ang + adjust);
      residue_addh(res, HE1, "HE1");
      Point3d HE2 = fix_atom(CG,CD,CE, hbdist, hbangle, tors_ang - adjust);
      residue_addh(res, HE2, "HE2");
      Point3d NZ1 = fix_atom(CD, CE, NZ, 1.20, torad(90.0), torad(180.0)+adjust);
      residue_addh(res, NZ1, "NZ1");
      Point3d NZ2 = fix_atom(CD, CE, NZ, 1.20, torad(90.0), torad(180.0)+torad(90.0));
      residue_addh(res, NZ2, "NZ2");
      Point3d NZ3 = fix_atom(CD, CE, NZ, 1.20, torad(90.0), torad(180.0)-adjust);
      residue_addh(res, NZ3, "NZ3");
}*/