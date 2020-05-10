/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Ofek Brazani (Azrieli college of engineering, ofek1b@gmail.com)
   
   Foursets code addition:
   Michal Gabay (Azrieli college of engineering, michalg552@gmail.com)
   Shira Yerushalmi (Azrieli college of engineering, shirushalmi@gmail.com)
------------------------------------------------------------------------- */

#include <cstdlib>
#include <cstring>
#include "fix_ave_atom.h"
#include "fix_reaxc_checkFourset.h"
#include "atom.h"
#include "update.h"
#include "pair_reaxc.h"
#include "modify.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "comm.h"
#include "force.h"
#include "compute.h"
#include "input.h"
#include "variable.h"
#include "memory.h"
#include "error.h"
#include "reaxc_list.h"
#include "reaxc_types.h"
#include "reaxc_defs.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixReaxCCheckFourset::FixReaxCCheckFourset(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  //check if the command in the "in" file is legal
  if (narg != 7) error->all(FLERR,"Illegal fix reax/c/checkFourset command");
  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);
  ntypes = atom->ntypes;
  nmax = atom->nmax;

  nevery_cond_check=10; //DEFAULT VALUE
  nevery = force->inumeric(FLERR,arg[3]);

  if (nevery <= 0 )
    error->all(FLERR,"Illegal fix reax/c/checkFourset command, illigal nevery");//**********
  nevery_cond_check=nevery;
  printf("\n\nnevery_cond_check %d\n\n",nevery_cond_check);

//for dists fp that follow the distance between atoms
  if (me == 0) {
      char *suffix = strrchr(arg[4],'.');
      if (suffix && strcmp(suffix,".gz") == 0) {
  #ifdef LAMMPS_GZIP
        char gzip[128];
        snprintf(gzip,128,"gzip -6 > %s",arg[4]);
  #ifdef _WIN32
        fp = _popen(gzip,"wb");
  #else
        fp = popen(gzip,"w");
  #endif
  #else
        error->one(FLERR,"Cannot open gzipped file");
  #endif
      } else fp = fopen(arg[4],"w");

      if (fp == NULL) {
        char str[128];
        snprintf(str,128,"Cannot open fix reax/c/checkFourset file %s",arg[4]);
        error->one(FLERR,str);
      }
    }

//follow distances and document them every this many steps
   nevery_dists_follow = force->inumeric(FLERR,arg[5]);
  if (nevery_dists_follow <= 0 )
    error->all(FLERR,"Illegal fix reax/c/checkFourset command, illigal follow dists nevery");//**********


    //seperate the dists files into many files any this number of timesteps
   nevery_file_dists = force->inumeric(FLERR,arg[6]);
  if (nevery_file_dists <= 0 )
    error->all(FLERR,"Illegal fix reax/c/checkFourset command, illigal dists file nevery");//**********


  if (atom->tag_consecutive() == 0)
    error->all(FLERR,"Atom IDs must be consecutive for fix reax/c/checkFourset");

  if(nevery_cond_check > nevery_dists_follow || nevery_dists_follow > nevery_file_dists || nevery_cond_check > nevery_file_dists)
    error->all(FLERR,"Illegal fix reax/c/checkFourset command, nevery_cond_check < nevery_dists_follow < nevery_file_dists");//**********

  fourset = NULL;
  o_c_pair_tags = NULL;
  n_tags=NULL;
  MAX_NUM_FOURSETS=atom->nlocal;
  timeout_timesteps_at_start_and_end=10000; //default value
  
  // Michal & Shira
  follow_selected_atoms = NULL;
  ATOMS_ARRAY_SIZE = atom->nlocal + 1;
  ////// for checking
  printf("\n -------- ATOMS_ARRAY_SIZE = %d\n", ATOMS_ARRAY_SIZE);

  allocate(); //allocate all the memory for each struct.
  strcpy(fp_suffix,arg[4]); //define the suffix for each dists file that follow dists between reactive atoms

  int _set_flag=set_mol_pattern();//set o_c_pair_tags, n_tags by the user input
  if(_set_flag!=0) error->all(FLERR,"Illegal \"Extra_Potential_Parameters\" file, illigal molecole file pattern");//**********
  
}

/* ---------------------------------------------------------------------- */

FixReaxCCheckFourset::~FixReaxCCheckFourset()
{
  MPI_Comm_rank(world,&me);
  destroy();
  if (me == 0) fclose(fp);
}

/* ---------------------------------------------------------------------- */

int FixReaxCCheckFourset::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixReaxCCheckFourset::setup(int /*vflag*/)
{
  end_of_step();
}

/* ---------------------------------------------------------------------- */

void FixReaxCCheckFourset::init()
{
  //reaxc singelton instance
  reaxc = (PairReaxC *) force->pair_match("reax/c",0);
  if (reaxc == NULL) error->all(FLERR,"Cannot use fix reax/c/checkFourset without "
                                "pair_style reax/c, reax/c/kk, or reax/c/omp");
}

/* ---------------------------------------------------------------------- */

void FixReaxCCheckFourset::end_of_step()
{
  Output_ReaxC_Bonds(update->ntimestep);
  if (me == 0) fflush(fp);
}

/* ---------------------------------------------------------------------- */

void FixReaxCCheckFourset::Output_ReaxC_Bonds(bigint /*ntimestep*/)
{
  FindNbr(lists);
}

/* ---------------------------------------------------------------------- */

/* this function update the two neigh lists of all the atoms:
    far_neigh_list- for non-bonded atoms with distance less then 10
    bond_list- for bonded atoms.
  using those lists to find legal foursets that stands with the distances
  conditions from the paper with only O-C pairs that the C atom bonded only
  to one O atom.*/
void FixReaxCCheckFourset::FindNbr(struct _reax_list * /*lists*/)
{
  
  int nlocal_tot = static_cast<int> (atom->nlocal);
  if (atom->nmax > nmax) {
    nmax = atom->nmax;
  }
  int num_of_epons=int(atom->nlocal/58.5);

  //do it only each nevery timestep
  if(update->ntimestep%nevery_cond_check!=0)
    return;

  //document distances every this many steps only
  int doc_dists_flag=0;
  if(update->ntimestep%nevery_dists_follow==0)  doc_dists_flag=1;

  if(fp!=NULL){
    //writing the header of the dists file
    if(update->ntimestep==0){
      fprintf(fp,"# totalTimesteps " BIGINT_FORMAT " \n",update->laststep);
      fprintf(fp,"# totalAtomNum %d \n",atom->nlocal);
      fprintf(fp,"# nevery_dists_follow %d\n",nevery_dists_follow);//document dists every this many steps
      fprintf(fp,"# nevery_file_dists %d",nevery_file_dists);//create new dist file every this many steps
    }
    //create new dists file to write to with the same suffix.
    if(update->ntimestep%nevery_file_dists==0 && update->ntimestep>0){
      fclose(fp);
      char fp_name[100];
      strcpy(fp_name,fp_suffix);
      strcat(fp_name,".");
      char ts[20];
      sprintf(ts, "%d", update->ntimestep);
      strcat(fp_name,ts);
      fp = fopen(fp_name,"w");
    }
    //write to the dists file the current time step as header to this timestep distances list
    if(fp!=NULL && doc_dists_flag==1)
      fprintf(fp,"\n# Timestep " BIGINT_FORMAT ,update->ntimestep);
  }

 //const known values.
  const int TYPE_C = 0;
  const int TYPE_H = 1;
  const int TYPE_O = 2;
  const int TYPE_N = 3;
  const int EPON_SIZE = 43;
  const int DETDA_SIZE = 31;
  
  
  //i1= O atom, i2= H atom, i3= N atom, i4= C atom
  int pi1, pi2, pi3, pi4; //the index value of each atom in the list
  int _o, _h, _n, _c; //the index of each atom in the main all atoms list
  int start_o, end_o, start_h, end_h, start_n, end_n, start_c, end_c, start_12, end_12; //pointers to the start and the end of each atom neigh list.
  int type_i1, type_i2, type_i3, type_i4, tag_i1, tag_i2; //tag&type value of each atom
  reax_list *far_nbrs, *bond_nbrs; //pointer to neigh lists
  far_neighbor_data *nbr_p_oh, *nbr_p_nc; //pointer to neigh in far neigh list
  bond_data *nbr_p_hn, *nbr_p_co, *nbr_p_12; //pointer to neigh in bond neigh list
  reax_atom *atom_i1, *atom_i2, *atom_i3, *atom_i4, *atom_i5; //pointer to atom struct

  far_nbrs = (reaxc->lists) + FAR_NBRS;
  bond_nbrs = (reaxc->lists) + BONDS;
  num_fourset = 0;

  //reset the neigh lists
  for(int nn=0; nn<MAX_NUM_FOURSETS; nn++){
    fourset[nn][0]=0;
    fourset[nn][1]=0;
    fourset[nn][2]=0;
    fourset[nn][3]=0;
  }

  for( _o = 0; _o < reaxc->system->N; ++_o ) {
    atom_i1 = &(reaxc->system->my_atoms[_o]);
    type_i1  = atom_i1->type;
    tag_i1 = atom_i1->orig_id;
    if(tag_i1<0) continue;

    if(doc_dists_flag==1){
      /* writing to the dists file distance only between O,C,N atoms known as legal candidate
        for foursets atoms to the rest of all atoms from their own far-neigh-list and bond-list*/
      if(fp!=NULL){
        int writing_flag=0;
        for(int wf=0; wf<4; wf++){
          if( (tag_i1-o_c_pair_tags[wf][0])%EPON_SIZE==0 || (tag_i1-o_c_pair_tags[wf][1])%EPON_SIZE==0 )
            writing_flag=1;
        }

        if( (tag_i1-(num_of_epons*EPON_SIZE)-n_tags[0])%DETDA_SIZE==0 || (tag_i1-(num_of_epons*EPON_SIZE)-n_tags[1])%DETDA_SIZE==0 )
          writing_flag=1;

        if(writing_flag==1)
          fprintf(fp,"\n# atom %d type %d ",tag_i1, type_i1+1);
        else continue;
      }
           
    //follow this atom neigh lists to write his distances from the rest of the atoms.
      start_12 = Start_Index(_o, bond_nbrs);
      end_12   = End_Index(_o, bond_nbrs);

      //follow and write bonds distances.
      for( pi1 = start_12; pi1 < end_12; ++pi1 ) {
        nbr_p_12 = &( bond_nbrs->select.bond_list[pi1] );
        _h = nbr_p_12->nbr;
        atom_i2 = &(reaxc->system->my_atoms[_h]);
        type_i2= atom_i2->type;
        tag_i2 = atom_i2->orig_id;

        if(fp!=NULL && tag_i2>0){
          int writing_flag=0;
          for(int wf=0; wf<4; wf++){
            if( (tag_i1-o_c_pair_tags[wf][0])%EPON_SIZE==0 || (tag_i1-o_c_pair_tags[wf][1])%EPON_SIZE==0 )
              writing_flag=1;
          }
            if( (tag_i1-(num_of_epons*EPON_SIZE)-n_tags[0])%DETDA_SIZE==0 || (tag_i1-(num_of_epons*EPON_SIZE)-n_tags[1])%DETDA_SIZE==0 )
              writing_flag=1;
            if(writing_flag==1)
              fprintf(fp,"%d %f ",tag_i2, nbr_p_12->d);
        }
      }
    }
      start_o = Start_Index(_o, far_nbrs);
      end_o   = End_Index(_o, far_nbrs);
    
    //follow and write far neigh distances.
    for( pi1 = start_o; pi1 < end_o; ++pi1 ) {
      nbr_p_oh = &( far_nbrs->select.far_nbr_list[pi1] );
      _h = nbr_p_oh->nbr;
      atom_i2 = &(reaxc->system->my_atoms[_h]);
      type_i2= atom_i2->type;
      tag_i2 = atom_i2->orig_id;

      if(fp!=NULL && tag_i2>0 && doc_dists_flag==1){
        int writing_flag=0;
        for(int wf=0; wf<4; wf++){
          if( (tag_i1-o_c_pair_tags[wf][0])%EPON_SIZE==0 || (tag_i1-o_c_pair_tags[wf][1])%EPON_SIZE==0 )
            writing_flag=1;
        }
          if( (tag_i1-(num_of_epons*EPON_SIZE)-n_tags[0])%DETDA_SIZE==0 || (tag_i1-(num_of_epons*EPON_SIZE)-n_tags[1])%DETDA_SIZE==0 )
            writing_flag=1;
          if(writing_flag==1)
            fprintf(fp,"%d %f ",tag_i2, nbr_p_oh->d);
      }
      
      /* if type_i1==O and type_i2==H look for fourset.
        else, continue write the distances file*/
      if(type_i1 != TYPE_O || type_i2 != TYPE_H) continue;
      //foursets list is full
      if(num_fourset>=MAX_NUM_FOURSETS) continue;
      //make sure to have free first&last timesteps from running the extra potential on fourset to let the system stable
      if(update->ntimestep<timeout_timesteps_at_start_and_end || (update->laststep - update->ntimestep) < timeout_timesteps_at_start_and_end) continue;

      //if O-H distance meets the paper condition, look for N atom
      if (1.3 <= nbr_p_oh->d && nbr_p_oh->d <= 8.0 ){
        start_h = Start_Index(_h, bond_nbrs);
        end_h = End_Index(_h, bond_nbrs);
        for( pi2 = start_h; pi2 < end_h; ++pi2 ){
          nbr_p_hn = &( bond_nbrs->select.bond_list[pi2] );
          _n=nbr_p_hn->nbr;
          atom_i3 = &(reaxc->system->my_atoms[_n]);
          type_i3= atom_i3->type;
          if(type_i3 != TYPE_N) continue;
          //if H-N distance meets the paper condition, look for C atom
          if (0.8 <= nbr_p_hn->d && nbr_p_hn->d <= 1.3 ){
            start_n = Start_Index(_n, far_nbrs);
            end_n = End_Index(_n, far_nbrs);
            for( pi3 = start_n; pi3 < end_n; ++pi3 ){
              nbr_p_nc = &( far_nbrs->select.far_nbr_list[pi3] );
              _c=nbr_p_nc->nbr;
              atom_i4 = &(reaxc->system->my_atoms[_c]);
              type_i4= atom_i4->type;
              if(type_i4 != TYPE_C) continue;
              //if N-C distance meets the paper condition, look for the founded O atom in C's neigh lists
              if(3.0 <= nbr_p_nc->d && nbr_p_nc->d <= 8.0 ){
                int _re=atom_i1->orig_id%EPON_SIZE;
                int _x=int( (atom_i1->orig_id - _re) / EPON_SIZE );
                int _optional_c_tag=0;
                //check if the O-C pait is legal (C is bonded to only this O atom)
                for(int oc=0; oc<4; oc++){
                  if(_re==o_c_pair_tags[oc][0])
                    _optional_c_tag=o_c_pair_tags[oc][1]+EPON_SIZE*_x;
                }
              
                if(atom_i4->orig_id != _optional_c_tag) continue;
                start_c = Start_Index(_c, bond_nbrs);
                end_c = End_Index(_c, bond_nbrs);
                //check if C and O are bonded
                for( pi4 = start_c; pi4 < end_c; ++pi4 ){
                  nbr_p_co = &( bond_nbrs->select.bond_list[pi4] );
                  atom_i5=&(reaxc->system->my_atoms[nbr_p_co->nbr]);
                  if(nbr_p_co->nbr != _o) continue;
                  //if C-O meets the paper condition, we found legal fourset.
                  if (0.9 <= nbr_p_co->d && nbr_p_co->d <= 2.2 ){
                    if(num_fourset<MAX_NUM_FOURSETS){
                      //add the fourset to the foursets list
                      fourset[num_fourset][0] = atom_i1->orig_id;
                      fourset[num_fourset][1] = atom_i2->orig_id;
                      fourset[num_fourset][2] = atom_i3->orig_id;
                      fourset[num_fourset][3] = atom_i4->orig_id;
                      num_fourset++;
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  //code that prints the foursets list any 100 timesteps.
  /*if(update->ntimestep%100==0){
    printf("\n\nat timestep %d\n",update->ntimestep);
    for(int nn=0; nn<num_fourset; nn++)
      printf("fourset #%d: %d %d %d %d\n",nn, fourset[nn][0], fourset[nn][1], fourset[nn][2], fourset[nn][3]);
    printf("\n");
  }*/

    // Michal & Shira
	if(num_fourset!=0){
		
		// Check which quartets meet the conditions and delete those that do not
		// Initialize the follow_selected_atoms array
		for(int i=0; i<ATOMS_ARRAY_SIZE; i++)
		{
			follow_selected_atoms[i] = false;
		}
		printf("\n in fix_reaxc: init array\n");
		
		
		// If an atom is already in another quartet, delete the current quartet in which it appears
		for(int i=0; i<num_fourset; i++)
		{
			int atom_id1 = fourset[i][0];
			int atom_id2 = fourset[i][1];
			int atom_id3 = fourset[i][2];
			int atom_id4 = fourset[i][3];
			int exist_atom = -1;
			if(follow_selected_atoms[atom_id1] == true){
				exist_atom = atom_id1;
			}
			else if(follow_selected_atoms[atom_id2] == true){
				exist_atom = atom_id2;
			}
			else if(follow_selected_atoms[atom_id3] == true){
				exist_atom = atom_id3;
			}
			else if(follow_selected_atoms[atom_id4] == true){
				exist_atom = atom_id4;
			}
			if(exist_atom != -1){
				// At this point, delete means to put zeros
				fourset[i][0] = 0;
				fourset[i][1] = 0;
				fourset[i][2] = 0;
				fourset[i][3] = 0;
				printf("\n delete fourset number %d, because atom id %d\n", i+1, exist_atom);
				continue;
			}
			else{	// update follow_selected_atoms with curr atoms
				follow_selected_atoms[atom_id1] = true;
				follow_selected_atoms[atom_id2] = true;
				follow_selected_atoms[atom_id3] = true;
				follow_selected_atoms[atom_id4] = true;
			}
		}
		
		
		// Check that remains quarters different from 0 (with no bugs - should remain at least 1)
		bool valid_quarter = false;
		if(fourset[0][0] != 0 && fourset[0][1] != 0 && fourset[0][2] != 0 && fourset[0][3] != 0){
			valid_quarter = true;
		}
		
		
		if(valid_quarter){
			// Try to apply the extra potential on the chosen fourset
			int apply_flag = reaxc->set_fourset(fourset, num_fourset);
			// if the apply succeeded
			if(apply_flag==1){
				// optional: print all the founded foursets
				//for(int nn=0; nn<num_fourset; nn++)
				//	printf("\n applied fourset #%d: %d %d %d %d\n",nn, fourset[nn][0], fourset[nn][1], fourset[nn][2], fourset[nn][3]);
				printf("\n\n start operate the potential\n");
				// write to the dists file the fourset we found and apply the extra potential on
				fprintf (fp,"\n# fourset O H N C at timestep " BIGINT_FORMAT " : ",update->ntimestep);
				fprintf(fp,"1/1- %d %d %d %d",fourset[0][0], fourset[0][1], fourset[0][2], fourset[0][3]);
			}
		}
		else{
			printf("\n ERROR: invalid zeros quarter \n");
		}
    }
  
}

/* ---------------------------------------------------------------------- */

int FixReaxCCheckFourset::nint(const double &r)
{
  int i = 0;
  if (r>0.0) i = static_cast<int>(r+0.5);
  else if (r<0.0) i = static_cast<int>(r-0.5);
  return i;
}

/* ---------------------------------------------------------------------- */

void FixReaxCCheckFourset::destroy()
{
  memory->destroy(fourset);
  memory->destroy(o_c_pair_tags);
  memory->destroy(n_tags);
  memory->destroy(fp_suffix);
}

/* ---------------------------------------------------------------------- */

void FixReaxCCheckFourset::allocate()
{
  memory->create(fourset,MAX_NUM_FOURSETS,4,"reax/c/checkFourset:fourset");
  memory->create(o_c_pair_tags,4,2,"reax/c/checkFourset:o_c_pair_tags");
  memory->create(n_tags,2,"reax/c/checkFourset:n_tags");
  memory->create(fp_suffix,100,"reax/c/checkFourset:fp_suffix");
  memory->create(follow_selected_atoms,ATOMS_ARRAY_SIZE,"reax/c/checkFourset:follow_selected_atoms");// Michal & Shira
}

/* ---------------------------------------------------------------------- */

double FixReaxCCheckFourset::memory_usage()
{
  double bytes;
  //bytes = 3.0*nmax*sizeof(double);
  bytes += 1.0*MAX_NUM_FOURSETS*4*sizeof(int);//fourset
  bytes += 1.0*4*2*sizeof(int);//o_c_pair_tags
  bytes += 1.0*2*sizeof(int);//n_tags
  bytes += 1.0*100*sizeof(char);//fp_suffix
  bytes += 1.0*ATOMS_ARRAY_SIZE*sizeof(bool);//follow_selected_atoms	// Michal & Shira

  return bytes;
}
/* ---------------------------------------------------------------------- */
//return 0 for success. else, return -1.
/*this function set the  o_c_pair_tags and n_tags parameters by the user input
  file "Extra_Potential_Parameters.txt" */
int FixReaxCCheckFourset::set_mol_pattern(){

//ofek
 FILE* parameters_fp = fopen("Extra_Potential_Parameters.txt","r");
  if (parameters_fp == NULL) {
    char str[128];
    snprintf(str,128,"Cannot open fix reax/c/checkFourset file Extra_Potential_Parameters.txt");
    error->one(FLERR,str);
    return -1;
  }
  char buff[1000];
  fread(buff, 1000, 1, parameters_fp);
  char *token = strtok(buff, "\n");
  int finish_flag=0;
  int temp;
  
  //reset structs
  for(int i=0; i<4; i++){
    o_c_pair_tags[i][0]=o_c_pair_tags[i][1]=0;
  }
  for(int i=0; i<2; i++){
    n_tags[i]=0;
  }

  int rtn_val;

  //set the o-c pair tags parameters
  while(token){
    //set time step to start and to end search for foursets to apply on the extra potential
    if(strcmp(token, "start_and_end_timeout_timesteps")==0){
      token = strtok(NULL, "\n");
      rtn_val=sscanf(token, "%d", &temp);
      if(rtn_val<=0){
        fclose(parameters_fp);
        return -1;
      }
      timeout_timesteps_at_start_and_end=temp;
      finish_flag++;
    }
    if(strcmp(token, "O-C pair tags")==0 || strcmp(token, "o-c pair tags")==0){
      for(int i=0; i<4; i++){
        for(int j=0; j<2; j++){
          if(j==0) token = strtok(NULL, " ");
          else token = strtok(NULL, "\n");
          rtn_val=sscanf(token, "%d", &temp);
          if(rtn_val<=0){
            fclose(parameters_fp);
            return -1;
          }
          o_c_pair_tags[i][j]=temp;
        }
      }
      finish_flag++;  
    }
    //set the N tags parameters
    if(strcmp(token, "n tags")==0 || strcmp(token, "N tags")==0){
      token = strtok(NULL, " ");
      rtn_val=sscanf(token, "%d", &temp);
      if(rtn_val<=0){
        fclose(parameters_fp);
        return -1;
      }
      n_tags[0]=temp;
      token = strtok(NULL, "\n");
      rtn_val= sscanf(token, "%d", &temp);
      if(rtn_val<=0){
        fclose(parameters_fp);
        return -1;
      }
      n_tags[1]=temp;
      finish_flag++;
    }

    //only if both of them defiened succefully
    if(finish_flag==3){
      for(int i=0; i<4; i++){
        printf("\nO-C pair %d %d\n",o_c_pair_tags[i][0],o_c_pair_tags[i][1]);
      }
      printf("\nn_tag %d %d\n",n_tags[0],n_tags[1]);
      fclose(parameters_fp);
      return 0;
    }
    token = strtok(NULL, "\n");
  }
  // if not both of them defiened succefully, ERROR
  if(finish_flag<3){
    fclose(parameters_fp);
    printf("\n not finish\n");
    return -1;
  }
}
