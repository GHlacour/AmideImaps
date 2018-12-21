/*
 * $Id: AmideImaps.c,v 1.18 2002/02/28 11:00:27 spoel Exp $
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 4.5.1
 * Copyright (c) 1991-2001, University of Groningen, The Netherlands
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 * 
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 * 
 * For more info, check our website at http://www.gromacs.org
 * 
 * And Hey:
 * GROtesk MACabre and Sinister
 */
static char *SRCID_g_potential_c = "$Id: g_potential.c,v 1.18 2002/02/28 11:00:27 spoel Exp $";
#include <math.h>
#include <ctype.h>
#include "sysstuff.h"
#include "string.h"
#include "typedefs.h"
#include "smalloc.h"
#include "macros.h"
#include "princ.h"
#include "rmpbc.h"
#include "vec.h"
#include "xvgr.h"
#include "pbc.h"
#include "copyrite.h"
#include "futil.h"
#include "statutil.h"
#include "rdgroup.h"
#include "gmx_fatal.h"
#include "AmideImaps_backbone_sidechain.h"
#include "bondf.h"
#include "string2.h"
#include "index.h"
#include "gstat.h"
#include "gmx_ana.h"

#define EPS0 8.85419E-12
#define ELC 1.60219E-19
#define Bohr 0.0529189379 // 1 Bohr=0.0529189379 nm
#define autorcm 219474.63068 // 1 Hartree = 219474.63068 cm-1
#define dtoau 0.393456 // 1 Debye = 0.393456 e Bohr
#define dmu 0.05296 // 0.031*2.84 // The change in dipole by 0->1 excitation in OH 
#define pi 3.1415927 // pi

// Index routine
int Sindex(int a,int b,int N){
  int ind;
  ind=a+N*b;
  return ind;
}

/* This program contruct the time dependent Hamiltonian and Dipoles
   The Hamiltonian and Dipoles are stored to disk in one of the following
   formats:
   Diagonal: only instantaneous eigenfrequencies, but full dipoles
   Local: Hamiltonain in local basis, full dipoles
   Eigen: Hamiltonain in eigen basis (including zeros), full dipoles

   Maps: Jansen Map, Skinner Map, Roy Map, Torri Map, Tasumi Map
*/

void calc_Ham(const char *fnTRX,t_topology *top,char input[256],output_env_t oenv)
{
  // GROMACS type variables
  rvec *x;
  matrix box;
  t_trxstatus *status;
  int nr_frames;
  int       ePBC=-1;
  int N_singles;
  int opls;
  int cyclic;
  int cModelTDC;
  int cModelTCC;
  int nModelTDC;
  int nModelTCC;
  int nModelGLDP;
  int nModelTasumi; 
  int TDCKrimm;
  int TDTorii;
  int nnTreatMap;
  real t;
  t_residue *residue;
  t_snap *snap;
  t_map *map;
  t_map *Promap;
  t_Skinner *SkinnerAmI;
  t_Tokmakoff *TokmakoffAmI;
  t_amide *amideGroups;
  //t_coupling *coupling;
  t_select *select;
  float *Hamiltonian;
  float *dipole;
  int i,numFormat,numMap;
  FILE *hamFile,*dipFile;
  t_pbc pbc;

  // Allocate memory
  snap=(t_snap *)calloc(1,sizeof(t_snap));
  //coupling=(t_coupling *)calloc(1,sizeof(t_coupling));
   // Read input files
  readInput(input,snap);
  map=(t_map *)calloc(1,sizeof(t_map));
  Promap=(t_map *)calloc(1,sizeof(t_map));
  SkinnerAmI=(t_Skinner *)calloc(1,sizeof(t_Skinner));
  TokmakoffAmI=(t_Tokmakoff *)calloc(1,sizeof(t_Tokmakoff));
  select=(t_select *)calloc(1,sizeof(t_select));
    
  /* Maping strings to number for different force fields, peptide types and electrostatic maps */
  opls=0;   
  if(!strcmp(snap->forcefield,"OPLSAA")) {
   opls=1; /* OPLSAA force field parameters used */
  }  
  cyclic=0;
  if(!strcmp(snap->peptide,"Cyclic")){
    cyclic=1; /* peptide type cyclic */
  } 
  numMap=0; 
  if(!strcmp(snap->ElectrostaticMap,"Skinner")) {
   readSkinnerAmI(snap,SkinnerAmI); /* Two Sites (CN) map */
   numMap=0;
  }
  if(!strcmp(snap->ElectrostaticMap,"Jansen")) { 
   readMap(snap,map); /* Four Sites (CON-D) map */
   numMap=1; 
  }
  if(!strcmp(snap->ElectrostaticMap,"Tokmakoff")) {
   readTokmakoffAmI(snap,TokmakoffAmI); /* One Site (O) map */
   numMap=2;
  }
  TDTorii=0;
  if(!strcmp(snap->TD,"Torii")) {
   TDTorii=1;
  }
  //printf("%d %d %d %d \n",opls,cyclic,numMap,TDTorii);

  readProMap(snap,Promap); /* Four Sites (CON-CD) proline map */
  //readCoupling(snap,coupling);
  /* Maping strings to number for different coupling models */
  TDCKrimm=0;
  if(!strcmp(snap->TDCtypeName,"Krimm")){
   TDCKrimm=1;
  }
  cModelTDC=0;
  if(!strcmp(snap->cModel,"TDC")){
   cModelTDC=1;
  }
  cModelTCC=0;
  if(!strcmp(snap->cModel,"TCC")){
   cModelTCC=1;
  }
  nModelTDC=0;
  if(!strcmp(snap->nModel,"TDC")){
   nModelTDC=1;
  }
  nModelTCC=0;
  if(!strcmp(snap->nModel,"TCC")){
   nModelTCC=1;
  }
  nModelTasumi=0; 
  if(!strcmp(snap->nModel,"Tasumi")){
   nModelTasumi=1;
  }
  nModelGLDP=0;
  if(!strcmp(snap->nModel,"GLDP")){
   nModelGLDP=1;
  }
  nnTreatMap=0;
  if (!strcmp(snap->nnTreat,"Map")){
   nnTreatMap=1;
  } 
  //printf("%d %d %d %d %d %d %d %d\n",cModelTDC,cModelTCC,nModelTDC,nModelTCC,TDCKrimm,nModelTasumi,nModelGLDP,nnTreatMap);
  // Read first GMX configuration and number of atoms
  if ((snap->Natoms = read_first_x(oenv,&status,fnTRX,&t,&x,box)) == 0) {
       gmx_fatal(FARGS,"Could not read coordinates from statusfile\n"); 
  }

  /* Identify backbone and sidechain amideI groups*/
  if((snap->NbackboneCO+snap->NsidechainCO)!=snap->Namidebonds) {
   printf("Fatal error!!!\n");
   printf("Backbone and sidechain chromophores do not add up to the total number of chromophores \n");
   exit(-1);
  }
  amideGroups=(t_amide *)calloc((snap->Namidebonds),sizeof(t_amide));
  //modifyTopology(snap,top);
  identifyBackboneCO(snap,top,cyclic,amideGroups);
  if(snap->NsidechainCO>0) {
    if(numMap!=0) {
     printf("Fatal error!!!\n");
     printf("Only Skinner map can be applicable for all units if sidechain chromophores are found. \n");
     exit(-1);
    }
    identifySidechainCO(snap,top,amideGroups);
  }

  // Check charges
  calcCharge(snap,top);

  // Construct Hamiltonian
  N_singles=snap->Namidebonds; 
  
  Hamiltonian=(float *)calloc((N_singles)*(N_singles),sizeof(float));
  dipole=(float *)calloc((N_singles)*3,sizeof(float));
  subsystem(snap,select);
  
  /*********** Start processing trajectory ***********/
   // Write output as a simple text file
  numFormat=0;
  if (!strcmp(snap->format,"TEXT")){
      hamFile=fopen(snap->energyFName,"w");
      dipFile=fopen(snap->dipoleFName,"w");
      numFormat=0;
  }

  // Write output as a binary file
  if (!strcmp(snap->format,"BIN")){
      hamFile=fopen(snap->energyFName,"wb");
      dipFile=fopen(snap->dipoleFName,"wb");
      numFormat=1; 
  }

  nr_frames=0;
  do {
        set_pbc(&pbc,ePBC,box);
        /* Use map to calculate diagonal and dipoles */
        calcDiag(&pbc,snap,map,Promap,SkinnerAmI,TokmakoffAmI,numMap,top,Hamiltonian,dipole,amideGroups,x,opls,cyclic,TDTorii,nnTreatMap,select,nr_frames);
        /* Use map to calculate nearest neighborh coupling 
           Use TDC/TCC to calculate remaining couplings */
        calcCoupling(&pbc,snap,cModelTDC,cModelTCC,nModelTDC,nModelTCC,TDCKrimm,nModelGLDP,nModelTasumi,top,x,Hamiltonian,cyclic,amideGroups);
        /* Save to disk */
        saveHam(snap,Hamiltonian,dipole,hamFile,dipFile,numFormat,nr_frames);
        nr_frames++;
  } while (read_next_x(oenv,status,&t,snap->Natoms,x,box));
  
  free(dipole),free(Hamiltonian);//free(coupling)
  free(snap),free(amideGroups),free(map);
  free(select),free(Promap),free(SkinnerAmI),free(TokmakoffAmI);
  fclose(hamFile),fclose(dipFile);
  /*********** done with status file **********/
  close_trj(status);
  sfree(x);  /* free memory used by coordinate array */
  
}

////////////////////////////////////////////////////////////////////////////////////////
//                                         Subroutines                               //
//////////////////////////////////////////////////////////////////////////////////////

/* Modify topology files */
//void modifyTopology(t_snap *snap,t_topology *top){
  //int i,j,k;
  //int res0;
  //int NatomsChain1;
  //int NatomsPerChain;
  //int NresiduePerChain;

  //NatomsChain1=0;
  //NatomsPerChain=0;
  //NresiduePerChain=0;
  
  //res0=top->atoms.atom[0].resind; 
  //for (i=0;i<snap->Natoms;i++){
      //if ((top->atoms.atom[i].resind)!=res0) {
         //NatomsChain1=i; //Number of atoms of the first residue
         //break;
      //}
  //} 
  //for (i=NatomsChain1;i<snap->Natoms;i++){
      //if ((top->atoms.atom[i].resind)==res0) {
         //NatomsPerChain=i; //Number of atoms of the first chain
         //NresiduePerChain=top->atoms.atom[i].resind+1; 
         //break;
      //}
  //}
  //k=0;
  //for(i=NatomsPerChain;i<snap->Natoms;i++) {
    //if((top->atoms.atom[i].resind)==res0) {
      //k++;
    //}
    //top->atoms.atom[i].resind+=k*NresiduePerChain;
  //}
  //printf("ResNum   AtomNum\n");
  //for(i=0;i<snap->Natoms;i++){
    //printf(" %d        %d\n",top->atoms.atom[i].resind,i);
  //}
//}

/* Select sub-system */
void subsystem(t_snap *snap,t_select *select){
  /* Initializing subsystems */
  select->all=0;
  select->water=0;
  select->protein=0;
  select->lipid=0;
  select->water_and_ion=0;
  select->na=0;
  select->cl=0;
  select->k=0;
  select->no_nnfs=0;

  if (!strcmp(snap->subSelect,"All")) {
    select->all=1;
  } else if (!strcmp(snap->subSelect,"Water")) {
    select->water=1;
  } else if (!strcmp(snap->subSelect,"Protein")){
    select->protein=1;
  } else if (!strcmp(snap->subSelect,"Lipid")) {
    select->lipid=1;
  } else if (!strcmp(snap->subSelect,"WaterAndIon")) {
    select->water_and_ion=1;
  } else if (!strcmp(snap->subSelect,"NA+")) {
    select->na=1;
  } else if (!strcmp(snap->subSelect,"CL-")) {
    select->cl=1;
  } else if (!strcmp(snap->subSelect,"K+")) {
    select->k=1;
  } else if (!strcmp(snap->subSelect,"NoNNFS")) {
    select->no_nnfs=1;
  }
  //printf("%d %d %d %d %d %d %d %d\n",select->all,select->water,select->protein,select->lipid,select->water_and_ion,select->na,select->cl,select->no_nnfs); 
  return;
}

/* Save Hamiltonian and dipole */
void saveHam(t_snap *snap,float *Hamiltonian,float *dipole,FILE *hamFile,FILE *dipFile,int numFormat,int nr_frames){
  int i,j,x;

  // Write as simple text file
  if (numFormat==0){
      fprintf(hamFile,"%d ",nr_frames);
     // Write Hamiltonian  
     for (i=0;i<snap->Namidebonds;i++){
       for(j=i;j<snap->Namidebonds;j++){
         fprintf(hamFile,"%f ",Hamiltonian[i+j*snap->Namidebonds]);
       }
      }
      fprintf(hamFile,"\n");
 
     // Write dipole 
     fprintf(dipFile,"%d ",nr_frames);
     for (x=0;x<3;x++){
      for (i=0;i<snap->Namidebonds;i++){
	fprintf(dipFile,"%f ",dipole[x+i*3]);
      }
     }
     fprintf(dipFile,"\n");
  }

  // Write as binary
  if (numFormat==1){
     // Write Hamiltonian 
     fwrite(&nr_frames,sizeof(int),1,hamFile);
     for (i=0;i<snap->Namidebonds;i++){
       for (j=i;j<snap->Namidebonds;j++){
         fwrite(&Hamiltonian[i+j*snap->Namidebonds],sizeof(float),1,hamFile);
       }
     }

    // Write Dipoles
    fwrite(&nr_frames,sizeof(int),1,dipFile);
    for (x=0;x<3;x++){
      for (i=0;i<snap->Namidebonds;i++){
        fwrite(&dipole[x+i*3],sizeof(float),1,dipFile);
      }
    }
  }
    
  return;
}

/* Identify the backbone amide I groups */
void identifyBackboneCO(t_snap *snap,t_topology *top,int cyclic,t_amide *amideGroups){
  int i,j;
  int t;
  float q;
  int orderN,orderC;
  int control;
  int amide;
  int at[6]={0,0,0,0,0,0};
  bool  proline;
  int resNum;  

  orderN=0,orderC=0;
  printf("\n");  
  amide=0;
  resNum=snap->COchains-1;
    
  // Find amidebonds
  for (j=0;j<snap->Nresidues;j++){ // Possible amidebonds
    proline=FALSE;
    control=0;
    for (i=0;i<snap->Natoms;i++){
      if (top->atoms.atom[i].resind==j+orderC){ // Test if C terminal atoms are present
        if (!strcmp(*top->atoms.resinfo[j].name,"FOR")){                // Residue FOR 
          if (!strcmp(*top->atoms.atomname[i],"C")) control+=8,at[0]=i;
          if (!strcmp(*top->atoms.atomname[i],"O")) control+=16,at[1]=i;
          if (!strcmp(*top->atoms.atomname[i],"H")) control+=32,at[5]=i;
          // printf("Residue%2d is %s %s\n",j,*top->atoms.resinfo[j].name,*top->atoms.atomname[i]);
        } else if (!strcmp(*top->atoms.resinfo[j].name,"ACE") || !strcmp(*top->atoms.resinfo[j].name,"ACG")){         // Residue ACE 
           if (!strcmp(*top->atoms.atomname[i],"C")) control+=8,at[0]=i;
           if (!strcmp(*top->atoms.atomname[i],"O")) control+=16,at[1]=i;
           if (!strcmp(*top->atoms.atomname[i],"CH3")) control+=32,at[5]=i;
           if (!strcmp(*top->atoms.atomname[i],"CA")) control+=32,at[5]=i; 
          // printf("Residue%2d is %s %s\n",j,*top->atoms.resinfo[j].name,*top->atoms.atomname[i]);
        } else {
            if (!strcmp(*top->atoms.atomname[i],"C")) control+=8,at[0]=i;
            if (!strcmp(*top->atoms.atomname[i],"O")) control+=16,at[1]=i;
            if (!strcmp(*top->atoms.atomname[i],"CA")) control+=32,at[5]=i;
            // printf("Residue%2d is %s %s\n",j,*top->atoms.resinfo[j].name,*top->atoms.atomname[i]);
        }
      }
      if(cyclic==1){
        if (j%snap->COchains==resNum){
         if (top->atoms.atom[i].resind==((j+1)-snap->COchains)){
          if (!strcmp(*top->atoms.resinfo[(j+1)-snap->COchains].name,"PRO")){  // Residue PRO
           proline=TRUE;
           if (!strcmp(*top->atoms.atomname[i],"N")) control++,at[2]=i;
           if (!strcmp(*top->atoms.atomname[i],"CD")) control+=2,at[3]=i;
           if (!strcmp(*top->atoms.atomname[i],"CA")) control+=4,at[4]=i;
           if (!strcmp(*top->atoms.atomname[i],"C2")) control+=4,at[4]=i;
           //printf("Residue%2d is %s %s\n",j,*top->atoms.resinfo[j].name,*top->atoms.atomname[i]);
          } else {
           if (!strcmp(*top->atoms.atomname[i],"N")) control++,at[2]=i;
           if (!strcmp(*top->atoms.atomname[i],"H")) control+=2,at[3]=i;
           if (!strcmp(*top->atoms.atomname[i],"HN")) control+=2,at[3]=i;
           if (!strcmp(*top->atoms.atomname[i],"CA")) control+=4,at[4]=i;
           if (!strcmp(*top->atoms.atomname[i],"C2")) control+=4,at[4]=i;
           //printf("Residue%2d is %s %s\n",j,*top->atoms.resinfo[j].name,*top->atoms.atomname[i]);
          }
         }
        } else { 
         if (top->atoms.atom[i].resind==j+1+orderN){ // Test if N terminal atoms are present
          if (!strcmp(*top->atoms.resinfo[j+1].name,"PRO")){  // Residue PRO
            proline=TRUE;
            if (!strcmp(*top->atoms.atomname[i],"N")) control++,at[2]=i;
            if (!strcmp(*top->atoms.atomname[i],"CD")) control+=2,at[3]=i;
            if (!strcmp(*top->atoms.atomname[i],"CA")) control+=4,at[4]=i;
            if (!strcmp(*top->atoms.atomname[i],"C2")) control+=4,at[4]=i;
            //printf("Residue%2d is %s %s\n",j,*top->atoms.resinfo[j].name,*top->atoms.atomname[i]);
          } else {
            if (!strcmp(*top->atoms.atomname[i],"N")) control++,at[2]=i;
            if (!strcmp(*top->atoms.atomname[i],"H")) control+=2,at[3]=i;
            if (!strcmp(*top->atoms.atomname[i],"HN")) control+=2,at[3]=i;
            if (!strcmp(*top->atoms.atomname[i],"CA")) control+=4,at[4]=i;
            if (!strcmp(*top->atoms.atomname[i],"C2")) control+=4,at[4]=i;
           //printf("Residue%2d is %s %s\n",j,*top->atoms.resinfo[j].name,*top->atoms.atomname[i]);
          }
         } 
        } 
        
      } else { 
        if (top->atoms.atom[i].resind==j+1+orderN){ // Test if N terminal atoms are present 
         if (!strcmp(*top->atoms.resinfo[j+1].name,"PRO")){  // Residue PRO 
           proline=TRUE;
           if (!strcmp(*top->atoms.atomname[i],"N")) control++,at[2]=i;
           if (!strcmp(*top->atoms.atomname[i],"CD")) control+=2,at[3]=i;
           if (!strcmp(*top->atoms.atomname[i],"CA")) control+=4,at[4]=i;
           if (!strcmp(*top->atoms.atomname[i],"C2")) control+=4,at[4]=i;
           //printf("Residue%2d is %s %s\n",j,*top->atoms.resinfo[j].name,*top->atoms.atomname[i]);
         } else if (!strcmp(*top->atoms.resinfo[j+1].name,"NAC")) {  // Residue NAC 
           if (!strcmp(*top->atoms.atomname[i],"N")) control++,at[2]=i;
           if (!strcmp(*top->atoms.atomname[i],"H")) control+=2,at[3]=i;
           if (!strcmp(*top->atoms.atomname[i],"HN")) control+=2,at[3]=i;
           if (!strcmp(*top->atoms.atomname[i],"CH3")) control+=4,at[4]=i;
           if (!strcmp(*top->atoms.atomname[i],"CA")) control+=4,at[4]=i;  
           //printf("Residue%2d is %s %s\n",j,*top->atoms.resinfo[j].name,*top->atoms.atomname[i]);
         } else { 
           if (!strcmp(*top->atoms.atomname[i],"N")) control++,at[2]=i;
           if (!strcmp(*top->atoms.atomname[i],"H")) control+=2,at[3]=i;
           if (!strcmp(*top->atoms.atomname[i],"HN")) control+=2,at[3]=i;
           if (!strcmp(*top->atoms.atomname[i],"CA")) control+=4,at[4]=i;
           if (!strcmp(*top->atoms.atomname[i],"C2")) control+=4,at[4]=i;
           //printf("Residue%2d is %s %s\n",j,*top->atoms.resinfo[j].name,*top->atoms.atomname[i]);
         }
        }
      }
    }

    if (control==63){ // We have an amide bond
      if (amide>snap->Namidebonds){
	printf("Fatal error!!!\n");
	printf("More amidebonds found in file than specified in input...\n");
	exit(-1);
      }
      amideGroups[amide].Nres=j+orderN;
      if((cyclic==1) && (j%snap->COchains==resNum)) {
        amideGroups[amide].Cres=((j+1)-snap->COchains);
      } else {
        amideGroups[amide].Cres=j+1+orderC;
      }
      amideGroups[amide].proline=proline;
      if (proline){
	printf("Amideunit number %d was identified as a proline amide bond.\n",amide);
      }
      for (i=0;i<6;i++) amideGroups[amide].a[i]=at[i];//printf(" %3d",at[i]+1);
      //printf(" Residue %d Amide unit %d Nres and Cres %d %d\n",j,amide,amideGroups[amide].Nres,amideGroups[amide].Cres); 
      amide++;
    }

  }
  // Test
  printf("\n");
  printf("Backbone Amide-I units: \n");
  printf(" Amide-I-units    Atom numbers \n");
  for (i=0;i<snap->NbackboneCO;i++){
    printf(" %d                %d %d %d %d %d %d\n",i,amideGroups[i].a[0],amideGroups[i].a[1],amideGroups[i].a[2],amideGroups[i].a[3],amideGroups[i].a[4],amideGroups[i].a[5]);
  }

  // Charge test
  //for (i=0;i<snap->NbackboneCO;i++){
   // q=0;
    //for (j=0;j<6;j++) q+=top->atoms.atom[amideGroups[i].a[j]].q;
    //printf("Amide group %2d charge: %6.3f\n",i,q);
  //}

  // Calculate charge difference
  q=0;
  for (j=0;j<4;j++) q+=top->atoms.atom[amideGroups[0].a[j]].q;
  snap->chargediff=-q;

  return; 
}

/* Identify the sidechain amide I groups */
void identifySidechainCO(t_snap *snap,t_topology *top,t_amide *amideGroups){
  int i,j;
  float q;
  int amide;
  int at[6]={0,0,0,0,0,0};

  printf("\n");
  amide=snap->NbackboneCO-1;

  // Find amidebonds
  for (j=0;j<snap->Nresidues;j++){ // Possible amidebonds
    if (!strcmp(*top->atoms.resinfo[j].name,"GLN")) {
      amide++;
      amideGroups[amide].Nres=amideGroups[amide-1].Nres+1;
      amideGroups[amide].Cres=amideGroups[amide-1].Cres+1;
      for (i=0;i<snap->Natoms;i++){
        if (top->atoms.atom[i].resind==j){
           if (!strcmp(*top->atoms.atomname[i],"CD"))  {
             amideGroups[amide].a[0]=i;
           } else if (!strcmp(*top->atoms.atomname[i],"OE1")) {
             amideGroups[amide].a[1]=i;
           } else if (!strcmp(*top->atoms.atomname[i],"CG"))  {
             amideGroups[amide].a[5]=i;
           } else if (!strcmp(*top->atoms.atomname[i],"NE2")) {
             amideGroups[amide].a[2]=i;
           } else if (!strcmp(*top->atoms.atomname[i],"HE21")) {
             amideGroups[amide].a[3]=i;
           } if (!strcmp(*top->atoms.atomname[i],"HE22")) {
             amideGroups[amide].a[4]=i;
           }
          // printf("Residue%2d is %s %s\n",j,*top->atoms.resinfo[j].name,*top->atoms.atomname[i]);
        }
      }
    } else if (!strcmp(*top->atoms.resinfo[j].name,"ASN")) {
      amide++;
      amideGroups[amide].Nres=amideGroups[amide-1].Nres+1;
      amideGroups[amide].Cres=amideGroups[amide-1].Cres+1;
      for (i=0;i<snap->Natoms;i++){
        if (top->atoms.atom[i].resind==j){
           if (!strcmp(*top->atoms.atomname[i],"CG"))  {
             amideGroups[amide].a[0]=i;           
           } else if (!strcmp(*top->atoms.atomname[i],"OD1")) {
             amideGroups[amide].a[1]=i;
           } else if (!strcmp(*top->atoms.atomname[i],"CB"))  {
             amideGroups[amide].a[5]=i;
           } else if (!strcmp(*top->atoms.atomname[i],"ND2")) {
             amideGroups[amide].a[2]=i;
           } else if (!strcmp(*top->atoms.atomname[i],"HD21")) {
             amideGroups[amide].a[3]=i;
           } if (!strcmp(*top->atoms.atomname[i],"HD22")) {
             amideGroups[amide].a[4]=i;
           }
          // printf("Residue%2d is %s %s\n",j,*top->atoms.resinfo[j].name,*top->atoms.atomname[i]);
        }
      }
    }
    
  }
  // Test
  printf("\n");
  printf("Sidechain Amide-I units: \n");
  printf(" Amide-I-units    Atom numbers \n");
  for (i=snap->NbackboneCO;i<snap->Namidebonds;i++){
    printf(" %d                %d %d %d %d %d %d\n",i,amideGroups[i].a[0],amideGroups[i].a[1],amideGroups[i].a[2],amideGroups[i].a[3],amideGroups[i].a[4],amideGroups[i].a[5]);
  }

  // Charge test
  //for (i=0;i<snap->Namidebonds;i++){
    //q=0;
    //for (j=0;j<6;j++) q+=top->atoms.atom[amideGroups[i].a[j]].q;
    //printf(" Sidechain amide group %2d charge: %6.3f\n",i,q);
  //} 
  return;
}

void calcCharge(t_snap *snap,t_topology *top){
  int N,a,b,i;
  float pol[3];
  float mag;
  float q;
  N=snap->Nresidues;
  q=0;
 
  // Loop over atoms
  for(a=0;a<snap->Natoms;a++){
    q+=top->atoms.atom[a].q;
  }
//  if (q>0.0001 || q<-0.0001){
//    printf("Cummulative charge NOT zero!\n");
//    exit(0);
//  } 
  return;
}

/* Diagonal Hamiltonian and Dipole moments for amide I units */
void calcDiag(t_pbc *pbc,t_snap *snap,t_map *map,t_map *Promap, t_Skinner *SkinnerAmI, t_Tokmakoff *TokmakoffAmI,int numMap,t_topology *top,float *Hamiltonian,float *dipole,t_amide *amideGroups,rvec *xG,int opls,int cyclic,int TDTorii,int nnTreatMap,t_select *select,int nr_frames) {

  float POT[7];
  float E[3][4],G[6][4];
  float Ecn[3][2];
  int i,j,jj,k,m,a,b,c,ip;
  int N;
  rvec x,y,z; // Unit vectors for molecular coordinate systems
  int inbond;
  rvec r,dip,mid;
  float d2,id,id3,id5;
  float aa2bohr=1.0/0.529177; 
  float bohr2aa=0.529177;
  float q;
  float qsum;
  float qcor,qcorC,qcorN;
  float nndamp;
  float nq;
  int residues; // Number of residues including solvent
  int nchargeGroups; // Number of charge groups in the system
  rvec *residue;
  real *mass;
  real *NormFac;
  real ma;
  rvec *chargeGroup;
  rvec ra,r1,r2,rb;
  rvec rcgs,rcm;
  int chainNR,ntermchain,ctermchain;
  int *AtomsInCgs;
  int Donotcontinue;

  N=snap->Namidebonds;
  residues=top->atoms.nres;
  nchargeGroups=top->cgs.nr;
  //printf("NumChrgeGroup %d\n",top->cgs.nr);
  residue=(rvec *)calloc(residues,sizeof(rvec));
  chargeGroup=(rvec *)calloc(nchargeGroups,sizeof(rvec));
  mass=(real *)calloc(residues,sizeof(real));
  NormFac=(real *)calloc(nchargeGroups,sizeof(real));
  AtomsInCgs=(int *)calloc(snap->Natoms,sizeof(int));

  // Clear residue and charge group position array
  for (i=0;i<residues;i++) residue[i][XX]=residue[i][YY]=residue[i][ZZ]=0;
  for (i=0;i<nchargeGroups;i++) chargeGroup[i][XX]=chargeGroup[i][YY]=chargeGroup[i][ZZ]=0;

  // Find residue position
  for (i=0;i<snap->Natoms;i++){    
    ma=top->atoms.atom[i].m;
    mass[top->atoms.atom[i].resind]+=ma;
    residue[top->atoms.atom[i].resind][XX]+=ma*xG[i][XX];
    residue[top->atoms.atom[i].resind][YY]+=ma*xG[i][YY];
    residue[top->atoms.atom[i].resind][ZZ]+=ma*xG[i][ZZ];
  }
  // Normalize
  for (i=0;i<residues;i++){
    if (mass[i]>0) residue[i][XX]/=mass[i],residue[i][YY]/=mass[i],residue[i][ZZ]/=mass[i];
  }

  // Find charge group position
  j=0;   
  for (i=0;i<snap->Natoms;i++){ 
    NormFac[j]+=1.0;
    chargeGroup[j][XX]+=xG[i][XX];
    chargeGroup[j][YY]+=xG[i][YY];
    chargeGroup[j][ZZ]+=xG[i][ZZ];
    AtomsInCgs[i]=j;
    //printf("atom %d charge group %d\n",i,AtomsInCgs[i]); 
    if((i+1)==(top->cgs.index[j+1])) j+=1;
  }
  // Normalize
  for (i=0;i<nchargeGroups;i++){
    if (NormFac[i]>0) chargeGroup[i][XX]/=NormFac[i],chargeGroup[i][YY]/=NormFac[i],chargeGroup[i][ZZ]/=NormFac[i];
  }
 
  // Loop over units
 for (i=0;i<N;i++){
  // Define coordinate system
   pbc_dx(pbc,xG[amideGroups[i].a[1]],xG[amideGroups[i].a[0]],x); // C->O
   unitv(x,x);
   pbc_dx(pbc,xG[amideGroups[i].a[2]],xG[amideGroups[i].a[0]],y); // C->N
   project(x,y);
   unitv(y,y);
   cprod(x,y,z);
   
   // Identify chain number
   chainNR=i/snap->COchains;
   ntermchain=chainNR*snap->COchains;
   ctermchain=(chainNR+1)*snap->COchains-1;     
   // Calculate Potential, Field and Gradient 
   for (j=0;j<4;j++){
     Donotcontinue=1;
     POT[j]=0.0; 
     for (k=0;k<3;k++) E[k][j]=0.0;
     for (k=0;k<6;k++) G[k][j]=0.0;
     if(amideGroups[i].proline!=TRUE){
      if((numMap==0) && ((j==1)||(j==3))){
        Donotcontinue=0;
      } else if((numMap==2) && ((j==0)||(j==2)||(j==3))){
        Donotcontinue=0;
      }
     } else {
       if((numMap==2) && ((j==0)||(j==2)||(j==3))){
        Donotcontinue=0;
       } 
     } 
     /* Proceed to calculate electric fields, gradients etc */
     if(Donotcontinue==1){
       qsum=0.0;
      // Loop over atoms 
      for (a=0;a<snap->Natoms;a++){
	  // Exclude atoms in the amide bond (including C(a))
	  inbond=0;
          if(i>=snap->NbackboneCO){
           for (b=0;b<5;b++) if (amideGroups[i].a[b]==a) inbond=1;
          } else { 
	    for (b=0;b<6;b++) if (amideGroups[i].a[b]==a) inbond=1;
            if((amideGroups[i].proline==TRUE) && (amideGroups[i].a[3]==a) && (numMap==2)) inbond=0; /* Include CD atom for Tokmakoff map */
          } 
	  if((numMap!=2) && (i<snap->NbackboneCO)){
	   // Exclude hydrogens on C(a) or on methyl carbon
	   if (top->atoms.atom[a].resind==amideGroups[i].Cres || top->atoms.atom[a].resind==amideGroups[i].Nres){
	    if (!strcmp(*top->atoms.atomname[a],"HA")) inbond=1;
	    if (!strcmp(*top->atoms.atomname[a],"HA1")) inbond=1;
	    if (!strcmp(*top->atoms.atomname[a],"HA2")) inbond=1;
	    if (!strcmp(*top->atoms.atomname[a],"HA3")) inbond=1;
            if (!strcmp(*top->atoms.atomname[a],"HH31")) inbond=1;
	    if (!strcmp(*top->atoms.atomname[a],"HH32")) inbond=1;
	    if (!strcmp(*top->atoms.atomname[a],"HH33")) inbond=1;
           } 

           if (top->atoms.atom[a].resind==amideGroups[i].Cres){
            if(amideGroups[i].proline==TRUE) {
              if (!strcmp(*top->atoms.atomname[a],"HD1")) inbond=1;
              if (!strcmp(*top->atoms.atomname[a],"HD2")) inbond=1;
            }
	   }

	   //Setting charges of Nearest Neighbor to zero 
	   if ((nnTreatMap==1) && (N>1)){ // nearest neighbor maps are used if there are more than one amide I units
            /************************
	    * Nterminal neighbour   * 
            ************************/
            if((i==ntermchain) && (cyclic==1)){ //For the first amide unit of a cyclic peptide chain
              if (top->atoms.atom[a].resind==amideGroups[ctermchain].Nres){ //exclude atoms from the last unit of a peptide chain 
               if (!strcmp(*top->atoms.atomname[a],"C")){
                 if (amideGroups[ctermchain].a[0]==a) inbond=1;
               }
               if (!strcmp(*top->atoms.atomname[a],"O")) {
                 if (amideGroups[ctermchain].a[1]==a) inbond=1;
               }
               if(opls!=1){
                 if (!strcmp(*top->atoms.atomname[a],"CA")) inbond=1; 
                 if (!strcmp(*top->atoms.atomname[a],"HA1")) inbond=1;
                 if (!strcmp(*top->atoms.atomname[a],"HA2")) inbond=1;
                 if (!strcmp(*top->atoms.atomname[a],"HA3")) inbond=1;
                 if (!strcmp(*top->atoms.atomname[a],"HA")) inbond=1; 
               }
              }
              if (top->atoms.atom[a].resind==amideGroups[ntermchain].Nres){ //exclude atoms from the first unit of a peptide chain
               if (!strcmp(*top->atoms.atomname[a],"N")){
                 if (amideGroups[ctermchain].a[2]==a) inbond=1;
               }
               if (!strcmp(*top->atoms.atomname[a],"H") || !strcmp(*top->atoms.atomname[a],"D")){
                 if (amideGroups[ctermchain].a[3]==a) inbond=1;
               }
               if(amideGroups[ctermchain].proline==TRUE){ // if the last unit of a chain is Proline
                 if (!strcmp(*top->atoms.atomname[a],"CD")) inbond=1;
                 if (!strcmp(*top->atoms.atomname[a],"HD1")) inbond=1;
                 if (!strcmp(*top->atoms.atomname[a],"HD2")) inbond=1;
                 if (!strcmp(*top->atoms.atomname[a],"CB")) inbond=1;
                 if (!strcmp(*top->atoms.atomname[a],"HB1")) inbond=1;
                 if (!strcmp(*top->atoms.atomname[a],"HB2")) inbond=1;
                 if (!strcmp(*top->atoms.atomname[a],"CG")) inbond=1;
                 if (!strcmp(*top->atoms.atomname[a],"HG1")) inbond=1;
                 if (!strcmp(*top->atoms.atomname[a],"HG2")) inbond=1;
               }
              }
            }
            if(i>ntermchain){ // units except the first one of a peptide chain
	     if (top->atoms.atom[a].resind==amideGroups[i].Nres-1){
	      if (!strcmp(*top->atoms.atomname[a],"C")){
		 if (amideGroups[i-1].a[0]==a) inbond=1;
	      }
	      if (!strcmp(*top->atoms.atomname[a],"O")) {
                   if (amideGroups[i-1].a[1]==a) inbond=1;
              } 
              if(opls!=1){ 
		 if (!strcmp(*top->atoms.atomname[a],"CA")) inbond=1;  
		 if (!strcmp(*top->atoms.atomname[a],"HA1")) inbond=1; 
		 if (!strcmp(*top->atoms.atomname[a],"HA2")) inbond=1; 
		 if (!strcmp(*top->atoms.atomname[a],"HA3")) inbond=1; 
		 if (!strcmp(*top->atoms.atomname[a],"HA")) inbond=1;  
                 if (!strcmp(*top->atoms.atomname[a],"CH3")) inbond=1; 
                 if (!strcmp(*top->atoms.atomname[a],"HH31")) inbond=1;
	         if (!strcmp(*top->atoms.atomname[a],"HH32")) inbond=1;
	         if (!strcmp(*top->atoms.atomname[a],"HH33")) inbond=1;
              }
	     } else if (top->atoms.atom[a].resind==amideGroups[i].Nres){
	       if (!strcmp(*top->atoms.atomname[a],"N")){
		 if (amideGroups[i-1].a[2]==a) inbond=1;
	       }	      
	       if (!strcmp(*top->atoms.atomname[a],"H") || !strcmp(*top->atoms.atomname[a],"D")){ 
		 if (amideGroups[i-1].a[3]==a) inbond=1;
               }
               if(amideGroups[i-1].proline==TRUE) {
                 if (!strcmp(*top->atoms.atomname[a],"CD")) inbond=1; 
                 if (!strcmp(*top->atoms.atomname[a],"HD1")) inbond=1;
                 if (!strcmp(*top->atoms.atomname[a],"HD2")) inbond=1;
                 if (!strcmp(*top->atoms.atomname[a],"CB")) inbond=1; 
                 if (!strcmp(*top->atoms.atomname[a],"HB1")) inbond=1;
                 if (!strcmp(*top->atoms.atomname[a],"HB2")) inbond=1;
                 if (!strcmp(*top->atoms.atomname[a],"CG")) inbond=1; 
                 if (!strcmp(*top->atoms.atomname[a],"HG1")) inbond=1;
                 if (!strcmp(*top->atoms.atomname[a],"HG2")) inbond=1;
               } 
	     }
            }
	    /*********************
	    * Cterminal neighbour *
            **********************/
            if((i==ctermchain) && (cyclic==1)){ // For the last amide unit of a cyclic peptide chain
              if (top->atoms.atom[a].resind==amideGroups[ntermchain].Nres){ //exclude atoms from N-term atoms of the first unit
               if (!strcmp(*top->atoms.atomname[a],"C")){
                 if (amideGroups[ntermchain].a[0]==a) inbond=1;
               }
               if (!strcmp(*top->atoms.atomname[a],"O")) {
                 if (amideGroups[ntermchain].a[1]==a) inbond=1;
               }
              }

              if (top->atoms.atom[a].resind==amideGroups[ntermchain].Cres){ //exclude atoms from C-term atoms of the first unit
               if (!strcmp(*top->atoms.atomname[a],"N")){
                 if (amideGroups[ntermchain].a[2]==a) inbond=1;
               }
               if (!strcmp(*top->atoms.atomname[a],"H") || !strcmp(*top->atoms.atomname[a],"D")){
                 if (amideGroups[ntermchain].a[3]==a) inbond=1;
               }
               if (!strcmp(*top->atoms.atomname[a],"CA")) inbond=1;
               if (!strcmp(*top->atoms.atomname[a],"HA1")) inbond=1;
               if (!strcmp(*top->atoms.atomname[a],"HA2")) inbond=1;
               if (!strcmp(*top->atoms.atomname[a],"HA3")) inbond=1;
               if (!strcmp(*top->atoms.atomname[a],"HA")) inbond=1;
               if(amideGroups[ntermchain].proline==TRUE){
                if (!strcmp(*top->atoms.atomname[a],"CD")) inbond=1;
                if (!strcmp(*top->atoms.atomname[a],"HD1")) inbond=1;
                if (!strcmp(*top->atoms.atomname[a],"HD2")) inbond=1;
                if (!strcmp(*top->atoms.atomname[a],"CB")) inbond=1;
                if (!strcmp(*top->atoms.atomname[a],"HB1")) inbond=1;
                if (!strcmp(*top->atoms.atomname[a],"HB2")) inbond=1;
                if (!strcmp(*top->atoms.atomname[a],"CG")) inbond=1;
                if (!strcmp(*top->atoms.atomname[a],"HG1")) inbond=1;
                if (!strcmp(*top->atoms.atomname[a],"HG2")) inbond=1;
               }

              }

            }
            /********************************************************************************************************************/
            if(amideGroups[i].proline==TRUE){ // exclude the atoms if they belong to proline   
             if((i==ctermchain) && (top->atoms.atom[a].resind==amideGroups[i].Cres)){ // which is the last unit of a peptide chain
              if (!strcmp(*top->atoms.atomname[a],"CB")) inbond=1;
              if (!strcmp(*top->atoms.atomname[a],"HB1")) inbond=1;
              if (!strcmp(*top->atoms.atomname[a],"HB2")) inbond=1;
              if (!strcmp(*top->atoms.atomname[a],"CG")) inbond=1;
              if (!strcmp(*top->atoms.atomname[a],"HG1")) inbond=1;
              if (!strcmp(*top->atoms.atomname[a],"HG2")) inbond=1;
             }
            }  
            /************************************************************************/
            if(i<ctermchain){ // For all units except the last one in a peptide chain
	     if (top->atoms.atom[a].resind==amideGroups[i].Cres){
	      if (!strcmp(*top->atoms.atomname[a],"C")){
	       if (amideGroups[i+1].a[0]==a) inbond=1;
	      }
              if (!strcmp(*top->atoms.atomname[a],"O")){
               if (amideGroups[i+1].a[1]==a) inbond=1;
              }
              if(amideGroups[i].proline==TRUE) {
                 if (!strcmp(*top->atoms.atomname[a],"CB")) inbond=1; 
                 if (!strcmp(*top->atoms.atomname[a],"HB1")) inbond=1;
                 if (!strcmp(*top->atoms.atomname[a],"HB2")) inbond=1;
                 if (!strcmp(*top->atoms.atomname[a],"CG")) inbond=1; 
                 if (!strcmp(*top->atoms.atomname[a],"HG1")) inbond=1;
                 if (!strcmp(*top->atoms.atomname[a],"HG2")) inbond=1;
              }
	     } else if (top->atoms.atom[a].resind==amideGroups[i].Cres+1){
	       if (!strcmp(*top->atoms.atomname[a],"N")){
		if (amideGroups[i+1].a[2]==a) inbond=1;
	       }
               if (!strcmp(*top->atoms.atomname[a],"H") || !strcmp(*top->atoms.atomname[a],"D")) {
                if (amideGroups[i+1].a[3]==a) inbond=1; 
               } 
               if(amideGroups[i+1].proline==TRUE) {  
                if (!strcmp(*top->atoms.atomname[a],"CD")) inbond=1; 
                if (!strcmp(*top->atoms.atomname[a],"HD1")) inbond=1; 
                if (!strcmp(*top->atoms.atomname[a],"HD2")) inbond=1;
                if (!strcmp(*top->atoms.atomname[a],"CB")) inbond=1; 
                if (!strcmp(*top->atoms.atomname[a],"HB1")) inbond=1;
                if (!strcmp(*top->atoms.atomname[a],"HB2")) inbond=1; 
                if (!strcmp(*top->atoms.atomname[a],"CG")) inbond=1; 
                if (!strcmp(*top->atoms.atomname[a],"HG1")) inbond=1;
                if (!strcmp(*top->atoms.atomname[a],"HG2")) inbond=1; 
              }
	      if (!strcmp(*top->atoms.atomname[a],"CA")) inbond=1;  
	      if (!strcmp(*top->atoms.atomname[a],"HA")) inbond=1;  
	      if (!strcmp(*top->atoms.atomname[a],"HA1")) inbond=1; 
	      if (!strcmp(*top->atoms.atomname[a],"HA2")) inbond=1; 
	      if (!strcmp(*top->atoms.atomname[a],"HA3")) inbond=1; 
              if (!strcmp(*top->atoms.atomname[a],"CH3")) inbond=1; 
              if (!strcmp(*top->atoms.atomname[a],"HH31")) inbond=1;
              if (!strcmp(*top->atoms.atomname[a],"HH32")) inbond=1;
              if (!strcmp(*top->atoms.atomname[a],"HH33")) inbond=1;
	     }
            }
            /***************************************************************************************/
            if((i==(ctermchain-1)) && (cyclic==1)){ //For the unit just before the last one of a cyclic peptide chain
              if (top->atoms.atom[a].resind==amideGroups[ntermchain].Nres){
               if (!strcmp(*top->atoms.atomname[a],"N")){
                 if (amideGroups[ctermchain].a[2]==a) inbond=1;
               }
               if (!strcmp(*top->atoms.atomname[a],"H") || !strcmp(*top->atoms.atomname[a],"D")){
                 if (amideGroups[ctermchain].a[3]==a) inbond=1;
               }
               if (!strcmp(*top->atoms.atomname[a],"CA")) inbond=1; 
               if (!strcmp(*top->atoms.atomname[a],"HA1")) inbond=1;
               if (!strcmp(*top->atoms.atomname[a],"HA2")) inbond=1;
               if (!strcmp(*top->atoms.atomname[a],"HA3")) inbond=1;
               if (!strcmp(*top->atoms.atomname[a],"HA")) inbond=1; 
               if(amideGroups[ctermchain].proline==TRUE){
                if (!strcmp(*top->atoms.atomname[a],"CD")) inbond=1; 
                if (!strcmp(*top->atoms.atomname[a],"HD1")) inbond=1;
                if (!strcmp(*top->atoms.atomname[a],"HD2")) inbond=1;
                if (!strcmp(*top->atoms.atomname[a],"CB")) inbond=1; 
                if (!strcmp(*top->atoms.atomname[a],"HB1")) inbond=1;
                if (!strcmp(*top->atoms.atomname[a],"HB2")) inbond=1;
                if (!strcmp(*top->atoms.atomname[a],"CG")) inbond=1; 
                if (!strcmp(*top->atoms.atomname[a],"HG1")) inbond=1;
                if (!strcmp(*top->atoms.atomname[a],"HG2")) inbond=1;
               }
              }
            }            
            	  
	   }
          } 
          
	  // Select particular subsystem
	  // include All molecules 
          if (select->all==1) {  }

          // Only Water
          else if (select->water==1){
            if (strcmp(*(top->atoms.resinfo[top->atoms.atom[a].resind].name),"SOL")) inbond=1;
          }

          // Only Water and Ion
          else if (select->water_and_ion==1){
            if (strcmp(*(top->atoms.resinfo[top->atoms.atom[a].resind].name),"SOL") ||
                strcmp(*(top->atoms.resinfo[top->atoms.atom[a].resind].name),"K+") ||
                strcmp(*(top->atoms.resinfo[top->atoms.atom[a].resind].name),"K") ||
                strcmp(*(top->atoms.resinfo[top->atoms.atom[a].resind].name),"NA+") ||
                strcmp(*(top->atoms.resinfo[top->atoms.atom[a].resind].name),"NA") ||
                strcmp(*(top->atoms.resinfo[top->atoms.atom[a].resind].name),"CL-") ||
                strcmp(*(top->atoms.resinfo[top->atoms.atom[a].resind].name),"CL")) inbond=1;
          }

          // Only Protein
          else if (select->protein==1){
	    if (!strcmp(*(top->atoms.resinfo[top->atoms.atom[a].resind].name),"SOL") || 
                !strcmp(*(top->atoms.resinfo[top->atoms.atom[a].resind].name),"DMPC") ||
                !strcmp(*(top->atoms.resinfo[top->atoms.atom[a].resind].name),"DOPC") ||  
                !strcmp(*(top->atoms.resinfo[top->atoms.atom[a].resind].name),"POPC") ||
                !strcmp(*(top->atoms.resinfo[top->atoms.atom[a].resind].name),"POPE") ||  
                !strcmp(*(top->atoms.resinfo[top->atoms.atom[a].resind].name),"DPPG") ||
                !strcmp(*(top->atoms.resinfo[top->atoms.atom[a].resind].name),"DPTAP") ||
                !strcmp(*(top->atoms.resinfo[top->atoms.atom[a].resind].name),"K+") ||
                !strcmp(*(top->atoms.resinfo[top->atoms.atom[a].resind].name),"K") ||
                !strcmp(*(top->atoms.resinfo[top->atoms.atom[a].resind].name),"CL-") ||
                !strcmp(*(top->atoms.resinfo[top->atoms.atom[a].resind].name),"CL") ||
                !strcmp(*(top->atoms.resinfo[top->atoms.atom[a].resind].name),"NA+") ||     
                !strcmp(*(top->atoms.resinfo[top->atoms.atom[a].resind].name),"NA")) inbond=1;
          }

          // Only Lipids
          else if (select->lipid==1){
            if (strcmp(*(top->atoms.resinfo[top->atoms.atom[a].resind].name),"DMPC") ||
                strcmp(*(top->atoms.resinfo[top->atoms.atom[a].resind].name),"POPC") ||
                strcmp(*(top->atoms.resinfo[top->atoms.atom[a].resind].name),"DOPC") ||
                strcmp(*(top->atoms.resinfo[top->atoms.atom[a].resind].name),"POPE") ||
                strcmp(*(top->atoms.resinfo[top->atoms.atom[a].resind].name),"DPPG") ||
                strcmp(*(top->atoms.resinfo[top->atoms.atom[a].resind].name),"DPTAP")) inbond=1; 
          }

          // Only CL- 
          else if (select->cl==1){
            if (strcmp(*(top->atoms.resinfo[top->atoms.atom[a].resind].name),"CL-") || 
                strcmp(*(top->atoms.resinfo[top->atoms.atom[a].resind].name),"CL")) inbond=1;
          }

          // Only NA+ 
          else if (select->na==1){
            if (strcmp(*(top->atoms.resinfo[top->atoms.atom[a].resind].name),"NA+") ||
                strcmp(*(top->atoms.resinfo[top->atoms.atom[a].resind].name),"NA")) inbond=1;
          }
          
         // Only K+ 
          else if (select->k==1){
            if (strcmp(*(top->atoms.resinfo[top->atoms.atom[a].resind].name),"K+") ||
                strcmp(*(top->atoms.resinfo[top->atoms.atom[a].resind].name),"K")) inbond=1;
          } 

          if(inbond==0){
            if((numMap==0) && (amideGroups[i].proline!=TRUE)){ /* Skinner map */
              pbc_dx(pbc,chargeGroup[AtomsInCgs[a]],chargeGroup[AtomsInCgs[amideGroups[i].a[j]]],rcgs);
              if(sqrt(iprod(rcgs,rcgs)) < (snap->cutoff)){ 
                // Find distance vector
                pbc_dx(pbc,xG[amideGroups[i].a[j]],chargeGroup[AtomsInCgs[a]],r1);
                pbc_dx(pbc,chargeGroup[AtomsInCgs[a]],xG[a],r2);
                rvec_add(r2,r1,r); 
              } else {
                inbond=1;
              }
            } else if((numMap==1) || (numMap==2) || (amideGroups[i].proline == TRUE)){ /* Jansen map */
                     pbc_dx(pbc,residue[top->atoms.atom[a].resind],residue[top->atoms.atom[amideGroups[i].a[0]].resind],rcm);
                     if (iprod(rcm,rcm)<snap->cutoff*snap->cutoff){
                       // Find distance vector
                       pbc_dx(pbc,xG[amideGroups[i].a[j]],xG[a],r); 
                     } else {
                       inbond=1; 
                     }
            } 
          }
          if(inbond==0){
          // Find length in Bohr
             svmul(aa2bohr*10,r,r);
             d2=iprod(r,r);
             id=1.0/sqrt(d2);
             id3=id*id*id;
              q=top->atoms.atom[a].q;
              qsum+=q;
              POT[j]+=id*q;
              electricField(E,x,y,z,r,q,id3,j);
              if((numMap==1) || ((amideGroups[i].proline == TRUE) && (numMap==0))){  
                 electricFieldGradient(G,x,y,z,r,d2,q,id3,j);
              } 
          } 
	
      }

      //if (fabs(qsum)>0.00001) {
	//if (j==0 && nr_frames==0){
    	  //printf("WARNING!!! ");
 	  //printf("TOTAL CHARGE NONZERO ");
	  //printf("For amide unit %d\n",i);
	  //printf(":Calculated sum %f\n",qsum);
	//} 
      //}
      
     }
   }
   //printf("%d %f %f\n",i,E[0][0],E[0][2]);
   // Calculate frequency
   if(amideGroups[i].proline==TRUE){
      if(numMap==2) { /* If Tokmakoff map is used */
        Hamiltonian[i*N+i]=TokmakoffAmI->omegaPro;
        // Field
        Hamiltonian[i*N+i]+=E[0][1]*TokmakoffAmI->cEPro[0][1];
      } else {
        Hamiltonian[i*N+i]=Promap->omega;
        // Field
        for (j=0;j<4;j++){
          for (k=0;k<3;k++){
            Hamiltonian[i*N+i]+=E[k][j]*Promap->cE[k][j];
          }
          // Field Gradients
          for (k=0;k<6;k++){
            Hamiltonian[i*N+i]+=G[k][j]*Promap->cG[k][j];
          }
        }
      }

   } else if (numMap==0) { /* If Skinner map is used */
     if(i>=snap->NbackboneCO) {
       // Gas phase 
       Hamiltonian[i*N+i]=SkinnerAmI->omega_S;
       // Field
       Hamiltonian[i*N+i]+=(E[0][0]*SkinnerAmI->cE[0][4])+(E[0][2]*SkinnerAmI->cE[0][6]);
     } else {
       // Gas phase
       Hamiltonian[i*N+i]=SkinnerAmI->omega;
       // Field
       Hamiltonian[i*N+i]+=(E[0][0]*SkinnerAmI->cE[0][0])+(E[0][2]*SkinnerAmI->cE[0][2]);
     }
     
   } else if (numMap==2) { /* If Tokmakoff map is used */
     Hamiltonian[i*N+i]=TokmakoffAmI->omega;
     // Field
     Hamiltonian[i*N+i]+=E[0][1]*TokmakoffAmI->cE[0][1];

   } else { /* If Jansen map is used */
     Hamiltonian[i*N+i]=map->omega;
     // Potential
     for (j=0;j<map->maptype;j++){
       Hamiltonian[i*N+i]+=POT[j]*map->p[j];
     }
     // Field
     for (j=0;j<4;j++){
       for (k=0;k<3;k++){
	 Hamiltonian[i*N+i]+=E[k][j]*map->cE[k][j];
       }
       // Gradients
       for (k=0;k<6;k++){
 	 Hamiltonian[i*N+i]+=G[k][j]*map->cG[k][j];
       }
     }

   }

   // Add shifts due to nearest neighbour using map
   if (nnTreatMap==1){
      // Don't do if we only account for solvent, sidechain or lipid
      if ((i>=snap->NbackboneCO) || (numMap==2) || (select->water==1) || (select->water_and_ion==1) || (select->lipid==1) || (select->cl==1) || (select->na==1) || (select->k==1) || (select->no_nnfs==1)){
      } else {
        nnMap(pbc,snap,top,Hamiltonian,amideGroups,i,cyclic,xG);
      }
   }
   
   //printf("%d %d %d %d %d %f\n",i,amideGroups[i].Nres,amideGroups[i].Cres,ntermchain,ctermchain,Hamiltonian[i*N+i]); 
   // Torii and Tasumi Transition Dipole Moments
   if(TDTorii==1){
     dip[0]=0.0;dip[1]=0.0;dip[2]=0.0; 
     TDM_Torii(pbc,snap,top,xG,amideGroups,i,dip);
     for (m=0;m<3;m++){
       dipole[3*i+m]=dip[m];
     }
   } else {
     if(numMap==0){
      printf("Exit Error: Please Use Torii map for the transition dipoles if the Skinner map is used \n");
      printf("In the input file insert the option: TransitionDipole Torii \n");
      exit(-1); 
     }
     if(amideGroups[i].proline==TRUE) {
        // Calculate dipoles in local coordinates
        for (m=0;m<3;m++){
         dip[m]=Promap->mu[m];
          for (j=0;j<4;j++){
            for (k=0;k<3;k++){
              dip[m]+=E[k][j]*Promap->mucE[m][k][j];
            }
            for (k=0;k<6;k++) dip[m]+=G[k][j]*Promap->mucG[m][k][j];
          }
        }
        // Transform to global coordinates
        for (m=0;m<3;m++){
          dipole[3*i+m]=dip[0]*x[m]+dip[1]*y[m]+dip[2]*z[m];
        }  
     } else {
       // Calculate dipole in local coordinates
       for (m=0;m<3;m++){
         dip[m]=map->mu[m];
         for (j=0;j<4;j++) dip[m]+=POT[j]*map->mup[m][j];
         for (j=0;j<4;j++){
           for (k=0;k<3;k++){
             dip[m]+=E[k][j]*map->mucE[m][k][j];
           }
           for (k=0;k<6;k++) dip[m]+=G[k][j]*map->mucG[m][k][j];
         }
       }

       // Transform to global coordinates
       for (m=0;m<3;m++){
         dipole[3*i+m]=dip[0]*x[m]+dip[1]*y[m]+dip[2]*z[m];
       }
     } 
   }
 }
  free(residue);
  free(mass);
  free(AtomsInCgs);
  free(chargeGroup);
  free(NormFac);
  return;
}

// Calculate the electric field for Jansen map
void electricField(float E[3][4],rvec x,rvec y,rvec z,rvec r,real q,real id3,int j){
  E[0][j]+=q*id3*iprod(r,x);
  E[1][j]+=q*id3*iprod(r,y);
  E[2][j]+=q*id3*iprod(r,z);
  return;
}

// Calculate the electric field gradient
void electricFieldGradient(float G[6][4],rvec x,rvec y,rvec z,rvec r,real d2,real q,real id3,int j){
  G[0][j]+=(1.0-3.0*iprod(x,r)*iprod(x,r)/d2)*q*id3;
  G[1][j]+=(1.0-3.0*iprod(y,r)*iprod(y,r)/d2)*q*id3;
  G[2][j]+=(1.0-3.0*iprod(z,r)*iprod(z,r)/d2)*q*id3;
  G[3][j]+=(-3.0*iprod(x,r)*iprod(y,r)/d2)*q*id3;
  G[4][j]+=(-3.0*iprod(x,r)*iprod(z,r)/d2)*q*id3;
  G[5][j]+=(-3.0*iprod(y,r)*iprod(z,r)/d2)*q*id3;
  return;
}

void nnMap(t_pbc *pbc,t_snap *snap,t_topology *top,float *Hamiltonian,t_amide *amideGroups,int i,int cyclic,rvec *xG){
  int N;
  int Nphi,Npsi;
  float t,u,y1,y2,y3,y4,x1l,x1h,x2l,x2h;
  int tem;
  float delta;
  int dim;
  float space;
  int j;
  int map; 
  int chainNR;
  int ntermchain;
  int ctermchain;
  int CO_per_chain;  // Carbonyl group per Chain
  t_rama rama;

  N=snap->Namidebonds;

  // Identify chain number
  chainNR=i/snap->COchains;
  ntermchain=chainNR*snap->COchains;
  ctermchain=(chainNR+1)*snap->COchains-1;

  dim=13;
  space=30;
  CO_per_chain=snap->COchains; 
  /**************************
  * Shift due to Nterminal  *
  ***************************/ 
  if (i==ntermchain){
   if(cyclic==1){
    //Check which NNFS map to be used
    map=cis_or_trans(pbc,snap,top,amideGroups,i,cyclic,xG);
    rama=Ramachandran(pbc,snap,top,amideGroups,i,ctermchain,xG);
    //printf("Ramachandran angles: Phi=%6.1f, Psi=%6.1f, map=%d, i=%d, j=%d \n",rama.phi,rama.psi,map,i,ctermchain);
    Nphi=(int) (rama.phi+180.0)/space;
    Npsi=(int) (rama.psi+180.0)/space;
    if (Nphi==dim-1) Nphi=dim-2;
    if (Npsi==dim-1) Npsi=dim-2;
    if (Nphi>=0 && Nphi<dim-1 && Npsi>=0 && Nphi<dim-1){
      x1l=Nphi*space-180;;
      x1h=x1l+space;
      x2l=Npsi*space-180;
      x2h=x2l+space;
      //Nterm
      y1=NtermShiftS(map,(Npsi)*dim+(Nphi));
      y2=NtermShiftS(map,(Npsi+1)*dim+(Nphi));
      y3=NtermShiftS(map,(Npsi+1)*dim+(Nphi+1));
      y4=NtermShiftS(map,(Npsi)*dim+(Nphi+1));
      //printf("%f %f %f %f\n",y1,y2,y3,y4);
      u=(rama.phi-x1l)/(x1h-x1l);
      t=(rama.psi-x2l)/(x2h-x2l);
      //printf("%f %f\n",t,u);
      delta=(1-t)*(1-u)*y1+t*(1-u)*y2+t*u*y3+(1-t)*u*y4;
      Hamiltonian[i*N+i]+=delta;
      //printf("NShift %f\n",delta);
    } else {
      printf("Ill defined Ramachandran angles between unit %d and %d!\n",i,i-1);
      printf("The nearest neighbor shift will be set to 0.\n");
    }
   }
  } else {
    map=cis_or_trans(pbc,snap,top,amideGroups,i,cyclic,xG);
    rama=Ramachandran(pbc,snap,top,amideGroups,i,i-1,xG);
    //printf("Ramachandran angles: Phi=%6.1f, Psi=%6.1f, map=%d, i=%d, j=%d \n",rama.phi,rama.psi,map,i,i-1);
    Nphi=(int) (rama.phi+180.0)/space;
    Npsi=(int) (rama.psi+180.0)/space;
    if (Nphi==dim-1) Nphi=dim-2;
    if (Npsi==dim-1) Npsi=dim-2;
    if (Nphi>=0 && Nphi<dim-1 && Npsi>=0 && Nphi<dim-1){
      x1l=Nphi*space-180;;
      x1h=x1l+space;
      x2l=Npsi*space-180;
      x2h=x2l+space;
      //Nterm
      y1=NtermShiftS(map,(Npsi)*dim+(Nphi));
      y2=NtermShiftS(map,(Npsi+1)*dim+(Nphi));
      y3=NtermShiftS(map,(Npsi+1)*dim+(Nphi+1));
      y4=NtermShiftS(map,(Npsi)*dim+(Nphi+1));
      //printf("%f %f %f %f\n",y1,y2,y3,y4);
      u=(rama.phi-x1l)/(x1h-x1l);
      t=(rama.psi-x2l)/(x2h-x2l);
      //printf("%f %f\n",t,u);
      delta=(1-t)*(1-u)*y1+t*(1-u)*y2+t*u*y3+(1-t)*u*y4;
      Hamiltonian[i*N+i]+=delta;
      //printf("NShift %f\n",delta);
    } else {
      printf("Ill defined Ramachandran angles between unit %d and %d!\n",i,i-1);
      printf("The nearest neighbor shift will be set to 0.\n");
    }
  }
  /**************************
  * Shift due to Cterminal  *
  ***************************/
  if (i==ctermchain){
   //including Ctermshift for the last unit of a cyclic peptide chain 
   if(cyclic==1){
    map=cis_or_trans(pbc,snap,top,amideGroups,ntermchain,cyclic,xG);
    rama=Ramachandran(pbc,snap,top,amideGroups,ntermchain,ctermchain,xG);
    //printf("Ramachandran angles: Phi=%6.1f, Psi=%6.1f, map=%d, i=%d, j=%d \n",rama.phi,rama.psi,map,ntermchain,ctermchain);
    Nphi=(int) (rama.phi+180.0)/space;
    Npsi=(int) (rama.psi+180.0)/space;
    //printf("N %d %d %d\n",Npsi,Nphi,Npsi*7+Nphi);
    if (Nphi>=0 && Nphi<dim && Npsi>=0 && Nphi<dim){
      x1l=Nphi*space-180;;
      x1h=x1l+space;
      x2l=Npsi*space-180;
      x2h=x2l+space;
      //Cterm
      y1=CtermShiftS(map,(Npsi)*dim+(Nphi));
      y2=CtermShiftS(map,(Npsi+1)*dim+(Nphi));
      y3=CtermShiftS(map,(Npsi+1)*dim+(Nphi+1));
      y4=CtermShiftS(map,(Npsi)*dim+(Nphi+1));
      //printf("%f %f %f %f\n",y1,y2,y3,y4);
      u=(rama.phi-x1l)/(x1h-x1l);
      t=(rama.psi-x2l)/(x2h-x2l);
      //printf("%f %f\n",t,u);
      delta=(1-t)*(1-u)*y1+t*(1-u)*y2+t*u*y3+(1-t)*u*y4;
      Hamiltonian[(i)*N+(i)]+=delta;
      //printf("CShift %f\n",delta);
    } else {
      printf("Ill defined Ramachandran angles between unit %d and %d!\n",i,0);
      printf("The nearest neighbor shift will be set to 0.\n");
    }
   }  
  } else { //Ctermshift for all units except the last one of a non-cyclic peptide chain 
    map=cis_or_trans(pbc,snap,top,amideGroups,i+1,cyclic,xG); 
    rama=Ramachandran(pbc,snap,top,amideGroups,i+1,i,xG);
    //printf("Ramachandran angles: Phi=%6.1f, Psi=%6.1f, map=%d, i=%d, j=%d \n",rama.phi,rama.psi,map,i+1,i); 
    Nphi=(int) (rama.phi+180.0)/space;
    Npsi=(int) (rama.psi+180.0)/space;
    //printf("N %d %d %d\n",Npsi,Nphi,Npsi*7+Nphi);
    if (Nphi>=0 && Nphi<dim && Npsi>=0 && Nphi<dim){
      x1l=Nphi*space-180;;
      x1h=x1l+space;
      x2l=Npsi*space-180;
      x2h=x2l+space;
      //Cterm
      y1=CtermShiftS(map,(Npsi)*dim+(Nphi));
      y2=CtermShiftS(map,(Npsi+1)*dim+(Nphi));
      y3=CtermShiftS(map,(Npsi+1)*dim+(Nphi+1));
      y4=CtermShiftS(map,(Npsi)*dim+(Nphi+1));
      //printf("%f %f %f %f\n",y1,y2,y3,y4);
      u=(rama.phi-x1l)/(x1h-x1l);
      t=(rama.psi-x2l)/(x2h-x2l);
      //printf("%f %f\n",t,u);
      delta=(1-t)*(1-u)*y1+t*(1-u)*y2+t*u*y3+(1-t)*u*y4;
      Hamiltonian[(i)*N+(i)]+=delta;
     // printf("CShift %f\n",delta);
    } else {
      printf("Ill defined Ramachandran angles between unit %d and %d!\n",i,i+1);
      printf("The nearest neighbor shift will be set to 0.\n");
    }
  } 
  
  return;
}

/* amide I - amide I  couplings */
void calcCoupling(t_pbc *pbc,t_snap *snap,int cModelTDC,int cModelTCC,int nModelTDC,int nModelTCC,int TDCKrimm,int nModelGLDP,int nModelTasumi,t_topology *top,rvec *xG,float *Hamiltonian,int cyclic,t_amide *amideGroups){
  int i,j,k,a,b,c;
  int N;
  float J;
  int Jcalculated;
  int NN; // Nearest neighbor map dimension
  int CO_per_chain;
  int chainNR,ntermchain,ctermchain;

  J=0;
  CO_per_chain=snap->COchains;
  N=snap->Namidebonds;

  //Loop over peptide unit pairs (j<i)
  for (i=0;i<N;i++){
    //Identify chain number
    chainNR=i/snap->COchains;
    ntermchain=chainNR*snap->COchains;
    ctermchain=(chainNR+1)*snap->COchains-1;

    for (j=0;j<i;j++){
      Jcalculated=0;

      if(j>=snap->NbackboneCO) { //couplings only between the sidechain chromophores
       if(TDCKrimm==1){
            J=calcTDCKrimm(pbc,snap,top,xG,amideGroups,i,j);
            Jcalculated=1;
          } else {
            J=calcTDCTasumi(pbc,snap,top,xG,amideGroups,i,j); //default Torii and Tasumi
            Jcalculated=1;
          }
      } else { //backbone-backbone and backbone-sidechain couplings
        //TDC model
        if (((cModelTDC==1) && (j!=i-1)) || ((j==i-1) && (nModelTDC==1))){
          if(cyclic==1){
           if((j==ntermchain) && (i==ctermchain)){
            J=calcGLDP(pbc,snap,top,xG,amideGroups,cyclic,i,j);
            Jcalculated=1;
           }
          }
          if (!Jcalculated) {
	   if(TDCKrimm==1){
            J=calcTDCKrimm(pbc,snap,top,xG,amideGroups,i,j); //Krimm TDC model
            Jcalculated=1;
           } else {
             J=calcTDCTasumi(pbc,snap,top,xG,amideGroups,i,j);//default Torii and Tasumi
             Jcalculated=1;
           }
          }
        }    
        
        //TCC model
        if (((cModelTCC==1) && (j!=i-1)) || ((j==i-1) && (nModelTCC==1))){
         if(cyclic==1){
           if((j==ntermchain) && (i==ctermchain)){ 
            J=calcGLDP(pbc,snap,top,xG,amideGroups,cyclic,i,j);
            Jcalculated=1; 
            //printf("J %d %d %f cm-1\n",i,j,J);
           } else {
             J=calcTCC2(pbc,snap,top,xG,amideGroups,i,j); //both old and new (for proline) TCC model
             Jcalculated=1;
           }
          } else {
            J=calcTCC2(pbc,snap,top,xG,amideGroups,i,j); //both old and new (for proline) TCC model
  	    Jcalculated=1;
          }
        }

        //Tasumi nnMAP 
        if ((j==i-1) && (nModelTasumi==1)){
          if ((i==ntermchain) && (cModelTDC==1)){  // Use long distance coupling  
           if(TDCKrimm==1){
             J=calcTDCKrimm(pbc,snap,top,xG,amideGroups,i,j);
             Jcalculated=1;
           } else {
             J=calcTDCTasumi(pbc,snap,top,xG,amideGroups,i,j);//default Torii and Tasumi
             Jcalculated=1;
           }   
          } else if ((i==ntermchain) && (cModelTCC==1)){
            J=calcTCC2(pbc,snap,top,xG,amideGroups,i,j);
            Jcalculated=1;
          } else {	
            J=calcTasumi(pbc,snap,top,xG,amideGroups,i,j);
	    Jcalculated=1;
          }
        }

        // GLDP nnMap 
        if ((j==i-1) && (nModelGLDP==1)){
          if ((i==ntermchain) && (cModelTDC==1)){  // Use long distance coupling   
            if(TDCKrimm==1){
              J=calcTDCKrimm(pbc,snap,top,xG,amideGroups,i,j);
              Jcalculated=1;
            } else {
              J=calcTDCTasumi(pbc,snap,top,xG,amideGroups,i,j);//default Torii and Tasumi
              Jcalculated=1;
            }
          } else if ((i==ntermchain) && (cModelTCC==1)){
            J=calcTCC2(pbc,snap,top,xG,amideGroups,i,j);
            Jcalculated=1;
          } else {
            J=calcGLDP(pbc,snap,top,xG,amideGroups,cyclic,i,j);
	    Jcalculated=1;
          }
        }

      } 
      if (Jcalculated==0){
	printf("The coupling between unit %d and %d could not be calculated.\n",i,j);
	printf("Please, specify a valid coupling model!\n");
	printf("Present model %s for nonnearestneighbors\n",snap->cModel);
	printf("and %s for nearestneighbors.\n",snap->nModel);
	exit(-1);
      }

      //      if (J==NAN) J=0; // Discard unphysical couplings
      
      // Update Hamiltonian
      Hamiltonian[i*N+j]=J;
      Hamiltonian[j*N+i]=J;
      //printf("J %d %d %f cm-1\n",i,j,J); 
    }
  }
  return;
}

/* TDC coupling (Krimm dipole) */
float calcTDCKrimm(t_pbc *pbc,t_snap *snap,t_topology *top,rvec *xG,t_amide *amideGroups,int i,int j){
  float J;
  rvec ri,rj,mi,mj,d;
  rvec COi,COj,CNi,CNj;
  float fourPiEpsilon; //219476*0.529177249*0.208194*0.208194; // In cm-1A^3/D^2
  float angle; // Angle with CO bond
  float displace; // Displacement of position
  float magnitude;
  float bond,r,ir3,ir5;
  t_rama rama;

  if (fabs(i-j)<1.5) rama=Ramachandran(pbc,snap,top,amideGroups,i,j,xG);
  // Krimm
  fourPiEpsilon=219476*0.529177249*0.208194*0.208194; // In cm-1A^3/D^2
  angle=20.0/180.0*3.14159265359;
  displace=0.868;
  displace=0.0868; // nm
  magnitude=0.348;
  // Knoester
  fourPiEpsilon=580;
  angle=20.0/180.0*3.14159265359;
  displace=0.868;
  displace=0.0868; // nm
  magnitude=1.0;

  // Locate dipoles
  //  dist(COi,atom[amideGroups[i].a[1]].x,atom[amideGroups[i].a[0]].x);
  pbc_dx(pbc,xG[amideGroups[i].a[1]],xG[amideGroups[i].a[0]],COi);
  bond=sqrt(iprod(COi,COi));
  svmul(displace/bond,COi,COi);
  rvec_add(COi,xG[amideGroups[i].a[0]],ri);
  //  dist(COj,atom[amideGroups[j].a[1]].x,atom[amideGroups[j].a[0]].x);
  pbc_dx(pbc,xG[amideGroups[j].a[1]],xG[amideGroups[j].a[0]],COj);
  bond=sqrt(iprod(COj,COj));
  svmul(displace/bond,COj,COj);
  rvec_add(COj,xG[amideGroups[j].a[0]],rj);
  //  printf("r(%d) %f %f %f  r(%d) %f %f %f\n",i,ri[0],ri[1],ri[2],j,rj[0],rj[1],rj[2]); 
  // Find directions
  unitv(COi,COi),unitv(COj,COj);
  //  dist(CNi,atom[amideGroups[i].a[2]].x,atom[amideGroups[i].a[0]].x);
  pbc_dx(pbc,xG[amideGroups[i].a[2]],xG[amideGroups[i].a[0]],CNi);
  //  dist(CNj,atom[amideGroups[j].a[2]].x,atom[amideGroups[j].a[0]].x);
  pbc_dx(pbc,xG[amideGroups[j].a[2]],xG[amideGroups[j].a[0]],CNj);
  project(COi,CNi),project(COj,CNj);
  //  printf("xy %f %f\n",iprod(COi,CNi),iprod(COj,CNj)); 
  unitv(CNi,CNi),unitv(CNj,CNj);
  bond=tan(angle);
  //  printf("Bond %f %f\n",bond,angle);
  svmul(-bond,CNi,CNi),svmul(-bond,CNj,CNj);
  rvec_add(COi,CNi,mi),rvec_add(COj,CNj,mj);
  unitv(mi,mi),unitv(mj,mj);
  svmul(-magnitude,mi,mi),svmul(-magnitude,mj,mj);
  // Calculate distance
  //  dist(d,rj,ri);
  pbc_dx(pbc,rj,ri,d);
  r=sqrt(iprod(d,d));
  r=r*10;
  d[0]*=10,d[1]*=10,d[2]*=10;
  //  printf("%d %d mu1 %f %f %f mu2 %f %f %f d %f %f %f %f\n",i,j,mi[0],mi[1],mi[2],mj[0],mj[1],mj[2],d[0],d[1],d[2],r);
  ir3=1.0/(r*r*r);
  ir5=ir3/(r*r);
  //  printf("r3 %f r5 %f\n",fourPiEpsilon*(iprod(mi,mj)*ir3),fourPiEpsilon*3.0*iprod(mi,d)*iprod(mj,d)*ir5);
  J=fourPiEpsilon*(iprod(mi,mj)*ir3-3.0*iprod(mi,d)*iprod(mj,d)*ir5);
  return J;
}

/* TDM (Torii and Tasumi dipole) */
void TDM_Torii(t_pbc *pbc,t_snap *snap,t_topology *top,rvec *xG,t_amide *amideGroups,int i,rvec dip){
  rvec mi,dri;
  rvec COi,CNi;
  rvec CO_i;
  rvec nCOi,nCNi;
  rvec dCOi,dCNi;
  int  k;
  float AoverEpsilon; 
  float magnitude;
  float theta; // Angle between the dipole and CO bond 
  float phi,magCO_1,magCO_2,magCO;
  float perpendicular;

  AoverEpsilon=0.1*(848619.0/1650.0); // A/epsilon conversion factor (JPCB 2011, 115, 3713)
  magnitude=2.73; //D A^-1 u^-1/2
  //theta=10.0/180.0*3.14159265359;
  theta=10*((4.0*atan(1.0))/180.0);
  // Locate CO and CN bonds
  pbc_dx(pbc,xG[amideGroups[i].a[1]],xG[amideGroups[i].a[0]],COi);
  pbc_dx(pbc,xG[amideGroups[i].a[2]],xG[amideGroups[i].a[0]],CNi);
  unitv(COi,nCOi); //Unit vector along CO
  unitv(CNi,nCNi); //Unit vector along CN
  svmul(0.665,nCOi,dCOi);//contribution of CO bond to dipole vector d
  svmul(0.258,nCNi,dCNi);//contribution of CN bond to dipole vector d
  rvec_add(dCOi,dCNi,dri);//location of dipole vector dri w.r.t atom C
  //Dipole mi vector
  phi=acos(iprod(nCOi,dri)/(sqrt(iprod(dri,dri))));
  magCO_1=(sqrt(iprod(dri,dri)))*cos(phi);
  perpendicular=(sqrt(iprod(dri,dri)))*sin(phi);
  magCO_2=perpendicular/tan(theta);
  magCO=magCO_1+magCO_2;
  svmul(-magCO,nCOi,CO_i); //reversing the direction of nCO to nOC
  rvec_add(CO_i,dri,mi);
  unitv(mi,mi);
  svmul(magnitude,mi,mi);
  for(k=0; k<3; k++) dip[k]=mi[k]; 
  return;
}

/* TDC coupling (Torii and Tasumi dipole) */
float calcTDCTasumi(t_pbc *pbc,t_snap *snap,t_topology *top,rvec *xG,t_amide *amideGroups,int i,int j){
  float J;
  rvec ri,rj,mi,mj,d,dri,drj;
  rvec COi,COj,CNi,CNj;
  rvec CO_i,CO_j;
  rvec nCOi,nCOj,nCNi,nCNj;
  rvec dCOi,dCOj,dCNi,dCNj;
  float AoverEpsilon; 
  float magnitude;
  float r,ir3,ir5;
  float theta; // Angle between the dipole and CO bond 
  float phi,magCO_1,magCO_2,magCO;
  float perpendicular;

  AoverEpsilon=0.1*(848619.0/1650.0); // A/epsilon conversion factor
  magnitude=2.73; //D A^-1 u^-1/2
  //theta=10.0/180.0*3.14159265359;
  theta=10*((4.0*atan(1.0))/180.0);
  // Locate CO and CN bonds
  pbc_dx(pbc,xG[amideGroups[i].a[1]],xG[amideGroups[i].a[0]],COi);
  pbc_dx(pbc,xG[amideGroups[j].a[1]],xG[amideGroups[j].a[0]],COj);
  pbc_dx(pbc,xG[amideGroups[i].a[2]],xG[amideGroups[i].a[0]],CNi);
  pbc_dx(pbc,xG[amideGroups[j].a[2]],xG[amideGroups[j].a[0]],CNj);
   
  unitv(COi,nCOi),unitv(COj,nCOj); //Unit vector along CO
  unitv(CNi,nCNi),unitv(CNj,nCNj); //Unit vector along CN
  svmul(0.0665,nCOi,dCOi);//contribution of CO bond to dipole vector d
  svmul(0.0258,nCNi,dCNi);//contribution of CN bond to dipole vector d
  svmul(0.0665,nCOj,dCOj);
  svmul(0.0258,nCNj,dCNj);
  rvec_add(dCOi,dCNi,dri);//location of dipole vector dri w.r.t atom C 
  rvec_add(xG[amideGroups[i].a[0]],dri,ri);//location of dipole vector in global coordinates
  rvec_add(dCOj,dCNj,drj);//location of dipole vector drj w.r.t atom C
  rvec_add(xG[amideGroups[j].a[0]],drj,rj);//location of dipole vector in global coordinates
  pbc_dx(pbc,rj,ri,d);//vector d that connects between ri and rj
   
  r=sqrt(iprod(d,d));//distance between the dipoles ri and rj 
  r=r*10;
  d[0]*=10,d[1]*=10,d[2]*=10;  
  ir3=1.0/(r*r*r);
  ir5=ir3/(r*r);
  //Dipole mi vector
  phi=acos(iprod(nCOi,dri)/(sqrt(iprod(dri,dri))));
  magCO_1=(sqrt(iprod(dri,dri)))*cos(phi);
  perpendicular=(sqrt(iprod(dri,dri)))*sin(phi);  
  magCO_2=perpendicular/tan(theta);
  magCO=magCO_1+magCO_2;
  svmul(-magCO,nCOi,CO_i); //reversing the direction of nCO to nOC
  rvec_add(CO_i,dri,mi); 
  unitv(mi,mi);
  svmul(magnitude,mi,mi);
  //Dipole mj vector
  phi=acos(iprod(nCOj,drj)/(sqrt(iprod(drj,drj))));
  magCO_1=(sqrt(iprod(drj,drj)))*cos(phi);
  perpendicular=(sqrt(iprod(drj,drj)))*sin(phi);  
  magCO_2=perpendicular/tan(theta);
  magCO=magCO_1+magCO_2;
  svmul(-magCO,nCOj,CO_j);//reversing the direction of nCO to nOC
  rvec_add(CO_j,drj,mj);
  unitv(mj,mj);
  svmul(magnitude,mj,mj);

  J=AoverEpsilon*(iprod(mi,mj)*ir3-3.0*iprod(mi,d)*iprod(mj,d)*ir5);
  return J;
}

// TCC coupling (Mulliken and MDC charges for proline and other units, respectively) 
float calcTCC2(t_pbc *pbc,t_snap *snap,t_topology *top,rvec *xG,t_amide *amideGroups,int i,int j){
  float JJ;
  int a,b,c;
  float q[6][2],dq[6][2];
  rvec ri,rj,d;
  rvec COi,COj,CNi,CNj,zi,zj;
  float v[3][6][2];
  rvec va,vb;
  float alpha[2];
  //  float fourPiEpsilon=219476.0*0.529177249; // In cm-1A^3
  float fourPiEpsilon=219476.0*0.529177249; // In cm-1 nm^3
  float r,ir,ir3,ir5;
  float bond;
  t_rama rama;
  int typei,typej;
  typei=typej=0;

  if(fabs(i-j)<1.5) rama=Ramachandran(pbc,snap,top,amideGroups,i,j,xG);
  //Check whether i/j th unit is proline
  if (amideGroups[i].proline==TRUE) typei=1;
  if (amideGroups[j].proline==TRUE) typej=1;
  // TCC parameters
  alpha[0]=0.028074,alpha[1]=0.0280896;
  q[0][0]=0.37173,q[1][0]=-0.53632,q[2][0]=-0.48418,q[3][0]=0.24278,q[4][0]=0.29527,q[5][0]=0.11072;//other units
  q[0][1]=0.159744,q[1][1]=-0.347870,q[2][1]=0.024241,q[3][1]=0.061036,q[4][1]=0.086851,q[5][1]=0.015998;//proline
  // Alpha carbons exchanged 26/3-2008 
  dq[0][0]=0.02845,dq[1][0]=0.01530,dq[2][0]=-0.01736,dq[3][0]=-0.00008,dq[4][0]=-0.00963,dq[5][0]=-0.01668;//other units
  dq[0][1]=0.010457,dq[1][1]=-0.044243,dq[2][1]=0.013641,dq[3][1]=0.003277,dq[4][1]=0.005295,dq[5][1]=0.011574;//proline
  //normal mode coordinates 
  v[0][0][0]=-0.831,v[1][0][0]=0.105,v[2][0][0]=0; // C
  v[0][1][0]=0.517,v[1][1][0]=-0.047,v[2][1][0]=0; // O
  v[0][2][0]=0.074,v[1][2][0]=-0.036,v[2][2][0]=0; // N
  v[0][3][0]=0.073,v[1][3][0]=-0.133,v[2][3][0]=0; // D
  v[0][4][0]=0,v[1][4][0]=0,v[2][4][0]=0; // C(a)
  v[0][5][0]=0,v[1][5][0]=0,v[2][5][0]=0; // C(a)-1
  v[0][0][1]=-0.836,v[1][0][1]=-0.004,v[2][0][1]=0; // C
  v[0][1][1]=0.537,v[1][1][1]=0.023,v[2][1][1]=0; // O
  v[0][2][1]=0.099,v[1][2][1]=-0.018,v[2][2][1]=0; // N
  v[0][3][1]=0,v[1][3][1]=0.0,v[2][3][1]=0; // D
  v[0][4][1]=0,v[1][4][1]=0,v[2][4][1]=0; // C(a)
  v[0][5][1]=0,v[1][5][1]=0,v[2][5][1]=0; // C(a)-1
   
  // Define molecular coordinate systems
  //  dist(COi,atom[amideGroups[i].a[1]].x,atom[amideGroups[i].a[0]].x);
  pbc_dx(pbc,xG[amideGroups[i].a[1]],xG[amideGroups[i].a[0]],COi);
  //  dist(COj,atom[amideGroups[j].a[1]].x,atom[amideGroups[j].a[0]].x);
  pbc_dx(pbc,xG[amideGroups[j].a[1]],xG[amideGroups[j].a[0]],COj);
  unitv(COi,COi),unitv(COj,COj);
  //  dist(CNi,atom[amideGroups[i].a[2]].x,atom[amideGroups[i].a[0]].x);
  pbc_dx(pbc,xG[amideGroups[i].a[2]],xG[amideGroups[i].a[0]],CNi);
  //  dist(CNj,atom[amideGroups[j].a[2]].x,atom[amideGroups[j].a[0]].x);
  pbc_dx(pbc,xG[amideGroups[j].a[2]],xG[amideGroups[j].a[0]],CNj);
  project(COi,CNi),project(COj,CNj);
  unitv(CNi,CNi),unitv(CNj,CNj);
  cprod(COi,CNi,zi),cprod(COj,CNj,zj);
  // CO - x, CN - y, z - z
  //  printf("%f %f %f  %f %f %f  %f %f %f\n",COi[0],COi[1],COi[2],CNi[0],CNi[1],CNi[2],zi[0],zi[1],zi[2]);

  JJ=0;
  // Loop over atoms
  for (a=0;a<6;a++){
    for (b=0;b<6;b++){
      // Find normal mode vectors in global coordinate system
      for (c=0;c<3;c++) va[c]=(COi[c]*v[0][a][typei]+CNi[c]*v[1][a][typei]+zi[c]*v[2][a][typei])*alpha[typei];
      for (c=0;c<3;c++) vb[c]=(COj[c]*v[0][b][typej]+CNj[c]*v[1][b][typej]+zj[c]*v[2][b][typej])*alpha[typej];
      // Find distance vector
      // dist(d,atom[amideGroups[i].a[a]].x,atom[amideGroups[j].a[b]].x);
      pbc_dx(pbc,xG[amideGroups[i].a[a]],xG[amideGroups[j].a[b]],d);
      r=sqrt(iprod(d,d));
      r=r*10;
      d[0]*=10,d[1]*=10,d[2]*=10;
      ir=1.0/r;
      ir3=ir*ir*ir;
      ir5=ir3*ir*ir;
      
      if (r<0.1) ir=1,ir3=1,ir5=1;
         
      JJ-=3.0*ir5*q[a][typei]*q[b][typej]*(iprod(vb,d)*iprod(va,d));
      //      if (i==11 && j==0) printf("r5 %d %d %d %d %f\n",i,j,a,b,JJ*fourPiEpsilon);
      JJ-=ir3*(dq[a][typei]*q[b][typej]*iprod(vb,d)-q[a][typei]*dq[b][typej]*iprod(va,d)-iprod(va,vb)*q[a][typei]*q[b][typej]);
      //      if (i==11 && j==0) printf("r3 %d %d %d %d %f\n",i,j,a,b,JJ*fourPiEpsilon);
      JJ+=ir*dq[a][typei]*dq[b][typej];
      //      if (i==11 && j==0) printf("r1 %d %d %d %d %f\n",i,j,a,b,JJ*fourPiEpsilon);     
     
    }
  }
  JJ*=fourPiEpsilon;

  return JJ;
}
// Calculate Ramachandran angles and extract coupling from Tasumi map
float calcTasumi(t_pbc *pbc,t_snap *snap,t_topology *top,rvec *xG,t_amide *amideGroups,int i,int j){
  float JJ;
  float phi,psi;
  float phic,phis;
  float psic,psis;
  //  float r1[3],r2[3],r3[3],r4[3];
  //  float n1[3],n2[3];
  int ri,rj;
  int Nphi,Npsi;
  float t,u,y1,y2,y3,y4,x1l,x1h,x2l,x2h;
  int tem;
  t_rama rama;

  JJ=0;
  rama=Ramachandran(pbc,snap,top,amideGroups,i,j,xG);
  phi=rama.phi,psi=rama.psi;
  //  printf("%f\n",Tasumi[4]);
  Nphi=(int) (phi+180.0)/30.0;
  Npsi=(int) (psi+180.0)/30.0;
  if (Nphi>=0 && Nphi<13 && Npsi>=0 && Nphi<13){
    x1l=Nphi*30.0-180;;
    x1h=x1l+30.0;
    x2l=Npsi*30.0-180;
    x2h=x2l+30.0;
    y1=Tasumi[(Nphi)*13+(Npsi)];
    y2=Tasumi[(Nphi+1)*13+(Npsi)];
    y3=Tasumi[(Nphi+1)*13+(Npsi+1)];
    y4=Tasumi[(Nphi)*13+(Npsi+1)];
    //    printf("%f %f %f %f\n",y1,y2,y3,y4);
    t=(phi-x1l)/(x1h-x1l);
    u=(psi-x2l)/(x2h-x2l);
    //    printf("%f %f\n",t,u);
    JJ=(1-t)*(1-u)*y1+t*(1-u)*y2+t*u*y3+(1-t)*u*y4;
    //    printf("%d %d\n",Nphi,Npsi);
  } else {
    printf("Ill defined Ramachandran angles between unit %d and %d!\n",i,j);
    printf("The nearest neighbor coupling will be set to 0.\n");
  }
  
  return JJ;
}

// Calculate Ramachandran angles and extract coupling from GLDP map
float calcGLDP(t_pbc *pbc,t_snap *snap,t_topology *top,rvec *xG,t_amide *amideGroups,int cyclic,int i,int j){
  float JJ;
  float phi,psi;
  float phic,phis;
  float psic,psis;
  int ri,rj;
  int Nphi,Npsi;
  float t,u,y1,y2,y3,y4,x1l,x1h,x2l,x2h;
  int tem;
  int dim;
  float space;
  int N;
  int map; 
  int chainNR,ntermchain,ctermchain; 
  t_rama rama;

  JJ=0;
  N=snap->Namidebonds;
  // Identify chain number
  chainNR=j/snap->COchains;
  ntermchain=chainNR*snap->COchains;
  ctermchain=(chainNR+1)*snap->COchains-1; 
 
  dim=13;
  space=30;
  if((i==ctermchain) && (j==ntermchain)  && (cyclic==1)){
    //Checking which NNC map should be used
    map=cis_or_trans(pbc,snap,top,amideGroups,j,cyclic,xG);
    //Find Ramachandran angles
    rama=Ramachandran(pbc,snap,top,amideGroups,j,i,xG);
    //printf("Ramachandran angles: Phi=%6.1f, Psi=%6.1f, map=%d, i=%d, j=%d \n",rama.phi,rama.psi,map,i,j); 
  } else {
   //Checking which NNC map should be used
   map=cis_or_trans(pbc,snap,top,amideGroups,i,cyclic,xG);
   //Find Ramachandran angles
   rama=Ramachandran(pbc,snap,top,amideGroups,i,j,xG);
   //printf("Ramachandran angles: Phi=%6.1f, Psi=%6.1f, map=%d, i=%d, j=%d \n",rama.phi,rama.psi,map,i,j);
  }
  Nphi=(int) (rama.phi+180.0)/space;
  Npsi=(int) (rama.psi+180.0)/space;
  //printf("N %d %d %d\n",Npsi,Nphi,Npsi*7+Nphi);
  if (Nphi>=0 && Nphi<dim && Npsi>=0 && Nphi<dim){
    x1l=Nphi*space-180;;
    x1h=x1l+space;
    x2l=Npsi*space-180;
    x2h=x2l+space;
    //Cterm
    y1=CouplingS(map,(Npsi)*dim+(Nphi));
    y2=CouplingS(map,(Npsi+1)*dim+(Nphi));
    y3=CouplingS(map,(Npsi+1)*dim+(Nphi+1));
    y4=CouplingS(map,(Npsi)*dim+(Nphi+1));
    //printf("%f %f %f %f\n",y1,y2,y3,y4);
    u=(rama.phi-x1l)/(x1h-x1l);
    t=(rama.psi-x2l)/(x2h-x2l);
    //printf("%f %f\n",t,u);
    JJ=(1-t)*(1-u)*y1+t*(1-u)*y2+t*u*y3+(1-t)*u*y4;
    //printf("CShift %f\n",delta);
  } else {
    printf("Ill defined Ramachandran angles between unit %d and %d!\n",i,i-1);
    printf("The nearest neighbor coupling will be set to 0.\n");
  }
  
  return JJ;
}
/////////////////////////////////////////////////////////////
// Extract couplings from different nearest neaighbor maps//
///////////////////////////////////////////////////////////
float CouplingS(int map,int index){
  float value;
  value=-1000.0;
  if (map==0){
    value=Coupling[index];
  } 
  if (map==1) {
     value=Coupling_transPro_transGly[index];
  } else if (map==2) {
     value=Coupling_cisPro_transGly[index];
  } else if (map==3) {
    value=Coupling_transGly_transPro[index];
  } else if (map==4) {
    value=Coupling_cisGly_transPro[index];
  } else if (map==5) {
    value=Coupling_transDPro_transGly[index]; 
  } else if (map==6) {
    value=Coupling_cisDPro_transGly[index];   
  }
  if (value == -1000.0) {
    printf("failed to extract NN couplings!\n");
    exit(-1);
  }
  return value;
}
////////////////////////////////////
//Extract shifts due to N-terminal/
//////////////////////////////////
float NtermShiftS(int map,int index){
  float value;
  value = -1000.0;
  if (map==0){
    value=NtermShift[index];
  }
  if (map==1) {
    value=NtermShift_transPro_transGly[index];
  } else if (map==2) {
    value=NtermShift_cisPro_transGly[index];
  } else if (map==3) {
    value=NtermShift_transGly_transPro[index];
  } else if (map==4) {
    value=NtermShift_cisGly_transPro[index];
  } else if (map==5) {
    value=NtermShift_transDPro_transGly[index];
  } else if (map==6) {
    value=NtermShift_cisDPro_transGly[index];
  }
  if (value == -1000.0) {
    printf("failed to extract Nterminal shift!\n");
    exit(-1);
  } 
  return value;
}

/* Extract shifts due to C-terminal */

float CtermShiftS(int map,int index){
  float value;
  value=-1000.0;
  if (map==0){
    value=CtermShift[index];
  } 
  if (map==1) {
    value=CtermShift_transPro_transGly[index];
  } else if (map==2) {
    value=CtermShift_cisPro_transGly[index];
  } else if (map==3) {
    value=CtermShift_transGly_transPro[index];
  } else if (map==4) {
    value=CtermShift_cisGly_transPro[index];
  } else if (map==5) {
    value=CtermShift_transDPro_transGly[index];
  } else if (map==6) {
    value=CtermShift_cisDPro_transGly[index];
  }
  if (value == -1000.0) {
    printf("failed to extract Cterminal shift!\n");
    exit(-1);
  }
  return value;
}

/* determining whether the configuration is cis- or trans- */
 
/* Cis or trans */
int cis_or_trans(t_pbc *pbc,t_snap *snap,t_topology *top,t_amide *amideGroups,int i,int cyclic,rvec *xG){
  int cis_trans;
  int N;
  float r1[3],r2[3];
  int chainNR,ntermchain,ctermchain;

  // Identify chain number
  chainNR=i/snap->COchains;
  ntermchain=chainNR*snap->COchains;
  ctermchain=(chainNR+1)*snap->COchains-1;

  N=snap->Namidebonds;
  cis_trans=0;
  //For the first and the last terminal units of a cyclic peptide
  if(cyclic==1){
    if (i==ntermchain) {
     //Pro-Gly
     if (amideGroups[ctermchain].proline==TRUE){
       pbc_dx(pbc,xG[amideGroups[ctermchain].a[1]],xG[amideGroups[ctermchain].a[0]],r1);
       pbc_dx(pbc,xG[amideGroups[ctermchain].a[3]],xG[amideGroups[ctermchain].a[2]],r2);
       if (iprod(r1,r2)<0){
         if(!strcmp(snap->ProType,"DProline")) {
           cis_trans=5; // (Trans-Trans)
         } else {
           cis_trans=1; // (Trans-Trans)
         }   
       } else {
         if(!strcmp(snap->ProType,"DProline")){
           cis_trans=6; // (Cis-Trans) 
         } else { 
           cis_trans=2; // (Cis-Trans)
         }
       }
     }
     //Gly-Pro
     if (amideGroups[ntermchain].proline==TRUE){
       pbc_dx(pbc,xG[amideGroups[ctermchain].a[1]],xG[amideGroups[ctermchain].a[0]],r1);
       pbc_dx(pbc,xG[amideGroups[ctermchain].a[3]],xG[amideGroups[ctermchain].a[2]],r2);
       if (iprod(r1,r2)<0){
        cis_trans=3; // (Trans-Trans)
       } else {
        cis_trans=4; // (Cis-Trans)
       }
     } 
    }
  }

  //For all units of a (non)-cyclic peptide   
  if (i>ntermchain){
     //Pro-Gly 
     if (amideGroups[i-1].proline==TRUE){
       pbc_dx(pbc,xG[amideGroups[i-1].a[1]],xG[amideGroups[i-1].a[0]],r1);
       pbc_dx(pbc,xG[amideGroups[i-1].a[3]],xG[amideGroups[i-1].a[2]],r2);
        if (iprod(r1,r2)<0){
         if(!strcmp(snap->ProType,"DProline")) {
           cis_trans=5; // (Trans-Trans)
         } else {
           cis_trans=1; // (Trans-Trans)
         }   
       } else {
         if(!strcmp(snap->ProType,"DProline")){
           cis_trans=6; // (Cis-Trans) 
         } else { 
           cis_trans=2; // (Cis-Trans)
         }
       }
     }
     //Gly-Pro
     if (amideGroups[i].proline==TRUE){
       pbc_dx(pbc,xG[amideGroups[i-1].a[1]],xG[amideGroups[i-1].a[0]],r1);
       pbc_dx(pbc,xG[amideGroups[i-1].a[3]],xG[amideGroups[i-1].a[2]],r2);
       if (iprod(r1,r2)<0){
         cis_trans=3; // (Trans-Trans)
       } else {
        cis_trans=4; // (Cis-Trans)
       }
     }
  }
  return cis_trans;
}

/* Calculate Ramachandran angles */
t_rama Ramachandran(t_pbc *pbc,t_snap *snap,t_topology *top,t_amide *amideGroups,int i,int j,rvec *xG){
  float phi,psi;
  float phic,phis;
  float psic,psis;
  rvec r1,r2,r3,r4;
  rvec n1,n2;
  int ri,rj;
  int Nphi,Npsi;
  float t,u,y1,y2,y3,y4,x1l,x1h,x2l,x2h;
  int tem;
  t_rama rama;

  float phi_g_rama, psi_g_rama; 
  int t1,t2,t3;
  real sign,cos_phi;
 
  ri=i,rj=j,tem=4;
  if (!strcmp("CtoN",snap->order)) ri=j,rj=i,tem=5; // Reverse order

  // Method of g_rama
  phi_g_rama=180/3.1416*dih_angle(xG[amideGroups[rj].a[0]],
                       xG[amideGroups[rj].a[2]],
                       xG[amideGroups[rj].a[tem]],
                       xG[amideGroups[ri].a[0]],
                       pbc,r1,r2,r3,n1,n2,&sign,&t1,&t2,&t3);


  psi_g_rama=180/3.1416*dih_angle(xG[amideGroups[rj].a[2]],
                       xG[amideGroups[rj].a[tem]],
                       xG[amideGroups[ri].a[0]],
                       xG[amideGroups[ri].a[2]],
                       pbc,r1,r2,r3,n1,n2,&sign,&t1,&t2,&t3);

//  printf("             g_rama method: Phi=%6.1f, Psi=%6.1f\n",phi_g_rama,psi_g_rama);

  rama.phi=phi_g_rama,rama.psi=psi_g_rama;
  return rama;
}

/* Read the overall input file */
void readInput(char inputFName[256],t_snap *snap){
  FILE *inputFile;
  char *pStatus;
  char Buffer[256];
  size_t LabelLength;
  char *pValue;
  int control,i;

  printf("Using input file '%s'.\n",inputFName);  

  // Open input file
  inputFile=fopen(inputFName, "r");
  if (inputFile == NULL) {
    printf("File not found!\n");
    exit(-1);
  }

  control=0;

  // Defaults
  snap->nterm=1,snap->cterm=1;

  // Read input data
  do {
    pStatus = fgets(&Buffer[0],sizeof(Buffer),inputFile);
    if (pStatus == NULL) {
      break;
    }
    
    // Compute LabelLength
    LabelLength = strcspn(&Buffer[0], " ");

    // Forcefield Keyword
    if (keyWordS("Forcefield",Buffer,snap->forcefield,LabelLength)==1) continue; 

    // Peptide type file keyword if Cyclic
    if (keyWordS("Peptidetype",Buffer,snap->peptide,LabelLength)==1) continue;

    // Proline type keyword: if DProline
    if (keyWordS("Prolinetype",Buffer,snap->ProType,LabelLength)==1) continue; 

    // Atoms keyword
    //if (keyWordI("Atoms",Buffer,&(snap->Natoms),LabelLength)==1) continue;

    // Types keyword
    //if (keyWordI("Types",Buffer,&(snap->Ntypes),LabelLength)==1) continue;

    // Residues keyword
    if (keyWordI("Residues",Buffer,&(snap->Nresidues),LabelLength)==1) continue;

    // Sidechain chromophores keyword
    if (keyWordI("BackboneAmidebonds",Buffer,&(snap->NbackboneCO),LabelLength)==1) continue;

    // Sidechain chromophores keyword
    if (keyWordI("SidechainAmidebonds",Buffer,&(snap->NsidechainCO),LabelLength)==1) continue;

    // Amide bonds
    if (keyWordI("TotalAmidebonds",Buffer,&(snap->Namidebonds),LabelLength)==1) continue;
    
    // How many CO group per chain
    if (keyWordI("COperChain",Buffer,&(snap->COchains),LabelLength)==1) continue;
   
    // Terminals
    //if (keyWordI("Nterm",Buffer,&(snap->nterm),LabelLength)==1) continue;
    //if (keyWordI("Cterm",Buffer,&(snap->cterm),LabelLength)==1) continue;

    // Cutoff keyword
    if (keyWordF("Cutoff",Buffer,&(snap->cutoff),LabelLength)==1) continue;

    // Energy output file
    if (keyWordS("Hamiltonianfile",Buffer,snap->energyFName,LabelLength)==1) continue;
    
    // Dipole output file
    if (keyWordS("Dipolefile",Buffer,snap->dipoleFName,LabelLength)==1) continue;

    // Format keyword
    if (keyWordS("Format",Buffer,snap->format,LabelLength)==1) continue;

    // Order keyword
    //if (keyWordS("Order",Buffer,snap->order,LabelLength)==1) continue;
    
    // Jansen Map file keyword
    if (keyWordS("EstaticMap",Buffer,snap->ElectrostaticMap,LabelLength)==1) continue;

    // Jansen Map file keyword
    if (keyWordS("JansenMapfile",Buffer,snap->mapJansenFName,LabelLength)==1) continue;
    
    // Skinner Map file keyword
    if (keyWordS("SkinnerMapfile",Buffer,snap->mapSkinnerFName,LabelLength)==1) continue;

    // Tokmakoff Map file keyword
    if (keyWordS("TokmakoffMapfile",Buffer,snap->mapTokmakoffFName,LabelLength)==1) continue;

    // Roy Map file keyword
    if (keyWordS("RoyMapfile",Buffer,snap->mapRoyFName,LabelLength)==1) continue;
    
    // Coupling file keyword
    //if (keyWordS("Couplingfile",Buffer,snap->couplingFName,LabelLength)==1) continue;

     // Couplingmodel keyword
    if (keyWordS("LongrangeCoupling",Buffer,snap->cModel,LabelLength)==1) continue;

    // Neighbormodel keyword
    if (keyWordS("NearestNeighborCoupling",Buffer,snap->nModel,LabelLength)==1) continue; 

    // Transition Dipole type keyword
    if (keyWordS("TransitionDipole",Buffer,snap->TD,LabelLength)==1) continue;

    // TDC type keyword
    if (keyWordS("TDCtype",Buffer,snap->TDCtypeName,LabelLength)==1) continue;

    // Nearest neighbour treatment
    if (keyWordS("NNTreatment",Buffer,snap->nnTreat,LabelLength)==1) continue;   
 
    // Include only part of system in frequency calculation
    // Can be All, Water, Backbone, Sidechains, Protein or Lipid 
    if (keyWordS("Select",Buffer,snap->subSelect,LabelLength)==1) continue;

    // Skip Doubles
    //if (keyWordS("Skip",Buffer,snap->skipDoubles,LabelLength)==1) continue;

    // Proton ID keyword
    //if (keyWord2I("ProtonID",Buffer,&(snap->protonId1),&(snap->protonId2),LabelLength)==1) continue;

    // Proton Charge keyword
    //if (keyWordF("OCharge",Buffer,&(snap->oCharge),LabelLength)==1) continue;
    //if (keyWordF("HCharge",Buffer,&(snap->hCharge),LabelLength)==1) continue;

  } while (TRUE);
  fclose(inputFile);

  return;
}

/* Read the Jansen frequency map */
void readMap(t_snap *snap,t_map *map)
{
  FILE *inputFile;
  int i,j,k;
  float omega2;
  float dummy;

  if (!strcmp(snap->mapJansenFName,"Hirst")){
    // Use Hirst map
    map->p[0]=-0.00999,map->p[1]=-0.00251,map->p[2]=-0.05964,map->p[3]=-0.02414,map->p[4]=-0.00051,map->p[5]=0.02648,
      map->p[6]=0.07033;
    for (i=0;i<7;i++) map->p[i]*=0.529177*219476;
    map->omega=1717;
    map->mu[0]=0.939692620786*0.348; // 0.348 Debye = Krimm dipole
    map->mu[1]=-0.342020143326*0.348;
    map->mu[2]=0;
    map->maptype=7;
  } else {
    // Read from generic file
    inputFile=fopen(snap->mapJansenFName,"r");
    if (inputFile == NULL) {
      printf("Jansen Map file not found!\n");
      exit(-1);
    }
    
    fscanf(inputFile,"%f %f %f",&map->omega,&map->k,&omega2);
    for (i=0;i<4;i++){
      fscanf(inputFile,"%f %f %f %f %f %f %f %f %f %f",&map->p[i],&map->cE[0][i],&map->cE[1][i],&map->cE[2][i],&map->cG[0][i],&map->cG[1][i],&map->cG[2][i],&map->cG[3][i],&map->cG[4][i],&map->cG[5][i]);
      if (omega2==0){
	// Rewrite to form Dw=<c|E> from Dw=w*k*<c'|E>
	for (j=0;j<3;j++) map->cE[j][i]*=map->omega*map->k;
	for (j=0;j<6;j++) map->cG[j][i]*=map->omega*map->k;
	map->p[i]*=map->omega*map->k;
      } else {
	map->omega=omega2;
      }
      
    }
    
    // Overtone frequency (skip)
    fscanf(inputFile,"%f %f %f",&dummy,&dummy,&dummy);
    for (i=0;i<4;i++){
      fscanf(inputFile,"%f %f %f %f %f %f %f %f %f %f",&dummy,&dummy,&dummy,&dummy,&dummy,&dummy,&dummy,&dummy,&dummy,&dummy);   
    }
    
    // Dipole
    for (k=0;k<3;k++){
      fscanf(inputFile,"%f %f %f",&map->mu[k],&map->muk[k],&omega2);
      for (i=0;i<4;i++){
	fscanf(inputFile,"%f %f %f %f %f %f %f %f %f %f",&map->mup[k][i],&map->mucE[k][0][i],&map->mucE[k][1][i],&map->mucE[k][2][i],&map->mucG[k][0][i],&map->mucG[k][1][i],&map->mucG[k][2][i],&map->mucG[k][3][i],&map->mucG[k][4][i],&map->mucG[k][5][i]);
	if (omega2==0){
	  // Rewrite to form Dw=<c|E> from Dw=w*k*<c'|E>
	  for (j=0;j<3;j++) map->mucE[k][j][i]*=map->mu[k]*map->muk[k];
	  for (j=0;j<6;j++) map->mucG[k][j][i]*=map->mu[k]*map->muk[k];
	  map->mup[k][i]*=map->mu[k]*map->muk[k];
	} else {
	  map->mu[k]=omega2;
	}       
      } 
    }
    map->maptype=4; // 4 atom map (so far only type supported)
    fclose(inputFile);
  }
  printf("Jansen Electrostatic map read!\n");
  return;
}

/* Read proline eletrostatic map file */
void readProMap(t_snap *snap,t_map *Promap)
{
  FILE *inputFile;
  int i,j,k;
  float omega2;
  float dummy;

  // Read from generic file
    inputFile=fopen(snap->mapRoyFName,"r");
    if (inputFile == NULL) {
      printf("Roy Electrostatic Map file not found!\n");
      exit(-1);
    }
   
  fscanf(inputFile,"%f %f %f",&Promap->omega,&Promap->k,&omega2);
  for (i=0;i<4;i++){
    fscanf(inputFile,"%f %f %f %f %f %f %f %f %f %f",&Promap->p[i],&Promap->cE[0][i],&Promap->cE[1][i],&Promap->cE[2][i],&Promap->cG[0][i],&Promap->cG[1][i],&Promap->cG[2][i],&Promap->cG[3][i],&Promap->cG[4][i],&Promap->cG[5][i]);
    if (omega2==0){
	// Rewrite to form Dw=<c|E> from Dw=w*k*<c'|E>
      for (j=0;j<3;j++) Promap->cE[j][i]*=Promap->omega*Promap->k;
	for (j=0;j<6;j++) Promap->cG[j][i]*=Promap->omega*Promap->k;
	Promap->p[i]*=Promap->omega*Promap->k;
    } else {
	Promap->omega=omega2;
    }
      
  }
    
  // Overtone frequency (skip)
  fscanf(inputFile,"%f %f %f",&dummy,&dummy,&dummy);
  for (i=0;i<4;i++){
    fscanf(inputFile,"%f %f %f %f %f %f %f %f %f %f",&dummy,&dummy,&dummy,&dummy,&dummy,&dummy,&dummy,&dummy,&dummy,&dummy);   
  }
    
  // Dipole
  for (k=0;k<3;k++){
    fscanf(inputFile,"%f %f %f",&Promap->mu[k],&Promap->muk[k],&omega2);
    for (i=0;i<4;i++){
      fscanf(inputFile,"%f %f %f %f %f %f %f %f %f %f",&Promap->mup[k][i],&Promap->mucE[k][0][i],&Promap->mucE[k][1][i],&Promap->mucE[k][2][i],&Promap->mucG[k][0][i],&Promap->mucG[k][1][i],&Promap->mucG[k][2][i],&Promap->mucG[k][3][i],&Promap->mucG[k][4][i],&Promap->mucG[k][5][i]);
      if (omega2==0){
	  // Rewrite to form Dw=<c|E> from Dw=w*k*<c'|E>
        for (j=0;j<3;j++) Promap->mucE[k][j][i]*=Promap->mu[k]*Promap->muk[k];
        for (j=0;j<6;j++) Promap->mucG[k][j][i]*=Promap->mu[k]*Promap->muk[k];
        Promap->mup[k][i]*=Promap->mu[k]*Promap->muk[k];
      } else {
        Promap->mu[k]=omega2;
      }       
    } 
  }
 Promap->maptype=4; // 4 atom map (so far only type supported)
 fclose(inputFile);
 printf("Proline electrostatic map read!\n");
 return;
}

/* Read Skinner eletrostatic map file */
void readSkinnerAmI(t_snap *snap,t_Skinner *SkinnerAmI){
  FILE *inputFile;
  int i,j;
  float omega2;
  float omega2_S;
  
  // Read from generic file
  inputFile=fopen(snap->mapSkinnerFName,"r");
  if (inputFile == NULL) {
     printf("Skinner Electrostatic map file not found!\n");
     exit(-1);
  }
  fscanf(inputFile,"%f %f %f",&SkinnerAmI->omega,&SkinnerAmI->k,&omega2);
  for (i=0;i<4;i++){
    fscanf(inputFile,"%f %f %f %f %f %f %f %f %f %f",&SkinnerAmI->p[i],&SkinnerAmI->cE[0][i],&SkinnerAmI->cE[1][i],&SkinnerAmI->cE[2][i],&SkinnerAmI->cG[0][i],&SkinnerAmI->cG[1][i],&SkinnerAmI->cG[2][i],&SkinnerAmI->cG[3][i],&SkinnerAmI->cG[4][i],&SkinnerAmI->cG[5][i]);
  }
  fscanf(inputFile,"%f %f %f",&SkinnerAmI->omega_S,&SkinnerAmI->k_S,&omega2_S); // S for sidechain
  for (i=4;i<8;i++){
    fscanf(inputFile,"%f %f %f %f %f %f %f %f %f %f",&SkinnerAmI->p[i],&SkinnerAmI->cE[0][i],&SkinnerAmI->cE[1][i],&SkinnerAmI->cE[2][i],&SkinnerAmI->cG[0][i],&SkinnerAmI->cG[1][i],&SkinnerAmI->cG[2][i],&SkinnerAmI->cG[3][i],&SkinnerAmI->cG[4][i],&SkinnerAmI->cG[5][i]);
  }
   
  fclose(inputFile);
  printf("Skinner electrostatic map read!\n");
  return;
}

/* Read Tokmakoff eletrostatic map file */
void readTokmakoffAmI(t_snap *snap,t_Tokmakoff *TokmakoffAmI){
  FILE *inputFile;
  int i,j;
  float omega2;
    
  // Read from generic file
  inputFile=fopen(snap->mapTokmakoffFName,"r");
  if (inputFile == NULL) {
     printf("Tokmakoff Electrostatic map file not found!\n");
     exit(-1);
  }
  fscanf(inputFile,"%f %f %f",&TokmakoffAmI->omega,&TokmakoffAmI->k,&omega2);
  for (i=0;i<4;i++){
    fscanf(inputFile,"%f %f %f %f %f %f %f %f %f %f",&TokmakoffAmI->p[i],&TokmakoffAmI->cE[0][i],&TokmakoffAmI->cE[1][i],&TokmakoffAmI->cE[2][i],&TokmakoffAmI->cG[0][i],&TokmakoffAmI->cG[1][i],&TokmakoffAmI->cG[2][i],&TokmakoffAmI->cG[3][i],&TokmakoffAmI->cG[4][i],&TokmakoffAmI->cG[5][i]);
  }
  fscanf(inputFile,"%f %f %f",&TokmakoffAmI->omegaPro,&TokmakoffAmI->k,&omega2);
  for (i=0;i<4;i++){
    fscanf(inputFile,"%f %f %f %f %f %f %f %f %f %f",&TokmakoffAmI->pPro[i],&TokmakoffAmI->cEPro[0][i],&TokmakoffAmI->cEPro[1][i],&TokmakoffAmI->cEPro[2][i],&TokmakoffAmI->cGPro[0][i],&TokmakoffAmI->cGPro[1][i],&TokmakoffAmI->cGPro[2][i],&TokmakoffAmI->cGPro[3][i],&TokmakoffAmI->cGPro[4][i],&TokmakoffAmI->cGPro[5][i]);
  }
  fclose(inputFile);
  printf("Tokmakoff electrostatic map read!\n");
  return;
}

/* Read the file with coupling model information */
//void readCoupling(t_snap *snap,t_coupling *coupling){
  //FILE *inputFile;
  //char *pStatus;
  //char Buffer[256];
  //char dummy[256];
  //size_t LabelLength;
  //char *pValue;
  //int control;

  //control=0;
  //inputFile=fopen(snap->couplingFName,"r");
  //if (inputFile == NULL) {
    //printf("Coupling file not found!\n");
    //exit(-1);
  //}
  //do {
    //pStatus = fgets(&Buffer[0],sizeof(Buffer),inputFile);
    //if (pStatus == NULL) {
      //break;
    //}
    
    // Compute LabelLength
    //LabelLength = strcspn(&Buffer[0], " ");

    // Couplingmodel keyword
    //if (keyWordS("Couplingmodel",Buffer,coupling->cModel,LabelLength)==1) continue;

    // Neighbormodel keyword
    //if (keyWordS("Neighbormodel",Buffer,coupling->nModel,LabelLength)==1) continue;
  
  //} while (TRUE);
  //fclose(inputFile);
  //return;
//}

// Read string
int keyWordS(char *keyWord,char *Buffer,char *value,size_t LabelLength){
  char *pValue;
  char dummy[256];
  if (!strncmp(&Buffer[0],&keyWord[0],LabelLength)){
    printf("%s:",keyWord);
    pValue = &Buffer[LabelLength];
    while (*pValue == ' '){
      pValue++;
    }
      
    sscanf(Buffer,"%s %s",dummy,value);
    printf(" %s\n",value);
    return 1;
  }
  return 0;
}

// Read integer
int keyWordI(char *keyWord,char *Buffer,int *ivalue,size_t LabelLength){
  char *pValue;
  char dummy[256];
  if (!strncmp(&Buffer[0],&keyWord[0],LabelLength)){
    printf("%s:",keyWord);
    pValue = &Buffer[LabelLength];
    while (*pValue == ' '){
      pValue++;
    }
      
    *ivalue=atoi(pValue);
    printf(" %d\n",*ivalue);
    return 1;
  }
  return 0;
}

// Read two integer input
int keyWord2I(char *keyWord,char *Buffer,int *i1,int *i2,size_t LabelLength){
  char *pValue;
  char dummy[256];
  if (!strncmp(&Buffer[0],&keyWord[0],LabelLength)){
    printf("%s:",keyWord);
    pValue = &Buffer[LabelLength];
    while (*pValue == ' '){
      pValue++;
    }
    sscanf(Buffer,"%s %d %d",dummy,i1,i2);
    printf(" %d %d\n",*i1,*i2);
    return 1;
  }
  return 0;
}

// Read float
int keyWordF(char *keyWord,char *Buffer,float *ivalue,size_t LabelLength){
  char *pValue;
  char dummy[256];
  if (!strncmp(&Buffer[0],&keyWord[0],LabelLength)){
    printf("%s:",keyWord);
    pValue = &Buffer[LabelLength];
    while (*pValue == ' '){
      pValue++;
    }
      
    *ivalue=atof(pValue);
    printf(" %f\n",*ivalue);
    return 1;
  }
  return 0;
}

// b is made orthogonal on a
void project(rvec a,rvec b){
  int i;
  float ip;
  ip=iprod(a,b)/iprod(a,a);
  for (i=XX;i<ZZ+1;i++) b[i]-=ip*a[i];
  return;
}


// Initial routine reading topology and command line input
int main(int argc,char *argv[])
{
  const char *desc[] = {
    "This program maps an MD trajectory to trajectories of amide-I Hamiltonian and transition dipoles. The diagonal elements of the Hamiltonian are amide-I frequencies and the off-diagonal elecments are couplings between amide-I sites. Frequency maps used were developed by The Skinner group at UW-Madison, WI,USA, The Jansen group, University of Groningen, The Netherlands, and The Tokmakoff group, University of Chicago. Transition charge coupling and nearest neighbor couplings are treated with the parametrizations developed by The Jansen Group. Transition dipole couplings are maped using Torii/Krimm model. Check the manual for details."
  };

  static char *bugs[] = {
    "Discarding slices for integration should not be necessary."
  };

  static int help=0;
  static char *input="input_AmideImaps";
  //  strcpy(&input[0], "input");
  //  static int end=0;
  //  static int  axis = 2;                      /* normal to memb. default z  */
  //  static char *axtitle="Z"; 
    t_pargs pa [] = {
      { "-i",   FALSE, etSTR, {&input}, "Input file"},
      { "-h",   FALSE, etINT, {&help}, 
        "Help." }};
  //   { "-e",  FALSE, etINT, {&end},
  //     "Number of last frame." }
  //  };

  char      **grpname;            	    /* groupnames                 */
  int       ngrps = 0,                      /* nr. of groups              */
            *ngx;                           /* sizes of groups            */
  int  ePBC;
  t_topology *top;                	    /* topology 		  */ 
  t_trxstatus *status;
  //  atom_id   **index;             	    /* indices for all groups     */
  t_filenm  fnm[] = {             	    /* files for g_order 	  */
    { efTRX, "-f", NULL,  ffREAD },    	    /* trajectory file 	          */
    //    { efNDX, NULL, NULL,  ffREAD },    	    /* index file 		  */
    { efTPX, NULL, NULL,  ffREAD },    	    /* topology file           	  */
    //{ efXVG,"-of","field", ffWRITE }, 	    /* xvgr output file 	  */
  };

#define NFILE asize(fnm)
#define NPA asize(pa)
  const char *fnTPS,*fnNDX,*fnDAT=NULL;
  gmx_bool       bSQ,bRDF;
  output_env_t oenv;

  CopyRight(stderr,argv[0]);
  parse_common_args(&argc,argv,PCA_CAN_VIEW | PCA_CAN_TIME | PCA_BE_NICE,
                    NFILE,fnm,NPA,pa,asize(desc),desc,0,NULL,&oenv);
  
  top = read_top(ftp2fn(efTPX,NFILE,fnm),&ePBC); /* read topology file */
  calc_Ham(ftp2fn(efTRX,NFILE,fnm),top,input,oenv);
  thanx(stderr);
  return 0;
}

