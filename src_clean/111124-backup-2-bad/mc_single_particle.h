#include "mc_utilities.h"

////////////////////////////////////////////////
// Generalized function for single Body moves //
////////////////////////////////////////////////

//Zhao's note: decompose the single body move into different sections: preparation, calculation, and acceptance //
//For easier manipulation of moves, to use for, for example, MLP //
inline void SingleBody_Prepare(Variables& Vars, size_t systemId, size_t SelectedMolInComponent, size_t SelectedComponent, int MoveType)
{
  Components& SystemComponents = Vars.SystemComponents[systemId];
  Simulations& Sims            = Vars.Sims[systemId];
  ForceField& FF               = Vars.device_FF;
  RandomNumber& Random         = Vars.Random;
  WidomStruct& Widom           = Vars.Widom[systemId];

  bool& Do_New  = Vars.TempVal.Do_New;
  bool& Do_Old  = Vars.TempVal.Do_Old;
  Do_New = false;
  Do_Old = false;

  size_t Molsize = SystemComponents.Moleculesize[SelectedComponent]; //Get the size of the selected Molecule
  //Set up Old position and New position arrays
  if(Molsize >= 1024)
  {
    throw std::runtime_error("Molecule size is greater than allocated size, Why so big?\n");
  }
  size_t& start_position = Vars.TempVal.start_position; 
  start_position = SelectedMolInComponent*SystemComponents.Moleculesize[SelectedComponent];

  SystemComponents.Moves[SelectedComponent].Record_Move_Total(MoveType);
  double3 MaxChange = {0.0, 0.0, 0.0};
  switch (MoveType)
  {
    case TRANSLATION:
    {
      Do_New = true; Do_Old = true;
      MaxChange = SystemComponents.MaxTranslation[SelectedComponent];
      break;
    }
    case ROTATION:
    {
      Do_New = true; Do_Old = true;
      MaxChange = SystemComponents.MaxRotation[SelectedComponent];
      break;
    }
    case SPECIAL_ROTATION:
    {
      Do_New = true; Do_Old = true;
      MaxChange = SystemComponents.MaxSpecialRotation[SelectedComponent];
      //Zhao's note: if we separate framework components, there might be lots of overlaps between different species (node and linker overlaps), we can turn this Overlap flag off//
      //printf("Performing move on %zu comp, %zu mol\n", SelectedComponent, SelectedMolInComponent);
      break;
    }
    case SINGLE_INSERTION:
    {
      Do_New = true;
      start_position = 0;
      break;
    } 
    case SINGLE_DELETION:
    {
      Do_Old = true;
      break;
    }
  }
  if(!Do_New && !Do_Old) throw std::runtime_error("Doing Nothing For Single Particle Move?\n");

  //Zhao's note: possible bug, you may only need 3 instead of 3 * N random numbers//
  Random.Check(Molsize);
  get_new_position<<<1, Molsize>>>(Sims, FF, start_position, SelectedComponent, MaxChange, Random.device_random, Random.offset, MoveType);
  Random.Update(Molsize);
}

inline MoveEnergy SingleBody_Calculation(Variables& Vars, size_t systemId, size_t SelectedMolInComponent, size_t SelectedComponent, int MoveType)
{
  Components& SystemComponents = Vars.SystemComponents[systemId];
  Simulations& Sims            = Vars.Sims[systemId];
  ForceField& FF               = Vars.device_FF;
  WidomStruct& Widom           = Vars.Widom[systemId];

  bool& Overlap = Vars.TempVal.Overlap; Overlap = true;
  bool& Do_New  = Vars.TempVal.Do_New;
  bool& Do_Old  = Vars.TempVal.Do_Old;

  size_t Molsize = SystemComponents.Moleculesize[SelectedComponent]; //Get the size of the selected Molecule
  // Setup for the pairwise calculation //
  // New Features: divide the Blocks into two parts: Host-Guest + Guest-Guest //
  
  size_t NHostAtom = 0; size_t NGuestAtom = 0;
  for(size_t i = 0; i < SystemComponents.NComponents.y; i++)
    NHostAtom += SystemComponents.Moleculesize[i] * SystemComponents.NumberOfMolecule_for_Component[i];
  for(size_t i = SystemComponents.NComponents.y; i < SystemComponents.NComponents.x; i++)
    NGuestAtom+= SystemComponents.Moleculesize[i] * SystemComponents.NumberOfMolecule_for_Component[i];

  //Zhao's note: Cross term, if the selected species is host atom, the crossAtom = guest, vice versa//
  size_t NCrossAtom = NHostAtom;
  if(SelectedComponent < SystemComponents.NComponents.y) //Framework component//
    NCrossAtom = NGuestAtom;
  size_t HH_Nthread=0; size_t HH_Nblock=0; Setup_threadblock(NHostAtom *  Molsize, &HH_Nblock, &HH_Nthread);
  size_t HG_Nthread=0; size_t HG_Nblock=0; Setup_threadblock(NCrossAtom * Molsize, &HG_Nblock, &HG_Nthread);
  size_t GG_Nthread=0; size_t GG_Nblock=0; Setup_threadblock(NGuestAtom * Molsize, &GG_Nblock, &GG_Nthread);

  size_t SameTypeNthread = 0;
  if(SelectedComponent < SystemComponents.NComponents.y) //Framework-Framework + Framework-Adsorbate//
  {GG_Nthread = 0; GG_Nblock = 0; SameTypeNthread = HH_Nthread; }
  else //Framework-Adsorbate + Adsorbate-Adsorbate//
  {HH_Nthread = 0; HH_Nblock = 0; SameTypeNthread = GG_Nthread; }

  size_t Nthread = std::max(SameTypeNthread, HG_Nthread);
  size_t Total_Nblock  = HH_Nblock + HG_Nblock + GG_Nblock;

  int3 NBlocks = {(int) HH_Nblock, (int) HG_Nblock, (int) GG_Nblock}; //x: HH_Nblock, y: HG_Nblock, z: GG_Nblock;
  //printf("Total_Comp: %zu, Host Comp: %zu, Adsorbate Comp: %zu\n", SystemComponents.NComponents.x, SystemComponents.NComponents.y, SystemComponents.NComponents.z);
  //printf("NHostAtom: %zu, HH_Nblock: %zu, HG_Nblock: %zu, NGuestAtom: %zu, GG_Nblock: %zu\n", NHostAtom, HH_Nblock, HG_Nblock, NGuestAtom, GG_Nblock);
  size_t Atomsize = 0;
  for(size_t ijk = 0; ijk < SystemComponents.NComponents.x; ijk++)
    Atomsize += SystemComponents.Moleculesize[ijk] * SystemComponents.NumberOfMolecule_for_Component[ijk];

  if(Atomsize != 0)
  {
    Calculate_Single_Body_Energy_VDWReal<<<Total_Nblock, Nthread, Nthread * 2 * sizeof(double)>>>(Sims.Box, Sims.d_a, Sims.Old, Sims.New, FF, Sims.Blocksum, SelectedComponent, Atomsize, Molsize, Sims.device_flag, NBlocks, Do_New, Do_Old, SystemComponents.NComponents);

    cudaMemcpy(SystemComponents.flag, Sims.device_flag, sizeof(bool), cudaMemcpyDeviceToHost);
  }
  MoveEnergy tot; 
  if(!SystemComponents.flag[0] || !Overlap)
  {
    double BlockResult[Total_Nblock + Total_Nblock];
    cudaMemcpy(BlockResult, Sims.Blocksum, 2 * Total_Nblock * sizeof(double), cudaMemcpyDeviceToHost);
   
    //VDW Part and Real Part Coulomb//
    for(size_t i = 0; i < HH_Nblock; i++) 
    {
      tot.HHVDW += BlockResult[i];
      tot.HHReal+= BlockResult[i + Total_Nblock];
      //if(MoveType == SPECIAL_ROTATION) printf("HH Block %zu, VDW: %.5f, Real: %.5f\n", i, BlockResult[i], BlockResult[i + Total_Nblock]);
    }
    for(size_t i = HH_Nblock; i < HH_Nblock + HG_Nblock; i++) 
    {
      tot.HGVDW += BlockResult[i];
      tot.HGReal+= BlockResult[i + Total_Nblock];
      //printf("HG Block %zu, VDW: %.5f, Real: %.5f\n", i, BlockResult[i], BlockResult[i + Total_Nblock]);
    }
    for(size_t i = HH_Nblock + HG_Nblock; i < Total_Nblock; i++)
    {
      tot.GGVDW += BlockResult[i];
      tot.GGReal+= BlockResult[i + Total_Nblock];
      //printf("GG Block %zu, VDW: %.5f, Real: %.5f\n", i, BlockResult[i], BlockResult[i + Total_Nblock]);
    }

    /*
    printf("HG_NBlock: %zu\n", Total_Nblock);
    printf("Separated VDW : %.5f (HH), %.5f (HG), %.5f (GG)\n", tot.HHVDW,  tot.HGVDW , tot.GGVDW);
    printf("Separated Real: %.5f (HH), %.5f (HG), %.5f (GG)\n", tot.HHReal, tot.HGReal, tot.GGReal);
    */

    // Calculate Ewald //
    bool EwaldPerformed = false;
    if(!FF.noCharges && SystemComponents.hasPartialCharge[SelectedComponent])
    {
      double2 newScale  = SystemComponents.Lambda[SelectedComponent].SET_SCALE(1.0);
      double2 EwaldE = GPU_EwaldDifference_General(Sims.Box, Sims.d_a, Sims.New, Sims.Old, FF, Sims.Blocksum, SystemComponents, SelectedComponent, MoveType, 0, newScale);
      if(HH_Nblock == 0)
      {
        tot.GGEwaldE = EwaldE.x;
        tot.HGEwaldE = EwaldE.y;
      }
      else
      {
        tot.HHEwaldE = EwaldE.x;
        tot.HGEwaldE = EwaldE.y;
        //printf("HHEwald: %.5f, HGEwald: %.5f\n", tot.HHEwaldE, tot.HGEwaldE);
      }
      EwaldPerformed = true;
    }
    if(SystemComponents.UseDNNforHostGuest)
    {
      //Calculate DNN//
      if(!EwaldPerformed) Prepare_DNN_InitialPositions(Sims.d_a, Sims.New, Sims.Old, SystemComponents, SelectedComponent, MoveType, 0);
      tot.DNN_E = DNN_Prediction_Move(SystemComponents, Sims, SelectedComponent, MoveType);
      double correction = tot.DNN_Correction(); //If use DNN, HGVDWReal and HGEwaldE are zeroed//
      if(fabs(correction) > SystemComponents.DNNDrift) //If there is a huge drift in the energy correction between DNN and Classical HostGuest//
      {
        //printf("TRANSLATION/ROTATION: Bad Prediction, reject the move!!!\n");
        switch(MoveType)
        {
          case TRANSLATION: case ROTATION:
          {
            SystemComponents.TranslationRotationDNNReject ++; break;
          }
          case SINGLE_INSERTION: case SINGLE_DELETION:
          {
            SystemComponents.SingleSwapDNNReject ++; break;
          }
        }
        WriteOutliers(SystemComponents, Sims, NEW, tot, correction); //Print New Locations//
        WriteOutliers(SystemComponents, Sims, OLD, tot, correction); //Print Old Locations//
        tot.zero();
        return tot;
      }
      SystemComponents.SingleMoveDNNDrift += fabs(correction);
    }
  }
  return tot;
}

inline void SingleBody_Acceptance(Variables& Vars, size_t systemId, size_t SelectedMolInComponent, size_t SelectedComponent, int MoveType, MoveEnergy& tot)
{
  Components& SystemComponents = Vars.SystemComponents[systemId];
  Simulations& Sims            = Vars.Sims[systemId];
  ForceField& FF               = Vars.device_FF;
  WidomStruct& Widom           = Vars.Widom[systemId];

  size_t Molsize = SystemComponents.Moleculesize[SelectedComponent]; //Get the size of the selected Molecule
  size_t& start_position = Vars.TempVal.start_position;
  bool Accept = false; 

  //Get Number of Molecules for this component (For updating TMMC)//
  double NMol = SystemComponents.NumberOfMolecule_for_Component[SelectedComponent];
  if(SystemComponents.hasfractionalMolecule[SelectedComponent]) NMol--;

  double preFactor = GetPrefactor(SystemComponents, Sims, SelectedComponent, MoveType);
  double Pacc      = preFactor * std::exp(-SystemComponents.Beta * tot.total());

  //Apply the bias according to the macrostate//
  if(MoveType == SINGLE_INSERTION || MoveType == SINGLE_DELETION)
  {
    SystemComponents.Tmmc[SelectedComponent].ApplyWLBias(preFactor, NMol, MoveType);
    SystemComponents.Tmmc[SelectedComponent].ApplyTMBias(preFactor, NMol, MoveType);
  }
  //if(MoveType == SINGLE_INSERTION) printf("SINGLE INSERTION, tot: %.5f, preFactor: %.5f, Pacc: %.5f\n", tot.total(), preFactor, Pacc);
  //if(MoveType == SINGLE_DELETION)  printf("SINGLE DELETION,  tot: %.5f, preFactor: %.5f, Pacc: %.5f\n", tot.total(), preFactor, Pacc);
  //

  double Random = Get_Uniform_Random();
  if(Random < preFactor * std::exp(-SystemComponents.Beta * tot.total())) Accept = true;

  switch(MoveType)
  {
    case TRANSLATION: case ROTATION: case SPECIAL_ROTATION:
    { 
      if(Accept)
      {
        update_translation_position<<<1,Molsize>>>(Sims.d_a, Sims.New, start_position, SelectedComponent);
        SystemComponents.Moves[SelectedComponent].Record_Move_Accept(MoveType);
        if(!FF.noCharges && SystemComponents.hasPartialCharge[SelectedComponent])
        {
          Update_Vector_Ewald(Sims.Box, false, SystemComponents, SelectedComponent);
        }
      }
      else {tot.zero(); };
      SystemComponents.Tmmc[SelectedComponent].Update(1.0, NMol, MoveType);
      break;
    }
    case SINGLE_INSERTION:
    {
      SystemComponents.Tmmc[SelectedComponent].TreatAccOutofBound(Accept, NMol, MoveType);
      if(Accept)
      {
        SystemComponents.Moves[SelectedComponent].Record_Move_Accept(MoveType);
        AcceptInsertion(SystemComponents, Sims, SelectedComponent, 0, FF.noCharges, SINGLE_INSERTION); //0: selectedTrial//
      }
      else {tot.zero(); };
      SystemComponents.Tmmc[SelectedComponent].Update(Pacc, NMol, MoveType);
      break;
    }
    case SINGLE_DELETION:
    {
      SystemComponents.Tmmc[SelectedComponent].TreatAccOutofBound(Accept, NMol, MoveType);
      if(Accept)
      {
        SystemComponents.Moves[SelectedComponent].Record_Move_Accept(MoveType);
        size_t UpdateLocation = SelectedMolInComponent * SystemComponents.Moleculesize[SelectedComponent];
        AcceptDeletion(SystemComponents, Sims, SelectedComponent, UpdateLocation, SelectedMolInComponent, FF.noCharges);
      }
      else {tot.zero(); };
      SystemComponents.Tmmc[SelectedComponent].Update(Pacc, NMol, MoveType);
      break;
    }
  }
  //return Accept;
}

inline MoveEnergy SingleBodyMove(Variables& Vars, size_t systemId, size_t SelectedMolInComponent, size_t SelectedComponent, int MoveType)
{
  SingleBody_Prepare(Vars, systemId, SelectedMolInComponent, SelectedComponent, MoveType);
  MoveEnergy tot = SingleBody_Calculation(Vars, systemId, SelectedMolInComponent, SelectedComponent, MoveType);
  if(!Vars.TempVal.Overlap)
    SingleBody_Acceptance(Vars, systemId, SelectedMolInComponent, SelectedComponent, MoveType, tot);
  return tot;
}
