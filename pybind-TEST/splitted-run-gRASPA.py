import gRASPA as g
# run steps independently #
W = g.Initialize()
W.NumberOfInitializationCycles = 10000
# run the single-box MC moves through python #
box_index = 0 # run simulation box 0 #
g.InitializeMC(W, box_index)

start = g.omp_get_wtime()
for i in range(0, W.NumberOfInitializationCycles):
  Steps = g.Determine_Number_Of_Steps(W, box_index, i)
  for s in range(0, Steps):
    g.Select_Box_Component_Molecule(W, box_index)
    # run single MC move #
    r = W.MCMoveVariables.RandomNumber
    comp = W.MCMoveVariables.component
    mol  = W.MCMoveVariables.molecule

    MoveType = g.MoveTypes.TRANSLATION
    if(r < W.SystemComponents[box_index].Moves[comp].TranslationProb):
      MoveType = g.MoveTypes.TRANSLATION
      if(W.SystemComponents[box_index].NumberOfMolecules[comp] == 0): continue
    elif(r < W.SystemComponents[box_index].Moves[comp].RotationProb):
      MoveType = g.MoveTypes.ROTATION
      if(W.SystemComponents[box_index].NumberOfMolecules[comp] == 0): continue
    elif(r < W.SystemComponents[box_index].Moves[comp].SwapProb):
      if(g.Get_Uniform_Random() < 0.5):
        MoveType = g.MoveTypes.SINGLE_INSERTION
      else:
        MoveType = g.MoveTypes.SINGLE_DELETION
        if(W.SystemComponents[box_index].NumberOfMolecules[comp] == 0): continue
    W.MCMoveVariables.MoveType = int(MoveType)
    # PERFORM MOVE #
    g.SingleBody_Prepare(W, systemId = box_index)
    DeltaE = g.SingleBody_Calculation(W, systemId = box_index)
    g.SingleBody_Acceptance(W, box_index, DeltaE)
    g.MoveEnergy_Add(W.SystemComponents[box_index].deltaE, DeltaE) # add DeltaE (move delta) to deltaE (total delta) #
  # Gather statistics and averages per cycle #
  g.GatherStatisticsDuringSimulation(W, box_index, i)
g.MCEndOfPhaseSummary(W) # end of initialization phase summary #
end = g.omp_get_wtime()
print(f"simulation took {end - start} seconds\n")

g.finalize(W) # end of simulation #
'''
# save atom data for system 0, component 1 (co2 in this case)
DataA = g.GetAllAtoms(V, 0, 1)
# Run some MC moves #
# run multiple times of rotation moves, on the 10th CO2 molecule #
g.SingleBody_Prepare(V, 0, 10, 1, 1) # prepare the MC move
DeltaE = g.SingleBody_Calculation(V, 0, 10, 1, 1) # calculate the classical energies (vdW + electrostatics)
DeltaE.print() # check the deltaE for the rotation move
g.SingleBody_Acceptance(V, 0, 10, 1, 1, DeltaE) # determine whether to accept or reject the move

# run it multiple times ...

# calculate current total energy
TotalE = g.get_total_energy(V, 0)
# print total energy
TotalE.print()
'''
