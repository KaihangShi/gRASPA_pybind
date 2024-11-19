import gRASPA as g
import numpy as np

# run steps independently #
W = g.Initialize()
W.NumberOfInitializationCycles = 100
W.NumberOfProductionCycles = 0
g.RUN(W)
g.finalize(W)

# run the single-box MC moves through python #
box_index = 0 # run simulation box 0 #
g.InitializeMC(W, box_index)

g.Select_Box_Component_Molecule(W, box_index)
# run single MC move #
r = W.MCMoveVariables.RandomNumber
comp = W.MCMoveVariables.component
mol  = W.MCMoveVariables.molecule

#MoveType = g.MoveTypes.SINGLE_INSERTION
MoveType = g.MoveTypes.SINGLE_DELETION
#MoveType = g.MoveTypes.ROTATION

W.MCMoveVariables.MoveType = int(MoveType)

# PERFORM MOVE #
g.SingleBody_Prepare(W, systemId = box_index)
DeltaE = g.SingleBody_Calculation(W, systemId = box_index)

# Get Adsorbate data (old)
# Get framework data, if rigid, should just run it once!#
FrameworkData = g.GetAllAtoms(W, systemId = 0, component = 0)
# Get trial positions (new), original dict will get new entries: "Trial_pos", "Trial_charge", ...
# If WholeConfig = False, then only the moved part is here #
# else, the moved part will be merged so that all molecules + trial molecule (moved) is copied here #
OldData = g.GetAllAtoms(W, systemId = 0, component = 1)
TrialData = g.GetTrialConfig(W, systemId = 0, component = 1, WholeConfig = False)
Data = g.GetTrialConfig(W, systemId = 0, component = 1, WholeConfig = True)
print(f"OldData: {OldData}\n")
print(f"TrialData: {TrialData}\n")
print(f"Data: {Data}\n")
# combine old config with trial config
#Data["Trial_pos"]

#g.SingleBody_Acceptance(W, box_index, DeltaE)


