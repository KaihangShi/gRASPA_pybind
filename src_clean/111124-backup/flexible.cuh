double GenerateBondLength(Components& SystemComponents, size_t BondType, size_t comp, size_t BondID)
{
  double BondLength = 0.0; 
  switch(BondType)
  {
    case HARMONIC_BOND:
    {
      BOND RefBond = SystemComponents.RefBond[comp][BondID];
      double ran1=1.0/sqrt(SystemComponents.Beta*RefBond.ForceConstant);
      double ran2=1.0/pow((SystemComponents.RefLength+3.0*ran1), 2);
      while(Get_Uniform_Random() > SQR(BondLength) * ran2)
      {
        BondLength = SystemComponents.RefLength + ran1 * Get_Gaussian_Random();
      }
      break;
    }
    case FIXED_BOND: case RIGID_BOND:
    {
      BondLength = SystemComponents.RefBond[comp][BondID].RefLength;
      break;
    }
  }
  return BondLength;
}

double GenerateBendAngle(Components& SystemComponents, size_t BondType, size_t comp, size_t BendID)
{
  double theta = 0.0;
  switch(BendType)
  {
    case HARMONIC_BEND:
    {
      BEND RefBend = SystemComponents.RefBend[comp][BondID];
      theta = 0.0;
      double energy= 0.0;
      while(Get_Uniform_Random() > SQR(std::sin(theta))*exp(-SystemComponents.Beta*energy))
      {
        theta = M_PI * Get_Uniform_Random();
        energy= 0.5 * RefBend.ForceConstant * SQR(theta);
      }
      break;
    }
    case FIXED_BEND: case RIGID_BEND:
    {
      theta = SystemComponents.RefBend[comp][BendID].RefAngle;
      break;
    }
  }
  return theta;
}

void GenerateTrialOrientationsMCScheme(
