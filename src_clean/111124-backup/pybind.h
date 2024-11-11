//PYBIND11
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <typeinfo>

namespace py = pybind11;

/*
//PYBIND11 STUFF//
template<typename T>
py::list Convert_Pointer_To_PyList(size_t size, T* point)
{
  py::list temp;
  for(size_t i = 0; i < size; i++) temp.append(point[i]);
  return temp;
}
*/
/////////////////////////////////
//SOME FUNCTIONS USED BY PYBIND//
/////////////////////////////////
double* get_epsilon_pointer(ForceField& FF)
{
  return FF.epsilon;
};
double* get_sigma_pointer(ForceField& FF)
{
  return FF.sigma;
};

double get_val_ptr(double* ptr, size_t location)
{
  return ptr[location];
}


////////////////////
//POINTER WRAPPERS//
////////////////////

template <class T> class ptr_wrapper
{
    public:
        ptr_wrapper() : ptr(nullptr) {}
        ptr_wrapper(T* ptr) : ptr(ptr) {}
        ptr_wrapper(const ptr_wrapper& other) : ptr(other.ptr) {}
        T& operator* () const { return *ptr; }
        T* operator->() const { return  ptr; }
        T* get() const { return ptr; }
        void destroy() { delete ptr; }
        T& operator[](std::size_t idx) const { return ptr[idx]; }
    private:
        T* ptr;
};

//void change_val_ptr(T* ptr, size_t location, T newval)
template <typename T>
void change_val_ptr(ptr_wrapper<T>& ptr, size_t location, T newval)
{
  ptr[location] = newval;
}

ptr_wrapper<double> get_ptr_double(Variables& Vars, std::string INPUT)
{ 
  if(INPUT == "FF.epsilon")    {return Vars.FF.epsilon;}
  else if(INPUT == "FF.sigma") {return Vars.FF.sigma;}
  else if(INPUT == "FF.shift") {return Vars.FF.shift;}
}

ptr_wrapper<int> get_ptr_int(Variables& Vars, std::string INPUT)
{ 
  if(INPUT == "FF.FFType"){return Vars.FF.FFType;}
}

template<typename T>
void print_ptr(ptr_wrapper<T> ptr)
{
  int  element = 1;
  bool regular = true;
  
  //if(std::is_same<T, int2>::value) element = 2;
  //if(std::is_same<T, int3>::value) element = 3;
  //if(std::is_same<T, double2>::value) element = 2;
  //if(std::is_same<T, double3>::value) element = 3;
  //if(std::is_same<T, Complex>::value)
  //{
  //  regular = false;
  //}
  for (int i = 0; i < 3; ++i)
  {
    if(regular)  //double.int.double2.double3//
    {
      if(element == 1) std::cout << ptr[i] << " ";
    }
    else
    {
      //if(std::is_same<T, Complex>::value) std::cout << " " << ptr[i].real << " " << ptr[i].imag << " ";
      //if(element == 2) std::cout << " " << ptr[i].x << " " << ptr[i].y << " ";
      //if(element == 3) std::cout << " " << ptr[i].x << " " << ptr[i].y << " " << ptr[i].z << " ";
    }
  }
  std::cout << "\n";
}

template <typename T>
T square(T x)
{
  return x * x;
}

PYBIND11_MODULE(gRASPA, m)
{

  py::class_<double2>(m, "double2")
    .def(py::init<>())
    .def_readwrite("x", &double2::x)
    .def_readwrite("y", &double2::y);

  py::class_<double3>(m, "double3")
    .def(py::init<>())
    .def_readwrite("x", &double3::x)
    .def_readwrite("y", &double3::y)
    .def_readwrite("z", &double3::z);

  py::class_<Complex>(m, "Complex")
    .def(py::init<>())
    .def_readwrite("real", &Complex::real)
    .def_readwrite("imag", &Complex::imag);

  m.def("get_ptr_int",    &get_ptr_int,    "Vars", "INPUT");
  m.def("get_ptr_double", &get_ptr_double, "Vars", "INPUT");
  //m.def("wrap_get_ptr", &wrap_get_ptr<int>, "Vars", "INPUT");
  m.def("square", square<double>);
  m.def("square", square<int>);

  m.def("print_ptr", &print_ptr<int>);
  m.def("print_ptr", &print_ptr<double>);
  //m.def("print_ptr", &print_ptr<double3>);
  m.def("get_val_ptr", &get_val_ptr, "ptr", "location");
  m.def("change_val_ptr", &change_val_ptr<double>, "ptr", "location", "newvalue");
  m.def("change_val_ptr", &change_val_ptr<int>, "ptr", "location", "newvalue");
  m.def("change_val_ptr", &change_val_ptr<size_t>, "ptr", "location", "newvalue");
  m.def("change_val_ptr", &change_val_ptr<double2>, "ptr", "location", "newvalue");
  m.def("change_val_ptr", &change_val_ptr<double3>, "ptr", "location", "newvalue");

  py::class_<RandomNumber>(m, "RandomNumber")
    .def(py::init<>())
    .def_readwrite("host_random", &RandomNumber::host_random)
    .def("AllocateRandom", &RandomNumber::AllocateRandom, "Allocate Space for RandomNumber on CPU")
    .def("DeviceRandom", &RandomNumber::DeviceRandom, "Allocate Space for RandomNumber on GPU")
    .def_readwrite("randomsize", &RandomNumber::randomsize)
    .def_readwrite("offset", &RandomNumber::offset)
    .def_readwrite("Rounds", &RandomNumber::Rounds);

  py::class_<MoveEnergy>(m, "MoveEnergy")
    .def(py::init<>())
    .def_readwrite("HHVDW", &MoveEnergy::HHVDW)
    .def_readwrite("HHReal", &MoveEnergy::HHReal)
    .def_readwrite("HHEwaldE", &MoveEnergy::HHEwaldE)

    .def_readwrite("HGVDW", &MoveEnergy::HGVDW)
    .def_readwrite("HGReal", &MoveEnergy::HGReal)
    .def_readwrite("HGEwaldE", &MoveEnergy::HGEwaldE)
 
    .def_readwrite("GGVDW", &MoveEnergy::GGVDW)
    .def_readwrite("GGReal", &MoveEnergy::GGReal)
    .def_readwrite("GGEwaldE", &MoveEnergy::GGEwaldE) 

    .def_readwrite("TailE", &MoveEnergy::TailE)
    .def_readwrite("DNN_E", &MoveEnergy::DNN_E)

    .def("total", &MoveEnergy::total, "Calculate the total energy based on the component values")
    .def("take_negative", &MoveEnergy::take_negative, "Calculate the negated values")
    .def("zero", &MoveEnergy::zero, "zero the components of energies")
    .def("print", &MoveEnergy::print, "print the energies");
  //m.def("DeviceRandom", &DeviceRandom, "RandomNumber on GPU");
  //m.def("set_value", &set_value, "Set value");
  m.def("RUN", &pymain, "RUN SIMULATION");

  py::enum_<SIMULATION_MODE>(m, "SIMULATION_MODE")
    .value("CREATE_MOLECULE", CREATE_MOLECULE)
    .value("INITIALIZATION", INITIALIZATION)
    .value("EQUILIBRATION", EQUILIBRATION)
    .value("PRODUCTION", PRODUCTION)
    .export_values();

  py::class_<ForceField>(m, "ForceField")
    .def(py::init<>())
    .def_readwrite("VDWRealBias", &ForceField::VDWRealBias)
    .def_readwrite("noCharges", &ForceField::noCharges)
    .def_readwrite("size", &ForceField::size)
    .def_readwrite("OverlapCriteria", &ForceField::OverlapCriteria)
    .def_readwrite("CutOffVDW", &ForceField::CutOffVDW)
    .def_readwrite("CutOffCoul", &ForceField::CutOffCoul)
    .def_readonly("epsilon", &ForceField::epsilon);
  m.def("get_epsilon_pointer", &get_epsilon_pointer, py::return_value_policy::reference, "FF");
  m.def("get_sigma_pointer",   &get_sigma_pointer,   py::return_value_policy::reference, "FF");

  py::class_<Boxsize>(m, "Boxsize")
    .def(py::init<>())
    .def_readwrite("Cubic", &Boxsize::Cubic)
    .def_readwrite("Alpha", &Boxsize::Alpha);

  py::class_<Simulations>(m, "Simulations")
    .def(py::init<>())
    .def_readwrite("Nblocks", &Simulations::Nblocks)
    .def_readwrite("Box", &Simulations::Box);

  py::class_<ptr_wrapper<double>>(m,"p");
  py::class_<ptr_wrapper<Simulations>>(m,"pSim");
  //m.def("get_ptr", &get_ptr<Simulations>, "Vars", "INPUT");

  py::class_<Atom_FF>(m, "Atom_FF")
    .def(py::init<>())
    .def_readwrite("Name",    &Atom_FF::Name)
    .def_readwrite("epsilon", &Atom_FF::epsilon)
    .def_readwrite("sigma",   &Atom_FF::sigma)
    .def_readwrite("shift",   &Atom_FF::shift)
    .def_readwrite("tail",    &Atom_FF::tail);

  py::class_<Tail>(m, "Tail")
    .def(py::init<>())
    .def_readwrite("UseTail", &Tail::UseTail)
    .def_readwrite("Energy",  &Tail::Energy);

  py::class_<Input_Container>(m, "Input_Container")
    .def(py::init<>())
    .def_readwrite("AtomFF", &Input_Container::AtomFF)
    .def_readwrite("Mix_Epsilon", &Input_Container::Mix_Epsilon)
    .def_readwrite("Mix_Sigma", &Input_Container::Mix_Sigma)
    .def_readwrite("Mix_Shift", &Input_Container::Mix_Shift)
    .def_readwrite("Mix_Tail", &Input_Container::Mix_Tail)
    .def_readwrite("CutOffVDW", &Input_Container::CutOffVDW)
    .def_readwrite("CutOffCoul", &Input_Container::CutOffCoul)
    .def_readwrite("VDWRealBias", &Input_Container::VDWRealBias)
    .def_readwrite("OverlapCriteria", &Input_Container::OverlapCriteria)
    .def_readwrite("noCharges", &Input_Container::noCharges);
    

  py::class_<Variables>(m, "Variables")
    .def(py::init<>())
    .def("set_TEST", &Variables::set_TEST)
    .def("get_TEST", &Variables::get_TEST)
    .def_readwrite("Input", &Variables::Input)
    .def_readwrite("TEST", &Variables::TEST)
    .def_readwrite("Ttwo", &Variables::Ttwo)
    .def_readwrite("FF", &Variables::FF)
    .def_readwrite("Random", &Variables::Random)
    .def_readwrite("NumberOfInitializationCycles", &Variables::NumberOfInitializationCycles)
    .def_readwrite("NumberOfEquilibrationCycles", &Variables::NumberOfEquilibrationCycles)
    .def_readwrite("NumberOfProductionCycles", &Variables::NumberOfProductionCycles)
    .def_readwrite("SimulationMode", &Variables::SimulationMode);

  //m.def("get_total_energy", &check_energy_wrapper, "get total energy", "Var", "SimulationIndex"_a=0);
  m.def("get_total_energy", &check_energy_wrapper, "get total energy", "Var", "SimulationIndex");
  //RUN MOVES//
  m.def("Run_Simulation_ForOneBox", &Run_Simulation_ForOneBox, "run simulation for selected box", "Vars", "box_index", "SimulationMode");

  m.def("ForceField_Processing", &ForceField_Processing, "Input");
  m.def("Copy_InputLoader_Data", &Copy_InputLoader_Data, "Vars");
  //m.def("RunMoves", &RunMoves, "Run a single Monte Carlo Move");
}
