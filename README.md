# gRASPA pybind: pybind patch for gRASPA (test)
* enjoy gRASPA and GPU backend while on python and Jupyter Notebook!

## Advantages
* **Access**:  Get access to gRASPA's internal variables
* **Control**: See what happens during an MC move
* **Create**:  Create your own MC moves with python code and packaged MC move parts

## How to use
1. download gRASPA_pybind and [Zhaoli2042/gRASPA_fork](https://github.com/Zhaoli2042/gRASPA_fork)
    * support for [snurrlab/gRASPA](https://github.com/snurr-group/gRASPA) will be available soon!
3. copy the files in gRASPA's [src_clean](https://github.com/snurr-group/gRASPA/tree/main/src_clean) to gRASPA_pybind's [pybind_src](https://github.com/Zhaoli2042/gRASPA_pybind/tree/main/pybind_src)
4. Compile using `./BIND_NVC_COMPILE`
    * If successful, you can see a **shared library file (gRASPA.so)** being generated
4. Copy `gRASPA.so` file to the `JUPYTER_NOTEBOOKS` folder
5. Run `jupyter notebook` and open one of the example notebooks!
