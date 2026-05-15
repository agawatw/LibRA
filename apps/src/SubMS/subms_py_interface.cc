#include <pybind11/pybind11.h>
#include <subms.h>

namespace py = pybind11;

PYBIND11_MODULE(subms2py, m) {
    m.def("SubMS_func", &SubMS_func,
          py::arg("MSNBuf")      = "",
          py::arg("OutMSBuf")    = "",
          py::arg("WhichColStr") = "data",
          py::arg("fieldStr")    = "*",
          py::arg("timeStr")     = "",
          py::arg("spwStr")      = "*",
          py::arg("baselineStr") = "",
          py::arg("scanStr")     = "",
          py::arg("arrayStr")    = "",
          py::arg("uvdistStr")   = "",
          py::arg("taqlStr")     = "",
          py::arg("corrStr")     = "",
          py::arg("intentStr")   = "",
          py::arg("obsStr")      = "",
          py::arg("combineStr")  = "",
          py::arg("integ")       = -1.0,
          py::arg("chanStep")    = 1,
          py::arg("overwrite")   = false,
          py::arg("verbose")     = false);
}
