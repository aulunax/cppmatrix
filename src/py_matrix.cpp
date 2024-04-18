#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>
#include "Matrix.h"

namespace py = pybind11;

PYBIND11_MODULE(py_matrix, m) {
    py::register_exception<MatrixSizeDisparityException>(m, "MatrixSizeDisparityException");
    py::register_exception<NotSquareMatrixException>(m, "NotSquareMatrixException");
    py::register_exception<MatrixEquationNoUniqueSolutionException>(m, "MatrixEquationNoUniqueSolutionException");

    py::class_<Dimensions>(m, "Dimensions")
        .def(py::init<int, int>())
        .def_readwrite("n", &Dimensions::n)
        .def_readwrite("m", &Dimensions::m)
        .def("__eq__", &Dimensions::operator==)
        .def("__ne__", &Dimensions::operator!=);

    py::class_<Matrix<double>>(m, "Matrix")
        .def(py::init<>())
        .def(py::init<std::string>())
        .def(py::init<int, int>())
        .def(py::init<int, int, double>())
        .def(py::init<Dimensions>())
        .def(py::init<Dimensions, double>())
        .def(py::init<const Matrix<double>&>())
        //.def(py::init<Matrix<double>&&>(), py::return_value_policy::move)

        // copy and move assignment
        .def("__copy__",  [](const Matrix<double> &self) {
        return Matrix<double>(self);
        })
        .def("__deepcopy__", [](const Matrix<double> &self, py::dict) {
            return Matrix<double>(self);
        })

        // operators
        .def(py::self + py::self)
        .def(py::self - py::self)
        .def(py::self * py::self)
        //.def(double() * py::self)
        .def(py::self * double())
        .def(py::self / double())
        .def(-py::self)

        .def(py::self / py::self)
        .def(py::self | py::self)

        .def("__getitem__", [](Matrix<double> &obj, std::tuple<int, int> indices) {
            int row = std::get<0>(indices);
            int col = std::get<1>(indices);
            if (row >= obj.size.n || col >= obj.size.m) {
                throw py::index_error();
            }
            return obj[Dimensions{row,col}];
        })
        .def("__setitem__", [](Matrix<double> &obj, std::tuple<int, int> indices, int value) {
            int row = std::get<0>(indices);
            int col = std::get<1>(indices);
            if (row >= obj.size.n || col >= obj.size.m) {
                throw py::index_error();
            }
            obj[Dimensions{row,col}] = value;
        })
        .def(py::self == py::self)
        .def(py::self != py::self)

        // methods
        .def("transpose", &Matrix<double>::transpose)
        .def("triu", &Matrix<double>::triu)
        .def("tril", &Matrix<double>::tril)
        .def("diag", &Matrix<double>::diag)
        .def("fill", &Matrix<double>::fill)
        .def("reserve", &Matrix<double>::reserve)
        .def("print", &Matrix<double>::print)

        .def("threadedMatrixOperation", &Matrix<double>::threadedMatrixOperation)
        // static methods
        .def_static("fillWithValue", &Matrix<double>::fillWithValue)
        .def_static("addMatrices", &Matrix<double>::addMatrices)
        .def_static("substractMatrices", &Matrix<double>::substractMatrices)
        .def_static("multiplyMatrices", &Matrix<double>::multiplyMatrices)
        .def_static("multiplyByConstant", &Matrix<double>::multiplyByConstant)
        .def_static("divideByConstant", &Matrix<double>::divideByConstant)
        .def_static("getDiagonal", &Matrix<double>::getDiagonal)
        .def_static("getTriangle", &Matrix<double>::getTriangle)
        .def_static("getTransposed", &Matrix<double>::getTransposed)
        .def_static("gaussianElimination", &Matrix<double>::gaussianElimination)
        

;
}