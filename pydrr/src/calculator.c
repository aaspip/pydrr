/*
This file has been written with love by Burak Topal (lexygon) for giving a simple example for writing C extensions for Python 3.5.
Feel free to cloning, sharing, editing and committing some new examples.
I have tried to explain each part basicly as I can.
For communicating with me:
mail: hi@buraktopal.xyz
github: github.com/lexygon

For documentation for Python/C API please visit https://docs.python.org/3/c-api/
*/

// importing Python C API Header
#include <Python.h>

// creating functions that returning PyObject.
static PyObject *addition(PyObject *self, PyObject *args){
  // variables for our parameters. our parameters that are coming from python will be stored in theese variables.
  int number1;
  int number2;
  int result;

  // Parsing our Python parameters to C variables.
  // "ii" means we are taking 2 integer variables from Python.
  // if we were taking 2 integer and 1 string that would be "iis".
  // after parsing python variables, this is sending them to number1 and number2 variables. ORDER IS IMPORTANT!!
  if (!PyArg_ParseTuple(args, "ii", &number1, &number2))
         // if sending parameters are not fitting to types, it will return NULL
         return NULL;

  // after parsing, we are doing our job.
  result = number1 + number2;

  // like before, this part is parsing our variable to a python value and returning it.
  // in here i means we are returning an integer that comes from result variable.
  return Py_BuildValue("i", result);

}

//same structure with addition
static PyObject *substraction(PyObject *self, PyObject *args){
  int number1;
  int number2;
  int result;

  if (!PyArg_ParseTuple(args, "ii", &number1, &number2))
         return NULL;

  result = number1 - number2;

  return Py_BuildValue("i", result);

}

//same structure with addition
static PyObject *division(PyObject *self, PyObject *args){
  int number1;
  int number2;
  int result;

  if (!PyArg_ParseTuple(args, "ii", &number1, &number2))
         return NULL;

  result = number1 / number2;

  return Py_BuildValue("i", result);

}

//same structure with addition
static PyObject *multiplication(PyObject *self, PyObject *args){
  int number1;
  int number2;
  int result;

  if (!PyArg_ParseTuple(args, "ii", &number1, &number2))
         return NULL;

  result = number1 * number2;

  return Py_BuildValue("i", result);

}

// documentation for each functions.
static char addition_document[] = "Document stuff for addition...";
static char substraction_document[] = "Document stuff for substraction...";
static char division_document[] = "Document stuff for division...";
static char multiplication_document[] = "Document stuff for multiplication...";

// defining our functions like below:
// function_name, function, METH_VARARGS flag, function documents
static PyMethodDef functions[] = {
  {"addition", addition, METH_VARARGS, addition_document},
  {"substraction", substraction, METH_VARARGS, substraction_document},
  {"division", division, METH_VARARGS, division_document},
  {"multiplication", multiplication, METH_VARARGS, multiplication_document},
  {NULL, NULL, 0, NULL}
};

// initializing our module informations and settings in this structure
// for more informations, check head part of this file. there are some important links out there.
static struct PyModuleDef calculatorModule = {
  PyModuleDef_HEAD_INIT, // head informations for Python C API. It is needed to be first member in this struct !!
  "calculator",  // module name
  NULL, // means that the module does not support sub-interpreters, because it has global state.
  -1,
  functions  // our functions list
};

// runs while initializing and calls module creation function.
PyMODINIT_FUNC PyInit_calculator(void){
  return PyModule_Create(&calculatorModule);
}
