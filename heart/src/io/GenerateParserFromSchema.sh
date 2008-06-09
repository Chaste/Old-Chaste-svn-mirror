#!/bin/bash

xsd cxx-tree --hxx-suffix .hpp --cxx-suffix .cpp --cxx-prologue "#define COVERAGE_IGNORE" --cxx-epilogue "#undef COVERAGE_IGNORE" --hxx-prologue "#define COVERAGE_IGNORE" --hxx-epilogue "#undef COVERAGE_IGNORE" --root-element ChasteParameters ChasteParameters.xsd
