# simple_ed_cpp
A simple exact diagonalization code for Fuzhou Summer School (2023/08/07)

A simple exact diagonalization code solving 1d XXZ model, containing dense/sparse algorithms with/without U(1) symmetry (total Sz conservation).  

Requires Intel MKL as the dense matrix solver. In case of any MKL-related problems on your own computer, you can try the code on your clusters (or the free computing resources provided by Shuguang for this summer school), where the MKL is for sure correctly placed.

The code aims to be written in a direct/naive way, therefore, can be directly compared with the pseudocode in slides.
