#pragma once

#include "common.h"

class Basis
{
private:
	int L;
	int nop;
	l_int Dim;
	l_int* State;
	l_int* Index;

public:
	Basis(int _L);
    Basis(int _L, int _nop);
	~Basis();

	int get_L(void) { return L; }
	int get_nop(void) { return nop; }
	l_int get_Dim(void) { return Dim; }

	l_int get_state(const l_int& index);
	l_int get_index(const l_int& s);

	// debug
	void PrintStates();
};


