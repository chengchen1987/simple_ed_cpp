#include <cstdlib>
#include <cstdio>
#include <iostream> 
#include <fstream>
#include <iomanip>
//#include <cmath>
//#include <cstring>
//#include <vector>
using namespace std;
#include "basis.h" 

Basis::Basis(int _L) :
	L(_L),
	Dim(1 << l_int(L))
{
	State = new l_int[Dim];
	Index = new l_int[1 << (l_int)L];

	for (l_int s = 0; s < (1 << (l_int)L); s++) {
		State[s] = s;
		Index[s] = s;
	}
}

Basis::Basis(int _L, int _nop) :
	L(_L),
	nop(_nop),
	Dim(nchoosek(L, nop))
{
	State = new l_int[Dim];
	Index = new l_int[1 << l_int(L)];

	l_int j = 0;
	int totalcharge;
	for (l_int s = 0; s < (1 << (l_int)L); s++) {
		totalcharge = numOfBit1(s);
		if (totalcharge == nop) {
			State[j] = s;
			Index[s] = j;
			j++;
		}
	}
	// debug
	// PrintStates();
}

Basis::~Basis()
{
	delete[] State;
}

l_int Basis::get_state(const l_int& s)
{
	return State[s];
}

l_int Basis::get_index(const l_int& s)
{
	return Index[s];
}

/*
l_int Basis::get_index(const l_int& s)
{
	if (s < 0) return -1;
	l_int bmin = 0, bmax = Dim;
	l_int b;
	while (bmin <= bmax) {
		b = bmin + (bmax - bmin) / 2;
		l_int aux = State[b];
		if (s == aux) {
			return b;
		}
		else if (s > aux) {
			bmin = b + 1;
		}
		else if (s < aux) {
			bmax = b - 1;
		}
	}
	return -1;
}
*/

void Basis::PrintStates()
{
	ofstream ofd("states");
	for (l_int i = 0; i < Dim; i++) {
		ofd << setw(4) << i << setw(8) << State[i];
		ofd << "  ";
		for (int j = 0; j < L; j++) {
			ofd << ((State[i] >> (L - 1 - j)) & 1);
		}
		ofd << endl;
	}
	ofd.close();
}

