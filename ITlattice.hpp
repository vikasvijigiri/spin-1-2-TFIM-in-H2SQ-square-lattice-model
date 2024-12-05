#ifndef _ITLATTICE_HPP_DEFINED_
#define _ITLATTICE_HPP_DEFINED_
#include <random>
#include "ITvariables.hpp"
#pragma once
using std::exp;

class ITlattice
{
	private:
	int spin;
	bool in_cluster;
	int cluster;
	public:
		ITlattice()
		{
			spin = 1;
			in_cluster = false;
			cluster = -1;
		}
	ITlattice(int s)
	{
		spin = s;
		in_cluster = false;
	}
	void set_cl(int c)
	{
		in_cluster = true;
		cluster = c;
	}
	double J, Kx, size;
	bool bonded(void)
	{
		return in_cluster;
	}
	void reset_bond(void)
	{
		in_cluster = false;
		cluster = -1;
	}
	int S(void)
	{
		return spin;
	}
	void set_S(int s)
	{
		spin = s;
	}
	int cl(void)
	{
		return cluster;
	}
	void flip(void)
	{
		spin *= -1;
	}
	double bond_prob(ITlattice s2, double beta, double J)
	{
		return (1 - exp(-1 *beta *J)) *spin *s2.S();
	}
	double ice_rule(ITlattice s2, double beta, double J)
	{
		return ((1 - spin *s2.S()) / 2);
	}
};
#endif
