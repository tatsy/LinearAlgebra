#ifndef _TRIPLET_H_
#define _TRIPLET_H_

namespace funopt {
	class Triplet {
	public:
		int i, j;
		double val;

		Triplet() :
			i(0), j(0), val(0.0)
		{
		}

		Triplet(int i_, int j_, double val_) :
			i(i_), j(j_), val(val_)
		{
		}

		Triplet(const Triplet& t) :
			i(t.i), j(t.j), val(t.val)
		{
		}

		Triplet& operator=(const Triplet& t)
		{
			this->i = t.i;
			this->j = t.j;
			this->val = t.val;
			return *this;
		}

		bool operator<(const Triplet& t) const
		{
			if(j != t.j) return j < t.j;
			if(i != t.i) return i < t.i;
			return val < t.val;
		}

		bool operator>(const Triplet& t) const
		{
			if(j != t.j) return j > t.j;
			if(i != t.i) return i > t.i;
			return val > t.val;
		}
	};
}

#endif
