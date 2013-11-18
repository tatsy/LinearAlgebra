#ifndef _MT_RAND_H_
#define _MT_RAND_H_

class MTRand {
public:
	MTRand();
	int randInt(int n) const;
	double randReal() const;
};

#endif
