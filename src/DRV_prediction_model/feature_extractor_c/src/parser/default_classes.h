#ifndef __DEFAULT_CLASS__
#define __DEFAULT_CLASS__

#include<string>
using namespace std;

class point2F;
class point2I;
class int2I;
class float2F;


class point2F{
public:
	point2F(float,float);
	point2F();
	~point2F();
public:
	float xF;
	float yF;
	string containerS;

public:
	bool equalX (const point2F&) const;
	bool equalY (const point2F&) const;
	bool equal (const point2F&) const;
	bool equalX_p2 (const point2F&) const;
	bool equalY_p2 (const point2F&) const;
	bool equal_p2 (const point2F&) const;
	bool operator < (const point2F&) const;
	bool operator == (const point2F&) const;
	bool operator != (const point2F&) const;
	point2F operator = (const point2F&);
	point2F operator += (const point2F&);
	point2F operator -= (const point2F&);
	point2F operator /= (const float&);
	point2F operator *= (const float&);

	
};

class point2I{
public:
	point2I(int,int);
	point2I();
	~point2I();
public:
	int xI;
	int yI;
	string containerS;

public:
	/*
	bool equalX (const point2I&);
	bool equalY (const point2I&);
	bool equal (const point2I&);
	*/
	bool operator <  (const point2I&) const;
	bool operator == (const point2I&) const;
	bool operator != (const point2I&) const;
	point2I operator = (const point2I&);
	point2I operator += (const point2I&);
	point2I operator -= (const point2I&);

	
};

class int2I{
public:
	int2I(int,int);
	int2I();
	~int2I();
public:
	int xI;
	int yI;
};

class float2F{
public:
	float2F(float,float);
	float2F();
	~float2F();
public:
	float xF;
	float yF;
};

// OLD EDGE
/* 
class edge{
public:
	edge(bool,string,point2F*);
	~edge();
public:
	bool direction;
	string layerS; //layer or direction
	point2F endPointsP22[2];
};
*/


#endif
