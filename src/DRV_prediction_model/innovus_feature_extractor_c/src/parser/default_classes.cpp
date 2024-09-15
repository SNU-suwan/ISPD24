#include <cmath>
#include "global_variables.h"
#include "default_classes.h"

using namespace std;


//point2F
point2F::point2F(float _xF, float _yF){
	xF=_xF;
	yF=_yF;
	containerS.clear();
}
point2F::point2F(){
	xF=0;
	yF=0;
	containerS.clear();
}
point2F::~point2F(){
}
bool point2F::equalX(const point2F& _pointP2) const{
	if(abs(xF-_pointP2.xF)<PRECISION) return true;
	else return false;
	return false;
}
bool point2F::equalY(const point2F& _pointP2) const{
	if(abs(yF-_pointP2.yF)<PRECISION) return true;
	else return false;
	return false;
}
bool point2F::equal(const point2F& _pointP2) const {
	if (equalX(_pointP2) && equalY(_pointP2)) return true;
	else return false;
	return false;
}
bool point2F::equalX_p2(const point2F& _pointP2) const{
	if(abs(xF-_pointP2.xF)<0.009) return true;
	else return false;
	return false;
}
bool point2F::equalY_p2(const point2F& _pointP2) const{
	if(abs(yF-_pointP2.yF)<0.009) return true;
	else return false;
	return false;
}
bool point2F::equal_p2(const point2F& _pointP2) const {
	if (equalX_p2(_pointP2) && equalY_p2(_pointP2)) return true;
	else return false;
	return false;
}


bool point2F::operator < (const point2F& _pointP2) const{
	return xF < _pointP2.xF;	
}
bool point2F::operator == (const point2F& _pointP2) const{
	if(xF==_pointP2.xF && yF==_pointP2.yF){
		return true;
	}
	return false;
}
bool point2F::operator != (const point2F& _pointP2) const{
	if(xF==_pointP2.xF && yF==_pointP2.yF){
		return false;
	}
	return true;
}
point2F point2F::operator = (const point2F& _pointP2){
	xF=_pointP2.xF;
	yF=_pointP2.yF;
	return *this;
}
point2F point2F::operator += (const point2F& _pointP2){
	xF+=_pointP2.xF;
	yF+=_pointP2.yF;
	return *this;
}
point2F point2F::operator -= (const point2F& _pointP2){
	xF-=_pointP2.xF;
	yF-=_pointP2.yF;
	return *this;
}
point2F point2F::operator /= (const float& _divider){
	xF/=_divider;
	yF/=_divider;
	return *this;
}
point2F point2F::operator *= (const float& _multiplier){
	xF*=_multiplier;
	yF*=_multiplier;
	return *this;
}
//point2I
point2I::point2I(int _xI, int _yI){
	xI=_xI;
	yI=_yI;
}
point2I::point2I(){
	xI=0;
	yI=0;
}
point2I::~point2I(){
}
/*
bool point2I::equalX(const point2I& _pointP2){
	if(abs(xI-_pointP2.xI)<PRECISION) return true;
	else return false;
	return false;
}
bool point2I::equalY(const point2I& _pointP2){
	if(abs(yI-_pointP2.yI)<PRECISION) return true;
	else return false;
	return false;
}
bool point2I::equal(const point2I& _pointP2){
	if (equalX(_pointP2) && equalY(_pointP2)) return true;
	else return false;
	return false;
}
*/
bool point2I::operator <  (const point2I& _pointP2) const{
	if (xI==_pointP2.xI) return yI<_pointP2.yI;
	return xI < _pointP2.xI;	
}
bool point2I::operator == (const point2I& _pointP2) const{
	if(xI==_pointP2.xI && yI==_pointP2.yI){
		return true;
	}
	return false;
}
bool point2I::operator != (const point2I& _pointP2) const{
	if(xI==_pointP2.xI && yI==_pointP2.yI){
		return false;
	}
	return true;
}
point2I point2I::operator = (const point2I& _pointP2){
	xI=_pointP2.xI;
	yI=_pointP2.yI;
	return *this;
}
point2I point2I::operator += (const point2I& _pointP2){
	xI+=_pointP2.xI;
	yI+=_pointP2.yI;
	return *this;
}
point2I point2I::operator -= (const point2I& _pointP2){
	xI-=_pointP2.xI;
	yI-=_pointP2.yI;
	return *this;
}


//int2I
int2I::int2I(int _xI, int _yI){
	xI=_xI;
	yI=_yI;
}
int2I::int2I(){
	xI=0;
	yI=0;
}
int2I::~int2I(){
}

//float2F
float2F::float2F(float _xF, float _yF){
	xF=_xF;
	yF=_yF;
}
float2F::float2F(){
	xF=0;
	yF=0;
}
float2F::~float2F(){
}


// OLD edge
/*
edge::edge(bool _direction,string _layerS,point2F* _pointP22){
	direction=_direction;
	layerS=_layerS;
	endPointsP22[0]=_pointP22[0];
	endPointsP22[1]=_pointP22[1];
}
edge::~edge(){
}*/


