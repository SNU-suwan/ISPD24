#pragma once
#include <cassert>
#include <iostream>

template <typename T>
class Point {
public:
    Point(T xVal = 0, T yVal = 0) noexcept
    : x(xVal), y(yVal) {}

    // operators
    const T& operator[](const unsigned d) const {
        assert(d == 0 || d == 1);
        return (d == 0 ? x : y);
    }
    Point operator+(const Point& p) const { return Point(x + p.x, y + p.y); }
    Point operator+(const T t) const { return Point(x + t, y + t); }
    Point operator-(const Point& p) const { return Point(x - p.x, y - p.y); }
    Point operator-(const T t) const { return Point(x - t, y - t); }
    Point operator/(const T t) const { return Point(x / t, y / t); }
    Point& operator+=(const Point& p) {
        x += p.x; y += p.y;
        return *this;
    }
    Point& operator+=(const T t) {
        x += t; y += t;
        return *this;
    }
    Point& operator-=(const Point& p) {
        x -= p.x; y -= p.y;
        return *this;
    }
    Point& operator-=(const T t) {
        x -= t; y -= t;
        return *this;
    }
	Point& operator*=(const T t) {
        x *= t; y *= t;
        return *this;
    }
    bool operator==(const Point& p) const { return x == p.x && y == p.y; }
    bool operator!=(const Point& p) const { return !(*this == p); }

    friend inline std::ostream& operator<<(std::ostream& os, const Point& p) {
        os << "(" << p.x << ", " << p.y << ")";
        return os;
    }

public:
    T x, y;
};


template <typename T>
class BBox {
public:
    BBox() {}
    BBox(const Point<T>& lb_, const Point<T>& rt_) {
        lb = lb_; rt = rt_;
        if (lb.x > rt.x) std::swap(lb.x, rt.x);
        if (lb.y > rt.y) std::swap(lb.y, rt.y);
    }
    BBox(T lb_x, T lb_y, T rt_x, T rt_y) {
        lb.x = lb_x; lb.y = lb_y; rt.x = rt_x; rt.y = rt_y;
        if (lb.x > rt.x) std::swap(lb.x, rt.x);
        if (lb.y > rt.y) std::swap(lb.y, rt.y);
    }

    // operators
    const Point<T>& operator[](const unsigned d) const {
        assert(d == 0 || d == 1);
        return (d == 0 ? lb : rt);
    }
    BBox<T> operator+(const BBox<T>& bbox) {
        return BBox<T>(std::min(lb.x, bbox.lb.x), std::min(lb.y, bbox.lb.y), std::max(rt.x, bbox.rt.x), std::max(rt.y, bbox.rt.y));
    }
    BBox<T> operator+(const Point<T>& p) { return BBox<T>(lb + p, rt + p); }
    BBox<T> operator-(const Point<T>& p) { return BBox<T>(lb - p, rt - p); }
    BBox<T> operator&(const BBox<T>& bbox) {
        auto intersect = [](T b1_1, T b1_2, T b2_1, T b2_2) {
            if (b1_1 <= b2_1 && b2_1 < b1_2 && b1_2 <= b2_2) return Point<T>(b2_1, b1_2);
            else if (b2_1 <= b1_1 && b1_1 < b2_2 && b2_2 <= b1_2) return Point<T>(b1_1, b2_2);
            else if (b1_1 <= b2_1 && b2_1 < b2_2 && b2_2 <= b1_2) return Point<T>(b2_1, b2_2);
            else if (b2_1 <= b1_1 && b1_1 < b1_2 && b1_2 <= b2_2) return Point<T>(b1_1, b1_2);
            else return Point<T>(-1, -1);
        };
        auto xInter = intersect(lb.x, rt.x, bbox.lb.x, bbox.rt.x);
        auto yInter = intersect(lb.y, rt.y, bbox.lb.y, bbox.rt.y);
        return BBox<T>(xInter[0], yInter[0], xInter[1], yInter[1]);
    }
    
	BBox<T> operator*=(const T p) {
		lb*=p;
		rt*=p;
		return *this;
	}
	BBox<T> operator-=(const T p) {
		lb-=p;
		rt-=p;
		return *this;
	}
	bool operator==(const BBox<T>& bbox) { return lb == bbox.lb && rt == bbox.rt; }
    bool operator!=(const BBox<T>& bbox) { return !(*this == bbox); }

    //check overlap
    bool overlap(const BBox<T>& bbox) {
        if(lb.x >= bbox.rt.x || bbox.lb.x >= rt.x) return false;
        if(lb.y >= bbox.rt.y || bbox.lb.y >= rt.y) return false;
        return true;
    }

	bool is_intersect(const BBox<T>& bbox){
		if(lb.x>=bbox.rt.x || rt.x<=bbox.lb.x || rt.y<=bbox.lb.y || lb.y>=bbox.rt.y){
			return false;
		}
		return true;
	}
    T xWidth() const { return rt.x - lb.x; }
    T yWidth() const { return rt.y - lb.y; }
    T getArea() const { return xWidth() * yWidth(); }

    friend inline std::ostream& operator<<(std::ostream& os, const BBox& b) {
        os << "(" << b.lb << ", " << b.rt << ")";
        return os;
    }

public:
    Point<T> lb, rt;

};


