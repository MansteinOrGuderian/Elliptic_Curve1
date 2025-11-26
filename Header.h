#pragma once
#include <iostream>
#include <vector>
#include <set>
#include <iomanip>
#include <string>

// Field F_p, p -- prime, p \neq 2, 3
struct Point {
	unsigned int x;
	unsigned int y;
	Point() { x = static_cast<unsigned>(0), y = static_cast<unsigned>(0); }; // infinity point.. bad decision
	Point(unsigned x, unsigned y) { this->x = x, this->y = y; };
	Point(const Point& point) { x = point.x, y = point.y; };
	~Point() {};
	bool operator==(const Point P) { return (this->x == P.x ? true : false); };
	friend std::ostream& operator<<(std::ostream&, const Point&);
};

struct Elliptic_Curve { // y^2 = x^3 + ax + b; Smooth condition : 4 a^3 + 27 b^2 \neq 0 -> no singular points (cusps or self-intersections)
	
	// Constants to change
	unsigned PRIME_P = 631;
	unsigned a = 30;
	unsigned b = 34; // can't be 0, cause (0, 0) used as infinity point

	std::vector<Point> points_collection;

	Point sum_neq_points(Point, Point);
	Point double_Point(Point);
	Point minus_Point(Point);
	void get_all_curve_points();
	unsigned int Millers_algorithm(Point);

	Point straight_line(Point, Point);
	Point tangent_line(Point);
	Point vertical_line(Point);
};

