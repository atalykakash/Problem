//  Copyright © 2017 AAkash. All rights reserved.

#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

const double eps = 1e-9;

#define square(a) a*a
#define den(a, b, c, d) (a*d - b*c)

struct Point {
    double x, y;
    
    bool operator< (const Point &point) const {
        return x < point.x - eps || (abs(x - point.x) < eps && y < point.y - eps);
    }
};

struct Circle {
    double x, y, r;
};

struct Line {
    
    double a, b, c;
    Point p1, p2;
    
    Line() {}
    Line(Point point1, Point point2) {
        a = point1.y - point2.y;
        b = point2.x - point1.x;
        c = -a * point1.x - b * point1.y;
        p1 = point1;
        p2 = point2;
        normalize();
    }
    
    void normalize() {
        double z = sqrt(square(a) + square(b));
        if(abs(z) > eps) {
            a /= z, b /= z, c /= z;
        }
    }
    
    double dist(Point p) {
        return a * p.x + b * p.y + c;
    }
};

bool cmp(Point p1, Point p2) {
    return p1.x < p2.x || (p1.x == p2.x && p1.y < p2.y);
}

bool cw(Point p1, Point p2, Point p3) {
    return p1.x * (p2.y - p3.y) + p2.x * (p3.y - p1.y) + p3.x * (p1.y - p2.y) < 0;
}

bool ccw(Point p1, Point p2, Point p3) {
    return p1.x * (p2.y - p3.y) + p2.x * (p3.y - p1.y) + p3.x * (p1.y - p2.y) > 0;
}


class Intersection {
    
    vector<Point> intersects;
    vector<Point> convex_hull_points;
    
public:
    
    bool belong(double l, double r, double a) {
        return min(l, r) <= a + eps && a <= max(l, r) + eps;
    }
    
    bool intersect(double a, double b, double c, double d) {
        if(a > b) swap(a, b);
        if(c > d) swap(c, d);
        return max(a, c) <= min(b, d) + eps;
    }
    
    void find_intersects_line_to_line(Point point1, Point point2, Point point3, Point point4) {
        
        if(!intersect(point1.x, point2.x, point3.x, point4.x) || !intersect(point1.y, point2.y, point3.y, point4.y))
            return;
        
        Line line1(point1, point2);
        Line line2(point3, point4);
        
        double den = den(line1.a, line1.b, line2.a, line2.b);
        
        if(abs(den) < eps) return;
        
        Point point;
        
        point.x = -den(line1.c, line1.b, line2.c, line2.b) / den;
        point.y = -den(line1.a, line1.c, line2.a, line2.c) / den;
        
        if(belong(point1.x, point2.x, point.x) &&
           belong(point1.y, point2.y, point.y) &&
           belong(point3.x, point4.x, point.x)  &&
           belong(point3.y, point4.y, point.y)) {
            
            intersects.push_back(point);
        }
    }
    
    void find_intersects_circle_to_circle(Circle &circle1, Circle &circle2) {
        
        Point point;
        
        double d;
        double a;
        double h;
        
        d = sqrt(square((circle2.x-circle1.x)) + square((circle2.y - circle1.y)));
        
        if(d > circle1.r + circle2.r || d < abs(circle1.r - circle2.r) || (d == 0 && circle1.r == circle2.r))
            return;
        
        a = (square(circle1.r) - square(circle2.r) + square(d)) / (2 * d);
        h = sqrt(square(circle1.r) - square(a));
        
        point.x = circle1.x + a * (circle2.x - circle1.x) / d;
        point.y = circle1.y + a * (circle2.y - circle1.y) / d;
        
        Point p1, p2;
        
        p1.x = point.x + h * (circle2.y - circle1.y) / d;
        p1.y = point.y - h * (circle2.x - circle1.x) / d;
        
        intersects.push_back(p1);
        
        if(abs(a) == circle1.r)
            return;
        
        p2.x = point.x - h * (circle2.y - circle1.y) / d;
        p2.y = point.y + h * (circle2.x - circle1.x) / d;
        
        intersects.push_back(p2);
    }
    
    void find_intersects_circle_to_line(Circle &circle, Line &line) {
        
        Point point;
        
        line.c = line.c + line.a * circle.x + line.b * circle.y;
        
        point.x = (-line.a * line.c)/(square(line.a) + square(line.b));
        point.y = (-line.b * line.c)/(square(line.a) + square(line.b));
        
        if(square(line.c) > square(circle.r) * (square(line.a) + square(line.b)) + eps) {
            return;
        } else if(abs(square(line.c) - square(circle.r) * (square(line.a) + square(line.b))) < eps) {
            point.x = point.x + circle.x;
            point.y = point.y + circle.y;
            intersects.push_back(point);
        } else {
            double distance = square(circle.r) - square(line.c)/(square(line.a) + square(line.b));
            double mult = sqrt(distance/(square(line.a) + square(line.b)));
            Point p1, p2;
            p1.x = point.x + line.b * mult + circle.x;
            p1.y = point.y - line.a * mult + circle.y;
            p2.x = point.x - line.b * mult + circle.x;
            p2.y = point.y + line.a * mult + circle.y;
            if(line.p1.x - p1.x < eps && p1.x - line.p2.x < eps) {
                intersects.push_back(p1);
            }
            if(line.p1.x - p2.x < eps && p2.x - line.p2.x < eps) {
                intersects.push_back(p2);
            }
        }
    }
    
    void convex_hull(vector<Point> &points) {
        
        if(points.size() == 1)
            return;
        
        sort(points.begin(), points.end(), cmp);
        
        Point p1 = points[0];
        Point p2 = points.back();
        
        vector<Point> up, down;
        
        up.push_back(p1);
        down.push_back(p1);
        
        for(size_t i = 1; i<points.size(); i++) {
            if(i == points.size()-1 || cw(p1, points[i], p2)) {
                while (up.size()>=2 && !cw(up[up.size()-2], up[up.size()-1], points[i]))
                    up.pop_back();
                up.push_back(points[i]);
            }
            if(i == points.size()-1 || ccw(p1, points[i], p2)) {
                while (down.size()>=2 && !ccw(down[down.size()-2], down[down.size()-1], points[i]))
                    down.pop_back();
                down.push_back(points[i]);
            }
        }
    
        for(size_t i = 0; i<down.size(); i++) {
            convex_hull_points.push_back(down[i]);
        }
        
        for(size_t i = up.size()-2; i>0; i--) {
            convex_hull_points.push_back(up[i]);
        }
    }
    
    double convex_hull_area(vector<Point> &points) {
        
        double area = 0;
        
        for(size_t i = 0; i<points.size()-1; i++) {
            area += (points[i].x * points[i+1].y - points[i].y * points[i+1].x);
        }
        
        area += points[0].y * points[points.size()-1].x - points[0].x * points[points.size()-1].y;
        
        return area/2;
    }
    
    void print_result() {
        
        cout<<intersects.size()<<endl;
        
        for(int i = 0; i<intersects.size(); i++) {
            printf("%.4f %.4f\n", intersects[i].x, intersects[i].y);
        }
        
        convex_hull(intersects);
        
        cout<<convex_hull_points.size()<<endl;
        
        for(int i = 0; i<convex_hull_points.size(); i++) {
            printf("%.4f %.4f\n", convex_hull_points[i].x, convex_hull_points[i].y);
        }
        
        cout<<convex_hull_area(convex_hull_points);
    }
};


int main() {
    
    int n;
    
    cin>>n;
    cin.ignore();
    
    char letter;
    vector<Circle> circles;
    vector<Line> lines;
    
    Circle circle;
    Line line;
    Point p1, p2;
    
    for(int i = 0; i<n; i++) {
        cin>>letter;
        if(letter == 'C') {
            cin>>circle.x>>circle.y>>circle.r;
            circles.push_back(circle);
        } else{
            cin>>p1.x>>p1.y>>p2.x>>p2.y;
            Line line(p1, p2);
            lines.push_back(line);
        }
    }
    
    Intersection intersection;
    
    for(int i = 0; i<lines.size(); i++) {
        for(int j = 0; j<circles.size(); j++) {
            intersection.find_intersects_circle_to_line(circles[j], lines[i]);
        }
    }
    
    for(int i = 0; i<lines.size(); i++) {
        for(int j = i+1; j<lines.size(); j++) {
            intersection.find_intersects_line_to_line(lines[i].p1, lines[i].p2, lines[j].p1, lines[j].p2);
        }
    }
    
    for(int i = 0; i<circles.size(); i++) {
        for(int j = i+1; j<circles.size(); j++) {
            intersection.find_intersects_circle_to_circle(circles[i], circles[j]);
        }
    }
    
    intersection.print_result();
    
    return 0;
}

