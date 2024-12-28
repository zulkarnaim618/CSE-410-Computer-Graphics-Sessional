#include<bits/stdc++.h>
#ifdef __linux__
    #include<GL/glut.h>
#elif WIN32
    #include <windows.h>
    #include <glut.h>
#endif
using namespace std;

#ifndef _HELPER_
#define _HELPER_

class Vector {
public:
    GLdouble x;
    GLdouble y;
    GLdouble z;
    Vector() {
        this->x = 0;
        this->y = 0;
        this->z = 0;
    }
    Vector(GLdouble x, GLdouble y, GLdouble z) {
        this->x = x;
        this->y = y;
        this->z = z;
    }
    Vector(Vector* v) {
        this->x = v->x;
        this->y = v->y;
        this->z = v->z;
    }
    void crossProduct(Vector* b, Vector* ans) {
        ans->x = this->y*b->z-this->z*b->y;
        ans->y = this->z*b->x-this->x*b->z;
        ans->z = this->x*b->y-this->y*b->x;
    }
    void normalize() {
        double l = sqrt(x*x+y*y+z*z);
        x = x/l;
        y = y/l;
        z = z/l;
    }
};

ostream& operator<<(ostream& out, Vector& c) {
    out<<"("<<c.x<<","<<c.y<<","<<c.z<<")";
    return out;
}

class Camera {
public:
    GLdouble refDistance;
    Vector* eye;
    Vector* front;
    Vector* up;
    Vector* right;
    Camera() {
        this->eye = new Vector();
        this->front = new Vector();
        this->up = new Vector();
        this->right = new Vector();
        this->refDistance = 0;
    }
    Camera(Vector* eye, Vector* front, Vector* up, Vector* right, GLdouble refDistance) {
        this->eye = eye;
        this->front = front;
        this->up = up;
        this->right = right;
        this->refDistance = refDistance;
    }
    ~Camera() {
        delete eye;
        delete front;
        delete up;
        delete right;
    }
};

class Ball {
public:
    GLdouble radius;
    Vector* center;
    Vector* front;
    Vector* up;
    Vector* right;
    Vector* direction;
    Ball() {
        this->center = new Vector();
        this->front = new Vector();
        this->up = new Vector();
        this->right = new Vector();
        this->direction = new Vector();
        this->radius = 0;
    }
    Ball(Vector* center, Vector* front, Vector* up, Vector* right, Vector* direction, GLdouble radius) {
        this->center = center;
        this->front = front;
        this->up = up;
        this->right = right;
        this->direction = direction;
        this->radius = radius;
    }
    ~Ball() {
        delete center;
        delete front;
        delete up;
        delete right;
        delete direction;
    }
};

#endif