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
    double x;
    double y;
    double z;
    Vector() {
        this->x = 0;
        this->y = 0;
        this->z = 0;
    }
    Vector(double x, double y, double z) {
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
    Vector scalarProduct(double a) {
        Vector v;
        v.x = a*this->x;
        v.y = a*this->y;
        v.z = a*this->z;
        return v;
    }
    double dotProduct(Vector b) {
        double ans = 0;
        ans += b.x*this->x;
        ans += b.y*this->y;
        ans += b.z*this->z;
        return ans;
    }
    Vector crossProduct(Vector b) {
        Vector ans;
        ans.x = this->y*b.z-this->z*b.y;
        ans.y = this->z*b.x-this->x*b.z;
        ans.z = this->x*b.y-this->y*b.x;
        return ans;
    }
    Vector operator+(Vector b) {
        Vector ans;
        ans.x = this->x + b.x;
        ans.y = this->y + b.y;
        ans.z = this->z + b.z;
        return ans;
    }
    Vector operator-(Vector b) {
        Vector ans;
        ans.x = this->x - b.x;
        ans.y = this->y - b.y;
        ans.z = this->z - b.z;
        return ans;
    }
};

ostream& operator<<(ostream& out, Vector& c) {
    out<<"("<<c.x<<","<<c.y<<","<<c.z<<")";
    return out;
}

double det(double a[3][3]) {
    double ans = 0;
    ans += a[0][0]*(a[1][1]*a[2][2]-a[1][2]*a[2][1]);
    ans -= a[0][1]*(a[1][0]*a[2][2]-a[2][0]*a[1][2]);
    ans += a[0][2]*(a[1][0]*a[2][1]-a[2][0]*a[1][1]);
    return ans;
}

class Color {
public:
    double r,g,b;
    Color() {
        r = g = b = 0;
    }
    Color(double r, double g, double b) {
        this->r = r;
        this->g = g;
        this->b = b;
    }
};

class Camera {
public:
    double refDistance;
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
    Camera(Vector* eye, Vector* front, Vector* up, Vector* right, double refDistance) {
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

class Ray {
public:
    Vector start;
    Vector dir;     //normalize
    Ray(Vector start, Vector dir) {
        this->start = start;
        this->dir = dir;
        this->dir.normalize();
    }
};

class PointLight {
public:
    Vector light_pos;
    Color color;
    PointLight() {

    }
    PointLight(Vector light_pos, Color color) {
        this->light_pos = light_pos;
        this->color = color;
    }
    void draw() {
        glColor3f(color.r,color.g,color.b);
        glPointSize(4);
        glBegin(GL_POINTS);{
            glVertex3f(light_pos.x,light_pos.y,light_pos.z);
        }glEnd();
    }
};

class SpotLight {
public:
    PointLight point_light;
    Vector light_direction;
    double cutoff_angle;
    SpotLight(Vector light_pos, Vector light_direction, double cutoff_angle, Color color) {
        PointLight p(light_pos,color);
        this->point_light = p;
        this->light_direction = light_direction;
        this->cutoff_angle = cutoff_angle;
    }
    void draw() {
        glColor3f(point_light.color.r,point_light.color.g,point_light.color.b);
        glPointSize(7);
        glBegin(GL_POINTS);{
            glVertex3f(point_light.light_pos.x,point_light.light_pos.y,point_light.light_pos.z);
        }glEnd();
    }
};

double epsilon = 1e-6;

class Object;

extern vector<Object*> objects;
extern vector<PointLight*> pointLights;
extern vector<SpotLight*> spotLights;
extern int recursionLevel;

class Object {
public:
    Vector reference_point;
    double height, width, length;
    Color color;
    double coEfficients[4];     //ambient, diffuse, specular, reflection
    int shine;
    Object() {
        
    }
    Object(Vector ref, double height, double width, double length) {
        this->reference_point = ref;
        this->height = height;
        this->width = width;
        this->length = length;
    }
    virtual void draw() {

    }
    void setColor(double r, double g, double b) {
        Color c(r,g,b);
        this->color = c;
    }
    void setShine(int shine) {
        this->shine = shine;
    }
    void setCoEfficients(double ambient, double diffuse, double specular, double reflection) {
        coEfficients[0] = ambient;
        coEfficients[1] = diffuse;
        coEfficients[2] = specular;
        coEfficients[3] = reflection;
    }
    virtual Color getColorAt(Vector intersectionPoint) {
        return this->color;
    }
    virtual double intersect(Ray ray, Color& color, int level) {
        return -1;
    }
};

class Sphere: public Object {
public:
    Sphere(Vector center, double radius) {
        reference_point = center;
        length = radius;
    }
    void draw() {
        int stackCount,sectorCount;
        stackCount = 30;
        sectorCount = 30;
        double x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4;
        Vector front(1,0,0),right(0,1,0),up(0,0,1);
        for(int i=0;i<stackCount;i++) {
            for (int j=0;j<sectorCount;j++) {
                glColor3f(color.r,color.g,color.b);
                x1 = length*cos((90-i*(180.0/stackCount))*(M_PI/180.0))*cos(j*(360.0/sectorCount)*(M_PI/180.0));
                y1 = length*cos((90-i*(180.0/stackCount))*(M_PI/180.0))*sin(j*(360.0/sectorCount)*(M_PI/180.0));
                z1 = length*sin((90-i*(180.0/stackCount))*(M_PI/180.0));
                x2 = length*cos((90-i*(180.0/stackCount))*(M_PI/180.0))*cos((j+1)*(360.0/sectorCount)*(M_PI/180.0));
                y2 = length*cos((90-i*(180.0/stackCount))*(M_PI/180.0))*sin((j+1)*(360.0/sectorCount)*(M_PI/180.0));
                z2 = length*sin((90-i*(180.0/stackCount))*(M_PI/180.0));
                x3 = length*cos((90-(i+1)*(180.0/stackCount))*(M_PI/180.0))*cos(j*(360.0/sectorCount)*(M_PI/180.0));
                y3 = length*cos((90-(i+1)*(180.0/stackCount))*(M_PI/180.0))*sin(j*(360.0/sectorCount)*(M_PI/180.0));
                z3 = length*sin((90-(i+1)*(180.0/stackCount))*(M_PI/180.0));
                x4 = length*cos((90-(i+1)*(180.0/stackCount))*(M_PI/180.0))*cos((j+1)*(360.0/sectorCount)*(M_PI/180.0));
                y4 = length*cos((90-(i+1)*(180.0/stackCount))*(M_PI/180.0))*sin((j+1)*(360.0/sectorCount)*(M_PI/180.0));
                z4 = length*sin((90-(i+1)*(180.0/stackCount))*(M_PI/180.0));
                glBegin(GL_QUAD_STRIP);{
                    glVertex3f(reference_point.x+x1*front.x+z1*up.x+y1*right.x,reference_point.y+x1*front.y+z1*up.y+y1*right.y,reference_point.z+x1*front.z+z1*up.z+y1*right.z);
                    glVertex3f(reference_point.x+x2*front.x+z2*up.x+y2*right.x,reference_point.y+x2*front.y+z2*up.y+y2*right.y,reference_point.z+x2*front.z+z2*up.z+y2*right.z);
                    glVertex3f(reference_point.x+x3*front.x+z3*up.x+y3*right.x,reference_point.y+x3*front.y+z3*up.y+y3*right.y,reference_point.z+x3*front.z+z3*up.z+y3*right.z);
                    glVertex3f(reference_point.x+x4*front.x+z4*up.x+y4*right.x,reference_point.y+x4*front.y+z4*up.y+y4*right.y,reference_point.z+x4*front.z+z4*up.z+y4*right.z);
                }glEnd();
            }
        }
    }
    double intersect(Ray ray, Color& color, int level) {
        ray.start = ray.start - reference_point;
        double tp = -(ray.start.dotProduct(ray.dir));
        double t = -1;
        if (ray.start.dotProduct(ray.start)>length*length && tp<0) {
            t = -1;
        }
        else {
            double d = ray.start.dotProduct(ray.start)-tp*tp;
            if (d>length*length) {
                t = -1;
            }
            else {
                double tt = length*length - d;
                if (ray.start.dotProduct(ray.start)>length*length) {
                    t = tp - sqrt(tt);
                }
                else {
                    t = tp + sqrt(tt);
                }
            }
        }
        ray.start = ray.start + reference_point;
        if (level==0) return t;
        Vector intersectionPoint = ray.start + ray.dir.scalarProduct(t);
        Color intersectionPointColor = getColorAt(intersectionPoint);
        color.r = intersectionPointColor.r*coEfficients[0];
        color.g = intersectionPointColor.g*coEfficients[0];
        color.b = intersectionPointColor.b*coEfficients[0];
        Vector normal = intersectionPoint-reference_point;
        normal.normalize();
        if (ray.dir.dotProduct(normal)>0) {
            normal.x = -normal.x;
            normal.y = -normal.y;
            normal.z = -normal.z;
        }
        for (int i=0;i<pointLights.size();i++) {
            Ray lightRay(pointLights[i]->light_pos,intersectionPoint-pointLights[i]->light_pos);
            Color dummyColor;
            Vector d = intersectionPoint - pointLights[i]->light_pos;
            double dis = sqrt(d.x*d.x+d.y*d.y+d.z*d.z);
            bool obscured = false;
            for (int j=0;j<objects.size();j++) {
                double tVal = objects[j]->intersect(lightRay,dummyColor,0);
                if (tVal>0 && tVal+epsilon<dis) {
                    obscured = true;
                    break;
                }
            }
            if (!obscured) {
                double lambert = max(-lightRay.dir.dotProduct(normal),0.0);
                Vector reflectedRayDir = lightRay.dir-normal.scalarProduct(2*normal.dotProduct(lightRay.dir));
                double phong = max(-reflectedRayDir.dotProduct(ray.dir),0.0);
                color.r += pointLights[i]->color.r*coEfficients[1]*lambert*intersectionPointColor.r;
                color.g += pointLights[i]->color.g*coEfficients[1]*lambert*intersectionPointColor.g;
                color.b += pointLights[i]->color.b*coEfficients[1]*lambert*intersectionPointColor.b;
                color.r += pointLights[i]->color.r*coEfficients[2]*pow(phong,shine)*intersectionPointColor.r;
                color.g += pointLights[i]->color.g*coEfficients[2]*pow(phong,shine)*intersectionPointColor.g;
                color.b += pointLights[i]->color.b*coEfficients[2]*pow(phong,shine)*intersectionPointColor.b;
            }
        }
        for (int i=0;i<spotLights.size();i++) {
            Ray lightRay(spotLights[i]->point_light.light_pos,intersectionPoint-spotLights[i]->point_light.light_pos);
            Vector e = spotLights[i]->light_direction;
            e.normalize();
            double ans = e.dotProduct(lightRay.dir);
            if ((acos(ans)*(180.0/M_PI)) <= spotLights[i]->cutoff_angle) {
                Color dummyColor;
                Vector d = intersectionPoint - spotLights[i]->point_light.light_pos;
                double dis = sqrt(d.x*d.x+d.y*d.y+d.z*d.z);
                bool obscured = false;
                for (int j=0;j<objects.size();j++) {
                    double tVal = objects[j]->intersect(lightRay,dummyColor,0);
                    if (tVal>0 && tVal+epsilon<dis) {
                        obscured = true;
                        break;
                    }
                }
                if (!obscured) {
                    double lambert = max(-lightRay.dir.dotProduct(normal),0.0);
                    Vector reflectedRayDir = lightRay.dir-normal.scalarProduct(2*normal.dotProduct(lightRay.dir));
                    double phong = max(-reflectedRayDir.dotProduct(ray.dir),0.0);
                    color.r +=  spotLights[i]->point_light.color.r*coEfficients[1]*lambert*intersectionPointColor.r;
                    color.g +=  spotLights[i]->point_light.color.g*coEfficients[1]*lambert*intersectionPointColor.g;
                    color.b +=  spotLights[i]->point_light.color.b*coEfficients[1]*lambert*intersectionPointColor.b;
                    color.r +=  spotLights[i]->point_light.color.r*coEfficients[2]*pow(phong,shine)*intersectionPointColor.r;
                    color.g +=  spotLights[i]->point_light.color.g*coEfficients[2]*pow(phong,shine)*intersectionPointColor.g;
                    color.b +=  spotLights[i]->point_light.color.b*coEfficients[2]*pow(phong,shine)*intersectionPointColor.b;
                }
            }
        }
        if (level>=recursionLevel) {
            return t;
        }
        Ray reflectedRay(intersectionPoint,ray.dir-normal.scalarProduct(2*normal.dotProduct(ray.dir)));
        reflectedRay.start = reflectedRay.start + reflectedRay.dir.scalarProduct(epsilon);
        double t0=-1, tMin;
        Object* obj = nullptr;
        Color dummyColor;
        for (int k=0;k<objects.size();k++) {
            t0 = objects[k]->intersect(reflectedRay,dummyColor,0);
            if (t0>0) {
                if (t0<tMin || obj==nullptr) {
                    tMin = t0;
                    obj = objects[k];
                }
            }
        }
        if (obj!=nullptr) {
            Color c;
            obj->intersect(reflectedRay,c,level+1);
            color.r += c.r*coEfficients[3];
            color.g += c.g*coEfficients[3];
            color.b += c.b*coEfficients[3];
        }
        return t;
    }
};

class Triangle: public Object {
public:
    Vector a,b,c;
    Triangle(Vector a, Vector b, Vector c) {
        this->a = a;
        this->b = b;
        this->c = c;
    }
    void draw() {
        glColor3f(color.r,color.g,color.b);
        glBegin(GL_TRIANGLES);{
            glVertex3f(a.x, a.y, a.z);
            glVertex3f(b.x, b.y, b.z);
            glVertex3f(c.x, c.y, c.z);
        }glEnd();
    }
    double intersect(Ray ray, Color& color, int level) {
        double A[3][3] = {
            {a.x-b.x,a.x-c.x,ray.dir.x},
            {a.y-b.y,a.y-c.y,ray.dir.y},
            {a.z-b.z,a.z-c.z,ray.dir.z}
        };
        double Abeta[3][3] = {
            {a.x-ray.start.x,a.x-c.x,ray.dir.x},
            {a.y-ray.start.y,a.y-c.y,ray.dir.y},
            {a.z-ray.start.z,a.z-c.z,ray.dir.z}
        };
        double Agamma[3][3] = {
            {a.x-b.x,a.x-ray.start.x,ray.dir.x},
            {a.y-b.y,a.y-ray.start.y,ray.dir.y},
            {a.z-b.z,a.z-ray.start.z,ray.dir.z}
        };
        double At[3][3] = {
            {a.x-b.x,a.x-c.x,a.x-ray.start.x},
            {a.y-b.y,a.y-c.y,a.y-ray.start.y},
            {a.z-b.z,a.z-c.z,a.z-ray.start.z}
        };
        double detA = det(A);
        double t = -1;
        if (detA==0) {
            t = -1;
        }
        else {
            double beta = det(Abeta)/detA;
            double gamma = det(Agamma)/detA;
            t = det(At)/detA;
            if (!(beta+gamma<1 && beta>0 && gamma>0 && t>0)) {
                t = -1;
            }
        }
        if (level==0) return t;
        Vector intersectionPoint = ray.start + ray.dir.scalarProduct(t);
        Color intersectionPointColor = getColorAt(intersectionPoint);
        color.r = intersectionPointColor.r*coEfficients[0];
        color.g = intersectionPointColor.g*coEfficients[0];
        color.b = intersectionPointColor.b*coEfficients[0];
        Vector normal = (b-a).crossProduct(c-a);
        normal.normalize();
        if (ray.dir.dotProduct(normal)>0) {
            normal.x = -normal.x;
            normal.y = -normal.y;
            normal.z = -normal.z;
        }
        for (int i=0;i<pointLights.size();i++) {
            Ray lightRay(pointLights[i]->light_pos,intersectionPoint-pointLights[i]->light_pos);
            Color dummyColor;
            Vector d = intersectionPoint - pointLights[i]->light_pos;
            double dis = sqrt(d.x*d.x+d.y*d.y+d.z*d.z);
            bool obscured = false;
            for (int j=0;j<objects.size();j++) {
                double tVal = objects[j]->intersect(lightRay,dummyColor,0);
                if (tVal>0 && tVal+epsilon<dis) {
                    obscured = true;
                    break;
                }
            }
            if (!obscured) {
                double lambert = max(-lightRay.dir.dotProduct(normal),0.0);
                Vector reflectedRayDir = lightRay.dir-normal.scalarProduct(2*normal.dotProduct(lightRay.dir));
                double phong = max(-reflectedRayDir.dotProduct(ray.dir),0.0);
                color.r += pointLights[i]->color.r*coEfficients[1]*lambert*intersectionPointColor.r;
                color.g += pointLights[i]->color.g*coEfficients[1]*lambert*intersectionPointColor.g;
                color.b += pointLights[i]->color.b*coEfficients[1]*lambert*intersectionPointColor.b;
                color.r += pointLights[i]->color.r*coEfficients[2]*pow(phong,shine)*intersectionPointColor.r;
                color.g += pointLights[i]->color.g*coEfficients[2]*pow(phong,shine)*intersectionPointColor.g;
                color.b += pointLights[i]->color.b*coEfficients[2]*pow(phong,shine)*intersectionPointColor.b;
            }
        }
        for (int i=0;i<spotLights.size();i++) {
            Ray lightRay(spotLights[i]->point_light.light_pos,intersectionPoint-spotLights[i]->point_light.light_pos);
            Vector e = spotLights[i]->light_direction;
            e.normalize();
            double ans = e.dotProduct(lightRay.dir);
            if ((acos(ans)*(180.0/M_PI)) <= spotLights[i]->cutoff_angle) {
                Color dummyColor;
                Vector d = intersectionPoint - spotLights[i]->point_light.light_pos;
                double dis = sqrt(d.x*d.x+d.y*d.y+d.z*d.z);
                bool obscured = false;
                for (int j=0;j<objects.size();j++) {
                    double tVal = objects[j]->intersect(lightRay,dummyColor,0);
                    if (tVal>0 && tVal+epsilon<dis) {
                        obscured = true;
                        break;
                    }
                }
                if (!obscured) {
                    double lambert = max(-lightRay.dir.dotProduct(normal),0.0);
                    Vector reflectedRayDir = lightRay.dir-normal.scalarProduct(2*normal.dotProduct(lightRay.dir));
                    double phong = max(-reflectedRayDir.dotProduct(ray.dir),0.0);
                    color.r +=  spotLights[i]->point_light.color.r*coEfficients[1]*lambert*intersectionPointColor.r;
                    color.g +=  spotLights[i]->point_light.color.g*coEfficients[1]*lambert*intersectionPointColor.g;
                    color.b +=  spotLights[i]->point_light.color.b*coEfficients[1]*lambert*intersectionPointColor.b;
                    color.r +=  spotLights[i]->point_light.color.r*coEfficients[2]*pow(phong,shine)*intersectionPointColor.r;
                    color.g +=  spotLights[i]->point_light.color.g*coEfficients[2]*pow(phong,shine)*intersectionPointColor.g;
                    color.b +=  spotLights[i]->point_light.color.b*coEfficients[2]*pow(phong,shine)*intersectionPointColor.b;
                }
            }
        }
        if (level>=recursionLevel) {
            return t;
        }
        Ray reflectedRay(intersectionPoint,ray.dir-normal.scalarProduct(2*normal.dotProduct(ray.dir)));
        reflectedRay.start = reflectedRay.start + reflectedRay.dir.scalarProduct(epsilon);
        double t0=-1, tMin;
        Object* obj = nullptr;
        Color dummyColor;
        for (int k=0;k<objects.size();k++) {
            t0 = objects[k]->intersect(reflectedRay,dummyColor,0);
            if (t0>0) {
                if (t0<tMin || obj==nullptr) {
                    tMin = t0;
                    obj = objects[k];
                }
            }
        }
        if (obj!=nullptr) {
            Color c;
            obj->intersect(reflectedRay,c,level+1);
            color.r += c.r*coEfficients[3];
            color.g += c.g*coEfficients[3];
            color.b += c.b*coEfficients[3];
        }
        return t;
    }
};

class General: public Object {
public:
    double a,b,c,d,e,f,g,h,i,j;
    General(double a, double b, double c, double d, double e, double f, double g, double h, double i, double j, Vector reference_point, double length, double width, double height) {
        this->a = a;
        this->b = b;
        this->c = c;
        this->d = d;
        this->e = e;
        this->f = f;
        this->g = g;
        this->h = h;
        this->i = i;
        this->j = j;
        this->reference_point = reference_point;
        this->length = length;
        this->width = width;
        this->height = height;
    }
    void draw() {
        
    }
    double intersect(Ray ray, Color& color, int level) {
        double A = ray.start.x;
        double B = ray.dir.x;
        double C = ray.start.y;
        double D = ray.dir.y;
        double E = ray.start.z;
        double F = ray.dir.z;
        double m = a*B*B + b*D*D + c*F*F + d*B*D + e*B*F + f*D*F;
        double n = 2*A*B*a + 2*C*D*b + 2*E*F*c + d*A*D + d*B*C + A*F*e + E*B*e + C*F*f + D*E*f + g*B + h*D + i*F;
        double o = a*A*A + b*C*C + c*E*E + g*A + h*C + i*E + d*A*C + e*A*E + f*E*C + j;
        double t = -1;
        if (m==0) {
            t = -o/n;
        }
        else if (n*n-4*m*o<0) {
            t = -1;
        }
        else {
            t = -1;
            double t1 = (-n + sqrt(n*n-4*m*o))/(2*m);
            double t2 = (-n - sqrt(n*n-4*m*o))/(2*m);
            double tmin = min(t1,t2);
            double tmax = max(t1,t2);
            if (tmin>0) {
                Vector intersectPoint = ray.start + ray.dir.scalarProduct(tmin);
                if ((length==0 || (intersectPoint.x>=reference_point.x && intersectPoint.x<=reference_point.x+length)) && (width==0 || (intersectPoint.y>=reference_point.y && intersectPoint.y<=reference_point.y+width)) && (height==0 || (intersectPoint.z>=reference_point.z && intersectPoint.z<=reference_point.z+height))) {
                    t = tmin;
                }
            }
            if (t==-1 && tmax>0) {
                Vector intersectPoint = ray.start + ray.dir.scalarProduct(tmax);
                if ((length==0 || (intersectPoint.x>=reference_point.x && intersectPoint.x<=reference_point.x+length)) && (width==0 || (intersectPoint.y>=reference_point.y && intersectPoint.y<=reference_point.y+width)) && (height==0 || (intersectPoint.z>=reference_point.z && intersectPoint.z<=reference_point.z+height))) {
                    t = tmax;
                }
            }
        }
        if (level==0) return t;
        Vector intersectionPoint = ray.start + ray.dir.scalarProduct(t);
        Color intersectionPointColor = getColorAt(intersectionPoint);
        color.r = intersectionPointColor.r*coEfficients[0];
        color.g = intersectionPointColor.g*coEfficients[0];
        color.b = intersectionPointColor.b*coEfficients[0];
        Vector normal(2*a*intersectionPoint.x+d*intersectionPoint.y+e*intersectionPoint.z+g,2*b*intersectionPoint.y+d*intersectionPoint.x+f*intersectionPoint.z+h,2*c*intersectionPoint.z+e*intersectionPoint.x+f*intersectionPoint.y+i);
        normal.normalize();
        if (ray.dir.dotProduct(normal)>0) {
            normal.x = -normal.x;
            normal.y = -normal.y;
            normal.z = -normal.z;
        }
        for (int ii=0;ii<pointLights.size();ii++) {
            Ray lightRay(pointLights[ii]->light_pos,intersectionPoint-pointLights[ii]->light_pos);
            Color dummyColor;
            Vector d = intersectionPoint - pointLights[ii]->light_pos;
            double dis = sqrt(d.x*d.x+d.y*d.y+d.z*d.z);
            bool obscured = false;
            for (int jj=0;jj<objects.size();jj++) {
                double tVal = objects[jj]->intersect(lightRay,dummyColor,0);
                if (tVal>0 && tVal+epsilon<dis) {
                    obscured = true;
                    break;
                }
            }
            if (!obscured) {
                double lambert = max(-lightRay.dir.dotProduct(normal),0.0);
                Vector reflectedRayDir = lightRay.dir-normal.scalarProduct(2*normal.dotProduct(lightRay.dir));
                double phong = max(-reflectedRayDir.dotProduct(ray.dir),0.0);
                color.r += pointLights[ii]->color.r*coEfficients[1]*lambert*intersectionPointColor.r;
                color.g += pointLights[ii]->color.g*coEfficients[1]*lambert*intersectionPointColor.g;
                color.b += pointLights[ii]->color.b*coEfficients[1]*lambert*intersectionPointColor.b;
                color.r += pointLights[ii]->color.r*coEfficients[2]*pow(phong,shine)*intersectionPointColor.r;
                color.g += pointLights[ii]->color.g*coEfficients[2]*pow(phong,shine)*intersectionPointColor.g;
                color.b += pointLights[ii]->color.b*coEfficients[2]*pow(phong,shine)*intersectionPointColor.b;
            }
        }
        for (int ii=0;ii<spotLights.size();ii++) {
            Ray lightRay(spotLights[ii]->point_light.light_pos,intersectionPoint-spotLights[ii]->point_light.light_pos);
            Vector e = spotLights[ii]->light_direction;
            e.normalize();
            double ans = e.dotProduct(lightRay.dir);
            if ((acos(ans)*(180.0/M_PI)) <= spotLights[ii]->cutoff_angle) {
                Color dummyColor;
                Vector d = intersectionPoint - spotLights[ii]->point_light.light_pos;
                double dis = sqrt(d.x*d.x+d.y*d.y+d.z*d.z);
                bool obscured = false;
                for (int jj=0;jj<objects.size();jj++) {
                    double tVal = objects[jj]->intersect(lightRay,dummyColor,0);
                    if (tVal>0 && tVal+epsilon<dis) {
                        obscured = true;
                        break;
                    }
                }
                if (!obscured) {
                    double lambert = max(-lightRay.dir.dotProduct(normal),0.0);
                    Vector reflectedRayDir = lightRay.dir-normal.scalarProduct(2*normal.dotProduct(lightRay.dir));
                    double phong = max(-reflectedRayDir.dotProduct(ray.dir),0.0);
                    color.r +=  spotLights[ii]->point_light.color.r*coEfficients[1]*lambert*intersectionPointColor.r;
                    color.g +=  spotLights[ii]->point_light.color.g*coEfficients[1]*lambert*intersectionPointColor.g;
                    color.b +=  spotLights[ii]->point_light.color.b*coEfficients[1]*lambert*intersectionPointColor.b;
                    color.r +=  spotLights[ii]->point_light.color.r*coEfficients[2]*pow(phong,shine)*intersectionPointColor.r;
                    color.g +=  spotLights[ii]->point_light.color.g*coEfficients[2]*pow(phong,shine)*intersectionPointColor.g;
                    color.b +=  spotLights[ii]->point_light.color.b*coEfficients[2]*pow(phong,shine)*intersectionPointColor.b;
                }
            }
        }
        if (level>=recursionLevel) {
            return t;
        }
        Ray reflectedRay(intersectionPoint,ray.dir-normal.scalarProduct(2*normal.dotProduct(ray.dir)));
        reflectedRay.start = reflectedRay.start + reflectedRay.dir.scalarProduct(epsilon);
        double t0=-1, tMin;
        Object* obj = nullptr;
        Color dummyColor;
        for (int k=0;k<objects.size();k++) {
            t0 = objects[k]->intersect(reflectedRay,dummyColor,0);
            if (t0>0) {
                if (t0<tMin || obj==nullptr) {
                    tMin = t0;
                    obj = objects[k];
                }
            }
        }
        if (obj!=nullptr) {
            Color c;
            obj->intersect(reflectedRay,c,level+1);
            color.r += c.r*coEfficients[3];
            color.g += c.g*coEfficients[3];
            color.b += c.b*coEfficients[3];
        }
        return t;
    }
};

class Floor: public Object {
public:
    Floor(double floorWidth, double tileWidth) {
        Vector ref(-floorWidth/2.0,-floorWidth/2.0,0);
        this->reference_point = ref;
        width = floorWidth;             // width being used a floor width
        length = tileWidth;             // length being used a tileWidth (assumed square tile)
    }
    void draw() {
        int tileCount = width/length;   // length is tileWidth
        for (int i=0;i<tileCount;i++) {
            for (int j=0;j<tileCount;j++) {
                if ((i+j)&1) glColor3f(color.r,color.g,color.b);
                else glColor3f(1.0-color.r,1.0-color.g,1.0-color.b);
                glBegin(GL_QUADS);{
                    glVertex3f(reference_point.x+i*length,reference_point.y+j*length,0);
                    glVertex3f(reference_point.x+i*length,reference_point.y+(j+1)*length,0);
                    glVertex3f(reference_point.x+(i+1)*length,reference_point.y+(j+1)*length,0);
                    glVertex3f(reference_point.x+(i+1)*length,reference_point.y+j*length,0);
                }glEnd();
            }
        }
    }
    double intersect(Ray ray, Color& color, int level) {
        Vector normal(0,0,1);
        double t = -1;
        if (ray.dir.dotProduct(normal)==0) {
            t = -1;
        }
        else {
            t = -(normal.dotProduct(ray.start))/(normal.dotProduct(ray.dir));
            Vector intersectPoint = ray.start + ray.dir.scalarProduct(t);
            if (!((intersectPoint.x<=width/2.0 && intersectPoint.x>=-width/2.0) && (intersectPoint.y<=width/2.0 && intersectPoint.y>=-width/2.0))) {
                t = -1;
            }
        }
        if (level==0) return t;
        Vector intersectionPoint = ray.start + ray.dir.scalarProduct(t);
        Color intersectionPointColor = getColorAt(intersectionPoint);
        color.r = intersectionPointColor.r*coEfficients[0];
        color.g = intersectionPointColor.g*coEfficients[0];
        color.b = intersectionPointColor.b*coEfficients[0];
        //
        if (ray.dir.dotProduct(normal)>0) {
            normal.x = -normal.x;
            normal.y = -normal.y;
            normal.z = -normal.z;
        }
        for (int i=0;i<pointLights.size();i++) {
            Ray lightRay(pointLights[i]->light_pos,intersectionPoint-pointLights[i]->light_pos);
            Color dummyColor;
            Vector d = intersectionPoint - pointLights[i]->light_pos;
            double dis = sqrt(d.x*d.x+d.y*d.y+d.z*d.z);
            bool obscured = false;
            for (int j=0;j<objects.size();j++) {
                double tVal = objects[j]->intersect(lightRay,dummyColor,0);
                if (tVal>0 && tVal+epsilon<dis) {
                    obscured = true;
                    break;
                }
            }
            if (!obscured) {
                double lambert = max(-lightRay.dir.dotProduct(normal),0.0);
                Vector reflectedRayDir = lightRay.dir-normal.scalarProduct(2*normal.dotProduct(lightRay.dir));
                double phong = max(-reflectedRayDir.dotProduct(ray.dir),0.0);
                color.r += pointLights[i]->color.r*coEfficients[1]*lambert*intersectionPointColor.r;
                color.g += pointLights[i]->color.g*coEfficients[1]*lambert*intersectionPointColor.g;
                color.b += pointLights[i]->color.b*coEfficients[1]*lambert*intersectionPointColor.b;
                color.r += pointLights[i]->color.r*coEfficients[2]*pow(phong,shine)*intersectionPointColor.r;
                color.g += pointLights[i]->color.g*coEfficients[2]*pow(phong,shine)*intersectionPointColor.g;
                color.b += pointLights[i]->color.b*coEfficients[2]*pow(phong,shine)*intersectionPointColor.b;
            }
        }
        for (int i=0;i<spotLights.size();i++) {
            Ray lightRay(spotLights[i]->point_light.light_pos,intersectionPoint-spotLights[i]->point_light.light_pos);
            Vector e = spotLights[i]->light_direction;
            e.normalize();
            double ans = e.dotProduct(lightRay.dir);
            if ((acos(ans)*(180.0/M_PI)) <= spotLights[i]->cutoff_angle) {
                Color dummyColor;
                Vector d = intersectionPoint - spotLights[i]->point_light.light_pos;
                double dis = sqrt(d.x*d.x+d.y*d.y+d.z*d.z);
                bool obscured = false;
                for (int j=0;j<objects.size();j++) {
                    double tVal = objects[j]->intersect(lightRay,dummyColor,0);
                    if (tVal>0 && tVal+epsilon<dis) {
                        obscured = true;
                        break;
                    }
                }
                if (!obscured) {
                    double lambert = max(-lightRay.dir.dotProduct(normal),0.0);
                    Vector reflectedRayDir = lightRay.dir-normal.scalarProduct(2*normal.dotProduct(lightRay.dir));
                    double phong = max(-reflectedRayDir.dotProduct(ray.dir),0.0);
                    color.r +=  spotLights[i]->point_light.color.r*coEfficients[1]*lambert*intersectionPointColor.r;
                    color.g +=  spotLights[i]->point_light.color.g*coEfficients[1]*lambert*intersectionPointColor.g;
                    color.b +=  spotLights[i]->point_light.color.b*coEfficients[1]*lambert*intersectionPointColor.b;
                    color.r +=  spotLights[i]->point_light.color.r*coEfficients[2]*pow(phong,shine)*intersectionPointColor.r;
                    color.g +=  spotLights[i]->point_light.color.g*coEfficients[2]*pow(phong,shine)*intersectionPointColor.g;
                    color.b +=  spotLights[i]->point_light.color.b*coEfficients[2]*pow(phong,shine)*intersectionPointColor.b;
                }
            }
        }
        if (level>=recursionLevel) {
            return t;
        }
        Ray reflectedRay(intersectionPoint,ray.dir-normal.scalarProduct(2*normal.dotProduct(ray.dir)));
        reflectedRay.start = reflectedRay.start + reflectedRay.dir.scalarProduct(epsilon);
        double t0=-1, tMin;
        Object* obj = nullptr;
        Color dummyColor;
        for (int k=0;k<objects.size();k++) {
            t0 = objects[k]->intersect(reflectedRay,dummyColor,0);
            if (t0>0) {
                if (t0<tMin || obj==nullptr) {
                    tMin = t0;
                    obj = objects[k];
                }
            }
        }
        if (obj!=nullptr) {
            Color c;
            obj->intersect(reflectedRay,c,level+1);
            color.r += c.r*coEfficients[3];
            color.g += c.g*coEfficients[3];
            color.b += c.b*coEfficients[3];
        }
        return t;
    }
    Color getColorAt(Vector intersectionPoint) {
        int i = (intersectionPoint.x - reference_point.x)/length;
        int j = (intersectionPoint.y - reference_point.y)/length;
        if ((i+j)&1) return this->color;
        return Color(1.0-color.r,1.0-color.g,1.0-color.b);
    }
};


#endif