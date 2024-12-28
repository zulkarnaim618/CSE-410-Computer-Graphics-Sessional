#include <bits/stdc++.h>
#include "bitmap_image.hpp"
using namespace std;


static unsigned long int g_seed = 1;
inline int randomt()
{
    g_seed = (214013 * g_seed + 2531011);
    return (g_seed >> 16) & 0x7FFF;
}

class Color {
 public:
    int r;
    int g;
    int b;
    Color() {
        r = g = b = 0;
    }
    Color(int r, int g, int b) {
        this->r = r;
        this->g = g;
        this->b = b;
    }
    Color(const Color& c) {
        r = c.r;
        g = c.g;
        b = c.b;
    }
};

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
    Vector(const Vector& v) {
        x = v.x;
        y = v.y;
        z = v.z;
    }
    Vector(Vector* v) {
        x = v->x;
        y = v->y;
        z = v->z;
    }
    void normalize() {
        double t = sqrt(x*x+y*y+z*z);
        x = x/t;
        y = y/t;
        z = z/t;
    }
    Vector crossProduct(Vector& b){
        Vector ans;
        ans.x = this->y*b.z-this->z*b.y;
        ans.y = this->z*b.x-this->x*b.z;
        ans.z = this->x*b.y-this->y*b.x;
        return ans;
    }
    Vector rodriguesFormulae(Vector& axisNormal, double angle) {
        Vector cross = this->crossProduct(axisNormal);
        Vector ans;
        ans.x = axisNormal.x*cos(angle) + (1-cos(angle))*this->dotProduct(axisNormal)*this->x+sin(angle)*cross.x;
        ans.y = axisNormal.y*cos(angle) + (1-cos(angle))*this->dotProduct(axisNormal)*this->y+sin(angle)*cross.y;
        ans.z = axisNormal.z*cos(angle) + (1-cos(angle))*this->dotProduct(axisNormal)*this->z+sin(angle)*cross.z;
        return ans;
    }
    double dotProduct(Vector& b) {
        double ans;
        ans = x*b.x+y*b.y+z*b.z;
        return ans;
    }
    void scalarMultiply(double s) {
        x = x*s;
        y = y*s;
        z = z*s;
    }
};
ostream& operator<<(ostream& out, Vector c) {
    out<<"("<<c.x<<","<<c.y<<","<<c.z<<")";
    return out;
}
class Matrix {
public:
    int r,c;
    double** mat;
    Matrix() {
        this->r = 4;
        this->c = 4;
        mat = new double*[r];
        for (int i=0;i<r;i++) {
            mat[i] = new double[c];
        }
        makeZero();
    }
    Matrix(int r, int c) {
        this->r = r;
        this->c = c;
        mat = new double*[r];
        for (int i=0;i<r;i++) {
            mat[i] = new double[c];
        }
        makeZero();
    }
    Matrix(const Matrix& m) {
        this->r = m.r;
        this->c = m.c;
        mat = new double*[r];
        for (int i=0;i<r;i++) {
            mat[i] = new double[c];
        }
        for (int i=0;i<r;i++) {
            for (int j=0;j<c;j++) {
                mat[i][j] = m.mat[i][j];
            }
        }
    }
    Matrix(Vector& v) {
        this->r = 4;
        this->c = 1;
        mat = new double*[r];
        for (int i=0;i<r;i++) {
            mat[i] = new double[c];
        }
        mat[0][0] = v.x;
        mat[1][0] = v.y;
        mat[2][0] = v.z;
        mat[3][0] = 1;
    }
    ~Matrix() {
        free();
    }
    void free() {
        for (int i=0;i<r;i++) {
            delete[] mat[i];
        }
        delete[] mat;
    }
    void makeZero() {
        for (int i=0;i<r;i++) {
            for (int j=0;j<c;j++) {
                mat[i][j] = 0;
            }
        }
    }
    void makeIdentity() {
        if (r!=c) {
            cout<<"Not square matrix"<<endl;
            return;
        }
        makeZero();
        for (int i=0;i<r;i++) {
            mat[i][i] = 1;
        }
    }
    Matrix matMul(Matrix& b) {
        Matrix ans(this->r,b.c);
        if (this->c!=b.r) cout<<"Dimenstion mismatch"<<endl;
        else {
            for (int i=0;i<r;i++) {
                for (int j=0;j<b.c;j++) {
                    for (int k=0;k<c;k++) {
                        ans.mat[i][j]+= this->mat[i][k]*b.mat[k][j];
                    }
                }
            }
        }
        return ans;
    }
    Matrix operator=(Matrix m) {
        free();
        this->r = m.r;
        this->c = m.c;
        mat = new double*[r];
        for (int i=0;i<r;i++) {
            mat[i] = new double[c];
        }
        for (int i=0;i<r;i++) {
            for (int j=0;j<c;j++) {
                mat[i][j] = m.mat[i][j];
            }
        }
        return *this;
    }
};
ostream& operator<<(ostream& out, Matrix m) {
    for (int i=0;i<m.r;i++) {
        for (int j=0;j<m.c;j++) {
            out<<m.mat[i][j]<<" ";
        }
        out<<endl;
    }
    return out;
}
class Triangle {
public:
    Vector* a;
    Vector* b;
    Vector* c;
    Color* color;
    Triangle() {
        a = new Vector();
        b = new Vector();
        c = new Vector();
        color = new Color();
    }
    Triangle(Vector& a, Vector& b, Vector& c, Color& color) {
        this->a = new Vector(a);
        this->b = new Vector(b);
        this->c = new Vector(c);
        this->color = new Color(color);
    }
    Triangle(const Triangle& t) {
        a = new Vector(*t.a);
        b = new Vector(*t.b);
        c = new Vector(*t.c);
        color = new Color(*t.color);
    }
    ~Triangle() {
        delete a;
        delete b;
        delete c;
        delete color;
    }
};
void outputPoint(ofstream& out, Matrix& m) {
    out<<fixed<<setprecision(7);
    out<<m.mat[0][0]<<" "<<m.mat[1][0]<<" "<<m.mat[2][0]<<" "<<endl;
}

bool contains(Vector* a, Vector* b, double scanLineY) {
    double maxY = max(a->y,b->y);
    double minY = min(a->y,b->y);
    if (minY<=scanLineY && scanLineY<=maxY && minY!=maxY) return true;
    return false;
}
double getCordinateValue(double x1, double y1, double x2, double y2, double x) {
    double y;
    y = y1 + ((x-x1)*(y1-y2))/(x1-x2);
    return y;
}

int main ()
{
    string input = "IOs/5/";
    string output = "output/5/";
    ifstream in(input+"scene.txt");
    ofstream out(output+"stage1.txt");
    stack<Matrix> s;
    Matrix m;
    m.makeIdentity();
    string command;
    double a,b,c;
    int cnt = 0;
    in>>a>>b>>c;
    Vector eye(a,b,c);
    in>>a>>b>>c;
    Vector look(a,b,c);
    in>>a>>b>>c;
    Vector up(a,b,c);
    double fovY,aspectRatio,near,far;
    in>>fovY>>aspectRatio>>near>>far;
    while(true) {
        in>>command;
        if (command=="triangle") {
            cnt++;
            for (int i=0;i<3;i++) {
                in>>a>>b>>c;
                Vector p(a,b,c);
                Matrix mP(p);
                Matrix ans = m.matMul(mP);
                outputPoint(out,ans);
            }
            out<<endl;
        }
        else if (command=="translate") {
            in>>a>>b>>c;
            Matrix t;
            t.makeIdentity();
            t.mat[0][3] = a;
            t.mat[1][3] = b;
            t.mat[2][3] = c;
            m = m.matMul(t);
        }
        else if (command=="scale") {
            in>>a>>b>>c;
            Matrix s;
            s.makeIdentity();
            s.mat[0][0] = a;
            s.mat[1][1] = b;
            s.mat[2][2] = c;
            m = m.matMul(s);
        }
        else if (command=="rotate") {
            double angle;
            in>>angle>>a>>b>>c;
            angle = (angle*acos(-1))/180.0;
            Vector dir(a,b,c);
            dir.normalize();
            Vector xNorm(1,0,0);
            Vector yNorm(0,1,0);
            Vector zNorm(0,0,1);
            Vector c1 = dir.rodriguesFormulae(xNorm,angle);
            Vector c2 = dir.rodriguesFormulae(yNorm,angle);
            Vector c3 = dir.rodriguesFormulae(zNorm,angle);
            Matrix r;
            r.makeIdentity();
            r.mat[0][0] = c1.x; r.mat[0][1] = c2.x; r.mat[0][2] = c3.x;
            r.mat[1][0] = c1.y; r.mat[1][1] = c2.y; r.mat[1][2] = c3.y;
            r.mat[2][0] = c1.z; r.mat[2][1] = c2.z; r.mat[2][2] = c3.z;
            m = m.matMul(r);
        }
        else if (command=="push") {
            s.push(m);
        }
        else if (command=="pop") {
            m = s.top();
            s.pop();
        }
        else if (command=="end") {
            break;
        }
    }
    in.close();
    out.close();
    in.open(output+"stage1.txt");
    out.open(output+"stage2.txt");
    Vector l(look.x-eye.x,look.y-eye.y,look.z-eye.z);
    l.normalize();
    Vector r = l.crossProduct(up);
    r.normalize();
    Vector u = r.crossProduct(l);
    Matrix t;
    t.makeIdentity();
    t.mat[0][3] = -eye.x;
    t.mat[1][3] = -eye.y;
    t.mat[2][3] = -eye.z;
    Matrix ro;
    ro.makeIdentity();
    ro.mat[0][0] = r.x; ro.mat[0][1] = r.y; ro.mat[0][2] = r.z;
    ro.mat[1][0] = u.x; ro.mat[1][1] = u.y; ro.mat[1][2] = u.z;
    ro.mat[2][0] = -l.x; ro.mat[2][1] = -l.y; ro.mat[2][2] = -l.z;
    Matrix v = ro.matMul(t);
    for (int i=0;i<cnt;i++) {
        for (int j=0;j<3;j++) {
            in>>a>>b>>c;
            Vector p(a,b,c);
            Matrix mP(p);
            Matrix ans = v.matMul(mP);
            outputPoint(out,ans);
        }
        out<<endl;
    }
    in.close();
    out.close();
    in.open(output+"stage2.txt");
    out.open(output+"stage3.txt");
    double fovX = fovY*aspectRatio;
    double tt = near * tan((fovY/2.0)*(acos(-1)/180.0));
    double rr = near * tan((fovX/2.0)*(acos(-1)/180.0));
    Matrix p;
    p.mat[0][0] = near/rr;
    p.mat[1][1] = near/tt;
    p.mat[2][2] =  -(far+near)/(far-near);
    p.mat[3][2] = -1.0;
    p.mat[2][3] = -(2.0*far*near)/(far-near);
    for (int i=0;i<cnt;i++) {
        for (int j=0;j<3;j++) {
            in>>a>>b>>c;
            Vector pp(a,b,c);
            Matrix mP(pp);
            Matrix ans = p.matMul(mP);
            ans.mat[0][0] = ans.mat[0][0]/ans.mat[3][0];
            ans.mat[1][0] = ans.mat[1][0]/ans.mat[3][0];
            ans.mat[2][0] = ans.mat[2][0]/ans.mat[3][0];

            outputPoint(out,ans);
        }
        out<<endl;
    }
    in.close();
    out.close();

    ifstream config(input+"config.txt");
    in.open(output+"stage3.txt");
    out.open(output+"z_buffer.txt");
    double width,height;
    config>>width>>height;
    Triangle** triangles = new Triangle*[cnt];
    for (int i=0;i<cnt;i++) {
        in>>a>>b>>c;
        Vector A(a,b,c);
        in>>a>>b>>c;
        Vector B(a,b,c);
        in>>a>>b>>c;
        Vector C(a,b,c);
        Color color;
        color.r = randomt()%256;
        color.g = randomt()%256;
        color.b = randomt()%256;
        triangles[i] = new Triangle(A,B,C,color);
    }
    double dx,dy,rightLimitX,leftLimitX,topLimitY,bottomLimitY;
    rightLimitX = 1;
    leftLimitX = -1;
    topLimitY = 1;
    bottomLimitY = -1;
    dx = (rightLimitX-leftLimitX)/width;
    dy = (topLimitY-bottomLimitY)/height;
    double topY,leftX;
    topY = topLimitY - dy/2.0;
    leftX = leftLimitX + dx/2.0;
    double** zBuffer;
    zBuffer = new double*[(int)height];
    for (int i=0;i<height;i++) {
        zBuffer[i] = new double[(int)width];
    }
    double zMax = 1.0;
    double zMin = -1.0;
    for (int i=0;i<height;i++) {
        for (int j=0;j<width;j++) {
            zBuffer[i][j] = zMax;
        }
    }
    bitmap_image img(width,height);
    for (int i=0;i<width;i++) {
        for (int j=0;j<height;j++) {
            img.set_pixel(i,j,0,0,0);
        }
    }
    //algo
    for (int i=0;i<cnt;i++) {
        double maxY = max(max(triangles[i]->a->y,triangles[i]->b->y),triangles[i]->c->y);
        double minY = min(min(triangles[i]->a->y,triangles[i]->b->y),triangles[i]->c->y);
        int topScanLine, bottomScanLine;
        topScanLine = ceil((topY-min(topY,maxY))/dy);
        bottomScanLine = min(height-1,floor((topY-min(topY,minY))/dy));
        for (int j=topScanLine;j<=bottomScanLine;j++) {
            double scanLineY = topY-j*dy;
            Vector *p1,*p2,*p3;
            if (contains(triangles[i]->a,triangles[i]->b,scanLineY)) {
                if (contains(triangles[i]->a,triangles[i]->c,scanLineY)) {
                    p1 = new Vector(triangles[i]->a);
                    p2 = new Vector(triangles[i]->b);
                    p3 = new Vector(triangles[i]->c);
                }
                else {
                    p1 = new Vector(triangles[i]->b);
                    p2 = new Vector(triangles[i]->a);
                    p3 = new Vector(triangles[i]->c);
                }
            }
            else {
                p1 = new Vector(triangles[i]->c);
                p2 = new Vector(triangles[i]->a);
                p3 = new Vector(triangles[i]->b);
            }
            double xa,za,xb,zb;
            xa = getCordinateValue(p1->y,p1->x,p2->y,p2->x,scanLineY);
            xb = getCordinateValue(p1->y,p1->x,p3->y,p3->x,scanLineY);
            za = getCordinateValue(p1->y,p1->z,p2->y,p2->z,scanLineY);
            zb = getCordinateValue(p1->y,p1->z,p3->y,p3->z,scanLineY);
            double maxX = max(xa,xb);
            double minX = min(xa,xb);
            int leftIntColumn, rightIntColumn;
            leftIntColumn = round((max(leftX,minX)-leftX)/dx);
            rightIntColumn = min(width-1,round((max(leftX,maxX)-leftX)/dx));
            for (int k=leftIntColumn;k<=rightIntColumn;k++) {
                double scanLineX = leftX + k*dx;
                double zValue;
                zValue = getCordinateValue(xa,za,xb,zb,scanLineX);
                if (zMin<=zValue && zValue<zBuffer[j][k]) {
                    zBuffer[j][k] = zValue;
                    img.set_pixel(k,j,triangles[i]->color->r,triangles[i]->color->g,triangles[i]->color->b);
                }
            }
            delete p1;
            delete p2;
            delete p3;
        }
    }
    //save
    img.save_image(output+"out.bmp");
    out<<fixed<<setprecision(6);
    for (int i=0;i<height;i++) {
        for (int j=0;j<width;j++) {
            if (zBuffer[i][j]<zMax) {
                out<<zBuffer[i][j]<<"\t";
            }
        }
        out<<endl;
    }
    //free
    for (int i=0;i<width;i++) {
        delete[] zBuffer[i];
    }
    delete[] zBuffer;

    for (int i=0;i<cnt;i++) {
        delete triangles[i];
    }
    delete[] triangles;
    config.close();
    in.close();
    out.close();
    return 0;
}
