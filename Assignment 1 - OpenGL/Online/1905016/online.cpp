#include<bits/stdc++.h>
#include "helper.h"
#ifdef __linux__
    #include<GL/glut.h>
#elif WIN32
    #include <windows.h>
    #include <glut.h>
#endif
using namespace std;


Camera* camera;
double transformRate = 0.03;
double scale = 1;

void initGL() {
    camera = new Camera(new Vector(0,0,20),new Vector(0,0,-1),new Vector(0,1,0),new Vector(1,0,0),50);
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
    glEnable(GL_DEPTH_TEST);
}

void drawTriangleSegment(GLfloat red, GLfloat green, GLfloat blue) {
    glColor3f(red,green,blue);
    glBegin(GL_TRIANGLES);{
        glVertex3f( 1.0, 0.0, 0.0);
        glVertex3f(0.0, 1.0, 0.0);
        glVertex3f(0.0, 0.0, 1.0);
    }glEnd();
}

void drawOctaHedron() {
    for (int i=0;i<2;i++) {
        glPushMatrix();
        for (int j=0;j<2;j++) {
            glPushMatrix();
            for (int k=0;k<2;k++) {
                glPushMatrix();
                glTranslatef(1/3.0-scale/3.0,1/3.0-scale/3.0,1/3.0-scale/3.0);
                glScalef(scale,scale,scale);
                if (k&1) drawTriangleSegment(0.0,1.0,1.0);
                else drawTriangleSegment(1.0,0.0,1.0);
                glPopMatrix();
                glRotatef(90,0,1,0);
            }
            glPopMatrix();
            glRotatef(180,0,1,0);
        }
        glPopMatrix();
        glRotatef(180,0,0,1);
    }
}

void drawSphereSegment(GLfloat red, GLfloat green, GLfloat blue) {
    double centerX,centerY,centerZ;
    double radius = (1-scale)*(1.0/sqrt(3));
    centerX = 0;
    centerY = 1*scale;
    centerZ = 0;
    int steps = 20;
    for (int i=0;i<steps;i++) {
        for (int j=0;j<steps;j++) {
            Vector* n1 = new Vector(0,-sin((-45+i*(90.0/steps))*(M_PI/180.0)),cos((-45+i*(90.0/steps))*(M_PI/180.0)));
            Vector* n2 = new Vector(0,-sin((-45+(i+1)*(90.0/steps))*(M_PI/180.0)),cos((-45+(i+1)*(90.0/steps))*(M_PI/180.0)));
            Vector* n3= new Vector(cos((-45+j*(90.0/steps))*(M_PI/180.0)),-sin((-45+j*(90.0/steps))*(M_PI/180.0)),0);
            Vector* n4 = new Vector(cos((-45+(j+1)*(90.0/steps))*(M_PI/180.0)),-sin((-45+(j+1)*(90.0/steps))*(M_PI/180.0)),0);
            Vector *a,*b,*c,*d;
            a = new Vector();
            b = new Vector();
            c = new Vector();
            d = new Vector();
            n1->crossProduct(n3,a);
            n1->crossProduct(n4,b);
            n2->crossProduct(n3,c);
            n2->crossProduct(n4,d);
            a->normalize();
            b->normalize();
            c->normalize();
            d->normalize();
            glColor3f(red,green,blue);
            glBegin(GL_QUADS);{
                glVertex3f(centerX+radius*a->x,centerY+radius*a->y,centerZ+radius*a->z);
                glVertex3f(centerX+radius*b->x,centerY+radius*b->y,centerZ+radius*b->z);
                glVertex3f(centerX+radius*d->x,centerY+radius*d->y,centerZ+radius*d->z);
                glVertex3f(centerX+radius*c->x,centerY+radius*c->y,centerZ+radius*c->z);

            }glEnd();
            delete n1;
            delete n2;
            delete n3;
            delete n4;
            delete a;
            delete b;
            delete c;
            delete d;
        }
    }
}

void drawSphere() {
    glPushMatrix();
    drawSphereSegment(1,0,0);
    glRotatef(180,1,0,0);
    drawSphereSegment(1,0,0);
    glPopMatrix();
    glPushMatrix();
    glRotatef(90,1,0,0);
    drawSphereSegment(0,1,0);
    glRotatef(180,1,0,0);
    drawSphereSegment(0,1,0);
    glPopMatrix();
    glPushMatrix();
    glRotatef(90,0,0,1);
    drawSphereSegment(0,0,1);
    glRotatef(180,1,0,0);
    drawSphereSegment(0,0,1);
    glPopMatrix();
}

void drawCylinderSegment() {
    int segmentCount = 100;
    double height = sqrt(2*scale*scale);
    double cylinderAngle = 70.5287794;
    double radius = (1-scale)*(((2.0/3.0)*(sin(((180-cylinderAngle)/2.0)*(M_PI/180.0))))/(sin(cylinderAngle*(M_PI/180.0))));

    for (int i=0;i<segmentCount;i++) {
        double x1,z1,x2,z2;
        x1 = radius*cos((-(cylinderAngle/2.0)+i*(cylinderAngle/segmentCount))*(M_PI/180.0));
        z1 = radius*sin((-(cylinderAngle/2.0)+i*(cylinderAngle/segmentCount))*(M_PI/180.0));
        x2 = radius*cos((-(cylinderAngle/2.0)+(i+1)*(cylinderAngle/segmentCount))*(M_PI/180.0));
        z2 = radius*sin((-(cylinderAngle/2.0)+(i+1)*(cylinderAngle/segmentCount))*(M_PI/180.0));
        glColor3f(1,1,0);
        glBegin(GL_QUADS);{
                glVertex3f(x1,height,z1);
                glVertex3f(x2,height,z2);
                glVertex3f(x2,0,z2);
                glVertex3f(x1,0,z1);
        }glEnd();
    }
}

void drawAllCylinder() {
    double height = sqrt(2*scale*scale);
    double cylinderAngle = 70.5287794;
    double radius = (1-scale)*(((2.0/3.0)*(sin(((180-cylinderAngle)/2.0)*(M_PI/180.0))))/(sin(cylinderAngle*(M_PI/180.0))));
    glPushMatrix();
    for (int i=0;i<4;i++) {
        glPushMatrix();
        glTranslatef(scale,0,0);
        glRotatef(45,0,0,1);
        drawCylinderSegment();
        glPopMatrix();
        glRotatef(90,0,0,1);
    }
    glPopMatrix();
    glPushMatrix();
    glRotatef(90,0,1,0);
    for (int i=0;i<4;i++) {
        glPushMatrix();
        glTranslatef(scale,0,0);
        glRotatef(45,0,0,1);
        drawCylinderSegment();
        glPopMatrix();
        glRotatef(90,0,0,1);
    }
    glPopMatrix();
    glPushMatrix();
    glRotatef(90,1,0,0);
    for (int i=0;i<4;i++) {
        glPushMatrix();
        glTranslatef(scale,0,0);
        glRotatef(45,0,0,1);
        drawCylinderSegment();
        glPopMatrix();
        glRotatef(90,0,0,1);
    }
    glPopMatrix();
}

double objRotationAngleRate = 5;
double objRotationAngle = 0;


double rangle = 0;
vector<pair<double,double>> cor;
double translate = 0;

void drawcircle() {
    int division = 100;
    double angle = 0;
    double del = 360.0/division;
    double radius = 3.0;
    double centerX = -3;
    double centerY = 0;
    glColor3f(0,0,1.0);
    for (int i=0;i<division;i++) {
        double x,y,x1,y1;
        x = centerX+radius*cos(angle*(M_PI/180.0));
        y = centerY+radius*sin(angle*(M_PI/180.0));
        x1 = centerX+radius*cos((angle+del)*(M_PI/180.0));
        y1 = centerY+radius*sin((angle+del)*(M_PI/180.0));
        glBegin(GL_LINES);{
                glVertex3f(x,y,0);
                glVertex3f(x1,y1,0);
        }glEnd();
        angle += del;
    }
    double x,y;
    x = centerX+radius*cos(rangle*(M_PI/180.0));
    y = centerY+radius*sin(rangle*(M_PI/180.0));
    glColor3f(1,1,1.0);
    glBegin(GL_LINES);{
                glVertex3f(x,y,0);
                glVertex3f(centerX,centerY,0);
    }glEnd();
    glBegin(GL_LINES);{
                glVertex3f(x,y,0);
                glVertex3f(centerX+radius+1,y,0);
    }glEnd();
     //sine
     for (int i=0;i<cor.size();i++) {
         cor[i].first += translate;
     }
     cor.push_back({centerX+radius+1,y});
     for (int i=cor.size()-1;i>0;i--) {
         glBegin(GL_LINES);{
                 glVertex3f(cor[i].first,cor[i].second,0);
                 glVertex3f(cor[i-1].first,cor[i-1].second,0);
         }glEnd();
     }
}

void display() {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    gluLookAt(camera->eye->x,camera->eye->y,camera->eye->z,	camera->eye->x+camera->refDistance*camera->front->x,camera->eye->y+camera->refDistance*camera->front->y,camera->eye->z+camera->refDistance*camera->front->z,camera->up->x,camera->up->y,camera->up->z);

    drawcircle();
    // glColor3f(1.0, 0.0, 0.0);
    // glBegin(GL_LINES);{
    //     glVertex3f(100,0,0);
    //     glVertex3f(-100,0,0);


    // }glEnd();
    // glColor3f(0.0, 1.0, 0.0);
    // glBegin(GL_LINES);{
    //     glVertex3f(0,-100,0);
    //     glVertex3f(0, 100,0);


    // }glEnd();
    // glColor3f(0.0, 0.0, 1.0);
    // glBegin(GL_LINES);{
    //     glVertex3f(0,0, 100);
    //     glVertex3f(0,0,-100);


    // }glEnd();

    // glBegin(GL_LINES);{
    //     glVertex3f(5,0,0);
    //     glVertex3f(5,30,0);


    // }glEnd();

    glutSwapBuffers();
}

void reshape(GLsizei width, GLsizei height) {

   if (height == 0) height = 1;
   GLfloat aspect = (GLfloat)width / (GLfloat)height;
   glViewport(0, 0, width, height);
   glMatrixMode(GL_PROJECTION);
   glLoadIdentity();
   gluPerspective(45.0f, aspect, 0.1f, 1000.0f);
}

double rate = .5;
double angle = 3;

void keyboardListener(unsigned char key, int x, int y) {
    if (key=='1') {
        camera->right->x = camera->right->x*cos(angle*(M_PI/180.0)) + camera->front->x*sin(angle*(M_PI/180.0));
        camera->right->y = camera->right->y*cos(angle*(M_PI/180.0)) + camera->front->y*sin(angle*(M_PI/180.0));
        camera->right->z = camera->right->z*cos(angle*(M_PI/180.0)) + camera->front->z*sin(angle*(M_PI/180.0));
        camera->up->crossProduct(camera->right,camera->front);
    }
    else if (key=='2') {
        camera->right->x = camera->right->x*cos(-angle*(M_PI/180.0)) + camera->front->x*sin(-angle*(M_PI/180.0));
        camera->right->y = camera->right->y*cos(-angle*(M_PI/180.0)) + camera->front->y*sin(-angle*(M_PI/180.0));
        camera->right->z = camera->right->z*cos(-angle*(M_PI/180.0)) + camera->front->z*sin(-angle*(M_PI/180.0));
        camera->up->crossProduct(camera->right,camera->front);
    }
    else if (key=='3') {
        camera->front->x = camera->front->x*cos(angle*(M_PI/180.0)) + camera->up->x*sin(angle*(M_PI/180.0));
        camera->front->y = camera->front->y*cos(angle*(M_PI/180.0)) + camera->up->y*sin(angle*(M_PI/180.0));
        camera->front->z = camera->front->z*cos(angle*(M_PI/180.0)) + camera->up->z*sin(angle*(M_PI/180.0));
        camera->right->crossProduct(camera->front,camera->up);
    }
    else if (key=='4') {
        camera->front->x = camera->front->x*cos(-angle*(M_PI/180.0)) + camera->up->x*sin(-angle*(M_PI/180.0));
        camera->front->y = camera->front->y*cos(-angle*(M_PI/180.0)) + camera->up->y*sin(-angle*(M_PI/180.0));
        camera->front->z = camera->front->z*cos(-angle*(M_PI/180.0)) + camera->up->z*sin(-angle*(M_PI/180.0));
        camera->right->crossProduct(camera->front,camera->up);
    }
    else if (key=='6') {
        camera->right->x = camera->right->x*cos(angle*(M_PI/180.0)) + camera->up->x*sin(angle*(M_PI/180.0));
        camera->right->y = camera->right->y*cos(angle*(M_PI/180.0)) + camera->up->y*sin(angle*(M_PI/180.0));
        camera->right->z = camera->right->z*cos(angle*(M_PI/180.0)) + camera->up->z*sin(angle*(M_PI/180.0));
        camera->right->crossProduct(camera->front,camera->up);
    }
    else if (key=='5') {
        camera->right->x = camera->right->x*cos(-angle*(M_PI/180.0)) + camera->up->x*sin(-angle*(M_PI/180.0));
        camera->right->y = camera->right->y*cos(-angle*(M_PI/180.0)) + camera->up->y*sin(-angle*(M_PI/180.0));
        camera->right->z = camera->right->z*cos(-angle*(M_PI/180.0)) + camera->up->z*sin(-angle*(M_PI/180.0));
        camera->right->crossProduct(camera->front,camera->up);
    }
    else if (key=='w') {
        camera->eye->x += camera->up->x*rate;
        camera->eye->y += camera->up->y*rate;
        camera->eye->z += camera->up->z*rate;
        double rAngle = atan(rate/camera->refDistance);
        camera->front->x = camera->front->x*cos(-rAngle) + camera->up->x*sin(-rAngle);
        camera->front->y = camera->front->y*cos(-rAngle) + camera->up->y*sin(-rAngle);
        camera->front->z = camera->front->z*cos(-rAngle) + camera->up->z*sin(-rAngle);
        camera->right->crossProduct(camera->front,camera->up);
    }
    else if (key=='s') {
        camera->eye->x -= camera->up->x*rate;
        camera->eye->y -= camera->up->y*rate;
        camera->eye->z -= camera->up->z*rate;
        double rAngle = atan(rate/camera->refDistance);
        camera->front->x = camera->front->x*cos(rAngle) + camera->up->x*sin(rAngle);
        camera->front->y = camera->front->y*cos(rAngle) + camera->up->y*sin(rAngle);
        camera->front->z = camera->front->z*cos(rAngle) + camera->up->z*sin(rAngle);
        camera->right->crossProduct(camera->front,camera->up);
    }
    else if (key==',') {
        scale -= transformRate;
        if (scale<0) scale = 0;
    }
    else if (key=='.') {
        scale += transformRate;
        if (scale>1) scale = 1;
    }
    else if (key=='d') {
        objRotationAngle += objRotationAngleRate;
    }
    else if (key=='a') {
        objRotationAngle -= objRotationAngleRate;
    }
    glutPostRedisplay();
}

void specialKeyListener(int key, int x,int y) {
    if (key==GLUT_KEY_RIGHT) {
        camera->eye->x += camera->right->x*rate;
        camera->eye->y += camera->right->y*rate;
        camera->eye->z += camera->right->z*rate;
    }
    else if (key==GLUT_KEY_LEFT) {
        camera->eye->x -= camera->right->x*rate;
        camera->eye->y -= camera->right->y*rate;
        camera->eye->z -= camera->right->z*rate;
    }
    else if (key==GLUT_KEY_UP) {
        camera->eye->x += camera->front->x*rate;
        camera->eye->y += camera->front->y*rate;
        camera->eye->z += camera->front->z*rate;
    }
    else if (key==GLUT_KEY_DOWN) {
        camera->eye->x -= camera->front->x*rate;
        camera->eye->y -= camera->front->y*rate;
        camera->eye->z -= camera->front->z*rate;
    }
    else if (key==GLUT_KEY_PAGE_UP) {
        camera->eye->x += camera->up->x*rate;
        camera->eye->y += camera->up->y*rate;
        camera->eye->z += camera->up->z*rate;
    }
    else if (key==GLUT_KEY_PAGE_DOWN) {
        camera->eye->x -= camera->up->x*rate;
        camera->eye->y -= camera->up->y*rate;
        camera->eye->z -= camera->up->z*rate;
    }
    glutPostRedisplay();
}
void timer (int val) {
    rangle += 1.0;
    translate = 0.01;
    glutTimerFunc(50, timer,0);
    glutPostRedisplay();
}

int main(int argc, char** argv) {
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB);
    glutInitWindowSize(800,800);
    glutInitWindowPosition(50, 50);
    glutCreateWindow("Online");
    glutDisplayFunc(display);
    //glutKeyboardFunc(keyboardListener);
    //glutSpecialFunc(specialKeyListener);
    glutReshapeFunc(reshape);
    initGL();
    glutTimerFunc(500,timer,0);
    glutMainLoop();
    return 0;
}
