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
Ball*  ball;

bool isManual;
int dirChanged;

void initGL() {
    camera = new Camera(new Vector(0,20,70),new Vector(0,0,-1),new Vector(0,1,0),new Vector(1,0,0),50);
    ball = new Ball(new Vector(0,5,0),new Vector(0,0,-1),new Vector(0,1,0),new Vector(1,0,0),new Vector(1.0/sqrt(2),0,1.0/sqrt(2)),5);
    isManual = true;
    dirChanged = 0;
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);            
    glEnable(GL_DEPTH_TEST);   
}

void drawArrow(Vector* center, Vector* direction, double radius, double height,double red,double green, double blue) {
    int segmentCount = 50;
    double cylinderAngle = 360;
    Vector* y = new Vector(direction->x,direction->y,direction->z);
    y->normalize();
    Vector* x = new Vector(direction->y,-direction->x,0);
    x->normalize();
    Vector* z = new Vector();
    y->crossProduct(x,z);
    for (int i=0;i<segmentCount;i++) {
        double x1,z1,x2,z2;
        x1 = radius*cos((-(cylinderAngle/2.0)+i*(cylinderAngle/segmentCount))*(M_PI/180.0));
        z1 = radius*sin((-(cylinderAngle/2.0)+i*(cylinderAngle/segmentCount))*(M_PI/180.0));
        x2 = radius*cos((-(cylinderAngle/2.0)+(i+1)*(cylinderAngle/segmentCount))*(M_PI/180.0));
        z2 = radius*sin((-(cylinderAngle/2.0)+(i+1)*(cylinderAngle/segmentCount))*(M_PI/180.0));
        glColor3f(red,green,blue);  
        glBegin(GL_QUADS);{
                glVertex3f(center->x+x1*x->x+z1*z->x,center->y+x1*x->y+z1*z->y,center->z+x1*x->z+z1*z->z);
                glVertex3f(center->x+height*y->x+x1*x->x+z1*z->x,center->y+height*y->y+x1*x->y+z1*z->y,center->z+height*y->z+x1*x->z+z1*z->z);
                glVertex3f(center->x+height*y->x+x2*x->x+z2*z->x,center->y+height*y->y+x2*x->y+z2*z->y,center->z+height*y->z+x2*x->z+z2*z->z);
                glVertex3f(center->x+x2*x->x+z2*z->x,center->y+x2*x->y+z2*z->y,center->z+x2*x->z+z2*z->z);
        }glEnd();
        glBegin(GL_TRIANGLES);{
                glVertex3f(center->x+x1*x->x+z1*z->x,center->y+x1*x->y+z1*z->y,center->z+x1*x->z+z1*z->z);
                glVertex3f(center->x+x2*x->x+z2*z->x,center->y+x2*x->y+z2*z->y,center->z+x2*x->z+z2*z->z);
                glVertex3f(center->x,center->y,center->z);
        }glEnd();
        x1 = 2*radius*cos((-(cylinderAngle/2.0)+i*(cylinderAngle/segmentCount))*(M_PI/180.0));
        z1 = 2*radius*sin((-(cylinderAngle/2.0)+i*(cylinderAngle/segmentCount))*(M_PI/180.0));
        x2 = 2*radius*cos((-(cylinderAngle/2.0)+(i+1)*(cylinderAngle/segmentCount))*(M_PI/180.0));
        z2 = 2*radius*sin((-(cylinderAngle/2.0)+(i+1)*(cylinderAngle/segmentCount))*(M_PI/180.0));
        glBegin(GL_TRIANGLES);{
                glVertex3f(center->x+height*y->x+x1*x->x+z1*z->x,center->y+height*y->y+x1*x->y+z1*z->y,center->z+height*y->z+x1*x->z+z1*z->z);
                glVertex3f(center->x+height*y->x+x2*x->x+z2*z->x,center->y+height*y->y+x2*x->y+z2*z->y,center->z+height*y->z+x2*x->z+z2*z->z);
                glVertex3f(center->x+height*y->x,center->y+height*y->y,center->z+height*y->z);
        }glEnd();
        glBegin(GL_TRIANGLES);{
                glVertex3f(center->x+height*y->x+x1*x->x+z1*z->x,center->y+height*y->y+x1*x->y+z1*z->y,center->z+height*y->z+x1*x->z+z1*z->z);
                glVertex3f(center->x+height*y->x+x2*x->x+z2*z->x,center->y+height*y->y+x2*x->y+z2*z->y,center->z+height*y->z+x2*x->z+z2*z->z);
                glVertex3f(center->x+height*(4.0/3.0)*y->x,center->y+height*(4.0/3.0)*y->y,center->z+height*(4.0/3.0)*y->z);
        }glEnd();
    }
}

int stackCount = 20;
int sectorCount = 32;
int division = 4;

void drawBall() {
    GLdouble x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4;
    for(int i=0;i<stackCount;i++) {
        for (int j=0;j<sectorCount;j++) {
            if (i>=stackCount/2) {
                if (j%(sectorCount/division)*2<sectorCount/division) glColor3f(0.0, 1.0, 0.0);
                else glColor3f(1.0, 0.0, 0.0);
            }
            else {
                if (j%(sectorCount/division)*2<sectorCount/division) glColor3f(1.0, 0.0, 0.0);
                else glColor3f(0.0, 1.0, 0.0);
            }
            x1 = ball->radius*cos((90-i*(180.0/stackCount))*(M_PI/180.0))*cos(j*(360.0/sectorCount)*(M_PI/180.0));
            y1 = ball->radius*cos((90-i*(180.0/stackCount))*(M_PI/180.0))*sin(j*(360.0/sectorCount)*(M_PI/180.0));
            z1 = ball->radius*sin((90-i*(180.0/stackCount))*(M_PI/180.0));
            x2 = ball->radius*cos((90-i*(180.0/stackCount))*(M_PI/180.0))*cos((j+1)*(360.0/sectorCount)*(M_PI/180.0));
            y2 = ball->radius*cos((90-i*(180.0/stackCount))*(M_PI/180.0))*sin((j+1)*(360.0/sectorCount)*(M_PI/180.0));
            z2 = ball->radius*sin((90-i*(180.0/stackCount))*(M_PI/180.0));
            x3 = ball->radius*cos((90-(i+1)*(180.0/stackCount))*(M_PI/180.0))*cos(j*(360.0/sectorCount)*(M_PI/180.0));
            y3 = ball->radius*cos((90-(i+1)*(180.0/stackCount))*(M_PI/180.0))*sin(j*(360.0/sectorCount)*(M_PI/180.0));
            z3 = ball->radius*sin((90-(i+1)*(180.0/stackCount))*(M_PI/180.0));
            x4 = ball->radius*cos((90-(i+1)*(180.0/stackCount))*(M_PI/180.0))*cos((j+1)*(360.0/sectorCount)*(M_PI/180.0));
            y4 = ball->radius*cos((90-(i+1)*(180.0/stackCount))*(M_PI/180.0))*sin((j+1)*(360.0/sectorCount)*(M_PI/180.0));
            z4 = ball->radius*sin((90-(i+1)*(180.0/stackCount))*(M_PI/180.0));
            glBegin(GL_QUAD_STRIP);{
                glVertex3f(ball->center->x+x1*ball->front->x+z1*ball->up->x+y1*ball->right->x,ball->center->y+x1*ball->front->y+z1*ball->up->y+y1*ball->right->y,ball->center->z+x1*ball->front->z+z1*ball->up->z+y1*ball->right->z);
                glVertex3f(ball->center->x+x2*ball->front->x+z2*ball->up->x+y2*ball->right->x,ball->center->y+x2*ball->front->y+z2*ball->up->y+y2*ball->right->y,ball->center->z+x2*ball->front->z+z2*ball->up->z+y2*ball->right->z);
                glVertex3f(ball->center->x+x3*ball->front->x+z3*ball->up->x+y3*ball->right->x,ball->center->y+x3*ball->front->y+z3*ball->up->y+y3*ball->right->y,ball->center->z+x3*ball->front->z+z3*ball->up->z+y3*ball->right->z);
                glVertex3f(ball->center->x+x4*ball->front->x+z4*ball->up->x+y4*ball->right->x,ball->center->y+x4*ball->front->y+z4*ball->up->y+y4*ball->right->y,ball->center->z+x4*ball->front->z+z4*ball->up->z+y4*ball->right->z);
            }glEnd();
        }
    }
    //ball axis
    // drawArrow(ball->center,ball->front,.2,10,0,.5,.5);
    // drawArrow(ball->center,ball->up,.2,10,.5,0,.5);
    // drawArrow(ball->center,ball->right,.2,10,.5,.5,0);

    //ball direction
    drawArrow(ball->center,ball->direction,.2,10,0,0,1);
}

double floorTileSize = 20;
int floorTileCountSideways = 60;
double billiardSize = 120;
double billiardHeight = 10;

void drawFloorAndBilliard() {
    glPushMatrix();
    glTranslatef(-(floorTileCountSideways/2)*floorTileSize,0,-(floorTileCountSideways/2)*floorTileSize);
    for (int i=0;i<floorTileCountSideways;i++) {
        glPushMatrix();
        glTranslatef(0,0,i*floorTileSize);
        for (int j=0;j<floorTileCountSideways;j++) {
            if ((i+j)&1) glColor3f(0.0, 0.0, 0.0);
            else glColor3f(1.0, 1.0, 1.0);
            glBegin(GL_QUADS);{
                glVertex3f(0,0,0);
                glVertex3f(0,0,floorTileSize);
                glVertex3f(floorTileSize,0,floorTileSize);
                glVertex3f(floorTileSize,0,0);
            }glEnd();
            glTranslatef(floorTileSize,0,0);
        }
        glPopMatrix();
    }
    glPopMatrix();
    glPushMatrix();
    glColor3f(1.0, 0.0, 0.0);
    glTranslatef(-billiardSize/2,0,-billiardSize/2);
    for (int i=0;i<2;i++) {
        glBegin(GL_QUADS);{
            glVertex3f(0,0,0);
            glVertex3f(0,billiardHeight,0);
            glVertex3f(billiardSize,billiardHeight,0);
            glVertex3f(billiardSize,0,0);
        }glEnd();
        glTranslatef(0,0,billiardSize);
    }
    glPopMatrix();
    glPushMatrix();
    glTranslatef(-billiardSize/2,0,-billiardSize/2);
    for (int i=0;i<2;i++) {
        glBegin(GL_QUADS);{
            glVertex3f(0,0,0);
            glVertex3f(0,billiardHeight,0);
            glVertex3f(0,billiardHeight,billiardSize);
            glVertex3f(0,0,billiardSize);
        }glEnd();
        glTranslatef(billiardSize,0,0);
    }
    glPopMatrix();
}
 
void display() {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); 
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    gluLookAt(camera->eye->x,camera->eye->y,camera->eye->z,	camera->eye->x+camera->refDistance*camera->front->x,camera->eye->y+camera->refDistance*camera->front->y,camera->eye->z+camera->refDistance*camera->front->z,camera->up->x,camera->up->y,camera->up->z);

    drawFloorAndBilliard();
    drawBall();

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

double rate = 2;
double angle = 3;

void change(double angle,double distance,Vector* direction,Vector* front) {
    Vector* rollAxis = new Vector(-direction->z,0,direction->x);
    double rollAxisComponentValue = rollAxis->x*front->x+rollAxis->y*front->y+rollAxis->z*front->z;
    Vector* rollAxisCompVec = new Vector(rollAxisComponentValue*rollAxis->x,rollAxisComponentValue*rollAxis->y,rollAxisComponentValue*rollAxis->z);
    Vector* other = new Vector(front->x-rollAxisCompVec->x,front->y-rollAxisCompVec->y,front->z-rollAxisCompVec->z);
    Vector* ans = new Vector();
    rollAxis->crossProduct(other,ans);
    other->x = other->x*cos(angle)+ans->x*sin(angle);
    other->y = other->y*cos(angle)+ans->y*sin(angle);
    other->z = other->z*cos(angle)+ans->z*sin(angle);
    front->x = rollAxisCompVec->x+other->x;
    front->y = rollAxisCompVec->y+other->y;
    front->z = rollAxisCompVec->z+other->z;
    delete rollAxis;
    delete rollAxisCompVec;
    delete other;
    delete ans;
}

double timeIntervalMills = 1;
double velocityInMills = 0.25;
double timeTillCollisionMills = -1;


void updateTimeTillCollisionMills() {
    double xTime,zTime;
    if (ball->direction->x<0) {
        xTime = (-(-billiardSize/2.0)+(ball->center->x-ball->radius))/(velocityInMills*abs(ball->direction->x));
    }
    else if (ball->direction->x==0) {
        xTime = -1;
    }
    else {
        xTime = ((billiardSize/2.0)-(ball->center->x+ball->radius))/(velocityInMills*abs(ball->direction->x));
    }
    if (ball->direction->z<0) {
        zTime = (-(-billiardSize/2.0)+(ball->center->z-ball->radius))/(velocityInMills*abs(ball->direction->z));
    }
    else if (ball->direction->z==0) {
        zTime = -1;
    }
    else {
        zTime = ((billiardSize/2.0)-(ball->center->z+ball->radius))/(velocityInMills*abs(ball->direction->z));
    }
    if (xTime==-1 || zTime==-1) {
        timeTillCollisionMills = xTime+zTime+1;
    }
    else {
        timeTillCollisionMills = min(xTime,zTime);
    }
}

void Timer(int value) {
    if (!isManual && value==dirChanged) {
        double distance = 0;
        if (timeIntervalMills<=timeTillCollisionMills) {
            distance = timeIntervalMills*velocityInMills;
            timeTillCollisionMills -= timeIntervalMills;
        }
        else {
            distance = timeTillCollisionMills*velocityInMills;
            timeTillCollisionMills = 0;
        }
        double angle = distance/ball->radius;
        angle = -angle;
        change(angle,distance,ball->direction,ball->front);
        change(angle,distance,ball->direction,ball->up);
        change(angle,distance,ball->direction,ball->right);
        ball->center->x += distance*ball->direction->x;
        ball->center->y += distance*ball->direction->y;
        ball->center->z += distance*ball->direction->z;
        if (timeTillCollisionMills==0) {
            double error = 0;   //0.000001;
            if (((-(-billiardSize/2.0)+(ball->center->x-ball->radius))<=error || ((billiardSize/2.0)-(ball->center->x+ball->radius))<=error) && abs(ball->direction->x)>0) {
                ball->direction->x = -ball->direction->x;
            }
            if (((-(-billiardSize/2.0)+(ball->center->z-ball->radius))<=error || ((billiardSize/2.0)-(ball->center->z+ball->radius))<=error) && abs(ball->direction->z)>0) {
                ball->direction->z = -ball->direction->z;
            }
            updateTimeTillCollisionMills();
        }
        if (timeIntervalMills<=timeTillCollisionMills) {
            glutTimerFunc(timeIntervalMills, Timer, dirChanged);
        }
        else {
            glutTimerFunc(timeTillCollisionMills, Timer, dirChanged);
        }
        glutPostRedisplay();
    }
}


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
    else if (key=='l') {
        GLdouble x,z;
        double angle = 5;
        x = -ball->direction->z;
        z = ball->direction->x;
        ball->direction->x = ball->direction->x*cos(angle*(M_PI/180.0))+x*sin(angle*(M_PI/180.0));
        ball->direction->z = ball->direction->z*cos(angle*(M_PI/180.0))+z*sin(angle*(M_PI/180.0));
        if (!isManual) {
            dirChanged = (dirChanged+1)%10000;
            updateTimeTillCollisionMills();
            if (timeIntervalMills<=timeTillCollisionMills) {
                glutTimerFunc(timeIntervalMills, Timer, dirChanged);
            }
            else {
                glutTimerFunc(timeTillCollisionMills, Timer, dirChanged);
            }
        }
    }
    else if (key=='j') {
        GLdouble x,z;
        double angle = 5;
        angle = -angle;
        x = -ball->direction->z;
        z = ball->direction->x;
        ball->direction->x = ball->direction->x*cos(angle*(M_PI/180.0))+x*sin(angle*(M_PI/180.0));
        ball->direction->z = ball->direction->z*cos(angle*(M_PI/180.0))+z*sin(angle*(M_PI/180.0));
        if (!isManual) {
            dirChanged = (dirChanged+1)%10000;
            updateTimeTillCollisionMills();
            if (timeIntervalMills<=timeTillCollisionMills) {
                glutTimerFunc(timeIntervalMills, Timer, dirChanged);
            }
            else {
                glutTimerFunc(timeTillCollisionMills, Timer, dirChanged);
            }
        }
    }
    else if (isManual && key=='i') {
        double distance = 1;
        while (distance>0) {
            double xDis,zDis;
            if (ball->direction->x<0) {
                xDis = (-(-billiardSize/2.0)+(ball->center->x-ball->radius));
            }
            else {
                xDis = ((billiardSize/2.0)-(ball->center->x+ball->radius));
            }
            if (ball->direction->z<0) {
                zDis = (-(-billiardSize/2.0)+(ball->center->z-ball->radius));
            }
            else {
                zDis = ((billiardSize/2.0)-(ball->center->z+ball->radius));
            }
            int collision = 0;
            double prevDis = distance;
            if (xDis<=distance*abs(ball->direction->x) || zDis<=distance*abs(ball->direction->z)) {
                if (xDis/abs(ball->direction->x)<zDis/abs(ball->direction->z)) {
                    collision = 1;
                    distance = xDis/abs(ball->direction->x);
                }
                else if (xDis/abs(ball->direction->x)==zDis/abs(ball->direction->z)) {
                    collision = 3;
                    distance = xDis/abs(ball->direction->x);
                }
                else {
                    collision = 2;
                    distance = zDis/abs(ball->direction->z);
                }
            }
            double angle = distance/ball->radius;
            angle = -angle;
            change(angle,distance,ball->direction,ball->front);
            change(angle,distance,ball->direction,ball->up);
            change(angle,distance,ball->direction,ball->right);
            ball->center->x += distance*ball->direction->x;
            ball->center->y += distance*ball->direction->y;
            ball->center->z += distance*ball->direction->z;
            if (collision==1) {
                ball->direction->x = -ball->direction->x;
            }
            else if (collision==2) {
                ball->direction->z = -ball->direction->z;
            }
            else if (collision==3) {
                ball->direction->x = -ball->direction->x;
                ball->direction->z = -ball->direction->z;
            }
            distance = prevDis-distance;
        }
    }
    else if (isManual && key=='k') {
        double distance = 1;
        while (distance>0) {
            double xDis,zDis;
            if (ball->direction->x>0) {
                xDis = (-(-billiardSize/2.0)+(ball->center->x-ball->radius));
            }
            else {
                xDis = ((billiardSize/2.0)-(ball->center->x+ball->radius));
            }
            if (ball->direction->z>0) {
                zDis = (-(-billiardSize/2.0)+(ball->center->z-ball->radius));
            }
            else {
                zDis = ((billiardSize/2.0)-(ball->center->z+ball->radius));
            }
            int collision = 0;
            double prevDis = distance;
            if (xDis<=distance*abs(ball->direction->x) || zDis<=distance*abs(ball->direction->z)) {
                if (xDis/abs(ball->direction->x)<zDis/abs(ball->direction->z)) {
                    collision = 1;
                    distance = xDis/abs(ball->direction->x);
                }
                else if (xDis/abs(ball->direction->x)==zDis/abs(ball->direction->z)) {
                    collision = 3;
                    distance = xDis/abs(ball->direction->x);
                }
                else {
                    collision = 2;
                    distance = zDis/abs(ball->direction->z);
                }
            }
            double angle = distance/ball->radius;
            //angle = -angle;
            change(angle,distance,ball->direction,ball->front);
            change(angle,distance,ball->direction,ball->up);
            change(angle,distance,ball->direction,ball->right);
            ball->center->x -= distance*ball->direction->x;
            ball->center->y -= distance*ball->direction->y;
            ball->center->z -= distance*ball->direction->z;
            if (collision==1) {
                ball->direction->x = -ball->direction->x;
            }
            else if (collision==2) {
                ball->direction->z = -ball->direction->z;
            }
            else if (collision==3) {
                ball->direction->x = -ball->direction->x;
                ball->direction->z = -ball->direction->z;
            }
            distance = prevDis-distance;
        }
    }
    else if (key==' ') {
        isManual = !isManual;
        if (!isManual) {
            updateTimeTillCollisionMills();
            if (timeIntervalMills<=timeTillCollisionMills) {
                glutTimerFunc(timeIntervalMills, Timer, dirChanged);
            }
            else {
                glutTimerFunc(timeTillCollisionMills, Timer, dirChanged);
            }
        }
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

int main(int argc, char** argv) {
    glutInit(&argc, argv);         
    glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB);
    glutInitWindowSize(800,800);   
    glutInitWindowPosition(50, 50); 
    glutCreateWindow("Rolling Ball");   
    glutDisplayFunc(display);
    glutKeyboardFunc(keyboardListener);
    glutSpecialFunc(specialKeyListener);
    glutReshapeFunc(reshape);
    initGL();
    glutMainLoop();                 
    return 0;
}