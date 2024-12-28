#include<bits/stdc++.h>
#include "1905016_classes.h"
#include "bitmap_image.hpp"
#ifdef __linux__
    #include<GL/glut.h>
#elif WIN32
    #include <windows.h>
    #include <glut.h>
#endif
using namespace std;


Camera* camera;
int recursionLevel,imageWidth,imageHeight;
int captureCount = 1;

vector<Object*> objects;
vector<PointLight*> pointLights;
vector<SpotLight*> spotLights;


void initGL() {
    camera = new Camera(new Vector(0,0,500),new Vector(0,0,-1),new Vector(0,1,0),new Vector(1,0,0),50);
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
    glEnable(GL_DEPTH_TEST);
}

bool showIntersected = false;
bool showObject = true;
bool showRays = false;

void display() {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    gluLookAt(camera->eye->x,camera->eye->y,camera->eye->z,	camera->eye->x+camera->refDistance*camera->front->x,camera->eye->y+camera->refDistance*camera->front->y,camera->eye->z+camera->refDistance*camera->front->z,camera->up->x,camera->up->y,camera->up->z);

    // axis
    glColor3f(1.0, 0.0, 0.0);
    glBegin(GL_LINES);{
        glVertex3f(100,0,0);
        glVertex3f(-100,0,0);
    }glEnd();
    glColor3f(0.0, 1.0, 0.0);
    glBegin(GL_LINES);{
        glVertex3f(0,-100,0);
        glVertex3f(0, 100,0);
    }glEnd();
    glColor3f(0.0, 0.0, 1.0);
    glBegin(GL_LINES);{
        glVertex3f(0,0, 100);
        glVertex3f(0,0,-100);
    }glEnd();

    // draw objects
    for (int i=0;i<objects.size();i++) {
        glPushMatrix();
        objects[i]->draw();
        glPopMatrix();
    }
    // point lights
    for (int i=0;i<pointLights.size();i++) {
        glPushMatrix();
        pointLights[i]->draw();
        glPopMatrix();
    }
    // spot lights
    for (int i=0;i<spotLights.size();i++) {
        glPushMatrix();
        spotLights[i]->draw();
        glPopMatrix();
    }

    glutSwapBuffers();
}

void reshape(GLsizei width, GLsizei height) {
   if (height == 0) height = 1;
   GLfloat aspect = (GLfloat)width / (GLfloat)height;
   glViewport(0, 0, width, height);
   glMatrixMode(GL_PROJECTION);
   glLoadIdentity();
   gluPerspective(45.0f, aspect, 0.1f, 5000.0f);
}

unsigned char getColor(double x) {
    unsigned char c;
    if (x<0) x = 0;
    else if (x>1.0) x = 1;
    c = x*255;
    return c;
}

double windowHeight = 400, windowWidth = 400;
double viewAngle = 45;

void capture() {
    // initialize and set background color
    bitmap_image img(imageWidth,imageHeight);
    for (int i=0;i<imageWidth;i++) {
        for (int j=0;j<imageHeight;j++) {
            img.set_pixel(i,j,0,0,0);
        }
    }
    double planeDistance = (windowHeight/2.0)/tan((viewAngle/2.0)*(M_PI/180.0));
    Vector topleft;
    topleft = *(camera->eye) + camera->front->scalarProduct(planeDistance) - camera->right->scalarProduct(windowWidth/2.0) + camera->up->scalarProduct(windowHeight/2.0);
    double du = windowWidth/imageWidth;
    double dv = windowHeight/imageHeight;
    topleft = topleft + camera->right->scalarProduct(0.5*du) - camera->up->scalarProduct(0.5*dv);
    for (int i=0;i<imageWidth;i++) {
        for (int j=0;j<imageHeight;j++) {
            double t=-1, tMin;
            Object* obj = nullptr;
            Vector curPixel;
            curPixel = topleft + camera->right->scalarProduct(i*du) - camera->up->scalarProduct(j*dv);
            Ray ray(*(camera->eye),curPixel-*(camera->eye));
            Color dummyColor;
            for (int k=0;k<objects.size();k++) {
                t = objects[k]->intersect(ray,dummyColor,0);
                if (t>0) {
                    if (t<tMin || obj==nullptr) {
                        tMin = t;
                        obj = objects[k];
                    }
                }
            }
            if (obj!=nullptr) {
                Color color;
                obj->intersect(ray,color,1);
                img.set_pixel(i,j,getColor(color.r),getColor(color.g),getColor(color.b));
            }
        }
    }
    //save image
    img.save_image("Output_1"+to_string(captureCount)+".bmp");
    cout<<"image "<<captureCount<<" captured"<<endl;
    captureCount++;
}

double rate = 5;
double angle = 3;

void keyboardListener(unsigned char key, int x, int y) {
    if (key=='0') {
        capture();
    }
    else if (key=='1') {
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

void loadData() {
    ifstream in("scene.txt");
    in>>recursionLevel;
    in>>imageWidth;
    imageHeight = imageWidth;
    int count;
    //objects
    in>>count;
    for (int i=0;i<count;i++) {
        string type;
        double a,b,c,d;
        int s;
        Object *object;
        in>>type;
        if (type=="sphere") {
            in>>a>>b>>c;
            Vector center(a,b,c);
            double radius;
            in>>radius;
            object = new Sphere(center,radius);
        }
        else if (type=="triangle") {
            in>>a>>b>>c;
            Vector A(a,b,c);
            in>>a>>b>>c;
            Vector B(a,b,c);
            in>>a>>b>>c;
            Vector C(a,b,c);
            object = new Triangle(A,B,C);
        }
        else if (type=="general") {
            double A,B,C,D,E,F,G,H,I,J;
            in>>A>>B>>C>>D>>E>>F>>G>>H>>I>>J;
            double length, width, height;
            in>>a>>b>>c>>length>>width>>height;
            Vector ref(a,b,c);
            object = new General(A,B,C,D,E,F,G,H,I,J,ref,length,width,height);

        }
        else {
            cout<<"invalid type "<<i<<endl;
            object = new Object();
        }
        in>>a>>b>>c;
        object->setColor(a,b,c);
        in>>a>>b>>c>>d;
        object->setCoEfficients(a,b,c,d);
        in>>s;
        object->setShine(s);
        objects.push_back(object);
    }
    // point lights
    in>>count;
    for (int i=0;i<count;i++) {
        double a,b,c;
        in>>a>>b>>c;
        Vector pos(a,b,c);
        in>>a>>b>>c;
        Color color(a,b,c);
        PointLight *p = new PointLight(pos,color);
        pointLights.push_back(p);
    }
    // spot lights
    in>>count;
    for (int i=0;i<count;i++) {
        double a,b,c;
        in>>a>>b>>c;
        Vector pos(a,b,c);
        in>>a>>b>>c;
        Color color(a,b,c);
        in>>a>>b>>c;
        Vector dir(a,b,c);
        double angle;
        in>>angle;
        SpotLight *s = new SpotLight(pos,dir,angle,color);
        spotLights.push_back(s);
    }
    in.close();
    // floor
    {
        Object *floor = new Floor(1000.0,20.0);
        floor->setColor(1,1,1);
        floor->setCoEfficients(.4,0.2,0.1,0.5);
        floor->setShine(5);
        objects.push_back(floor);
    }
}

int main(int argc, char** argv) {
    loadData();
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB);
    glutInitWindowSize(800,800);
    glutInitWindowPosition(50, 50);
    glutCreateWindow("Raytracing");
    glutDisplayFunc(display);
    glutKeyboardFunc(keyboardListener);
    glutSpecialFunc(specialKeyListener);
    glutReshapeFunc(reshape);
    initGL();
    glutMainLoop();
    //free
    delete camera;
    //free objects
    for (int i=0;i<objects.size();i++) {
        delete objects[i];
    }
    objects.clear();
    //free point lights
    for (int i=0;i<pointLights.size();i++) {
        delete pointLights[i];
    }
    pointLights.clear();
    //free spot lights
    for (int i=0;i<spotLights.size();i++) {
        delete spotLights[i];
    }
    spotLights.clear();
    return 0;
}
