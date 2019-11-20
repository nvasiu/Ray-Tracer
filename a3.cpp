#include <chrono>
#include <thread>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <math.h>
#include <time.h>
#include <GL/glut.h>

using namespace std;

class Point {

	public:
	float x = 0;
	float y = 0;
	float z = 0;

	Point() {}

	Point(float xcoord, float ycoord, float zcoord) {
		x = xcoord;
		y = ycoord;
		z = zcoord;
	}

	static Point subtract(Point p1, Point p2) {
		float newx = p1.x - p2.x;
		float newy = p1.y - p2.y;
		float newz = p1.z - p2.z;
		return Point(newx,newy,newz);
	}

};

class Color {

	public:
	float r = 0;
	float g = 0;
	float b = 0;

	Color() {}

	Color(float red, float green, float blue) {
		r = red;
		g = green;
		b = blue;
	}

};

class Ray {

	public:
	Point p0;
	Point v;

	Ray() {}

	Ray(Point start, Point vector) {
		p0 = start;
		v = vector;
	}

};

class Sphere {

	public:
	Point center;
	float radius;
	Color kd;
	Color ks;
	int q;
	float kr;
	float kt;
	float refractionIndex;

	Sphere() {
		center = Point(0,0,0);
		radius = 1;
		kd = Color(0,0,0);
		ks = Color(0,0,0);
		q = 1;
		kr = 1;
		kt = 1;
		refractionIndex = 1;
	}

	Sphere(Point newcenter, float newr, Color newkd, Color newks, int newq, float newkr, float newkt, float newrefraction) {
		center = newcenter;
		radius = newr;
		kd = newkd;
		ks = newks;
		q = newq;
		kr = newkr;
		kt = newkt;
		refractionIndex = newrefraction;
	}

	float intersects(Ray testRay) {

		Point OC = Point::subtract(testRay.p0, center);
		
		float a = testRay.v.x*testRay.v.x + testRay.v.y*testRay.v.y + testRay.v.z*testRay.v.z;
		float b = 2*(OC.x*testRay.v.x + OC.y*testRay.v.y + OC.z*testRay.v.z);
		float c = (OC.x*OC.x + OC.y*OC.y + OC.z*OC.z) - (radius*radius);

		float discriminant = b*b-4*a*c;
		if (discriminant < 0) {
			return -1;
		}

		float t1 = (-b - sqrt(discriminant))/2;
		float t2 = (-b + sqrt(discriminant))/2;

		if (t1 < t2) {return t1;}
		else {return t2;}

	}

};

class Plane {

	public:
	float A;
	float B;
	float C;
	float D;
	Color kd;
	Color ks;
	int q;
	float kr;
	float kt;
	float refractionIndex;

	Plane() {
		A = 0;
		B = 0;
		C = 0;
		D = 0;
		kd = Color(0,0,0);
		ks = Color(0,0,0);
		q = 1;
		kr = 1;
		kt = 1;
		refractionIndex = 1;
	}

	Plane(float newA, float newB, float newC, float newD, Color newkd, Color newks, int newq, float newkr, float newkt, float newrefraction) {
		A = newA;
		B = newB;
		C = newC;
		D = newD;
		kd = newkd;
		ks = newks;
		q = newq;
		kr = newkr;
		kt = newkt;
		refractionIndex = newrefraction;
	}

};

class Light {

	public:
	Point location;
	Color intensity;

	Light() {
		location = Point(0,0,0);
		intensity = Color(0,0,0);
	}

	Light(Point newloc, Color newint) {
		location = newloc;
		intensity = newint;
	}

};

int windowDimension = 400;
int fovy = 45;
int rayTreeDepth = 4;
Sphere sphereList[3];
Plane planeList[1];
Light lightList[1];
Color backgroundColor(0.7,0.7,0.9);

Color trace(Ray r, int depth) {

	float firstT = 10000000;
	Sphere *firstSphere = nullptr;
	
	for (int i = 0; i < sizeof(sphereList)/sizeof(Sphere); i++) {
		float t = sphereList[i].intersects(r);
		if (t > 0 && t < firstT) {
			firstT = t;
			firstSphere = &sphereList[i];
		}
	}
	
	if (firstSphere == nullptr) {
		return backgroundColor;
	} else {
		return firstSphere->kd;
	}
	
}

void render() {
	
	for (int x = 0; x < windowDimension; x++) {
		for (int y = 0; y < windowDimension; y++) {
			Point start(0,0,0);
			
			double vx = tan(fovy)/x;
			double vy = tan(fovy)/y;

			Point vector(vx,vy,-1);

			Ray newray(start, vector);
			
			Color pixelColor = trace(newray, rayTreeDepth);
			
			glBegin(GL_POINTS);
				glColor3f(pixelColor.r, pixelColor.g, pixelColor.b);
				glVertex2i(x,y);
			glEnd();
			
		}
	}

}

void display(void) {
	
	glClear(GL_COLOR_BUFFER_BIT);

	render();

	glFlush();
	glutSwapBuffers();

}

int main(int argc, char** argv) {

	// Window
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
	glutInitWindowPosition(500, 200);
	glutInitWindowSize(windowDimension, windowDimension);
	glutCreateWindow("Assignment 3");

	glClearColor(1.0, 1.0, 1.0, 0.0);
	glMatrixMode(GL_PROJECTION);
	gluOrtho2D(0.0, windowDimension, 0.0, windowDimension);

	sphereList[0] = Sphere( Point(0,0,-50), 1, Color(1,0,0), Color(1,0,0), 1, 1, 1, 1);
	sphereList[1] = Sphere( Point(-3,0,-50), 1, Color(0,1,0), Color(1,0,0), 1, 1, 1, 1);
	sphereList[2] = Sphere( Point(3,0,-50), 1, Color(0,0,1), Color(1,0,0), 1, 1, 1, 1);

	glutDisplayFunc(display);
	glutMainLoop();	

}
