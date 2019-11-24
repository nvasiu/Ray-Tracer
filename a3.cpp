#include <chrono>
#include <thread>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <math.h>
#include <time.h>
#include <GL/glut.h>
#include <iostream>

#define inf 1e7

using namespace std;

class Vector {

	public:
	double x = 0;
	double y = 0;
	double z = 0;

	Vector() {}

	Vector(double xcoord, double ycoord, double zcoord) {
		x = xcoord;
		y = ycoord;
		z = zcoord;
	}

	Vector operator - (Vector p2) {
		return Vector(x-p2.x,y-p2.y,z-p2.z);
	}

	Vector operator + (Vector p2) {
		return Vector(x+p2.x,y+p2.y,z+p2.z);
	}

	Vector operator * (double f) {
		return Vector(x*f,y*f,z*f);
	}

	Vector operator * (Vector p2) {
		return Vector(x*p2.x,y*p2.y,z*p2.z);
	}

	static double dot(Vector p1, Vector p2) {
		return p1.x*p2.x + p1.y*p2.y + p1.z*p2.z;
	}

	static double getLength(Vector p1) {
		return sqrt(p1.x*p1.x + p1.y*p1.y + p1.z*p1.z);
	}

	static Vector normalize(Vector p1) {
		double length = sqrt(p1.x*p1.x + p1.y*p1.y + p1.z*p1.z);
		return Vector(p1.x/length, p1.y/length, p1.z/length);
	}

};

class Color {

	public:
	double r = 0;
	double g = 0;
	double b = 0;

	Color() {}

	Color(double red, double green, double blue) {
		r = red;
		g = green;
		b = blue;
	}

	Color operator * (double f) {
		double newr = r*f;
		double newg = g*f;
		double newb = b*f;
		return Color(newr, newg, newb);
	}

	Color operator * (Color f) {
		double newr = r*f.r;
		double newg = g*f.g;
		double newb = b*f.b;
		return Color(newr, newg, newb);
	}

	Color operator + (double f) {
		double newr = r+f;
		double newg = g+f;
		double newb = b+f;
		return Color(newr, newg, newb);
	}

	Color operator + (Color f) {
		double newr = r+f.r;
		double newg = g+f.g;
		double newb = b+f.b;
		return Color(newr, newg, newb);
	}

	bool operator == (Color f) {
		bool equal = true;
		if (r != f.r) equal = false;
		if (g != f.g) equal = false;
		if (b != f.b) equal = false;
		return equal;
	}

};

class Ray {

	public:
	Vector p0;
	Vector v;

	Ray() {}

	Ray(Vector start, Vector vector) {
		p0 = start;
		v = vector;
	}

};

class Sphere {

	public:
	Vector center;
	double radius;
	Color kd;
	Color ks;
	double q;
	double kr;
	double kt;
	double refractionIndex;

	Sphere() {
		center = Vector(0,0,0);
		radius = 1;
		kd = Color(0,0,0);
		ks = Color(0,0,0);
		q = 1;
		kr = 1;
		kt = 1;
		refractionIndex = 1;
	}

	Sphere(Vector newcenter, double newr, Color newkd, Color newks, int newq, double newkr, double newkt, double newrefraction) {
		center = newcenter;
		radius = newr;
		kd = newkd;
		ks = newks;
		q = newq;
		kr = newkr;
		kt = newkt;
		refractionIndex = newrefraction;
	}

	double intersects(Ray testRay) {

		Vector OC = testRay.p0-center;
		
		double a = Vector::dot(testRay.v, testRay.v);
		double b = 2*Vector::dot(OC, testRay.v);
		double c = Vector::dot(OC, OC) - (radius*radius);

		double discriminant = b*b-4*a*c;
		if (discriminant < 0) {
			return inf;
		}

		double t1 = (-b - sqrt(discriminant))/2;
		double t2 = (-b + sqrt(discriminant))/2;

		if (t1 < t2) {return t1;}
		else {return t2;}

	}

	Vector normal(Vector p) {
		return (p-center) * (-1/radius);
	}

};

class Plane {

	public:
	double A;
	double B;
	double C;
	double D;
	Color kd;
	Color ks;
	double q;
	double kr;
	double kt;
	double refractionIndex;

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

	Plane(double newA, double newB, double newC, double newD, Color newkd, Color newks, int newq, double newkr, double newkt, double newrefraction) {
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

	double intersects(Ray testRay) {

		Vector n(A,B,C);
		double denominator = Vector::dot(n, testRay.v);
	
		if (abs(denominator) <= 1e-6) return inf;

		double t = ( Vector::dot(n, testRay.p0) + D) / denominator;

		if (t <= 1e-6) return inf;

		return t;

	}

};

class Light {

	public:
	Vector location;
	Color intensity;

	Light() {
		location = Vector(0,0,0);
		intensity = Color(0,0,0);
	}

	Light(Vector newloc, Color newint) {
		location = newloc;
		intensity = newint;
	}

};

int windowDimension = 1000;
int fovy = 70;
int rayTreeDepth = 4;
int numSpheres;
Sphere sphereList[10];
int numPlanes;
Plane planeList[5];
int numLights;
Light lightList[5];
Color backgroundColor(0.2,0.2,0.2);

Color trace(Ray r, int depth);
Color phongLightingSphere(Ray r, Sphere s, double t);
Color phongLightingPlane(Ray r, Plane p, double t);
void render();

Color phongLightingSphere(Ray r, Sphere s, double t) {

	Color I = backgroundColor*s.kd;

	for (int i = 0; i < numLights; i++) {
		
		Vector P = r.p0 + r.v * t;
		Vector N = s.normal(P);
		Vector L = P - lightList[i].location;

		double NL = Vector::dot(N, L);

		double Si = 0;

		Vector revL = lightList[i].location - P;
		Ray shadowRay(P, revL);
		Color firstIntersect = trace(shadowRay, -1);

		if ( (firstIntersect == backgroundColor) ) Si = 1;

		Vector R = N * 2 * NL - L;
		double RV = Vector::dot(R, r.v);
		
		Color sum;

		if (NL > 0) {
			sum = sum + (s.kd*(NL/(Vector::getLength(N)*Vector::getLength(L))));
		}

		if (RV > 0) {
			sum = sum + (s.ks*pow((RV/(Vector::getLength(R)*Vector::getLength(r.v))),s.q));
		}

		I = I + lightList[i].intensity * Si * sum;

	}

	return I;

}

Color phongLightingPlane(Ray r, Plane p, double t) {

	Color I = backgroundColor*p.kd;

	for (int i = 0; i < numLights; i++) {
		
		Vector P = r.p0 + r.v * t;
		Vector N(-p.A, -p.B, -p.C);
		Vector L = P - lightList[i].location;

		double NL = Vector::dot(N, L);

		double Si = 0;

		Vector revL = lightList[i].location - P;
		Ray shadowRay(P, revL);
		Color firstIntersect = trace(shadowRay, -1);

		if ( (firstIntersect == backgroundColor) ) Si = 1;

		Vector R = N * 2 * NL - L;
		double RV = Vector::dot(R, r.v);
		
		Color sum;

		if (NL > 0) {
			sum = sum + (p.kd*(NL/(Vector::getLength(N)*Vector::getLength(L))));
		}

		if (RV > 0) {
			sum = sum + (p.ks*pow((RV/(Vector::getLength(R)*Vector::getLength(r.v))),p.q));
		}

		I = I + lightList[i].intensity * Si * sum;

	}

	return I;

}

Color trace(Ray r, int depth) {

	if (depth == 0) {
		return backgroundColor;
	}

	double firstT = inf;
	Sphere *firstSphere = nullptr;
	Plane *firstPlane = nullptr;
	int firstObject = 0;
	
	for (int i = 0; i < numSpheres; i++) {
		double t = sphereList[i].intersects(r);
		if (t > 0 && t < firstT) {
			firstT = t;
			firstSphere = &sphereList[i];
			firstObject = 1;
		}
	}

	for (int i = 0; i < numPlanes; i++) {
		double t = planeList[i].intersects(r);
		if (t > 1 && t < firstT) {
			firstT = t;
			firstPlane = &planeList[i];
			firstObject = 2;
		}
	}

	if (firstObject == 1) {
		if (depth == -1) {
			return firstSphere->kd;
		} else {
			return phongLightingSphere(r, *firstSphere, firstT);
		}
	} else if (firstObject == 2) {
		if (depth == -1) {
			return firstPlane->kd;
		} else {
			return phongLightingPlane(r, *firstPlane, firstT);
		}
	} else {
		return backgroundColor;
	}
	
}

void render() {
	
	for (int x = 0; x < windowDimension; x++) {
		for (int y = 0; y < windowDimension; y++) {
			Vector start(0,0,0);

			double vx = (2*(x+0.5)/windowDimension) - 1;
			double vy = (2*(y+0.5)/windowDimension) - 1;
			
			Vector vector(vx*tan(fovy/2),vy*tan(fovy/2),-1);

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

	numSpheres = 3;
	numPlanes = 1;
	numLights = 1;
	lightList[0] = Light(Vector(-100,100,100), Color(2,2,2));
	lightList[1] = Light(Vector(-5,10,-5), Color(0.2,0.2,0.2));
	planeList[0] = Plane(0.0, 1.0, 0.0, -1.0, Color(0.2,0.2,0.3), Color(0.3,0.3,0.3), 1, 0.5, 0.0, 1.5);
	//planeList[1] = Plane(0.0, 1.0, 0.0, 4.0, Color(0.5,0.5,0.7), Color(1,0,0), 32, 0.5, 0.0, 1.5);
	sphereList[0] = Sphere( Vector(1,0.5,-10), 1.5, Color(0.5,0,0), Color(1,1,1), 10, 0, 0, 0);
	sphereList[1] = Sphere( Vector(-1,0,-8), 1, Color(0,0.5,0), Color(1,1,1), 10, 0, 0.0, 1.5);
	sphereList[2] = Sphere( Vector(0,-0.75,-4), .25, Color(0,0,0.5), Color(1,1,1), 10, 0, 0.0, 1.5);

	glutDisplayFunc(display);
	glutMainLoop();	

}
