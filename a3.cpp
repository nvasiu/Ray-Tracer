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

int windowDimension = 400;
int fovy = 70;
int rayTreeDepth;
int numSpheres;
Sphere sphereList[10];
int numPlanes;
Plane planeList[5];
int numLights;
Light lightList[5];
Color backgroundColor;

Color trace(Ray r, int depth);
Color phongLightingSphere(Ray r, Sphere s, double t, int depth);
Color phongLightingPlane(Ray r, Plane p, double t, int depth);
void render();

Color phongLightingSphere(Ray r, Sphere s, double t, int depth) {

	Color I = backgroundColor*s.kd; // Ambient Light

	// Reflection
	Vector P = r.p0 + r.v * t;
	Vector N = s.normal(P);
	N = Vector::normalize(N);
	Vector R = N * 2 * Vector::dot(N, r.v) - r.v;
	R = Vector::normalize(R);
	Ray reflectedRay( P + (N*1e-6), R*-1);

	Color reflectedColor = trace(reflectedRay, depth-1);

	I = I + reflectedColor * s.kr;

	// Refraction
	double c1 = Vector::dot(N, r.v);

	double etai = 1, etat = s.refractionIndex;
	double eta = etai / etat;

	double c2 = sqrt( 1 - (eta*eta) * (1 - (c1*c1)) );

	Vector Rr = (r.v * eta) + (N * (eta*c1-c2));
	Rr = Vector::normalize(Rr);

	Ray refractedRay( P + (N*1e-6), Rr );

	Color refractedColor = trace(refractedRay, depth-1);

	I = I + refractedColor * s.kt;

	// For each light source
	for (int i = 0; i < numLights; i++) {
		
		Color sum;

		// Diffuse
		Vector L = P - lightList[i].location;
		L = Vector::normalize(L);
		double NL = Vector::dot(N, L);

		if (NL > 0) {
			sum = sum + (s.kd*(NL/(Vector::getLength(N)*Vector::getLength(L))));
		}

		// Specular
		Vector Ri = N * 2 * NL - L;
		Ri = Vector::normalize(Ri);
		double RiV = Vector::dot(Ri, r.v);

		if (RiV > 0) {
			sum = sum + (s.ks*pow((RiV/(Vector::getLength(Ri)*Vector::getLength(r.v))),s.q));
		}

		// Shadow
		double Si = 0;
		Vector revL = lightList[i].location - P;
		Ray shadowRay(P, revL);
		Color firstObjectToIntersect = trace(shadowRay, -1);
		if ( (firstObjectToIntersect == backgroundColor) ) Si = 1;

		I = I + lightList[i].intensity * Si * sum;

	}

	return I;

}

Color phongLightingPlane(Ray r, Plane p, double t, int depth) {

	Color I = backgroundColor*p.kd; // Ambient Light

	// Reflection
	Vector P = r.p0 + r.v * t;
	Vector N(-p.A, -p.B, -p.C);
	N = Vector::normalize(N);
	Vector R = N * 2 * Vector::dot(N, r.v) - r.v;
	R = Vector::normalize(R);
	Ray reflectedRay( P + (N*1e-6), R*-1);

	Color reflectedColor = trace(reflectedRay, depth-1);

	I = I + reflectedColor * p.kr;

	// Refraction
	double c1 = Vector::dot(N, r.v);

	double etai = 1, etat = p.refractionIndex;
	double eta = etai / etat;

	double c2 = sqrt( 1 - (eta*eta) * (1 - (c1*c1)) );

	Vector Rr = (r.v * eta) + (N * (eta*c1-c2));
	Rr = Vector::normalize(Rr);

	Ray refractedRay( P + (N*1e-6), Rr );

	Color refractedColor = trace(refractedRay, depth-1);

	I = I + refractedColor * p.kt;

	// For each light source
	for (int i = 0; i < numLights; i++) {
		
		Color sum;

		// Diffuse
		Vector L = P - lightList[i].location;
		L = Vector::normalize(L);
		double NL = Vector::dot(N, L);

		if (NL > 0) {
			sum = sum + (p.kd*(NL/(Vector::getLength(N)*Vector::getLength(L))));
		}

		// Specular
		Vector Ri = N * 2 * NL - L;
		Ri = Vector::normalize(Ri);
		double RiV = Vector::dot(Ri, r.v);

		if (RiV > 0) {
			sum = sum + (p.ks*pow((RiV/(Vector::getLength(Ri)*Vector::getLength(r.v))),p.q));
		}

		// Shadow
		double Si = 0;
		Vector revL = lightList[i].location - P;
		Ray shadowRay(P, revL);
		Color firstObjectToIntersect = trace(shadowRay, -1);
		if ( (firstObjectToIntersect == backgroundColor) ) Si = 1;

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
			return phongLightingSphere(r, *firstSphere, firstT, depth);
		}
	} else if (firstObject == 2) {
		if (depth == -1) {
			return firstPlane->kd;
		} else {
			return phongLightingPlane(r, *firstPlane, firstT, depth);
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
			vector = Vector::normalize(vector);

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

	int scene = atoi(argv[1]);

	switch (scene) {
	
		case 1:
			rayTreeDepth = 3;
			backgroundColor = Color(0.2,0.2,0.2);
			numSpheres = 3;
			numPlanes = 0;
			numLights = 1;

			lightList[0] = Light(Vector(-100,100,100), Color(2,2,2));

			sphereList[0] = Sphere( Vector(1.5,1,-6), 1, Color(0.5,0.0,0.0), Color(1,1,1), 100, 0.0, 0.0, 0.0);
			sphereList[2] = Sphere( Vector(0,-1,-6), 1, Color(0.0,0.5,0.0), Color(1,1,1), 30, 0.0, 0.0, 0.0);
			sphereList[1] = Sphere( Vector(-1.5,1,-6), 1, Color(0.0,0.0,0.5), Color(1,1,1), 5, 0.0, 0.0, 0.0);
		break;
		
		case 2:
			rayTreeDepth = 3;
			backgroundColor = Color(0.0,0.0,0.0);
			numSpheres = 3;
			numPlanes = 1;
			numLights = 2;
			
			lightList[0] = Light(Vector(-100,170,100), Color(1,0.5,0.5));
			lightList[1] = Light(Vector(100,170,100), Color(0.5,0.5,1));

			planeList[0] = Plane(0.0, 1.0, 0.0, -1.0, Color(0.2,0.2,0.3), Color(0.3,0.3,0.3), 1, 0.0, 0.0, 1.5);

			sphereList[0] = Sphere( Vector(0,0.5,-5), 1.5, Color(0.2,0.2,0.2), Color(1,1,1), 10, 1.0, 0, 0);
			sphereList[1] = Sphere( Vector(-0.75,-0.75,-2.5), .25, Color(0,0.5,0), Color(1,1,1), 10, 0.5, 0.0, 1.5);
			sphereList[2] = Sphere( Vector(0.75,-0.75,-2.5), .25, Color(0,0,0.5), Color(1,1,1), 10, 0.5, 0.0, 1.5);
		break;

		case 3:
			rayTreeDepth = 3;
			backgroundColor = Color(0.3,0.3,0.3);
			numSpheres = 3;
			numPlanes = 1;
			numLights = 1;
			
			lightList[0] = Light(Vector(-100,100,100), Color(1,1,0.9));

			planeList[0] = Plane(0.0, 1.0, 0.0, -1.0, Color(0.2,0.2,0.3), Color(0.3,0.3,0.3), 1, 0.0, 0.0, 1.5);

			sphereList[0] = Sphere( Vector(1,-0.25,-4), 0.75, Color(0.1,0.5,0.5), Color(0.9,0.9,0.9), 30, 0.5, 1.0, 1.2);
			sphereList[1] = Sphere( Vector(-1,-0.25,-4), 0.75, Color(0.05,0.1,0.0), Color(1,1,1), 10, 0.2, 1.0, 1.0);
			sphereList[2] = Sphere( Vector(0,0,-6), 1.0, Color(1,0,1), Color(1,1,1), 10, 0.0, 0.0, 1.0);
		break;

		case 4:
			rayTreeDepth = 3;
			backgroundColor = Color(0.0,0.0,0.0);
			numSpheres = 7;
			numPlanes = 0;
			numLights = 2;
			
			lightList[0] = Light(Vector(50,100,100), Color(1,1,0.9));

			sphereList[0] = Sphere( Vector(1,-5,-20), 5, Color(1,0.4,0.0), Color(1,1,1), 10, 1.0, 0, 0.0);
			sphereList[1] = Sphere( Vector(-5,-2,-30), 1, Color(1,1,1), Color(1,1,1), 100, 0.1, 0.0, 1.5);
			sphereList[2] = Sphere( Vector(-6,-1,-25), 1, Color(0,0.5,1), Color(1,1,1), 10, 0.5, 0.0, 1.5);
			sphereList[3] = Sphere( Vector(-6,0.5,-20), 1, Color(1,0,0), Color(1,1,1), 10, 0.0, 1.0, 1.5);
			sphereList[4] = Sphere( Vector(-4.5,2,-15), 1, Color(1,0.8,0), Color(1,1,1), 10, 0.0, 1.0, 1.0);
			sphereList[5] = Sphere( Vector(-3,3,-13), 1, Color(0,1,0), Color(1,1,1), 10, 1.0, 0.0, 1.5);
			sphereList[6] = Sphere( Vector(-1,3,-10), 1, Color(1,0.5,1), Color(1,1,1), 10, 0.5, 1.0, 1.0);
		break;

	}

	glutDisplayFunc(display);
	glutMainLoop();	

}
