#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include <windows.h>
#include<iostream>
#include <glut.h>
#include<iostream>
#include<sstream>
#include "classes.h"
#include "bitmap_image.hpp"
#include <vector>

#define pi (2*acos(0.0))
#define FOV 80
#define window_height 500
#define window_width 500

using namespace std;

double cameraHeight;
double cameraAngle;
double cameraRight;
double sphereRadius;
int drawgrid;
int drawaxes;
double angle;

int image_width;
int image_height;
int numObject;
int numLight;
struct point pos, u, r, l;

///externs
vector<Object*> objects;
vector<point> lights;
int recursion_level;

void drawAxes()
{
	if(drawaxes==1)
	{
		glBegin(GL_LINES);{
		    glColor3f(1, 0, 0);
			glVertex3f( 500,0,0);
			glVertex3f(-500,0,0);

			glColor3f(0, 1, 0);
			glVertex3f(0,-500,0);
			glVertex3f(0, 500,0);

            glColor3f(0, 0, 1);
			glVertex3f(0,0, 500);
			glVertex3f(0,0,-500);
		}glEnd();
	}
}


void drawGrid()
{
	int i;
	if(drawgrid==1)
	{
		glColor3f(0.6, 0.6, 0.6);	//grey
		glBegin(GL_LINES);{
			for(i=-8;i<=8;i++){

				if(i==0)
					continue;	//SKIP the MAIN axes

				//lines parallel to Y-axis
				glVertex3f(i*10, -90, 0);
				glVertex3f(i*10,  90, 0);

				//lines parallel to X-axis
				glVertex3f(-90, i*10, 0);
				glVertex3f( 90, i*10, 0);
			}
		}glEnd();
	}
}

void drawSquare(double a)
{
    glColor3f(1.0, 1.0, 1.0);
	glBegin(GL_QUADS);{
		glVertex3f( a, a,0);
		glVertex3f( a,-a,0);
		glVertex3f(-a,-a,0);
		glVertex3f(-a, a,0);
	}glEnd();
}


void drawCircle(double radius,int segments)
{
    int i;
    struct point points[100];
    glColor3f(0.7,0.7,0.7);
    //generate points
    for(i=0;i<=segments;i++)
    {
        points[i].x=radius*cos(((double)i/(double)segments)*2*pi);
        points[i].y=radius*sin(((double)i/(double)segments)*2*pi);
    }
    //draw segments using generated points
    for(i=0;i<segments;i++)
    {
        glBegin(GL_LINES);
        {
			glVertex3f(points[i].x,points[i].y,0);
			glVertex3f(points[i+1].x,points[i+1].y,0);
        }
        glEnd();
    }
}

void drawCone(double radius,double height,int segments)
{
    int i;
    double shade;
    struct point points[100];
    //generate points
    for(i=0;i<=segments;i++)
    {
        points[i].x=radius*cos(((double)i/(double)segments)*2*pi);
        points[i].y=radius*sin(((double)i/(double)segments)*2*pi);
    }
    //draw triangles using generated points
    for(i=0;i<segments;i++)
    {
        //create shading effect
        if(i<segments/2)shade=2*(double)i/(double)segments;
        else shade=2*(1.0-(double)i/(double)segments);
        glColor3f(shade,shade,shade);

        glBegin(GL_TRIANGLES);
        {
            glVertex3f(0,0,height);
			glVertex3f(points[i].x,points[i].y,0);
			glVertex3f(points[i+1].x,points[i+1].y,0);
        }
        glEnd();
    }
}

void drawPartCylinder(double radius, double height, int slices, int stacks)
{
    int i, j;
    double h, r;
    struct point points[100][100];
    //generate points
    for(i = 0; i <= stacks; i++)
    {
        r = radius;
        h = height * sin(((double)i / (double)stacks) * (pi/2.0));;
        for(j = 0; j <= slices; j++){
            points[i][j].x = r * cos(((double)j / (double)slices) * 0.5 * pi);
            points[i][j].y = r * sin(((double)j / (double)slices) * 0.5 * pi);
            points[i][j].z = h;
        }
    }
    //draw triangles using generated points
    glColor3f(0, 1, 0);
    for(i=0;i<stacks;i++)
    {
        for(j=0; j<slices; j++) {
            glBegin(GL_QUADS);
            {
                glVertex3f(points[i][j].x,points[i][j].y,points[i][j].z);
                glVertex3f(points[i][j+1].x,points[i][j+1].y,points[i][j+1].z);
                glVertex3f(points[i+1][j+1].x,points[i+1][j+1].y,points[i+1][j+1].z);
                glVertex3f(points[i+1][j].x,points[i+1][j].y,points[i+1][j].z);

                glVertex3f(points[i][j].x,points[i][j].y,-points[i][j].z);
                glVertex3f(points[i][j+1].x,points[i][j+1].y,-points[i][j+1].z);
                glVertex3f(points[i+1][j+1].x,points[i+1][j+1].y,-points[i+1][j+1].z);
                glVertex3f(points[i+1][j].x,points[i+1][j].y,-points[i+1][j].z);
            }
            glEnd();
        }

    }
}

void drawPartSphere(double radius,int slices,int stacks)
{
	struct point points[100][100];
	int i,j;
	double h,r;
	//generate points
	for(i=0;i<=stacks;i++)
	{
		h=radius*sin(((double)i/(double)stacks)*(pi/2.0));
		r=radius*cos(((double)i/(double)stacks)*(pi/2.0));
		for(j=0;j<=slices;j++)
		{
			points[i][j].x=r*cos(((double)j/(double)slices)*2*pi);
			points[i][j].y=r*sin(((double)j/(double)slices)*2*pi);
			points[i][j].z=h;
		}
	}
	//draw quads using generated points
	for(i=0;i<stacks;i++)
	{
        glColor3f(1, 0, 0);
		for(j=0;j<slices/4;j++)
		{
			glBegin(GL_QUADS);{
			    //upper hemisphere
				glVertex3f(points[i][j].x,points[i][j].y,points[i][j].z);
				glVertex3f(points[i][j+1].x,points[i][j+1].y,points[i][j+1].z);
				glVertex3f(points[i+1][j+1].x,points[i+1][j+1].y,points[i+1][j+1].z);
				glVertex3f(points[i+1][j].x,points[i+1][j].y,points[i+1][j].z);
			}glEnd();
		}
	}
}

void drawSphere(double radius,int slices,int stacks)
{
	struct point points[100][100];
	int i,j;
	double h,r;
	//generate points
	for(i=0;i<=stacks;i++)
	{
		h=radius*sin(((double)i/(double)stacks)*(pi/2.0));
		r=radius*cos(((double)i/(double)stacks)*(pi/2.0));
		for(j=0;j<=slices;j++)
		{
			points[i][j].x=r*cos(((double)j/(double)slices)*2*pi);
			points[i][j].y=r*sin(((double)j/(double)slices)*2*pi);
			points[i][j].z=h;
		}
	}
	//draw quads using generated points
	for(i=0;i<stacks;i++)
	{
        glColor3f((double)i/(double)stacks,(double)i/(double)stacks,(double)i/(double)stacks);
		for(j=0;j<slices;j++)
		{
			glBegin(GL_QUADS);{
			    //upper hemisphere
				glVertex3f(points[i][j].x,points[i][j].y,points[i][j].z);
				glVertex3f(points[i][j+1].x,points[i][j+1].y,points[i][j+1].z);
				glVertex3f(points[i+1][j+1].x,points[i+1][j+1].y,points[i+1][j+1].z);
				glVertex3f(points[i+1][j].x,points[i+1][j].y,points[i+1][j].z);
                //lower hemisphere
                glVertex3f(points[i][j].x,points[i][j].y,-points[i][j].z);
				glVertex3f(points[i][j+1].x,points[i][j+1].y,-points[i][j+1].z);
				glVertex3f(points[i+1][j+1].x,points[i+1][j+1].y,-points[i+1][j+1].z);
				glVertex3f(points[i+1][j].x,points[i+1][j].y,-points[i+1][j].z);
			}glEnd();
		}
	}
}

void drawcube2box()
{
    glPushMatrix();
    glTranslatef(36 - sphereRadius, 36 - sphereRadius,0);
    drawPartCylinder(sphereRadius, 36 - sphereRadius, 40, 20);
    glPopMatrix();

    glPushMatrix();
    glRotatef(90.0, 0, 0, 1);
    glTranslatef(36 - sphereRadius, 36 - sphereRadius,0);
    drawPartCylinder(sphereRadius, 36 - sphereRadius, 40, 20);
    glPopMatrix();

    glPushMatrix();
    glRotatef(-90.0, 0, 0, 1);
    glTranslatef(36 - sphereRadius, 36 - sphereRadius,0);
    drawPartCylinder(sphereRadius, 36 - sphereRadius, 40, 20);
    glPopMatrix();

    glPushMatrix();
    glRotatef(180.0, 0, 0, 1);
    glTranslatef(36 - sphereRadius, 36 - sphereRadius,0);
    drawPartCylinder(sphereRadius, 36 - sphereRadius, 40, 20);
    glPopMatrix();

    glPushMatrix();
    glRotatef(90.0, 0, 1, 0);
    glTranslatef(36 - sphereRadius, 36 - sphereRadius,0);
    drawPartCylinder(sphereRadius, 36 - sphereRadius, 40, 20);
    glPopMatrix();

    glPushMatrix();
    glRotatef(90.0, 0, 0, 1);
    glRotatef(90.0, 0, 1, 0);
    glTranslatef(36 - sphereRadius, 36 - sphereRadius,0);
    drawPartCylinder(sphereRadius, 36 - sphereRadius, 40, 20);
    glPopMatrix();

    glPushMatrix();
    glRotatef(-90.0, 0, 0, 1);
    glRotatef(90.0, 0, 1, 0);
    glTranslatef(36 - sphereRadius, 36 - sphereRadius,0);
    drawPartCylinder(sphereRadius, 36 - sphereRadius, 40, 20);
    glPopMatrix();

    glPushMatrix();
    glRotatef(180.0, 0, 0, 1);
    glRotatef(90.0, 0, 1, 0);
    glTranslatef(36 - sphereRadius, 36 - sphereRadius,0);
    drawPartCylinder(sphereRadius, 36 - sphereRadius, 40, 20);
    glPopMatrix();

    glPushMatrix();
    glRotatef(90.0, 1, 0, 0);
    glTranslatef(36 - sphereRadius, 36 - sphereRadius,0);
    drawPartCylinder(sphereRadius, 36 - sphereRadius, 40, 20);
    glPopMatrix();

    glPushMatrix();
    glRotatef(90.0, 0, 0, 1);
    glRotatef(90.0, 1, 0, 0);
    glTranslatef(36 - sphereRadius, 36 - sphereRadius,0);
    drawPartCylinder(sphereRadius, 36 - sphereRadius, 40, 20);
    glPopMatrix();

    glPushMatrix();
    glRotatef(-90.0, 0, 0, 1);
    glRotatef(90.0, 1, 0, 0);
    glTranslatef(36 - sphereRadius, 36 - sphereRadius,0);
    drawPartCylinder(sphereRadius, 36 - sphereRadius, 40, 20);
    glPopMatrix();

    glPushMatrix();
    glRotatef(180.0, 0, 0, 1);
    glRotatef(90.0, 1, 0, 0);
    glTranslatef(36 - sphereRadius, 36 - sphereRadius,0);
    drawPartCylinder(sphereRadius, 36 - sphereRadius, 40, 20);
    glPopMatrix();


    glPushMatrix();
    glTranslatef(36 - sphereRadius, 36 - sphereRadius, 36 - sphereRadius);
    drawPartSphere(sphereRadius, 40, 20);
    glPopMatrix();

    glPushMatrix();
    glRotatef(90.0, 0, 0, 1);
    glTranslatef(36 - sphereRadius, 36 - sphereRadius, 36 - sphereRadius);
    drawPartSphere(sphereRadius, 40, 20);
    glPopMatrix();

    glPushMatrix();
    glRotatef(-90.0, 0, 0, 1);
    glTranslatef(36 - sphereRadius, 36 - sphereRadius, 36 - sphereRadius);
    drawPartSphere(sphereRadius, 40, 20);
    glPopMatrix();

    glPushMatrix();
    glRotatef(180.0, 0, 0, 1);
    glTranslatef(36 - sphereRadius, 36 - sphereRadius, 36 - sphereRadius);
    drawPartSphere(sphereRadius, 40, 20);
    glPopMatrix();

    glPushMatrix();
    glRotatef(90.0, 0, 1, 0);
    glTranslatef(36 - sphereRadius, 36 - sphereRadius, 36 - sphereRadius);
    drawPartSphere(sphereRadius, 40, 20);
    glPopMatrix();


    glPushMatrix();
    glRotatef(90.0, 0, 0, 1);
    glRotatef(90.0, 0, 1, 0);
    glTranslatef(36 - sphereRadius, 36 - sphereRadius, 36 - sphereRadius);
    drawPartSphere(sphereRadius, 40, 20);
    glPopMatrix();

    glPushMatrix();
    glRotatef(-90.0, 0, 0, 1);
    glRotatef(90.0, 0, 1, 0);
    glTranslatef(36 - sphereRadius, 36 - sphereRadius, 36 - sphereRadius);
    drawPartSphere(sphereRadius, 40, 20);
    glPopMatrix();

    glPushMatrix();
    glRotatef(180.0, 0, 0, 1);
    glRotatef(90.0, 0, 1, 0);
    glTranslatef(36 - sphereRadius, 36 - sphereRadius, 36 - sphereRadius);
    drawPartSphere(sphereRadius, 40, 20);
    glPopMatrix();


    glPushMatrix();
    glRotatef(90.0, 1, 0, 0);
    glTranslatef(0, 0, -36);
    drawSquare(36 - sphereRadius);
    glPopMatrix();

    glPushMatrix();
    glRotatef(-90.0, 1, 0, 0);
    glTranslatef(0, 0, -36);
    drawSquare(36 - sphereRadius);
    glPopMatrix();

    glPushMatrix();
    glRotatef(90.0, 0, 1, 0);
    glTranslatef(0, 0, -36);
    drawSquare(36 - sphereRadius);
    glPopMatrix();

    glPushMatrix();
    glRotatef(-90.0, 0, 1, 0);
    glTranslatef(0, 0, -36);
    drawSquare(36 - sphereRadius);
    glPopMatrix();

    glPushMatrix();
    glRotatef(90.0, 0, 0, 1);
    glTranslatef(0, 0, -36);
    drawSquare(36 - sphereRadius);
    glPopMatrix();

    glPushMatrix();
    glRotatef(-90.0, 0, 0, 1);
    glTranslatef(0, 0, -36);
    drawSquare(36 - sphereRadius);
    glPopMatrix();

}
void drawSS()
{

    glColor3f(1,0,0);
    drawSquare(20);

    glRotatef(angle,0,0,1);
    glTranslatef(110,0,0);
    glRotatef(2*angle,0,0,1);
    glColor3f(0,1,0);
    drawSquare(15);

    glPushMatrix();
    {
        glRotatef(angle,0,0,1);
        glTranslatef(60,0,0);
        glRotatef(2*angle,0,0,1);
        glColor3f(0,0,1);
        drawSquare(10);
    }
    glPopMatrix();

    glRotatef(3*angle,0,0,1);
    glTranslatef(40,0,0);
    glRotatef(4*angle,0,0,1);
    glColor3f(1,1,0);
    drawSquare(5);
}

void drawLine(){
    glLineWidth(0);
    glColor3f(1.0, 0.0, 0.0);
    glBegin(GL_LINES);
    glVertex3f(10.0, 0.0, 0.0);
    glVertex3f(25, 100, 0);
    glEnd();
}

void capture()
{
    bitmap_image bitmap_image(image_width, image_height);
    double plane_distance = (window_height/2.0)/tan(FOV * pi / 360);

    struct point topleft;
    topleft.x = pos.x + l.x * plane_distance - r.x * (window_width/2.0) + u.x * (window_height/2.0);
    topleft.y = pos.y + l.y * plane_distance - r.y * (window_width/2.0) + u.y * (window_height/2.0);
    topleft.z = pos.z + l.z * plane_distance - r.z * (window_width/2.0) + u.z * (window_height/2.0);

    //printf("topleft is %f, %f, %f\n", topleft.x, topleft.y, topleft.z);

    double du = (1.0*window_width)/image_width;
    double dv = (1.0*window_height)/image_height;

    //cout<<du<<" "<<dv<<endl;

    //ofstream fout("corners.txt");
    for(int i = 0; i < image_height; i++){
        for(int j = 0; j < image_width; j++){
            //printf("%d ", i*j);
            struct point corner;
            corner.x = topleft.x + r.x * j * du - u.x * i * dv;
            corner.y = topleft.y + r.y * j * du - u.y * i * dv;
            corner.z = topleft.z + r.z * j * du - u.z * i * dv;

            //fout<<corner.x<<", "<<corner.y<<", "<<corner.z<<endl;
            struct point diff;
            diff.x = corner.x - pos.x;
            diff.y = corner.y - pos.y;
            diff.z = corner.z - pos.z;

            diff = normalize(diff);
            Ray newRay(pos, diff);

            //printf("%f %f %f\n", diff2.x, diff2.y, diff2.z);
            int nearest = -1;

            struct point dummyColor;
            double t_min = 9999;
            for(int k = 0; k < objects.size(); k++){
                double t = objects[k]->intersect(&newRay, &dummyColor, 0);

                if(t <= 0) continue;
                else if(t < t_min) {
                    t_min = t;
                    nearest = k;
                }
            }

            if(nearest != -1){
                //cout<<"not inside yet"<<endl;
                double t = objects[nearest]->intersect(&newRay, &dummyColor, 1);

                //if(dummyColor.x == 1) printf("                          yes");
                bitmap_image.set_pixel(j, i, dummyColor.x * 255, dummyColor.y * 255, dummyColor.z * 255);
            }
        }
    }

    bitmap_image.save_image("out.bmp");
    exit(0);
}

void keyboardListener(unsigned char key, int x,int y){
	switch(key){
        struct point cross;
        case '0': //call capture
            capture();
            break;
		case '1': //look left
			cross.x = u.y * r.z - u.z * r.y;
			cross.y = u.z * r.x - u.x * r.z;
			cross.z = u.x * r.y - u.y * r.x;
			r.x = r.x * cos(3 * pi / 180) + cross.x * sin(3 * pi / 180);
			r.y = r.y * cos(3 * pi / 180) + cross.y * sin(3 * pi / 180);
			r.z = r.z * cos(3 * pi / 180) + cross.z * sin(3 * pi / 180);
            r = normalize(r);

			cross.x = u.y * l.z - u.z * l.y;
			cross.y = u.z * l.x - u.x * l.z;
			cross.z = u.x * l.y - u.y * l.x;
			l.x = l.x * cos(3 * pi / 180) + cross.x * sin(3 * pi / 180);
			l.y = l.y * cos(3 * pi / 180) + cross.y * sin(3 * pi / 180);
			l.z = l.z * cos(3 * pi / 180) + cross.z * sin(3 * pi / 180);
			l = normalize(l);
			break;

		case '2': //look right
			cross.x = r.y * u.z - r.z * u.y;
			cross.y = r.z * u.x - r.x * u.z;
			cross.z = r.x * u.y - r.y * u.x;
			r.x = r.x * cos(3 * pi / 180) + cross.x * sin(3 * pi / 180);
			r.y = r.y * cos(3 * pi / 180) + cross.y * sin(3 * pi / 180);
			r.z = r.z * cos(3 * pi / 180) + cross.z * sin(3 * pi / 180);
            r = normalize(r);

			cross.x = l.y * u.z - l.z * u.y;
			cross.y = l.z * u.x - l.x * u.z;
			cross.z = l.x * u.y - l.y * u.x;
			l.x = l.x * cos(3 * pi / 180) + cross.x * sin(3 * pi / 180);
			l.y = l.y * cos(3 * pi / 180) + cross.y * sin(3 * pi / 180);
			l.z = l.z * cos(3 * pi / 180) + cross.z * sin(3 * pi / 180);
			l = normalize(l);
			break;

        case '3': //look up
			cross.x = r.y * l.z - r.z * l.y;
			cross.y = r.z * l.x - r.x * l.z;
			cross.z = r.x * l.y - r.y * l.x;
			l.x = l.x * cos(3 * pi / 180) + cross.x * sin(3 * pi / 180);
			l.y = l.y * cos(3 * pi / 180) + cross.y * sin(3 * pi / 180);
			l.z = l.z * cos(3 * pi / 180) + cross.z * sin(3 * pi / 180);
            l = normalize(l);

			cross.x = r.y * u.z - r.z * u.y;
			cross.y = r.z * u.x - r.x * u.z;
			cross.z = r.x * u.y - r.y * u.x;
			u.x = u.x * cos(3 * pi / 180) + cross.x * sin(3 * pi / 180);
			u.y = u.y * cos(3 * pi / 180) + cross.y * sin(3 * pi / 180);
			u.z = u.z * cos(3 * pi / 180) + cross.z * sin(3 * pi / 180);
			u = normalize(u);
			break;

        case '4': //look down
            cross.x = l.y * r.z - l.z * r.y;
			cross.y = l.z * r.x - l.x * r.z;
			cross.z = l.x * r.y - l.y * r.x;
			l.x = l.x * cos(3 * pi / 180) + cross.x * sin(3 * pi / 180);
			l.y = l.y * cos(3 * pi / 180) + cross.y * sin(3 * pi / 180);
			l.z = l.z * cos(3 * pi / 180) + cross.z * sin(3 * pi / 180);
            l = normalize(l);

            cross.x = u.y * r.z - u.z * r.y;
            cross.y = u.z * r.x - u.x * r.z;
            cross.z = u.x * r.y - u.y * r.x;
            u.x = u.x * cos(3 * pi / 180) + cross.x * sin(3 * pi / 180);
            u.y = u.y * cos(3 * pi / 180) + cross.y * sin(3 * pi / 180);
            u.z = u.z * cos(3 * pi / 180) + cross.z * sin(3 * pi / 180);
            u = normalize(u);
            break;

        case '5': //tilt clockwise
            cross.x = l.y * r.z - l.z * r.y;
            cross.y = l.z * r.x - l.x * r.z;
            cross.z = l.x * r.y - l.y * r.x;
            r.x = r.x * cos(3 * pi / 180) + cross.x * sin(3 * pi / 180);
            r.y = r.y * cos(3 * pi / 180) + cross.y * sin(3 * pi / 180);
            r.z = r.z * cos(3 * pi / 180) + cross.z * sin(3 * pi / 180);
            r = normalize(r);

            cross.x = l.y * u.z - l.z * u.y;
            cross.y = l.z * u.x - l.x * u.z;
            cross.z = l.x * u.y - l.y * u.x;
            u.x = u.x * cos(3 * pi / 180) + cross.x * sin(3 * pi / 180);
            u.y = u.y * cos(3 * pi / 180) + cross.y * sin(3 * pi / 180);
            u.z = u.z * cos(3 * pi / 180) + cross.z * sin(3 * pi / 180);
            u = normalize(u);
            break;
        case '6': //tilt counterclockwise
            cross.x = r.y * l.z - r.z * l.y;
            cross.y = r.z * l.x - r.x * l.z;
            cross.z = r.x * l.y - r.y * l.x;
            r.x = r.x * cos(3 * pi / 180) + cross.x * sin(3 * pi / 180);
            r.y = r.y * cos(3 * pi / 180) + cross.y * sin(3 * pi / 180);
            r.z = r.z * cos(3 * pi / 180) + cross.z * sin(3 * pi / 180);
            r = normalize(r);

            cross.x = u.y * l.z - u.z * l.y;
            cross.y = u.z * l.x - u.x * l.z;
            cross.z = u.x * l.y - u.y * l.x;
            u.x = u.x * cos(3 * pi / 180) + cross.x * sin(3 * pi / 180);
            u.y = u.y * cos(3 * pi / 180) + cross.y * sin(3 * pi / 180);
            u.z = u.z * cos(3 * pi / 180) + cross.z * sin(3 * pi / 180);
            u = normalize(u);
            break;
		default:
			break;
	}
}


void specialKeyListener(int key, int x,int y){
	switch(key){
		case GLUT_KEY_DOWN:		//down arrow key
			pos.x -= 3*l.x;
			pos.y -= 3*l.y;
			pos.z -= 3*l.z;
			break;
		case GLUT_KEY_UP:		// up arrow key
		    pos.x += l.x;
		    pos.y += l.y;
		    pos.z += l.z;
			break;

		case GLUT_KEY_RIGHT:
			pos.x += r.x;
			pos.y += r.y;
			pos.z += r.z;
			break;
		case GLUT_KEY_LEFT:
			pos.x -= r.x;
			pos.y -= r.y;
			pos.z -= r.z;
			break;

		case GLUT_KEY_PAGE_UP:
		    pos.x += u.x;
		    pos.y += u.y;
		    pos.z += u.z;
			break;
		case GLUT_KEY_PAGE_DOWN:
		    pos.x -= u.x;
		    pos.y -= u.y;
		    pos.z -= u.z;
			break;

		case GLUT_KEY_INSERT:
			break;

		case GLUT_KEY_HOME:
		    if(sphereRadius < 36.0)  sphereRadius += 2.0;
			break;
		case GLUT_KEY_END:
		    if(sphereRadius > 0.0)   sphereRadius -= 2.0;
			break;

		default:
			break;

	}
}


void mouseListener(int button, int state, int x, int y){	//x, y is the x-y of the screen (2D)
	switch(button){
		case GLUT_LEFT_BUTTON:
			if(state == GLUT_DOWN){		// 2 times?? in ONE click? -- solution is checking DOWN or UP
				//drawaxes=1-drawaxes;
			}
			break;

		case GLUT_RIGHT_BUTTON:
			//........
			break;

		case GLUT_MIDDLE_BUTTON:
			//........
			break;

		default:
			break;
	}
}



void display(){

	//clear the display
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glClearColor(0,0,0,0);	//color black
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	/********************
	/ set-up camera here
	********************/
	//load the correct matrix -- MODEL-VIEW matrix
	glMatrixMode(GL_MODELVIEW);

	//initialize the matrix
	glLoadIdentity();

	//now give three info
	//1. where is the camera (viewer)?
	//2. where is the camera looking?
	//3. Which direction is the camera's UP direction?

	//gluLookAt(100,100,100,	0,0,0,	0,0,1);
	//gluLookAt(200*cos(cameraAngle), 200*sin(cameraAngle), cameraHeight,		0,0,0,		0,0,1);

	//gluLookAt(0,0,200,	0,0,0,	0,1,0);

    gluLookAt(pos.x, pos.y, pos.z,		pos.x + l.x, pos.y + l.y, pos.z + l.z,		u.x, u.y, u.z);
    //gluLookAt(0, -200, 10,		0, 0, 0,		0, 0, 1);

	//again select MODEL-VIEW
	glMatrixMode(GL_MODELVIEW);


	/****************************
	/ Add your objects from here
	****************************/
	//add objects

	drawAxes();

	for(int i = 0; i < objects.size(); i++){
        objects[i]->draw();
	}

	for(int i = 0; i < lights.size(); i++){
        glBegin(GL_POINTS);
        glColor3f(1, 1, 1);
        glVertex3f(lights[i].x, lights[i].y, lights[i].z);
        //glVertex3f(lights[i].x+1, lights[i].y+1, lights[i].z+1);
        //glVertex3f(lights[i].x-1, lights[i].y+1, lights[i].z+1);
        //glVertex3f(lights[i].x+1, lights[i].y-1, lights[i].z+1);
        glEnd();
	}
	//drawGrid();

    //glColor3f(1,0,0);
    //drawSquare(10);

    //drawSS();
    //drawcube2box();
    //drawCircle(30,24);

    //drawCone(20,50,24);

	//drawSphere(30,24,20);

    //drawLine();


	//ADD this line in the end --- if you use double buffer (i.e. GL_DOUBLE)
	glutSwapBuffers();
}


void animate(){
	angle+=0.05;
	//codes for any changes in Models, Camera
	glutPostRedisplay();
}

void init(){
	//codes for initialization
	drawgrid=0;
	drawaxes=1;
	cameraHeight=150.0;
	sphereRadius = 0.0;
	cameraAngle=1.0;
	angle=0;

	//clear the screen
	glClearColor(0,0,0,0);

	/************************
	/ set-up projection here
	************************/
	//load the PROJECTION matrix
	glMatrixMode(GL_PROJECTION);

	//initialize the matrix
	glLoadIdentity();

	//give PERSPECTIVE parameters
	gluPerspective(FOV,	1,	1,	1000.0);
	//field of view in the Y (vertically)
	//aspect ratio that determines the field of view in the X direction (horizontally)
	//near distance
	//far distance
	pos.x = 0;
	pos.y = -200;
	pos.z = 70;
	u.x = 0;
	u.y = 0;
	u.z = 1;
	//r.x = -1/sqrt(2);
	//r.y = 1/sqrt(2);
	//r.z = 0;
	//l.x = -1/sqrt(2);
	//l.y = -1/sqrt(2);
	//l.z = 0;
	l.x = 0;
	l.y = 1;
	l.z = 0;
	r.x = 1;
	r.y = 0;
	r.z = 0;
}

void loadActualData()
{
    ifstream fin("scene.txt");
    fin>>recursion_level;
    //cout<<recursion_level<<endl;

    fin>>image_width;
    image_height = image_width;
    //cout<<image_width <<" "<<image_height<<endl;

    fin>>numObject;
    //cout<<numObject<<endl;

    Object *temp;
    for(int i = 0; i < numObject;){
        string str;
        fin>>str;
        //cout<<str<<endl;
        if(str.compare("sphere") == 0 || str.compare("Sphere") == 0){
            struct point center;
            double radius;
            fin>>center.x>>center.y>>center.z;
            //cout<<center.x<<" "<<center.y<<" "<<center.z<<endl;
            fin>>radius;
            //cout<<radius<<endl;
            temp = new Sphere(center, radius);
            double r, g, b, a, d, s, re;
            fin>>r>>g>>b;
            //cout<<r<<" "<<g<<" "<<b<<endl;
            temp->setColor(r, g, b);
            fin>>a>>d>>s>>re;
            //cout<<a<<" "<<d<<" "<<s<<" "<<re<<endl;
            temp->setCoEfficients(a, d, s, re);
            int shine;
            fin>>shine;
            //cout<<shine<<endl;
            temp->setShine(shine);
            objects.push_back(temp);
            i++;
        }
        else if(str.compare("triangle") == 0 || str.compare("Triangle") == 0){
            struct point aa, bb, c;
            fin>>aa.x>>aa.y>>aa.z;

            fin>>bb.x>>bb.y>>bb.z;

            fin>>c.x>>c.y>>c.z;

            //cout<<aa.x<<" "<<aa.y<<" "<<aa.z<<endl;
            //cout<<bb.x<<" "<<bb.y<<" "<<bb.z<<endl;
            //cout<<c.x<<" "<<c.y<<" "<<c.z<<endl;
            temp = new Triangle(aa, bb, c);
            double r, g, b, a, d, s, re;
            fin>>r>>g>>b;
            //cout<<r<<" "<<g<<" "<<b<<endl;
            temp->setColor(r, g, b);
            fin>>a>>d>>s>>re;
            //cout<<a<<" "<<d<<" "<<s<<" "<<re<<endl;
            temp->setCoEfficients(a, d, s, re);
            int shine;
            fin>>shine;
            //cout<<shine<<endl;
            temp->setShine(shine);
            objects.push_back(temp);
            i++;
        }
        else if(str.compare("general") == 0 || str.compare("General") == 0){
            double A, B, C, D, E, F, G, H, I, J, lens, wids, heis;
            fin>>A>>B>>C>>D>>E>>F>>G>>H>>I>>J;
            //cout<<A<<" "<<B<<" "<<C<<" "<<D<<" "<<E<<" "<<F<<" "<<G<<" "<<H<<" "<<I<<" "<<J<<endl;
            struct point ref_point;
            fin>>ref_point.x>>ref_point.y>>ref_point.z>>lens>>wids>>heis;
            //cout<<ref_point.x<<" "<<ref_point.y<<" "<<ref_point.z<<" "<<lens<<" "<<wids<<" "<<heis<<endl;;

            temp = new generalQuadratic(A, B, C, D, E, F, G, H, I, J, ref_point, lens, wids, heis);
            double r, g, b, a, d, s, re;
            fin>>r>>g>>b;
            //cout<<r<<" "<<g<<" "<<b<<endl;
            temp->setColor(r, g, b);
            fin>>a>>d>>s>>re;
            //cout<<a<<" "<<d<<" "<<s<<" "<<re<<endl;
            temp->setCoEfficients(a, d, s, re);
            int shine;
            fin>>shine;
            //cout<<shine<<endl;
            temp->setShine(shine);
            objects.push_back(temp);
            i++;
        }
    }
    fin>>numLight;
    //cout<<numLight<<endl;
    for(int i = 0; i < numLight; i++){
        struct point light;
        fin>>light.x>>light.y>>light.z;
        //cout<<light.x<<" "<<light.y<<" "<<light.z<<endl;
        lights.push_back(light);
    }

    temp = new Floor(1000, 20);
    temp->setCoEfficients(0.4, 0.2, 0.2, 0.2);
    temp->setShine(1);
    objects.push_back(temp);
}

void loadTestData()
{
    image_width = image_height = 768;
    recursion_level = 4;

    Object *temp;
    struct point C;
    double R;
    C.x = 10;
    C.y = 50;
    C.z = 10;
    R = 10.0;
    temp = new Sphere(C, R);
    temp->setColor(1, 0, 0);
    temp->setCoEfficients(0.4, 0.2, 0.2, 0.2);
    temp->setShine(1);

    objects.push_back(temp);

    struct point light1;
    light1.x = -30;
    light1.y = 20;
    light1.z = 20;
    lights.push_back(light1);

    /*light1.x = -60;
    light1.y = -20;
    light1.z = 10;
    lights.push_back(light1);*/

    C.x = -60;
    C.y = 15;
    C.z = 25;
    R = 25.0;
    temp = new Sphere(C, R);
    temp->setColor(1, 1, 0);
    temp->setCoEfficients(0.4, 0.2, 0.2, 0.2);
    temp->setShine(1);

    objects.push_back(temp);

    struct point a, b, c;
    a.x = -40, a.y = -40, a.z = 0;
    b.x = 40, b.y = -40, b.z = 0;
    c.x = 0, c.y = -40, c.z = 40;

    temp = new Triangle(a, b, c);
    temp->setColor(1, 0, 0);
    temp->setCoEfficients(0.4, 0.2, 0.2, 0.2);
    temp->setShine(1);
    //objects.push_back(temp);

    temp = new Floor(1000, 20);
    temp->setCoEfficients(0.4, 0.2, 0.2, 0.2);
    temp->setShine(1);
    objects.push_back(temp);

    struct point r;
    r.x = 0;
    r.y = 0;
    r.z = 0;

    /*temp = new generalQuadratic(1, 1, 1, 0, 0, 0, 0, 0, 0, -100, r, 0, 0, 20);
    temp->setColor(0, 1, 0);
    temp->setCoEfficients(0.4, 0.2, 0.1, 0.3);
    temp->setShine(10);
    objects.push_back(temp);*/

    temp = new generalQuadratic(0.0625, 0.04, 0.04, 0, 0, 0, 0, 0, 0, -36, r, 0, 0, 15);
    temp->setColor(1, 0, 0);
    temp->setCoEfficients(0.4, 0.2, 0.1, 0.3);
    temp->setShine(15);
    objects.push_back(temp);

}

int main(int argc, char **argv){
	//loadTestData();
	loadActualData();
	glutInit(&argc,argv);
	glutInitWindowSize(window_width, window_height);
	glutInitWindowPosition(0, 0);
	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB);	//Depth, Double buffer, RGB color

	glutCreateWindow("My OpenGL Program");

	init();

	glEnable(GL_DEPTH_TEST);	//enable Depth Testing

	glutDisplayFunc(display);	//display callback function
	glutIdleFunc(animate);		//what you want to do in the idle time (when no drawing is occuring)

	glutKeyboardFunc(keyboardListener);
	glutSpecialFunc(specialKeyListener);
	glutMouseFunc(mouseListener);

	glutMainLoop();		//The main loop of OpenGL

	return 0;
}
