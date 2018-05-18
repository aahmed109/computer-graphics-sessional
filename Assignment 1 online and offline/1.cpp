#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#include <windows.h>
#include <glut.h>

#define pi (2*acos(0.0))

using namespace std;

double cameraHeight;
double cameraAngle;
double cameraRight;
double sphereRadius;
int drawgrid;
int drawaxes;
double angle;

struct point
{
	double x,y,z;
};

struct point pos, u, r, l;

void drawAxes()
{
	if(drawaxes==1)
	{
		glColor3f(1.0, 1.0, 1.0);
		glBegin(GL_LINES);{
			glVertex3f( 100,0,0);
			glVertex3f(-100,0,0);

			glVertex3f(0,-100,0);
			glVertex3f(0, 100,0);

			glVertex3f(0,0, 100);
			glVertex3f(0,0,-100);
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
        h = height * sin(((double)i / (double)stacks) * (pi/2));;
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
		h=radius*sin(((double)i/(double)stacks)*(pi/2));
		r=radius*cos(((double)i/(double)stacks)*(pi/2));
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
		h=radius*sin(((double)i/(double)stacks)*(pi/2));
		r=radius*cos(((double)i/(double)stacks)*(pi/2));
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

void keyboardListener(unsigned char key, int x,int y){
	switch(key){
        struct point cross;
		case '1': //look left
			cross.x = u.y * r.z - u.z * r.y;
			cross.y = u.z * r.x - u.x * r.z;
			cross.z = u.x * r.y - u.y * r.x;
			r.x = r.x * cos(3 * pi / 180) + cross.x * sin(3 * pi / 180);
			r.y = r.y * cos(3 * pi / 180) + cross.y * sin(3 * pi / 180);
			r.z = r.z * cos(3 * pi / 180) + cross.z * sin(3 * pi / 180);

			cross.x = u.y * l.z - u.z * l.y;
			cross.y = u.z * l.x - u.x * l.z;
			cross.z = u.x * l.y - u.y * l.x;
			l.x = l.x * cos(3 * pi / 180) + cross.x * sin(3 * pi / 180);
			l.y = l.y * cos(3 * pi / 180) + cross.y * sin(3 * pi / 180);
			l.z = l.z * cos(3 * pi / 180) + cross.z * sin(3 * pi / 180);
			break;

		case '2': //look right
			cross.x = r.y * u.z - r.z * u.y;
			cross.y = r.z * u.x - r.x * u.z;
			cross.z = r.x * u.y - r.y * u.x;
			r.x = r.x * cos(3 * pi / 180) + cross.x * sin(3 * pi / 180);
			r.y = r.y * cos(3 * pi / 180) + cross.y * sin(3 * pi / 180);
			r.z = r.z * cos(3 * pi / 180) + cross.z * sin(3 * pi / 180);

			cross.x = l.y * u.z - l.z * u.y;
			cross.y = l.z * u.x - l.x * u.z;
			cross.z = l.x * u.y - l.y * u.x;
			l.x = l.x * cos(3 * pi / 180) + cross.x * sin(3 * pi / 180);
			l.y = l.y * cos(3 * pi / 180) + cross.y * sin(3 * pi / 180);
			l.z = l.z * cos(3 * pi / 180) + cross.z * sin(3 * pi / 180);
			break;

        case '3': //look up
			cross.x = r.y * l.z - r.z * l.y;
			cross.y = r.z * l.x - r.x * l.z;
			cross.z = r.x * l.y - r.y * l.x;
			l.x = l.x * cos(3 * pi / 180) + cross.x * sin(3 * pi / 180);
			l.y = l.y * cos(3 * pi / 180) + cross.y * sin(3 * pi / 180);
			l.z = l.z * cos(3 * pi / 180) + cross.z * sin(3 * pi / 180);

			cross.x = r.y * u.z - r.z * u.y;
			cross.y = r.z * u.x - r.x * u.z;
			cross.z = r.x * u.y - r.y * u.x;
			u.x = u.x * cos(3 * pi / 180) + cross.x * sin(3 * pi / 180);
			u.y = u.y * cos(3 * pi / 180) + cross.y * sin(3 * pi / 180);
			u.z = u.z * cos(3 * pi / 180) + cross.z * sin(3 * pi / 180);
			break;

        case '4': //look down
            cross.x = l.y * r.z - l.z * r.y;
			cross.y = l.z * r.x - l.x * r.z;
			cross.z = l.x * r.y - l.y * r.x;
			l.x = l.x * cos(3 * pi / 180) + cross.x * sin(3 * pi / 180);
			l.y = l.y * cos(3 * pi / 180) + cross.y * sin(3 * pi / 180);
			l.z = l.z * cos(3 * pi / 180) + cross.z * sin(3 * pi / 180);

            cross.x = u.y * r.z - u.z * r.y;
            cross.y = u.z * r.x - u.x * r.z;
            cross.z = u.x * r.y - u.y * r.x;
            u.x = u.x * cos(3 * pi / 180) + cross.x * sin(3 * pi / 180);
            u.y = u.y * cos(3 * pi / 180) + cross.y * sin(3 * pi / 180);
            u.z = u.z * cos(3 * pi / 180) + cross.z * sin(3 * pi / 180);
            break;

        case '5': //tilt clockwise
            cross.x = l.y * r.z - l.z * r.y;
            cross.y = l.z * r.x - l.x * r.z;
            cross.z = l.x * r.y - l.y * r.x;
            r.x = r.x * cos(3 * pi / 180) + cross.x * sin(3 * pi / 180);
            r.y = r.y * cos(3 * pi / 180) + cross.y * sin(3 * pi / 180);
            r.z = r.z * cos(3 * pi / 180) + cross.z * sin(3 * pi / 180);

            cross.x = l.y * u.z - l.z * u.y;
            cross.y = l.z * u.x - l.x * u.z;
            cross.z = l.x * u.y - l.y * u.x;
            u.x = u.x * cos(3 * pi / 180) + cross.x * sin(3 * pi / 180);
            u.y = u.y * cos(3 * pi / 180) + cross.y * sin(3 * pi / 180);
            u.z = u.z * cos(3 * pi / 180) + cross.z * sin(3 * pi / 180);
            break;
        case '6': //tilt counterclockwise
            cross.x = r.y * l.z - r.z * l.y;
            cross.y = r.z * l.x - r.x * l.z;
            cross.z = r.x * l.y - r.y * l.x;
            r.x = r.x * cos(3 * pi / 180) + cross.x * sin(3 * pi / 180);
            r.y = r.y * cos(3 * pi / 180) + cross.y * sin(3 * pi / 180);
            r.z = r.z * cos(3 * pi / 180) + cross.z * sin(3 * pi / 180);

            cross.x = u.y * l.z - u.z * l.y;
            cross.y = u.z * l.x - u.x * l.z;
            cross.z = u.x * l.y - u.y * l.x;
            u.x = u.x * cos(3 * pi / 180) + cross.x * sin(3 * pi / 180);
            u.y = u.y * cos(3 * pi / 180) + cross.y * sin(3 * pi / 180);
            u.z = u.z * cos(3 * pi / 180) + cross.z * sin(3 * pi / 180);
            break;
		default:
			break;
	}
}


void specialKeyListener(int key, int x,int y){
	switch(key){
		case GLUT_KEY_DOWN:		//down arrow key
			pos.x -= l.x;
			pos.y -= l.y;
			break;
		case GLUT_KEY_UP:		// up arrow key
		    pos.x += l.x;
		    pos.y += l.y;
			break;

		case GLUT_KEY_RIGHT:
			pos.x += r.x;
			pos.y += r.y;
			break;
		case GLUT_KEY_LEFT:
			pos.x -= r.x;
			pos.y -= r.y;
			break;

		case GLUT_KEY_PAGE_UP:
		    pos.z += u.z;
			break;
		case GLUT_KEY_PAGE_DOWN:
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
	gluLookAt(pos.x, pos.y, pos.z,		pos.x + l.x, pos.y + l.y, pos.z + l.z,		u.x, u.y, u.z);
	//gluLookAt(0,0,200,	0,0,0,	0,1,0);


	//again select MODEL-VIEW
	glMatrixMode(GL_MODELVIEW);


	/****************************
	/ Add your objects from here
	****************************/
	//add objects

	drawAxes();
	drawGrid();

    //glColor3f(1,0,0);
    //drawSquare(10);

    //drawSS();
    drawcube2box();
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
	gluPerspective(80,	1,	1,	1000.0);
	//field of view in the Y (vertically)
	//aspect ratio that determines the field of view in the X direction (horizontally)
	//near distance
	//far distance
	pos.x = 150;
	pos.y = 150;
	pos.z = 36;
	u.x = 0;
	u.y = 0;
	u.z = 1.5;
	r.x = -sqrt(1.5);
	r.y = sqrt(1.5);
	r.z = 0;
	l.x = -sqrt(1.5);
	l.y = -sqrt(1.5);
	l.z = 0;
}

int main(int argc, char **argv){
	glutInit(&argc,argv);
	glutInitWindowSize(500, 500);
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
