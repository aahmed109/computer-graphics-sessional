#include<stdio.h>
#include<algorithm>
#include<vector>
#include<glut.h>
#include "bitmap_image.hpp"
#ifndef _classes_h
#define _classes_h

#define pi (2*acos(0.0))
#define AMBIENT 0
#define DIFFUSE 1
#define SPECULAR 2
#define REFLECTION 3

#define REFINDEX 1.7

#define EPSILON 0.0000001

using namespace std;
///important and internal classes will be written here

struct point
{
	double x,y,z;
};

extern int recursion_level;
extern vector<point> lights;

class Ray
{
public:
    struct point start;
    struct point dir;

    Ray(struct point s, struct point d){
        start = s;
        dir = d;
    }

    Ray(){ ///in case more constructor needed

    }
};

struct point normalize(struct point vec)
{
    double mod = sqrt(pow(vec.x,2) + pow(vec.y,2) + pow(vec.z,2));
    vec.x /= mod;
    vec.y /= mod;
    vec.z /= mod;

    return vec;
}

class Object
{
public:
    double dotProduct(struct point a, struct point b)
    {
        return a.x * b.x + a.y * b.y + a.z * b.z;
    }

    struct point crossProduct(struct point a, struct point b)
    {
        struct point cross;
        cross.x = a.y * b.z - a.z * b.y;
        cross.y = a.z * b.x - a.x * b.z;
        cross.z = a.x * b.y - a.y * b.x;

        return cross;
    }

    struct point reference_point;
	double height, width, length, source_factor = 1.0;
	int Shine;
	double color[3];
	double co_efficients[4];    ///what's the purrrpose of this?

	Object(){ }
	virtual struct point getNormal(struct point pointed){
        return normalize(pointed);
	}

	virtual void draw(){}
	virtual double getT(Ray* r){
        return -1;
	}
	virtual double intersect(Ray *r, struct point* current_color, int level){
        return -1;
	}

    struct point getReflection(Ray *r, struct point pointed)
    {
        struct point reflection;
        reflection.x = r->dir.x - pointed.x * 2.0 * dotProduct(r->dir, pointed);
        reflection.y = r->dir.y - pointed.y * 2.0 * dotProduct(r->dir, pointed);
        reflection.z = r->dir.z - pointed.z * 2.0 * dotProduct(r->dir, pointed);

        return normalize(reflection);
    }

    struct point getRefraction(Ray *r, struct point pointed)
    {
        struct point refraction;
        double product = dotProduct(pointed, r->dir);
        double k = 1 - pow(REFINDEX,2) * (1 - pow(product,2));
        if(k < 0){
            refraction.x = 0;
            refraction.y = 0;
            refraction.z = 0;
        }
        else{
            refraction.x = REFINDEX * r->dir.x - (REFINDEX * product + sqrt(k)) * pointed.x;
            refraction.y = REFINDEX * r->dir.y - (REFINDEX * product + sqrt(k)) * pointed.y;
            refraction.z = REFINDEX * r->dir.z - (REFINDEX * product + sqrt(k)) * pointed.z;
        }

        return refraction;
    }

    void setColor(double r, double g, double b) {this->color[0] = r; this->color[1] = g; this->color[2] = b;};
    void setShine(int shine)    {this->Shine = shine;};
    void setCoEfficients(double a, double b, double c, double d)    {this->co_efficients[0] = a; this->co_efficients[1] = b; this->co_efficients[2] = c; this->co_efficients[3] = d;}

};
extern vector<Object*> objects;

class Sphere: public Object{
public:
    Sphere(struct point Center, double Radius){
        reference_point = Center;
        length = Radius;
        height = Radius;
    }

    void draw(){
        glColor3f(color[0], color[1], color[2]);

        struct point points[100][100];
        int i,j;
        double h,r;
        int stacks = 20;
        int slices = 24;
        //generate points
        for(i=0;i<=stacks;i++)
        {
            h=length*sin(((double)i/(double)stacks)*(pi/2.0));
            r=length*cos(((double)i/(double)stacks)*(pi/2.0));
            for(j=0;j<=slices;j++)
            {
                points[i][j].x=r*cos(((double)j/(double)slices)*2.0*pi);
                points[i][j].y=r*sin(((double)j/(double)slices)*2.0*pi);
                points[i][j].z=h;
            }
        }
        //draw quads using generated points
        for(i=0;i<stacks;i++)
        {
            for(j=0;j<slices;j++)
            {
                glBegin(GL_QUADS);{
                    //upper hemisphere
                    glVertex3f(points[i][j].x + reference_point.x, points[i][j].y + reference_point.y, points[i][j].z + reference_point.z);
                    glVertex3f(points[i][j+1].x + reference_point.x, points[i][j+1].y + reference_point.y, points[i][j+1].z + reference_point.z);
                    glVertex3f(points[i+1][j+1].x + reference_point.x, points[i+1][j+1].y + reference_point.y, points[i+1][j+1].z + reference_point.z);
                    glVertex3f(points[i+1][j].x + reference_point.x, points[i+1][j].y + reference_point.y, points[i+1][j].z + reference_point.z);
                    //lower hemisphere
                    glVertex3f(points[i][j].x + reference_point.x,points[i][j].y + reference_point.y,-points[i][j].z + reference_point.z);
                    glVertex3f(points[i][j+1].x + reference_point.x,points[i][j+1].y + reference_point.y,-points[i][j+1].z + reference_point.z);
                    glVertex3f(points[i+1][j+1].x + reference_point.x,points[i+1][j+1].y + reference_point.y,-points[i+1][j+1].z + reference_point.z);
                    glVertex3f(points[i+1][j].x + reference_point.x,points[i+1][j].y + reference_point.y,-points[i+1][j].z + reference_point.z);
                }glEnd();
            }
        }
    }

    double getT(Ray* r)
    {
        struct point s;
        s.x = r->start.x - reference_point.x;
        s.y = r->start.y - reference_point.y;
        s.z = r->start.z - reference_point.z;

        double A = dotProduct(r->dir, r->dir);
        double B = 2.0 * dotProduct(r->dir, s);
        double C = dotProduct(s, s) - pow(length,2);
        double d = pow(B,2) - 4 * A * C;

        if(d < 0) return -1;

        double t1 = (-B + sqrt(d)) / (2 * A);
        double t2 = (-B - sqrt(d)) / (2 * A);

        return t1<t2?t1:t2;
    }

    struct point getNormal(struct point pointed)
    {
        struct point normal;
        normal.x = pointed.x - reference_point.x;
        normal.y = pointed.y - reference_point.y;
        normal.z = pointed.z - reference_point.z;

        return normalize(normal);
    }



    double intersect(Ray *r, struct point* current_color, int level)
    {
        //printf("%f %f %f ", r->dir.x, r->dir.y, r->dir.z);

        double t = getT(r);

        //printf("here0\n");
        if(t <= 0) return -1;

        //printf("here1\n");
        if(level == 0) return t;

        //printf("here2\n");
        struct point intersection;
        intersection.x = r->start.x + r->dir.x * t;
        intersection.y = r->start.y + r->dir.y * t;
        intersection.z = r->start.z + r->dir.z * t;

        current_color->x = color[0] * co_efficients[AMBIENT];
        current_color->y = color[1] * co_efficients[AMBIENT];
        current_color->z = color[2] * co_efficients[AMBIENT];

        struct point normal = getNormal(intersection);
        struct point reflection = getReflection(r, normal);
        struct point refraction = getRefraction(r, normal);

        for(int i = 0; i < lights.size(); i++){
            struct point ldir;
            ldir.x = lights[i].x - intersection.x;
            ldir.y = lights[i].y - intersection.y;
            ldir.z = lights[i].z - intersection.z;

            double len = sqrt(pow(ldir.x,2) + pow(ldir.y,2) + pow(ldir.z,2));
            ldir = normalize(ldir);

            struct point lstart;
            lstart.x = intersection.x + ldir.x * 1.0;
            lstart.y = intersection.y + ldir.y * 1.0;
            lstart.z = intersection.z + ldir.z * 1.0;


            Ray L(lstart, ldir);

            bool obscured = false;

            for(int j = 0; j < objects.size(); j++){
                double possibleobscure = objects[j]->getT(&L);
                //printf(",%f %d, ", possibleobscure, j);
                if(possibleobscure > 0 && abs(possibleobscure) < len){
                    obscured = true;
                    break;
                }
            }
            //printf("\n");

            if(!obscured){
                //printf("direct");
                //L.dir.x = -L.dir.x;
                //L.dir.y = -L.dir.y;
                //L.dir.z = -L.dir.z;
                double lambert = dotProduct(L.dir, normal);
                double phong = dotProduct(reflection, r->dir)/2;

                lambert = lambert > 0?lambert:0;
                phong = phong > 0?phong:0;

                current_color->x += source_factor * color[0] * (lambert * co_efficients[DIFFUSE] + pow(phong, Shine) * co_efficients[SPECULAR]);
                current_color->y += source_factor * color[1] * (lambert * co_efficients[DIFFUSE] + pow(phong, Shine) * co_efficients[SPECULAR]);
                current_color->z += source_factor * color[2] * (lambert * co_efficients[DIFFUSE] + pow(phong, Shine) * co_efficients[SPECULAR]);
            }

            if(level < recursion_level){
                struct point rs;
                rs.x = intersection.x + reflection.x * 1.0;
                rs.y = intersection.y + reflection.y * 1.0;
                rs.z = intersection.z + reflection.z * 1.0;

                //reflection = normalize(reflection);

                Ray reflectionRay(rs, reflection);
                int nearest = -1;

                struct point reflectColor;
                double t_min = 9999;
                for(int k = 0; k < objects.size(); k++){
                    double t = objects[k]->getT(&reflectionRay);

                    if(t <= 0) continue;
                    else if(t < t_min) {
                        t_min = t;
                        nearest = k;
                    }
                }

                if(nearest != -1){
                    double t = objects[nearest]->intersect(&reflectionRay, &reflectColor, level + 1);
                    if(t != -1){
                        current_color->x += reflectColor.x*co_efficients[REFLECTION];
                        current_color->y += reflectColor.y*co_efficients[REFLECTION];
                        current_color->z += reflectColor.z*co_efficients[REFLECTION];
                    }

                }

                rs.x = intersection.x + refraction.x * 1.0;
                rs.y = intersection.y + refraction.y * 1.0;
                rs.z = intersection.z + refraction.z * 1.0;

                Ray refractionRay(rs, refraction);

                nearest = -1;

                struct point refractColor;
                t_min = 9999;
                for(int k = 0; k < objects.size(); k++){
                    double t = objects[k]->getT(&refractionRay);

                    if(t <= 0) continue;
                    else if(t < t_min) {
                        t_min = t;
                        nearest = k;
                    }
                }

                if(nearest != -1){
                    double t = objects[nearest]->intersect(&refractionRay, &refractColor, level + 1);
                    if(t != -1){
                        current_color->x += refractColor.x*REFINDEX;
                        current_color->y += refractColor.y*REFINDEX;
                        current_color->z += refractColor.z*REFINDEX;
                    }

                }
            }

            //Check whether all current_color pixel value is within 1 or 0 if not set it
        if(current_color->x > 1)    current_color->x = 1;
        if(current_color->x < 0)    current_color->x = 0;
        if(current_color->y > 1)    current_color->y = 1;
        if(current_color->y < 0)    current_color->y = 0;
        if(current_color->z > 1)    current_color->z = 1;
        if(current_color->z < 0)    current_color->z = 0;
        }

        return t;
    }
};

class Floor: public Object
{
public:
    bitmap_image texture;
    double texture_height, texture_width;
    Floor(double floor_width, double tile_width){
        reference_point.x = -floor_width/2.0;
        reference_point.y = -floor_width/2.0;
        reference_point.z = 0;
        length = tile_width;
        height = 0;
        texture = bitmap_image("texture.bmp");
        texture_height = texture.height()/1000.0;
        texture_width = texture.width()/1000.0;
    }

    struct point getNormal(struct point pointed){
        struct point normal;
        normal.x = 0;
        normal.y = 0;
        normal.z = 1;
        return normal;
    }

    double getT(Ray *r){
        struct point normal = getNormal(reference_point);
        return ((-1.0) * dotProduct(normal, r->start) / dotProduct(normal, r->dir));
    }

    double intersect(Ray *r, struct point* current_color, int level){
        double t = getT(r);

        struct point intersection;
        intersection.x = r->start.x + r->dir.x * t;
        intersection.y = r->start.y + r->dir.y * t;
        intersection.z = r->start.z + r->dir.z * t;

        if(intersection.x < reference_point.x || intersection.x > -reference_point.x || intersection.y < reference_point.y || intersection.y > -reference_point.y){
            return -1;
        }
        if(level == 0) return t;

        int x = (intersection.x - reference_point.x) / length;
        int y = (intersection.y - reference_point.y) / length;

        if((x + y) % 2 == 0){
            color[0] = color[1] = color[2] = 0;
        }
        else {
            color[0] = color[1] = color[2] = 1;
        }

        unsigned char red, green, blue;
        int xx = (intersection.x + abs(reference_point.x)) * texture_width;
        int yy = (intersection.y + abs(reference_point.y)) * texture_height;

        texture.get_pixel(xx, yy, red, green, blue);


        current_color->x = color[0] * co_efficients[AMBIENT] * red / 255.0;
        current_color->y = color[1] * co_efficients[AMBIENT] * green / 255.0;
        current_color->z = color[2] * co_efficients[AMBIENT] * blue /255.0;

        struct point normal = getNormal(intersection);
        struct point reflection = getReflection(r, normal);

        for(int i = 0; i < lights.size(); i++){
            struct point ldir;
            ldir.x = lights[i].x - intersection.x;
            ldir.y = lights[i].y - intersection.y;
            ldir.z = lights[i].z - intersection.z;

            double len = sqrt(pow(ldir.x,2) + pow(ldir.y,2) + pow(ldir.z,2));

            ldir = normalize(ldir);

            struct point lstart;
            lstart.x = intersection.x + ldir.x * 1.0;
            lstart.y = intersection.y + ldir.y * 1.0;
            lstart.z = intersection.z + ldir.z * 1.0;

            Ray L(lstart, ldir);

            bool obscured = false;

            for(int j = 0; j < objects.size(); j++){
                double possibleobscure = objects[j]->getT(&L);
                //printf("%f ", possibleobscure);
                if(possibleobscure > 0 && abs(possibleobscure) < len){
                    obscured = true;
                    break;
                }
            }
            //printf("\n");

            if(!obscured){
                //printf("direct");
                //L.dir.x = -L.dir.x;
                //L.dir.y = -L.dir.y;
                //L.dir.z = -L.dir.z;
                double lambert = dotProduct(L.dir, normal);
                //double phong = dotProduct(reflection, r->dir);
                struct point oppo;
                oppo.x = -r->dir.x;
                oppo.y = -r->dir.y;
                oppo.z = -r->dir.z;
                double phong = dotProduct(getReflection(&L, normal), oppo);
                //if(lambert > 0 || phong > 0) printf("emon to hobar kotha na\n");
                lambert = lambert > 0?lambert:0;
                phong = phong > 0?phong:0;

                current_color->x += source_factor * color[0] * (lambert * co_efficients[DIFFUSE] + pow(phong, Shine) * co_efficients[SPECULAR]);
                current_color->y += source_factor * color[1] * (lambert * co_efficients[DIFFUSE] + pow(phong, Shine) * co_efficients[SPECULAR]);
                current_color->z += source_factor * color[2] * (lambert * co_efficients[DIFFUSE] + pow(phong, Shine) * co_efficients[SPECULAR]);
            }

            if(level < recursion_level){
                struct point rs;
                rs.x = intersection.x + reflection.x * 1.0;
                rs.y = intersection.y + reflection.y * 1.0;
                rs.z = intersection.z + reflection.z * 1.0;

                //reflection = normalize(reflection);

                Ray reflectionRay(rs, reflection);
                int nearest = -1;

                struct point reflectColor;
                double t_min = 9999;
                for(int k = 0; k < objects.size(); k++){
                    double t = objects[k]->getT(&reflectionRay);

                    if(t <= 0) continue;
                    else if(t < t_min) {
                        t_min = t;
                        nearest = k;
                    }
                }

                if(nearest != -1){
                    double t = objects[nearest]->intersect(&reflectionRay, &reflectColor, level + 1);
                    if(t!=-1){
                        current_color->x += reflectColor.x*co_efficients[REFLECTION];
                        current_color->y += reflectColor.y*co_efficients[REFLECTION];
                        current_color->z += reflectColor.z*co_efficients[REFLECTION];
                    }
                }
            }

            //Check whether all current_color pixel value is within 1 or 0 if not set it
        if(current_color->x > 1)    current_color->x = 1;
        if(current_color->x < 0)    current_color->x = 0;
        if(current_color->y > 1)    current_color->y = 1;
        if(current_color->y < 0)    current_color->y = 0;
        if(current_color->z > 1)    current_color->z = 1;
        if(current_color->z < 0)    current_color->z = 0;
        }

        return t;
    }


    void draw(){
        struct point temp = reference_point;

        int black = 1;
        while(temp.y < -reference_point.y){
            while(temp.x < -reference_point.x){
                if(black == 1){
                    glColor3f(0, 0, 0);
                }
                else{
                    glColor3f(1, 1, 1);
                }

                black = 1 - black;
                glBegin(GL_QUADS);
                    glVertex3f(temp.x, temp.y, temp.z);
                    glVertex3f(temp.x, temp.y + length, temp.z);
                    glVertex3f(temp.x + length, temp.y + length, temp.z);
                    glVertex3f(temp.x + length, temp.y, temp.z);
                glEnd();

                temp.x += length;
            }
            temp.x = reference_point.x;
            temp.y += length;
            black = 1 - black;
        }
    }
};

class Triangle: public Object
{
public:
    struct point a, b, c;
    Triangle(struct point a, struct point b, struct point c){
        this->a = a;
        this->b = b;
        this->c = c;
    }

    void draw(){
        glColor3f(color[0], color[1], color[2]);
        glBegin(GL_TRIANGLES);
        {
            glVertex3f(a.x, a.y, a.z);
            glVertex3f(b.x, b.y, b.z);
            glVertex3f(c.x, c.y, c.z);
        }
        glEnd();
    }

    struct point getNormal(struct point pointed)
    {
        struct point u, v;
        u.x = b.x - a.x;
        u.y = b.y - a.y;
        u.z = b.z - a.z;

        v.x = c.x - a.x;
        v.y = c.y - a.y;
        v.z = c.z - a.z;
        struct point normal = crossProduct(u, v);
        normal = normalize(normal);

        if(dotProduct(pointed, normal) > 0){
            normal.x = -normal.x;
            normal.y = -normal.y;
            normal.z = -normal.z;
        }

        return normal;
    }

    double getT(Ray *r)
    {
        struct point edge1, edge2;
        edge1.x = b.x - a.x;
        edge1.y = b.y - a.y;
        edge1.z = b.z - a.z;

        edge2.x = c.x - a.x;
        edge2.y = c.y - a.y;
        edge2.z = c.z - a.z;

        struct point tempostart = r->start;
        struct point tempodir = r->dir;

        struct point h = crossProduct(tempodir, edge2);
        double aa = dotProduct(edge1, h);
        if(aa > -EPSILON && aa < EPSILON){
            return -1;
        }
        double f = 1/aa;

        struct point s;
        s.x = tempostart.x - a.x;
        s.y = tempostart.y - a.y;
        s.z = tempostart.z - a.z;

        double u = f * dotProduct(s, h);
        if(u < 0.0 || u > 1.0) return -1;

        struct point q = crossProduct(s, edge1);
        double v = f * dotProduct(tempodir, q);
        if(v < 0.0 || u + v > 1.0) return -1;

        double t = f * dotProduct(edge2, q);
        if(t > EPSILON) return t;

        return -1;

    }

    double intersect(Ray *r, struct point* current_color, int level)
    {
        double t = getT(r);
        if(t<=0){
            return -1;
        }
        if(level == 0){
            return t;
        }

        struct point intersection;
        intersection.x = r->start.x + r->dir.x * t;
        intersection.y = r->start.y + r->dir.y * t;
        intersection.z = r->start.z + r->dir.z * t;

        current_color->x = color[0] * co_efficients[AMBIENT];
        current_color->y = color[1] * co_efficients[AMBIENT];
        current_color->z = color[2] * co_efficients[AMBIENT];


        struct point normal = getNormal(r->dir);
        struct point reflection = getReflection(r, normal);

        for(int i = 0; i < lights.size(); i++){
            struct point ldir;
            ldir.x = lights[i].x - intersection.x;
            ldir.y = lights[i].y - intersection.y;
            ldir.z = lights[i].z - intersection.z;

            double len = sqrt(pow(ldir.x,2) + pow(ldir.y,2) + pow(ldir.z,2));
            ldir = normalize(ldir);

            struct point lstart;
            lstart.x = intersection.x + ldir.x * 1.0;
            lstart.y = intersection.y + ldir.y * 1.0;
            lstart.z = intersection.z + ldir.z * 1.0;


            Ray L(lstart, ldir);

            bool obscured = false;

            for(int j = 0; j < objects.size(); j++){
                double possibleobscure = objects[j]->getT(&L);
                //printf(",%f %d, ", possibleobscure, j);
                if(possibleobscure > 0 && abs(possibleobscure) <= len){
                    obscured = true;
                    break;
                }
            }
            //printf("\n");

            if(!obscured){
                //printf("direct");

                double lambert = dotProduct(L.dir, normal);
                double phong = dotProduct(reflection, r->dir);

                //if(lambert > 0 || phong > 0) printf("emon to hobar kotha na\n");
                lambert = lambert > 0?lambert:0;
                phong = phong > 0?phong:0;

                current_color->x += source_factor * color[0] * (lambert * co_efficients[DIFFUSE] + pow(phong, Shine) * co_efficients[SPECULAR]);
                current_color->y += source_factor * color[1] * (lambert * co_efficients[DIFFUSE] + pow(phong, Shine) * co_efficients[SPECULAR]);
                current_color->z += source_factor * color[2] * (lambert * co_efficients[DIFFUSE] + pow(phong, Shine) * co_efficients[SPECULAR]);
            }

            if(level < recursion_level){
                struct point rs;
                rs.x = intersection.x + reflection.x * 1.0;
                rs.y = intersection.y + reflection.y * 1.0;
                rs.z = intersection.z + reflection.z * 1.0;

                //reflection = normalize(reflection);

                Ray reflectionRay(rs, reflection);
                int nearest = -1;

                struct point reflectColor;
                double t_min = 9999;
                for(int k = 0; k < objects.size(); k++){
                    double t = objects[k]->getT(&reflectionRay);

                    if(t <= 0) continue;
                    else if(t < t_min) {
                        t_min = t;
                        nearest = k;
                    }
                }

                if(nearest != -1){
                    double t = objects[nearest]->intersect(&reflectionRay, &reflectColor, level + 1);
                    if(t != -1){
                        current_color->x += reflectColor.x*co_efficients[REFLECTION];
                        current_color->y += reflectColor.y*co_efficients[REFLECTION];
                        current_color->z += reflectColor.z*co_efficients[REFLECTION];
                    }


                }
            }

            //Check whether all current_color pixel value is within 1 or 0 if not set it
        if(current_color->x > 1)    current_color->x = 1;
        if(current_color->x < 0)    current_color->x = 0;
        if(current_color->y > 1)    current_color->y = 1;
        if(current_color->y < 0)    current_color->y = 0;
        if(current_color->z > 1)    current_color->z = 1;
        if(current_color->z < 0)    current_color->z = 0;
        }

        return t;
    }
};

class generalQuadratic: public Object
{
public:
    double A, B, C, D, E, F, G, H, I, J;
    generalQuadratic(double A, double B, double C, double D, double E, double F, double G, double H, double I, double J, struct point reference_point, double length, double width, double height){
        this->A = A;
        this->B = B;
        this->C = C;
        this->D = D;
        this->E = E;
        this->F = F;
        this->G = G;
        this->H = H;
        this->I = I;
        this->J = J;

        this->reference_point = reference_point;
        this->length = length;
        this->width = width;
        this->height = height;
    }

    void draw(){

    }

    struct point getNormal(struct point intersect)
    {
        struct point normal;
        normal.x = 2 * A * intersect.x + D * intersect.y + E * intersect.z + G;
        normal.y = 2 * B * intersect.y + D * intersect.x + F * intersect.z + H;
        normal.z = 2 * C * intersect.z + E * intersect.x + F * intersect.y + I;

        return normalize(normal);
    }

    double getT(Ray *r){
        //printf("%f %f %f %f %f %f %f %f %f %f\n", A, B, C, D, E, F, G, H, I, J);
        double xo = r->start.x, xd = r->dir.x, yo = r->start.y, yd = r->dir.y, zo = r->start.z, zd = r->dir.z;
        double a = A * pow(xd,2) + B * pow(yd,2) + C * pow(zd,2) + D * xd * yd + E * xd * zd + F * yd * zd;
        double b = 2 * A * xo * xd + 2 * B * yo * yd + 2 * C * zo * zd + D * (xo * yd + yo * xd) + E * (xo * zd + zo * xd) + F * (yo * zd + yd * zo) + G * xd + H * yd + I * zd;
        double c = A * pow(xo,2) + B * pow(yo,2) + C * pow(zo,2) + D * xo * yo + E * xo * zo + F * yo * zo + G * xo + H * yo + I * zo + J;

        if(a == 0){
            //printf("%f %f %f\n\n", a, b, c);
            return -c/b;
        }
        double d = pow(b,2) - 4 * a * c;
        //printf("%f\n",d);
        //if(d>0) printf("\n\n\nohhhh\n\n\n\n");
        if(d < 0) return -1;

        //printf("no?\n");
        double t0 = (-b - sqrt(d))/(2 * a);
        double t1 = (-b + sqrt(d))/(2 * a);

        struct point intersection1, intersection2;
        intersection1.x = r->start.x + r->dir.x * t0;
        intersection1.y = r->start.y + r->dir.y * t0;
        intersection1.z = r->start.z + r->dir.z * t0;

        intersection2.x = r->start.x + r->dir.x * t1;
        intersection2.y = r->start.y + r->dir.y * t1;
        intersection2.z = r->start.z + r->dir.z * t1;

        bool true1 = true, true2 = true;

        if(length != 0){
            if((intersection1.x < reference_point.x) || (intersection1.x > reference_point.x + length)) true1 = false;
            if((intersection2.x < reference_point.x) || (intersection2.x > reference_point.x + length)) true2 = false;
        }
        if(width != 0){
            if((intersection1.y < reference_point.y) || (intersection1.y > reference_point.y + width)) true1 = false;
            if((intersection2.y < reference_point.y) || (intersection2.y > reference_point.y + width)) true2 = false;

        }
        if(height != 0){
            if((intersection1.z < reference_point.z) || (intersection1.z > reference_point.z + height)) true1 = false;
            if((intersection2.z < reference_point.z) || (intersection2.z > reference_point.z + height)) true2 = false;
        }
        //printf("%f %f\n", t0, t1);
        if(true1 && true2) return t0<t1?t0:t1;
        else if(true1) return t0;
        else if(true2) return t1;
        return -1;
    }

    double intersect(Ray *r, struct point* current_color, int level)
    {
        double t = getT(r);
        //printf("%f ", t);
        if(t<=0) return -1;

        if(level == 0) return t;
        //printf("next\n");
        struct point intersection;
        intersection.x = r->start.x + r->dir.x * t;
        intersection.y = r->start.y + r->dir.y * t;
        intersection.z = r->start.z + r->dir.z * t;

        current_color->x = color[0] * co_efficients[AMBIENT];
        current_color->y = color[1] * co_efficients[AMBIENT];
        current_color->z = color[2] * co_efficients[AMBIENT];

        struct point normal = getNormal(intersection);
        struct point reflection = getReflection(r, normal);

        for(int i = 0; i < lights.size(); i++){
            struct point ldir;
            ldir.x = lights[i].x - intersection.x;
            ldir.y = lights[i].y - intersection.y;
            ldir.z = lights[i].z - intersection.z;

            double len = sqrt(pow(ldir.x,2) + pow(ldir.y,2) + pow(ldir.z,2));
            ldir = normalize(ldir);

            struct point lstart;
            lstart.x = intersection.x + ldir.x * 1.0;
            lstart.y = intersection.y + ldir.y * 1.0;
            lstart.z = intersection.z + ldir.z * 1.0;


            Ray L(lstart, ldir);

            bool obscured = false;

            for(int j = 0; j < objects.size(); j++){
                double possibleobscure = objects[j]->getT(&L);
                //printf(",%f %d, ", possibleobscure, j);
                if(possibleobscure > 0 && abs(possibleobscure) < len){
                    obscured = true;
                    break;
                }
            }
            //printf("\n");

            if(!obscured){
                //printf("direct");
                //L.dir.x = -L.dir.x;
                //L.dir.y = -L.dir.y;
                //L.dir.z = -L.dir.z;
                double lambert = dotProduct(L.dir, normal);
                double phong = dotProduct(reflection, r->dir);


                lambert = lambert > 0?lambert:0;
                phong = phong > 0?phong:0;

                current_color->x += source_factor * color[0] * (lambert * co_efficients[DIFFUSE] + pow(phong, Shine) * co_efficients[SPECULAR]);
                current_color->y += source_factor * color[1] * (lambert * co_efficients[DIFFUSE] + pow(phong, Shine) * co_efficients[SPECULAR]);
                current_color->z += source_factor * color[2] * (lambert * co_efficients[DIFFUSE] + pow(phong, Shine) * co_efficients[SPECULAR]);
            }

            if(level < recursion_level){
                struct point rs;
                rs.x = intersection.x + reflection.x * 1.0;
                rs.y = intersection.y + reflection.y * 1.0;
                rs.z = intersection.z + reflection.z * 1.0;

                //reflection = normalize(reflection);

                Ray reflectionRay(rs, reflection);
                int nearest = -1;

                struct point reflectColor;
                double t_min = 9999;
                for(int k = 0; k < objects.size(); k++){
                    double t = objects[k]->getT(&reflectionRay);

                    if(t <= 0) continue;
                    else if(t < t_min) {
                        t_min = t;
                        nearest = k;
                    }
                }

                if(nearest != -1){
                    double t = objects[nearest]->intersect(&reflectionRay, &reflectColor, level + 1);
                    if(t != -1){
                        current_color->x += reflectColor.x*co_efficients[REFLECTION];
                        current_color->y += reflectColor.y*co_efficients[REFLECTION];
                        current_color->z += reflectColor.z*co_efficients[REFLECTION];
                    }
                }
            }

            //Check whether all current_color pixel value is within 1 or 0 if not set it
        if(current_color->x > 1)    current_color->x = 1;
        if(current_color->x < 0)    current_color->x = 0;
        if(current_color->y > 1)    current_color->y = 1;
        if(current_color->y < 0)    current_color->y = 0;
        if(current_color->z > 1)    current_color->z = 1;
        if(current_color->z < 0)    current_color->z = 0;
        }

        return t;
    }
};

#endif // _classes_h
