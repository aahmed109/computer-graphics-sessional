#include<stdio.h>
#include<fstream>
#include<string>
#include<iostream>
#include<sstream>
#include<vector>
#include<algorithm>
#include<cmath>
#include<iomanip>
#include<random>
#include "bitmap_image.hpp"

#define PI 3.14159265

using namespace std;

typedef vector<double> Point;
typedef vector<int> Color;

int Screen_Width;
int Screen_Height;
double x_left_limit;
double x_right_limit;
double y_top_limit;
double y_bottom_limit;
double z_front_limit;
double z_rear_limit;

double z_max;

double dx;
double dy;

double Top_y;
double Left_x;

double** Z_buffer;

class Triangle{
private:
    Point vector_from_points(Point A, Point B)
    {
        Point V_1;
        V_1.push_back(B[0] - A[0]);
        V_1.push_back(B[1] - A[1]);
        V_1.push_back(B[2] - A[2]);

        return V_1;
    }

    Point cross_product(Point A, Point B)
    {
        Point cross;

        cross.push_back(A[1] * B[2] - A[2] * B[1]);
        cross.push_back(A[2] * B[0] - A[0] * B[2]);
        cross.push_back(A[0] * B[1] - A[1] * B[0]);

        return cross;
    }

public:
    double plane_A, plane_B, plane_C, plane_D;
    vector<Point> theTriangle;
    Color color;
    Triangle(){
    }
    void generate_color(){
        for(int i = 0; i < 3; i++)  color.push_back((rand() % 256));
    }

    void plane_equation(){
        Point A_not = theTriangle.front();
        Point V_1, V_2;

        V_1 = vector_from_points(theTriangle.front(), theTriangle[1]);
        V_2 = vector_from_points(theTriangle.front(), theTriangle[2]);

        Point product = cross_product(V_1, V_2);

        plane_A = product[0];
        plane_B = product[1];
        plane_C = product[2];
        plane_D = -(product[0] * A_not[0] + product[1] * A_not[1] + product[2] * A_not[2]);
    }

    ~Triangle(){

    }
};



vector<double> get_intersect(Triangle& t, double y)
{
    vector<double> v;
    for(int k = 0; k < 3; k++)
    {
        double ax = t.theTriangle[k][0];
        double ay = t.theTriangle[k][1];

        double bx = t.theTriangle[(k + 1) % 3][0];
        double by = t.theTriangle[(k + 1) % 3][1];

        if(ay == by) continue;

        if((y >= ay && y <= by) || (y >= by && y <= ay))
        {
            double x = (bx - ax) * (y - ay) / (by - ay) + ax;
            v.push_back(x);
        }
    }
    sort(v.begin(), v.end());
    return v;
}


void read_config_file(string configs[])
{
    ifstream fin("config.txt");
    string str;

    int lc = 0;
    while(getline(fin, str)) configs[lc++] = str;
    fin.close();

    for(int i = 0; i < 4; i++){
        stringstream ss(configs[i]);

        if(i == 0) {
            ss>>Screen_Width;
            ss>>Screen_Height;
        }
        else if(i == 1) {
            ss>>x_left_limit;
            x_right_limit = -x_left_limit;
        }
        else if(i == 2) {
            ss>>y_bottom_limit;
            y_top_limit = -y_bottom_limit;
        }
        else{
            ss>>z_front_limit;
            ss>>z_rear_limit;
        }
    }
}

void set_Z_buffer()
{
    Z_buffer = new double*[Screen_Height];
    for(int i = 0; i < Screen_Height; i++){
        Z_buffer[i] = new double[Screen_Width];
    }

    for(int i = 0; i < Screen_Height; i++){
        for(int j = 0; j < Screen_Width; j++){
            Z_buffer[i][j] = z_max;
        }
    }
}

void read_stage_3_file(string triangles[])
{
    ifstream fin("stage3.txt");
    string str;

    int lc = 0;
    while(getline(fin, str)){
        if(str.compare("")) triangles[lc++] = str;
    }

    fin.close();
}

void find_triangles(Triangle tri[], string triangles[], int n)
{
    Point point;
    for(int i = 0; i < n; i++){
        for(int j = 0; j < 3; j++){
            stringstream ss(triangles[i * 3 + j]);
            double s;
            for(int k = 0; k < 3; k++){
                ss>>s;
                point.push_back(s);
            }
            tri[i].theTriangle.push_back(point);
            point.clear();
        }
        tri[i].generate_color();

        tri[i].plane_equation();
    }
}

void print_triangle(Triangle tri[], int n)
{
        for(int i = 0; i < n; i++){
        for(int j = 0; j < 3; j++){
            for(int k = 0; k < 3; k++){
                printf("%f ", tri[i].theTriangle[j][k]);
            }
            cout<<endl;
        }
        cout<<endl;

        for(int j = 0; j < 3; j++){
            printf("%u ", tri[i].color[j]);
        }
        cout<<endl;
        cout<<endl;
    }
}

void apply_Z_buffer_algo(bitmap_image *image, Triangle tri[], int n)
{
    for(int i = 0; i < n; i++){

        Triangle one = tri[i];

        double top_scan_line = -9999.0;
        double bottom_scan_line = 9999.0;
        for(int j = 0; j < 3; j++){
            if(one.theTriangle[j][1] > top_scan_line) top_scan_line = one.theTriangle[j][1];
            if(one.theTriangle[j][1] < bottom_scan_line) bottom_scan_line = one.theTriangle[j][1];
        }

        if(top_scan_line > y_top_limit) top_scan_line = Top_y;

        if(bottom_scan_line < y_bottom_limit) bottom_scan_line = y_bottom_limit;

        int top_pixel = (int)((Top_y - top_scan_line) / dy);
        int bottom_pixel = (int)((Top_y - bottom_scan_line) / dy);

        int jj = top_pixel;

        double left_X = 9999.0, right_X = -9999.0;

        while(jj <= bottom_pixel){
            double j = Top_y - jj*dy;

            vector<double> v = get_intersect(one, j);
            if (v.size() < 2)
            {
                jj++;
                continue;
            }
            left_X = v[0];
            right_X = v[1];

            if(left_X < x_left_limit){
                left_X = Left_x;
            }
            if(right_X > x_right_limit){
                right_X = x_right_limit;
            }

            int left_pixel = (left_X - Left_x) / dx;
            int right_pixel = (right_X - Left_x) / dx;

            int kk = left_pixel;

            while(kk <= right_pixel){
                double k = Left_x + kk*dx;

                double z = (-one.plane_D - one.plane_A * k - one.plane_B * j) / one.plane_C;

                if(kk > Screen_Width - 1 || jj > Screen_Height -1) {
                    break;
                }

                if(z >= z_front_limit && z <= z_rear_limit && Z_buffer[jj][kk] > z){
                    Z_buffer[jj][kk] = z;
                    image->set_pixel(kk, jj, (int) one.color[0], (int) one.color[1], (int) one.color[2]);
                }

                kk++;
            }
            jj++;
        }
    }
}

void write_Z_buffer()
{
    ofstream fout("z_buffer.txt");
    for(int i = 0; i < Screen_Height; i++){
        for(int j = 0; j < Screen_Width; j++){
            if(Z_buffer[i][j] < z_max){
                fout <<std::fixed <<std::setprecision(6)<< Z_buffer[i][j]<<"\t";
            }
        }
        fout<<"\n";
    }
    fout.close();
}

int main()
{
    string configs[4];
    read_config_file(configs);

    dx = (x_right_limit - x_left_limit) / Screen_Width;
    dy = (y_top_limit - y_bottom_limit) / Screen_Height;


    Top_y = y_top_limit - (dy / 2);
    Left_x = x_left_limit + (dx / 2);

    z_max = z_front_limit>z_rear_limit?z_front_limit:z_rear_limit;

    set_Z_buffer();

    bitmap_image image(Screen_Width, Screen_Height);

    ifstream fin("stage3.txt");
    string str;

    int lines = 0;

    while(getline(fin, str)) {
        if(str.compare("")) lines++;
    }
    fin.close();

    string triangles[lines];
    read_stage_3_file(triangles);

    Triangle tri[lines / 3];

    find_triangles(tri, triangles, lines / 3);

    apply_Z_buffer_algo(&image, tri, lines/3);

    image.save_image("1.bmp");

    write_Z_buffer();
}
