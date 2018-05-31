#include<stdio.h>
#include<fstream>
#include<string>
#include<iostream>
#include<sstream>
#include<vector>
#include<cmath>
#include<iomanip>
#include<random>
#include "bitmap_image.hpp"

#define PI 3.14159265

using namespace std;

typedef vector<double> Point;
typedef vector<unsigned char> Color;

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
public:
    vector<Point> theTriangle;
    Color color;
    Triangle(){
    }
    void generate_color(){
        for(int i = 0; i < 3; i++)  color.push_back((unsigned char)(rand() % 255));
    }

    ~Triangle(){
        //delete theTriangle;
    }
};

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
            tri[i].generate_color();

            point.clear();
        }
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

void apply_Z_buffer_algo(Triangle tri[], int n)
{
    for(int i = 0; i < n; i++){
        Triangle one = tri[i];
        double top_scan_line = -9999.0;
        double bottom_scan_line = 9999.0;
        for(int j = 0; j < 3; j++){
            if(one.theTriangle[j][1] > top_scan_line) top_scan_line = one.theTriangle[j][1];
            if(one.theTriangle[j][1] < bottom_scan_line) bottom_scan_line = one.theTriangle[j][1];
        }
        printf("%f %f\n", top_scan_line, bottom_scan_line);
    }
}

void write_Z_buffer()
{
    ofstream fout("z_buffer.txt");
    for(int i = 0; i < Screen_Height; i++){
        for(int j = 0; j < Screen_Width; j++){
            if(Z_buffer[i][j] < z_max){
                fout <<std::fixed <<std::setprecision(6)<< Z_buffer[i][j]<<" ";
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
    //printf("%f %f\n", dx, dy);

    Top_y = 1 - (dy / 2);
    Left_x = -1 + (dx / 2);
    //printf("%f %f\n", Top_y, Left_x);

    z_max = z_front_limit>z_rear_limit?z_front_limit:z_rear_limit;

    set_Z_buffer();

    bitmap_image image(Screen_Width, Screen_Height);
    ///not sure if the following nested loop is necessary, i.e., if we
    ///really need to set the background color black, since the image
    ///looks black just after being initialized. will see later
    /*for(int i = 0; i < Screen_Height; i++){
        for(int j = 0; j < Screen_Width; j++){
            image.set_pixel(i, j, 0, 0, 0);
        }
    }*/

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

    print_triangle(tri, lines / 3);

    apply_Z_buffer_algo(tri, lines/3);

    image.save_image("1.bmp");

    write_Z_buffer();
}
