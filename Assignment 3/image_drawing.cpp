#include "bitmap_image.hpp"

using namespace std;

typedef vector<unsigned char> Color;
int main(){

    bitmap_image image(200,100);
    Color color;
    color.push_back(100);
    color.push_back(100);
    color.push_back(100);
    for(int i=0;i<200;i++){
        for(int j=0;j<100;j++){
            image.set_pixel(i,j,(int) color[0], (int) color[1], (int) color[2]);
        }
    }

    image.save_image("test.bmp");;

    return 0;
}
