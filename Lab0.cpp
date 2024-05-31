#include <GL/glut.h>
#include <GL/gl.h>
#include <cmath> 
#include <cstdlib>
#include <cstdio>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>


#define ROWS 3
#define COLUMNS 3

using namespace std;
string filename;
string in_path;

int checkpoint = 0;
vector<float> vertex_x_arr;
vector<float> vertex_y_arr;
vector<float> vertex_Red_arr;
vector<float> vertex_Green_arr;
vector<float> vertex_Blue_arr;
vector<float> mode_arr;

vector<float> vertex_x_arr_view1;
vector<float> vertex_y_arr_view1;
vector<float> vertex_Red_arr_view1;
vector<float> vertex_Green_arr_view1;
vector<float> vertex_Blue_arr_view1;
vector<float> mode_arr_view1;

vector<float> vertex_x_arr_view2;
vector<float> vertex_y_arr_view2;
vector<float> vertex_Red_arr_view2;
vector<float> vertex_Green_arr_view2;
vector<float> vertex_Blue_arr_view2;
vector<float> mode_arr_view2;

vector<float> vertex_x_arr_view3;
vector<float> vertex_y_arr_view3;
vector<float> vertex_Red_arr_view3;
vector<float> vertex_Green_arr_view3;
vector<float> vertex_Blue_arr_view3;
vector<float> mode_arr_view3;

vector<float> vertex_x_arr_view;
vector<float> vertex_y_arr_view;
vector<float> vertex_Red_arr_view;
vector<float> vertex_Green_arr_view;
vector<float> vertex_Blue_arr_view;
vector<float> mode_arr_view;

vector<int> x_point_arr;
vector<int> y_point_arr;
vector<float> Red_arr;
vector<float> Green_arr;
vector<float> Blue_arr;
vector<float> view_mode_arr;
//correct
void arrange(string str);
void Inner_product(float Operator[][COLUMNS], float Transform[][COLUMNS]);
void SCALE(float scale_x, float scale_y);
void ROTATE(float angle);
void TRANSLATE(float translate_x, float translate_y);
float random_func();
void randomRGB();
void drawALine(int x_start, int x_end, int y_start, int y_end, float Red, float Green, float Blue);
void SaveVertex(float x_point, float y_point, float SaveRed, float SaveGreen, float SaveBlue);
void TRIANGLE();
void SQUARE();
void Line_and_Point(float x_start, float x_end, float y_start, float y_end, float x_point, float y_point);

float tmp[8] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
float tmp_scale[2];
float view_boundary_x[5];
float view_boundary_y[5];

float Transform_matrix[3][3] = {
    {1.0, 0.0, 0.0},
    {0.0, 1.0, 0.0},
    {0.0, 0.0, 1.0}
};
float Unit_matrix[3][3] = {
    {1.0, 0.0, 0.0},
    {0.0, 1.0, 0.0},
    {0.0, 0.0, 1.0}
};

float Point_matrix[3] = { 0.0, 0.0, 0.0 };
string InOrOut = "in";
float Red, Green, Blue;
int i;



void Multiply_point(float Operator[][3], float Point[])
{
    float Zero_matrix[3] = { 0.0, 0.0, 0.0 };
    for (int i = 0; i < ROWS; i++) {
        for (int k = 0; k < 3; k++) {
            Zero_matrix[i] = Zero_matrix[i] + Operator[i][k] * Point[k];
            Point_matrix[i] = Zero_matrix[i];
        }
    }
    for (int i = 0; i < ROWS; i++) {
        //cout << "origin point: " << Point[i] << endl;
        //cout << "transform point: " << Point_matrix[i] << endl;
    }

}

void CLEARDATA() {
    vertex_x_arr.clear();
    vertex_y_arr.clear();
    vertex_Red_arr.clear();
    vertex_Green_arr.clear();
    vertex_Blue_arr.clear();
    mode_arr.clear();
}
int mode_add_minus = 0;
void Show() {
    int point = 0;
    glPointSize(1);
    glBegin(GL_POINTS);
    drawALine(tmp[4], tmp[5], tmp[6], tmp[6], 1.0, 1.0, 1.0);
    drawALine(tmp[4], tmp[5], tmp[7], tmp[7], 1.0, 1.0, 1.0);
    drawALine(tmp[4], tmp[4], tmp[6], tmp[7], 1.0, 1.0, 1.0);
    drawALine(tmp[5], tmp[5], tmp[6], tmp[7], 1.0, 1.0, 1.0);
    
    for (int i = 0; i < mode_arr.size(); i++) {
        point++;
        for (int j = 0; j < mode_arr_view[i]; j++) {
            drawALine((int)round(vertex_x_arr_view[point]), (int)round(vertex_x_arr_view[point - 1]), (int)round(vertex_y_arr_view[point]), (int)round(vertex_y_arr_view[point - 1]), vertex_Red_arr_view[point], vertex_Green_arr_view[point], vertex_Blue_arr_view[point]);
            point++;
        }
    }

    int left = view_boundary_x[1];
    int right = view_boundary_x[2];
    int top = view_boundary_y[0];
    int bottom = view_boundary_y[1];
    int INorOUT = 0;
    float RR = 1.0, GG = 1.0, BB = 1.0;
    for (int yy = bottom; yy < top; yy++) {
        for (int xx = left; xx < right; xx++) {
            point = 0;
            //cout << "(" << xx << ", " << yy << " )" << endl;
            for (int i = 0; i < mode_arr.size(); i++) {
                point++;
                INorOUT = 0;
                for (int j = 0; j < mode_arr_view[i]; j++) {
                    //drawALine((int)round(vertex_x_arr_view[point]), (int)round(vertex_x_arr_view[point - 1]), (int)round(vertex_y_arr_view[point]), (int)round(vertex_y_arr_view[point - 1]), vertex_Red_arr_view[point], vertex_Green_arr_view[point], vertex_Blue_arr_view[point]);
                    if (j == 0) {
                        RR = vertex_Red_arr_view[point];
                        GG = vertex_Green_arr_view[point];
                        BB = vertex_Blue_arr_view[point];
                    }
                    Line_and_Point(vertex_x_arr_view[point], vertex_x_arr_view[point - 1], vertex_y_arr_view[point], vertex_y_arr_view[point - 1], xx, yy);
                    if (InOrOut == "out") {
                        //cout << "out" << endl;
                        INorOUT++;
                    }
                    point++;
                }
                if (INorOUT == mode_arr_view[i] || INorOUT == 0) {
                    for (int j = 0; j < mode_arr_view[i]; j++) {
                        //drawALine((int)round(vertex_x_arr_view[point]), (int)round(vertex_x_arr_view[point - 1]), (int)round(vertex_y_arr_view[point]), (int)round(vertex_y_arr_view[point - 1]), vertex_Red_arr_view[point], vertex_Green_arr_view[point], vertex_Blue_arr_view[point]);
                        glColor3f(RR, GG, BB);
                        glVertex2i(xx, yy);
                    }
                }
            }
        }
    }

    glEnd();
    glutSwapBuffers();
    vertex_x_arr_view1.clear();
    vertex_y_arr_view1.clear();
    vertex_Red_arr_view1.clear();
    vertex_Green_arr_view1.clear();
    vertex_Blue_arr_view1.clear();
    mode_arr_view1.clear();
    vertex_x_arr_view2.clear();
    vertex_y_arr_view2.clear();
    vertex_Red_arr_view2.clear();
    vertex_Green_arr_view2.clear();
    vertex_Blue_arr_view2.clear();
    mode_arr_view2.clear();
    vertex_x_arr_view3.clear();
    vertex_y_arr_view3.clear();
    vertex_Red_arr_view3.clear();
    vertex_Green_arr_view3.clear();
    vertex_Blue_arr_view3.clear();
    mode_arr_view3.clear();
    vertex_x_arr_view.clear();
    vertex_y_arr_view.clear();
    vertex_Red_arr_view.clear();
    vertex_Green_arr_view.clear();
    vertex_Blue_arr_view.clear();
    mode_arr_view.clear();
}


void Line_and_Point(float x_start, float x_end, float y_start, float y_end, float x_point, float y_point) {
    float dx1 = x_end - x_start;
    float dy1 = y_end - y_start;
    float dx2 = x_point - x_start;
    float dy2 = y_point - y_start;
    float z_direction = dx1 * dy2 - dx2 * dy1;
    
    if (z_direction < 0) {
        InOrOut = "out";
    }
    else if (z_direction == 0) {
        InOrOut = "out";
    }
    else if (z_direction > 0) {
        InOrOut = "in";
    }
}

float interaction_x = 0.0;
float interaction_y = 0.0;
void Interaction(float x1, float x2, float x3, float x4, float y1, float y2, float y3, float y4) {
    float m1 = (y1 - y2) / (x1 - x2);
    float m2 = (y3 - y4) / (x3 - x4);
    if (x1 == x2) {
        if (x3 == x4) { //case 5
            interaction_x = x4;
            interaction_y = y4;
        }
        else if (y3 == y4) { //case 3
            interaction_x = x1;
            interaction_y = y3;
        }
        else { //case 2
            interaction_x = x1;
            interaction_y = m2 * (interaction_x - x3) + y3;
        }
    }
    else if (x3 == x4) {
        if (y1 == y2) { //case 7
            interaction_x = x3;
            interaction_y = y1;
        }
        else { //case 4
            interaction_x = x3;
            interaction_y = m1 * (interaction_x - x1) + y1;
        }
    }
    else if (y1 == y2) {
        if (y3 == y4) { //case 9
            interaction_x = x4;
            interaction_y = y4;
        }
        else{ //case6
            interaction_y = y1;
            interaction_x = (interaction_y - y3) / m2 + x3;
        }
    }
    else if (y3 == y4) { //case 8
        interaction_y = y3;
        interaction_x = (interaction_y - y1) / m1 + x1;
    }
    else{ //case 1
        interaction_x = (m1 * x1 - m2 * x3 - y1 + y3) / (m1 - m2);
        interaction_y = m1 * (interaction_x - x1) + y1;
    }
}


void VIEW(float wxl, float wxr, float wyb, float wyt, float vxl, float vxr, float vyb, float vyt) {
    float scale_x = (vxr - vxl) / (wxr - wxl);
    float scale_y = (vyt - vyb) / (wyt - wyb);
    float center_formar_x = 0.5 * (wxr + wxl);
    float center_formar_y = 0.5 * (wyt + wyb);
    float center_latter_x = 0.5 * (vxr + vxl);
    float center_latter_y = 0.5 * (vyt + vyb);
    float translate_x = center_latter_x - center_formar_x;
    float translate_y = center_latter_y - center_formar_y;
    view_boundary_x[0] = vxl;
    view_boundary_x[1] = vxl;
    view_boundary_x[2] = vxr;
    view_boundary_x[3] = vxr;
    view_boundary_x[4] = vxl;
    view_boundary_y[0] = vyt;
    view_boundary_y[1] = vyb;
    view_boundary_y[2] = vyb;
    view_boundary_y[3] = vyt;
    view_boundary_y[4] = vyt;

    float view_matrix[3][3] = {
        {1.0, 0.0, 0.0},
        {0.0, 1.0, 0.0},
        {0.0, 0.0, 1.0}
    };
    float view_inverse_translate_matrix[3][3] = {
        {1.0, 0.0, -center_formar_x},
        {0.0, 1.0, -center_formar_y},
        {0.0, 0.0,              1.0}
    };
    float view_scale_matrix[3][3] = {
        {scale_x, 0.0, 0.0},
        {0.0, scale_y, 0.0},
        {0.0, 0.0, 1.0}
    };
    float view_translate_matrix[3][3] = {
        {1.0, 0.0, center_latter_x},
        {0.0, 1.0, center_latter_y},
        {0.0, 0.0,             1.0}
    };

    for (int i = 0; i < ROWS; i++) {
        for (int j = 0; j < COLUMNS; j++) {
            float a = 0.0;
            for (int k = 0; k < 3; k++) {
                a = a + view_inverse_translate_matrix[i][k] * view_matrix[k][j];
            }
            view_matrix[i][j] = a;
        }
    }
    for (int i = 0; i < ROWS; i++) {
        for (int j = 0; j < COLUMNS; j++) {
            float a = 0.0;
            for (int k = 0; k < 3; k++) {
                a = a + view_scale_matrix[i][k] * view_matrix[k][j];
            }
            view_matrix[i][j] = a;
        }
    }
    for (int i = 0; i < ROWS; i++) {
        for (int j = 0; j < COLUMNS; j++) {
            float a = 0.0;
            for (int k = 0; k < 3; k++) {
                a = a + view_translate_matrix[i][k] * view_matrix[k][j];
            }
            view_matrix[i][j] = a;
        }
    }

    float vertex_point[3] = { 0.0, 0.0, 1.0 };
    //cout << vertex_x_arr.size() << " " << endl;
    //cout << mode_arr.size() << " " << endl;
    int point = 0;
    float tmp_x = 0.0;
    float tmp_y = 0.0;
    int mode_num = 0;
    /*
        for (int i = 0; i < mode_arr.size(); i++) {
        for (int j = 0; j < mode_arr[i] + 1; j++) {
            vertex_point[0] = vertex_x_arr[point];
            vertex_point[1] = vertex_y_arr[point];
            Multiply_point(view_matrix, vertex_point);
            vertex_x_arr_view.push_back(Point_matrix[0]);
            vertex_y_arr_view.push_back(Point_matrix[1]);
            point++;
        }
    }*/
    
    //cout << "view " << "(" << vxl << ", " << vyt << ")" << "(" << vxl << ", " << vyb << ")" << endl;
    
    // left boundary
    point = 0;
    int k = 0;
    for (int i = 0; i < mode_arr.size(); i++) { 
        mode_num = mode_arr[i];
        int tt = 0;
        if (mode_arr[i] == 0) {
            vertex_x_arr_view1.push_back(vertex_x_arr[point]);
            vertex_y_arr_view1.push_back(vertex_y_arr[point]);
            vertex_Red_arr_view1.push_back(vertex_Red_arr[point]);
            vertex_Green_arr_view1.push_back(vertex_Green_arr[point]);
            vertex_Blue_arr_view1.push_back(vertex_Blue_arr[point]);
        }
        for (int j = 0; j < mode_arr[i]; j++) {

            vertex_point[0] = vertex_x_arr[point];
            vertex_point[1] = vertex_y_arr[point];
            Multiply_point(view_matrix, vertex_point);
            float point1_x = Point_matrix[0];
            float point1_y = Point_matrix[1];
            Line_and_Point(view_boundary_x[k], view_boundary_x[k + 1], view_boundary_y[k], view_boundary_y[k + 1], point1_x, point1_y);
            string location1 = InOrOut;

            vertex_point[0] = vertex_x_arr[point + 1];
            vertex_point[1] = vertex_y_arr[point + 1];
            Multiply_point(view_matrix, vertex_point);
            float point2_x = Point_matrix[0];
            float point2_y = Point_matrix[1];
            Line_and_Point(view_boundary_x[k], view_boundary_x[k + 1], view_boundary_y[k], view_boundary_y[k + 1], point2_x, point2_y);
            string location2 = InOrOut;

            //cout << location1 << " " << location2 << endl;
            if (location1 == "in" && location2 == "out") {
                //cout << "point1 = (" << point1_x << ", " << point1_y << ") ";
                //cout << "point2 = (" << point2_x << ", " << point2_y << ")" << endl;
                Interaction(view_boundary_x[k], view_boundary_x[k + 1], point1_x, point2_x, view_boundary_y[k], view_boundary_y[k + 1], point1_y, point2_y);
                vertex_x_arr_view1.push_back(interaction_x);
                vertex_y_arr_view1.push_back(interaction_y);
                vertex_Red_arr_view1.push_back(vertex_Red_arr[point]);
                vertex_Green_arr_view1.push_back(vertex_Green_arr[point]);
                vertex_Blue_arr_view1.push_back(vertex_Blue_arr[point]);
                //cout << " => interaction = (" << interaction_x << ", " << interaction_y << ")" << endl;
                if (tt == 0) {
                    tmp_x = interaction_x;
                    tmp_y = interaction_y;
                    tt = 1;
                }
            }
            else if (location1 == "out" && location2 == "out") {
                mode_num--;
            }
            else if (location1 == "out" && location2 == "in") {
                //cout << "point1 = (" << point1_x << ", " << point1_y << ") ";
                //cout << "point2 = (" << point2_x << ", " << point2_y << ")" << endl;
                Interaction(view_boundary_x[k], view_boundary_x[k + 1], point1_x, point2_x, view_boundary_y[k], view_boundary_y[k + 1], point1_y, point2_y);
                vertex_x_arr_view1.push_back(interaction_x);
                vertex_y_arr_view1.push_back(interaction_y);                    
                vertex_Red_arr_view1.push_back(vertex_Red_arr[point]);
                vertex_Green_arr_view1.push_back(vertex_Green_arr[point]);
                vertex_Blue_arr_view1.push_back(vertex_Blue_arr[point]);
                //cout << " => interaction = (" << interaction_x << ", " << interaction_y << ") ";
                if (tt == 0) {
                    tmp_x = interaction_x;
                    tmp_y = interaction_y;
                    tt = 1;
                }
                vertex_x_arr_view1.push_back(point2_x);
                vertex_y_arr_view1.push_back(point2_y);
                vertex_Red_arr_view1.push_back(vertex_Red_arr[point]);
                vertex_Green_arr_view1.push_back(vertex_Green_arr[point]);
                vertex_Blue_arr_view1.push_back(vertex_Blue_arr[point]);
                //cout << " => point2 = (" << point2_x << ", " << point2_y << ")" << endl;
                mode_num++;
                    
            }
            else if (location1 == "in" && location2 == "in") {
                //cout << "point1 = (" << point1_x << ", " << point1_y << ") ";
                //cout << "point2 = (" << point2_x << ", " << point2_y << ")" << endl;
                vertex_x_arr_view1.push_back(point2_x);
                vertex_y_arr_view1.push_back(point2_y);
                vertex_Red_arr_view1.push_back(vertex_Red_arr[point]);
                vertex_Green_arr_view1.push_back(vertex_Green_arr[point]);
                vertex_Blue_arr_view1.push_back(vertex_Blue_arr[point]);
                //cout << " => point2 = (" << point2_x << ", " << point2_y << ")" << endl;
                if (tt == 0) {
                    tmp_x = point2_x;
                    tmp_y = point2_y;
                    tt = 1;
                }
            }
                
            if (j == mode_arr[i] - 1) {
                vertex_x_arr_view1.push_back(tmp_x);
                vertex_y_arr_view1.push_back(tmp_y);
                vertex_Red_arr_view1.push_back(vertex_Red_arr[point]);
                vertex_Green_arr_view1.push_back(vertex_Green_arr[point]);
                vertex_Blue_arr_view1.push_back(vertex_Blue_arr[point]);
                //cout << " => origin = (" << tmp_x << ", " << tmp_y << ")" << endl;
            }
            point++;
        }
        point++;
        mode_arr_view1.push_back(mode_num);
        //cout << "mode num: " << mode_num << endl;
    }

    // bottom boundary
    point = 0;
    k = 1;
    for (int i = 0; i < mode_arr.size(); i++) {
        mode_num = mode_arr_view1[i];
        int tt = 0;
        if (mode_arr_view1[i] == 0) {
            vertex_x_arr_view2.push_back(vertex_x_arr_view1[point]);
            vertex_y_arr_view2.push_back(vertex_y_arr_view1[point]);
            vertex_Red_arr_view2.push_back(vertex_Red_arr_view1[point]);
            vertex_Green_arr_view2.push_back(vertex_Green_arr_view1[point]);
            vertex_Blue_arr_view2.push_back(vertex_Blue_arr_view1[point]);
        }
        for (int j = 0; j < mode_arr_view1[i]; j++) {
            float point1_x = vertex_x_arr_view1[point];
            float point1_y = vertex_y_arr_view1[point];
            Line_and_Point(view_boundary_x[k], view_boundary_x[k + 1], view_boundary_y[k], view_boundary_y[k + 1], point1_x, point1_y);
            string location1 = InOrOut;

            float point2_x = vertex_x_arr_view1[point + 1];
            float point2_y = vertex_y_arr_view1[point + 1];
            Line_and_Point(view_boundary_x[k], view_boundary_x[k + 1], view_boundary_y[k], view_boundary_y[k + 1], point2_x, point2_y);
            string location2 = InOrOut;

            //cout << location1 << " " << location2 << endl;
            if (location1 == "in" && location2 == "out") {
                //cout << "point1 = (" << point1_x << ", " << point1_y << ") ";
                //cout << "point2 = (" << point2_x << ", " << point2_y << ")" << endl;
                Interaction(view_boundary_x[k], view_boundary_x[k + 1], point1_x, point2_x, view_boundary_y[k], view_boundary_y[k + 1], point1_y, point2_y);
                vertex_x_arr_view2.push_back(interaction_x);
                vertex_y_arr_view2.push_back(interaction_y);
                vertex_Red_arr_view2.push_back(vertex_Red_arr_view1[point]);
                vertex_Green_arr_view2.push_back(vertex_Green_arr_view1[point]);
                vertex_Blue_arr_view2.push_back(vertex_Blue_arr_view1[point]);
                //cout << " => interaction = (" << interaction_x << ", " << interaction_y << ")" << endl;
                if (tt == 0) {
                    tmp_x = interaction_x;
                    tmp_y = interaction_y;
                    tt = 1;
                }
            }
            else if (location1 == "out" && location2 == "out") {
                mode_num--;
            }
            else if (location1 == "out" && location2 == "in") {
                //cout << "point1 = (" << point1_x << ", " << point1_y << ") ";
                //cout << "point2 = (" << point2_x << ", " << point2_y << ")" << endl;
                Interaction(view_boundary_x[k], view_boundary_x[k + 1], point1_x, point2_x, view_boundary_y[k], view_boundary_y[k + 1], point1_y, point2_y);
                vertex_x_arr_view2.push_back(interaction_x);
                vertex_y_arr_view2.push_back(interaction_y);
                vertex_Red_arr_view2.push_back(vertex_Red_arr_view1[point]);
                vertex_Green_arr_view2.push_back(vertex_Green_arr_view1[point]);
                vertex_Blue_arr_view2.push_back(vertex_Blue_arr_view1[point]);
                //cout << " => interaction = (" << interaction_x << ", " << interaction_y << ") ";
                if (tt == 0) {
                    tmp_x = interaction_x;
                    tmp_y = interaction_y;
                    tt = 1;
                }
                vertex_x_arr_view2.push_back(point2_x);
                vertex_y_arr_view2.push_back(point2_y);
                vertex_Red_arr_view2.push_back(vertex_Red_arr_view1[point]);
                vertex_Green_arr_view2.push_back(vertex_Green_arr_view1[point]);
                vertex_Blue_arr_view2.push_back(vertex_Blue_arr_view1[point]);
                //cout << " => point2 = (" << point2_x << ", " << point2_y << ")" << endl;
                mode_num++;

            }
            else if (location1 == "in" && location2 == "in") {
                //cout << "point1 = (" << point1_x << ", " << point1_y << ") ";
                //cout << "point2 = (" << point2_x << ", " << point2_y << ")" << endl;
                vertex_x_arr_view2.push_back(point2_x);
                vertex_y_arr_view2.push_back(point2_y);
                vertex_Red_arr_view2.push_back(vertex_Red_arr_view1[point]);
                vertex_Green_arr_view2.push_back(vertex_Green_arr_view1[point]);
                vertex_Blue_arr_view2.push_back(vertex_Blue_arr_view1[point]);

                //cout << " => point2 = (" << point2_x << ", " << point2_y << ")" << endl;
                if (tt == 0) {
                    tmp_x = point2_x;
                    tmp_y = point2_y;
                    tt = 1;
                }
            }

            if (j == mode_arr_view1[i] - 1) {
                vertex_x_arr_view2.push_back(tmp_x);
                vertex_y_arr_view2.push_back(tmp_y);
                vertex_Red_arr_view2.push_back(vertex_Red_arr_view1[point]);
                vertex_Green_arr_view2.push_back(vertex_Green_arr_view1[point]);
                vertex_Blue_arr_view2.push_back(vertex_Blue_arr_view1[point]);
                //cout << " => origin = (" << tmp_x << ", " << tmp_y << ")" << endl;
            }
            point++;
        }
        point++;

        mode_arr_view2.push_back(mode_num);
        //cout << "mode num: " << mode_num << endl;
    }

    // right boundary
    point = 0;
    k = 2;
    for (int i = 0; i < mode_arr.size(); i++) { // 第幾個多邊形
        mode_num = mode_arr_view2[i];
        int tt = 0;
        if (mode_arr_view2[i] == 0) {
            vertex_x_arr_view3.push_back(vertex_x_arr_view2[point]);
            vertex_y_arr_view3.push_back(vertex_y_arr_view2[point]);
            vertex_Red_arr_view3.push_back(vertex_Red_arr_view2[point]);
            vertex_Green_arr_view3.push_back(vertex_Green_arr_view2[point]);
            vertex_Blue_arr_view3.push_back(vertex_Blue_arr_view2[point]);
        }
        for (int j = 0; j < mode_arr_view2[i]; j++) {
            float point1_x = vertex_x_arr_view2[point];
            float point1_y = vertex_y_arr_view2[point];
            Line_and_Point(view_boundary_x[k], view_boundary_x[k + 1], view_boundary_y[k], view_boundary_y[k + 1], point1_x, point1_y);
            string location1 = InOrOut;

            float point2_x = vertex_x_arr_view2[point + 1];
            float point2_y = vertex_y_arr_view2[point + 1];
            Line_and_Point(view_boundary_x[k], view_boundary_x[k + 1], view_boundary_y[k], view_boundary_y[k + 1], point2_x, point2_y);
            string location2 = InOrOut;

            //cout << location1 << " " << location2 << endl;
            if (location1 == "in" && location2 == "out") {
                //cout << "point1 = (" << point1_x << ", " << point1_y << ") ";
                //cout << "point2 = (" << point2_x << ", " << point2_y << ")" << endl;
                Interaction(view_boundary_x[k], view_boundary_x[k + 1], point1_x, point2_x, view_boundary_y[k], view_boundary_y[k + 1], point1_y, point2_y);
                vertex_x_arr_view3.push_back(interaction_x);
                vertex_y_arr_view3.push_back(interaction_y);
                vertex_Red_arr_view3.push_back(vertex_Red_arr_view2[point]);
                vertex_Green_arr_view3.push_back(vertex_Green_arr_view2[point]);
                vertex_Blue_arr_view3.push_back(vertex_Blue_arr_view2[point]);
                //cout << " => interaction = (" << interaction_x << ", " << interaction_y << ")" << endl;
                if (tt == 0) {
                    tmp_x = interaction_x;
                    tmp_y = interaction_y;
                    tt = 1;
                }
            }
            else if (location1 == "out" && location2 == "out") {
                mode_num--;
            }
            else if (location1 == "out" && location2 == "in") {
                //cout << "point1 = (" << point1_x << ", " << point1_y << ") ";
                //cout << "point2 = (" << point2_x << ", " << point2_y << ")" << endl;
                Interaction(view_boundary_x[k], view_boundary_x[k + 1], point1_x, point2_x, view_boundary_y[k], view_boundary_y[k + 1], point1_y, point2_y);
                vertex_x_arr_view3.push_back(interaction_x);
                vertex_y_arr_view3.push_back(interaction_y);
                vertex_Red_arr_view3.push_back(vertex_Red_arr_view2[point]);
                vertex_Green_arr_view3.push_back(vertex_Green_arr_view2[point]);
                vertex_Blue_arr_view3.push_back(vertex_Blue_arr_view2[point]);
                //cout << " => interaction = (" << interaction_x << ", " << interaction_y << ") ";
                if (tt == 0) {
                    tmp_x = interaction_x;
                    tmp_y = interaction_y;
                    tt = 1;
                }
                vertex_x_arr_view3.push_back(point2_x);
                vertex_y_arr_view3.push_back(point2_y);
                vertex_Red_arr_view3.push_back(vertex_Red_arr_view2[point]);
                vertex_Green_arr_view3.push_back(vertex_Green_arr_view2[point]);
                vertex_Blue_arr_view3.push_back(vertex_Blue_arr_view2[point]);
                //cout << " => point2 = (" << point2_x << ", " << point2_y << ")" << endl;
                mode_num++;

            }
            else if (location1 == "in" && location2 == "in") {
                //cout << "point1 = (" << point1_x << ", " << point1_y << ") ";
                //cout << "point2 = (" << point2_x << ", " << point2_y << ")" << endl;
                vertex_x_arr_view3.push_back(point2_x);
                vertex_y_arr_view3.push_back(point2_y);
                vertex_Red_arr_view3.push_back(vertex_Red_arr_view2[point]);
                vertex_Green_arr_view3.push_back(vertex_Green_arr_view2[point]);
                vertex_Blue_arr_view3.push_back(vertex_Blue_arr_view2[point]);

                //cout << " => point2 = (" << point2_x << ", " << point2_y << ")" << endl;
                if (tt == 0) {
                    tmp_x = point2_x;
                    tmp_y = point2_y;
                    tt = 1;
                }
            }

            if (j == mode_arr_view2[i] - 1) {
                vertex_x_arr_view3.push_back(tmp_x);
                vertex_y_arr_view3.push_back(tmp_y);
                vertex_Red_arr_view3.push_back(vertex_Red_arr_view2[point]);
                vertex_Green_arr_view3.push_back(vertex_Green_arr_view2[point]);
                vertex_Blue_arr_view3.push_back(vertex_Blue_arr_view2[point]);
                //cout << " => origin = (" << tmp_x << ", " << tmp_y << ")" << endl;
            }
            point++;
        }

        point++;

        mode_arr_view3.push_back(mode_num);
        //cout << "mode num: " << mode_num << endl;
    }

    // up boundary
    point = 0;
    k = 3;
    for (int i = 0; i < mode_arr.size(); i++) {
        mode_num = mode_arr_view3[i];
        int tt = 0;
        if (mode_arr_view3[i] == 0) {
            vertex_x_arr_view.push_back(vertex_x_arr_view3[point]);
            vertex_y_arr_view.push_back(vertex_y_arr_view3[point]);
            vertex_Red_arr_view.push_back(vertex_Red_arr_view3[point]);
            vertex_Green_arr_view.push_back(vertex_Green_arr_view3[point]);
            vertex_Blue_arr_view.push_back(vertex_Blue_arr_view3[point]);
        }
        for (int j = 0; j < mode_arr_view3[i]; j++) {
            float point1_x = vertex_x_arr_view3[point];
            float point1_y = vertex_y_arr_view3[point];
            Line_and_Point(view_boundary_x[k], view_boundary_x[k + 1], view_boundary_y[k], view_boundary_y[k + 1], point1_x, point1_y);
            string location1 = InOrOut;

            float point2_x = vertex_x_arr_view3[point + 1];
            float point2_y = vertex_y_arr_view3[point + 1];
            Line_and_Point(view_boundary_x[k], view_boundary_x[k + 1], view_boundary_y[k], view_boundary_y[k + 1], point2_x, point2_y);
            string location2 = InOrOut;
            //cout << location1 << " " << location2 << endl;
            if (location1 == "in" && location2 == "out") {
                //cout << "point1 = (" << point1_x << ", " << point1_y << ") ";
                //cout << "point2 = (" << point2_x << ", " << point2_y << ")" << endl;
                Interaction(view_boundary_x[k], view_boundary_x[k + 1], point1_x, point2_x, view_boundary_y[k], view_boundary_y[k + 1], point1_y, point2_y);
                vertex_x_arr_view.push_back(interaction_x);
                vertex_y_arr_view.push_back(interaction_y);
                vertex_Red_arr_view.push_back(vertex_Red_arr_view3[point]);
                vertex_Green_arr_view.push_back(vertex_Green_arr_view3[point]);
                vertex_Blue_arr_view.push_back(vertex_Blue_arr_view3[point]);
                //cout << " => interaction = (" << interaction_x << ", " << interaction_y << ")" << endl;
                if (tt == 0) {
                    tmp_x = interaction_x;
                    tmp_y = interaction_y;
                    tt = 1;
                }
            }
            else if (location1 == "out" && location2 == "out") {
                mode_num--;
            }
            else if (location1 == "out" && location2 == "in") {
                //cout << "point1 = (" << point1_x << ", " << point1_y << ") ";
                //cout << "point2 = (" << point2_x << ", " << point2_y << ")" << endl;
                Interaction(view_boundary_x[k], view_boundary_x[k + 1], point1_x, point2_x, view_boundary_y[k], view_boundary_y[k + 1], point1_y, point2_y);
                vertex_x_arr_view.push_back(interaction_x);
                vertex_y_arr_view.push_back(interaction_y);
                vertex_Red_arr_view.push_back(vertex_Red_arr_view3[point]);
                vertex_Green_arr_view.push_back(vertex_Green_arr_view3[point]);
                vertex_Blue_arr_view.push_back(vertex_Blue_arr_view3[point]);
                //cout << " => interaction = (" << interaction_x << ", " << interaction_y << ") ";
                if (tt == 0) {
                    tmp_x = interaction_x;
                    tmp_y = interaction_y;
                    tt = 1;
                }
                vertex_x_arr_view.push_back(point2_x);
                vertex_y_arr_view.push_back(point2_y);
                vertex_Red_arr_view.push_back(vertex_Red_arr_view3[point]);
                vertex_Green_arr_view.push_back(vertex_Green_arr_view3[point]);
                vertex_Blue_arr_view.push_back(vertex_Blue_arr_view3[point]);
                //cout << " => point2 = (" << point2_x << ", " << point2_y << ")" << endl;
                mode_num++;

            }
            else if (location1 == "in" && location2 == "in") {
                //cout << "point1 = (" << point1_x << ", " << point1_y << ") ";
                //cout << "point2 = (" << point2_x << ", " << point2_y << ")" << endl;
                vertex_x_arr_view.push_back(point2_x);
                vertex_y_arr_view.push_back(point2_y);
                vertex_Red_arr_view.push_back(vertex_Red_arr_view3[point]);
                vertex_Green_arr_view.push_back(vertex_Green_arr_view3[point]);
                vertex_Blue_arr_view.push_back(vertex_Blue_arr_view3[point]);

                //cout << " => point2 = (" << point2_x << ", " << point2_y << ")" << endl;
                if (tt == 0) {
                    tmp_x = point2_x;
                    tmp_y = point2_y;
                    tt = 1;
                }
            }

            if (j == mode_arr_view3[i] - 1) {
                vertex_x_arr_view.push_back(tmp_x);
                vertex_y_arr_view.push_back(tmp_y);
                vertex_Red_arr_view.push_back(vertex_Red_arr_view3[point]);
                vertex_Green_arr_view.push_back(vertex_Green_arr_view3[point]);
                vertex_Blue_arr_view.push_back(vertex_Blue_arr_view3[point]);
                //cout << " => origin = (" << tmp_x << ", " << tmp_y << ")" << endl;
            }
            point++;
        }
        point++;
        mode_arr_view.push_back(mode_num);
        //cout << "mode num: " << mode_num << endl;
    }


    /*
    // without clipping
    for (int i = 0; i < vertex_x_arr.size(); i++) {
        vertex_point[0] = vertex_x_arr[i];
        vertex_point[1] = vertex_y_arr[i];
        Multiply_point(view_matrix, vertex_point);
        vertex_x_arr_view.push_back(Point_matrix[0]);
        vertex_y_arr_view.push_back(Point_matrix[1]);
    }
    */
    Show();

}

void display() {

    fstream fp;
    fp.open(in_path, ios::in);
    if (!fp) {
        //cout << "Failed to open the file.\n";
        exit(0);
    }
    string InputLine, Instruction;
    int a = 0;
    
    while (getline(fp, InputLine)) {
        a = 0;

        while (a < 1) {
            a++;
            Instruction = InputLine.substr(0, InputLine.find(" "));
            if (a == 1) {
                if (Instruction == "reset") {
                    for (int i = 0; i < ROWS; i++) {
                        for (int j = 0; j < COLUMNS; j++) {
                            Transform_matrix[i][j] = Unit_matrix[i][j];
                        }
                    }
                }
                else if (Instruction == "translate") {
                    InputLine = InputLine.substr(InputLine.find(" ") + 1, InputLine.length());
                    arrange(InputLine);
                    TRANSLATE(tmp[0], tmp[1]);
                }
                else if (Instruction == "scale") {
                    InputLine = InputLine.substr(InputLine.find(" ") + 1, InputLine.length());
                    arrange(InputLine);
                    SCALE(tmp[0], tmp[1]);
                }
                else if (Instruction == "rotate") {
                    InputLine = InputLine.substr(InputLine.find(" ") + 1, InputLine.length());
                    arrange(InputLine);
                    ROTATE(tmp[0]);
                }
                else if (Instruction == "view") {
                    InputLine = InputLine.substr(InputLine.find(" ") + 1, InputLine.length());
                    arrange(InputLine);
                    VIEW(tmp[0], tmp[1], tmp[2], tmp[3], tmp[4], tmp[5], tmp[6], tmp[7]);
                    fgetc(stdin);
                }
                else if (Instruction == "clearData") {
                    CLEARDATA();
                }
                else if (Instruction == "clearScreen") {
                    glClear(GL_COLOR_BUFFER_BIT);
                    glutSwapBuffers();
                }
                else if (Instruction == "square") {
                    SQUARE();
                }
                else if (Instruction == "triangle") {
                    TRIANGLE();
                }
                else if (Instruction == "end") {
                    exit(0);
                }
                else if (Instruction == "#") {
                    break;
                }
                else if (Instruction == "") {
                    break;
                }
                else {
                    break;
                }
            }
            InputLine = InputLine.substr(InputLine.find(" ") + 1, InputLine.length());
            if (InputLine.find(" ") == -1) {
                break;
            }
        }
    }
    fp.close();
    
}

int main(int argc, char* argv[]) {

    filename = argv[1];
    in_path.append("../Data/").append(filename);
    //cout << in_path << endl;
    glutInit(&argc, argv);
    fgetc(stdin);
    glutInitDisplayMode(GLUT_RGBA | GLUT_SINGLE);
    glutInitWindowSize(800, 800);
    glutInitWindowPosition(10, 10);
    gluOrtho2D(0, 800, 0, 800);

    glutCreateWindow("Your First GLUT  Window!");
    glutDisplayFunc(display);
    gluOrtho2D(0, 800, 0, 800);
    glutMainLoop();

    return 0;
}

float random_func() {
    float randomValue = (float)rand() / RAND_MAX;
    return randomValue;
}

void randomRGB() {
    Red = random_func();
    Green = random_func();
    Blue = random_func();
    //glColor3f(Red, Green, Blue);
}

void arrange(string str) {
    int a = 0;
    string w = "";
    float w_f;
    for (auto x : str)
    {
        if (x == ' ')
        {
            if (w == "") {
            }
            else {
                w_f = strtof(w.c_str(), nullptr);
                tmp[a] = w_f;
                a++;
            }
            w = "";
        }
        else {
            w = w + x;
        }
    }
    w_f = strtof(w.c_str(), nullptr);
    tmp[a] = w_f;
}

void Inner_product(float Operator[][COLUMNS], float Transform[][COLUMNS]) {
    float Zero_matrix[3][3] = {
        {0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0}
    };
    float a = 0.0;
    for (int i = 0; i < ROWS; i++) {
        for (int j = 0; j < COLUMNS; j++) {
            a = 0.0;
            for (int k = 0; k < 3; k++) {
                a = a + Operator[i][k] * Transform[k][j];
            }
            Zero_matrix[i][j] = a;
        }
    }
    for (int i = 0; i < ROWS; i++) {
        for (int j = 0; j < COLUMNS; j++) {
            Transform_matrix[i][j] = Zero_matrix[i][j];
                  // cout << Transform_matrix[i][j];
        }
    }
}

void SCALE(float scale_x, float scale_y) {
    float scale_matrix[3][3] = {
        {scale_x,     0.0, 0.0},
        {    0.0, scale_y, 0.0},
        {    0.0,     0.0, 1.0}
    };
    Inner_product(scale_matrix, Transform_matrix);
}

void ROTATE(float angle) {
    float Degree = angle * 3.14159 / 180;
    float rotate_matrix[3][3] = {
        {cos(Degree), -sin(Degree), 0.0},
        {sin(Degree),  cos(Degree), 0.0},
        {        0.0,          0.0, 1.0}
    };
    Inner_product(rotate_matrix, Transform_matrix);
}

void TRANSLATE(float translate_x, float translate_y) {
    float translate_matrix[3][3] = {
        {1.0, 0.0, translate_x},
        {0.0, 1.0, translate_y},
        {0.0, 0.0,         1.0}
    };
    Inner_product(translate_matrix, Transform_matrix);
}

void SaveVertex(float x_point, float y_point, float SaveRed, float SaveGreen, float SaveBlue) {
    /*
    cout << "x = " << x_point << endl;
    cout << "y = " << y_point << endl;
    cout << "R = " << SaveRed << endl;
    cout << "G = " << SaveGreen << endl;
    cout << "B = " << SaveBlue << endl;
    */
    vertex_x_arr.push_back(x_point);
    vertex_y_arr.push_back(y_point);
    vertex_Red_arr.push_back(SaveRed);
    vertex_Green_arr.push_back(SaveGreen);
    vertex_Blue_arr.push_back(SaveBlue);
}

void TRIANGLE() {
    randomRGB();
    mode_arr.push_back(3);
    float Triangle_point1[3] = { 0.0,  1.0, 1.0 };
    float Triangle_point2[3] = { -1.0, -1.0, 1.0 };
    float Triangle_point3[3] = { 1.0, -1.0, 1.0 };
    float Triangle_point1_tranform[3];
    float Triangle_point2_tranform[3];
    float Triangle_point3_tranform[3];

    Multiply_point(Transform_matrix, Triangle_point1);
    for (int i = 0; i < ROWS; i++) {
        Triangle_point1_tranform[i] = Point_matrix[i];
    }
    Multiply_point(Transform_matrix, Triangle_point2);
    for (int i = 0; i < ROWS; i++) {
        Triangle_point2_tranform[i] = Point_matrix[i];
    }
    Multiply_point(Transform_matrix, Triangle_point3);
    for (int i = 0; i < ROWS; i++) {
        Triangle_point3_tranform[i] = Point_matrix[i];
    }

    SaveVertex(Triangle_point1_tranform[0], Triangle_point1_tranform[1], Red, Green, Blue);
    SaveVertex(Triangle_point2_tranform[0], Triangle_point2_tranform[1], Red, Green, Blue);
    SaveVertex(Triangle_point3_tranform[0], Triangle_point3_tranform[1], Red, Green, Blue);
    SaveVertex(Triangle_point1_tranform[0], Triangle_point1_tranform[1], Red, Green, Blue);
}

void SQUARE() {
    randomRGB();
    mode_arr.push_back(4);
    float Square_point1[3] = { 1.0,  1.0, 1.0 };
    float Square_point2[3] = { -1.0,  1.0, 1.0 };
    float Square_point3[3] = { -1.0, -1.0, 1.0 };
    float Square_point4[3] = { 1.0, -1.0, 1.0 };
    float Square_point1_tranform[3];
    float Square_point2_tranform[3];
    float Square_point3_tranform[3];
    float Square_point4_tranform[3];
    Multiply_point(Transform_matrix, Square_point1);
    for (int i = 0; i < ROWS; i++) {
        Square_point1_tranform[i] = Point_matrix[i];
    }
    Multiply_point(Transform_matrix, Square_point2);
    for (int i = 0; i < ROWS; i++) {
        Square_point2_tranform[i] = Point_matrix[i];
    }
    Multiply_point(Transform_matrix, Square_point3);
    for (int i = 0; i < ROWS; i++) {
        Square_point3_tranform[i] = Point_matrix[i];
    }
    Multiply_point(Transform_matrix, Square_point4);
    for (int i = 0; i < ROWS; i++) {
        Square_point4_tranform[i] = Point_matrix[i];
    }

    SaveVertex(Square_point1_tranform[0], Square_point1_tranform[1], Red, Green, Blue);
    SaveVertex(Square_point2_tranform[0], Square_point2_tranform[1], Red, Green, Blue);
    SaveVertex(Square_point3_tranform[0], Square_point3_tranform[1], Red, Green, Blue);
    SaveVertex(Square_point4_tranform[0], Square_point4_tranform[1], Red, Green, Blue);
    SaveVertex(Square_point1_tranform[0], Square_point1_tranform[1], Red, Green, Blue);

}

void drawALine(int x_start, int x_end, int y_start, int y_end, float R, float G, float B) {
    //cout << "drawALine" << endl;
    int dx = x_end - x_start;
    int dy = y_end - y_start;
    if (abs(dy) < abs(dx)) {
        if (x_end > x_start) {
            for (i = x_start; i <= x_end; i++) {
                glColor3f(R, G, B);
                glVertex2i(i, round(y_start + dy * (i - x_start) / dx));
                //savePixel(i, (int)round(y_start + dy * (i - x_start) / dx), R, G, B);
            }
        }
        else {
            for (i = x_start; i >= x_end; i--) {
                glColor3f(R, G, B);
                glVertex2i(i, round(y_start + dy * (i - x_start) / dx));
                //savePixel(i, (int)round(y_start + dy * (i - x_start) / dx), R, G, B);
            }
        }
    }
    else if (abs(dx) == 0 && abs(dy) == 0) {
        glColor3f(R, G, B);
        glVertex2i(x_start, y_start);
        //savePixel(x_start, y_start, , G, B);
    }
    else if (abs(dx) <= abs(dy)) {
        if (y_end > y_start) {
            for (i = y_start; i <= y_end; i++) {
                glColor3f(R, G, B);
                glVertex2i(round(x_start + dx * (i - y_start) / dy), i);
                //savePixel((int)round(x_start + dx * (i - y_start) / dy), i, R, G, B);
            }
        }
        else {
            for (i = y_start; i >= y_end; i--) {
                glColor3f(R, G, B);
                glVertex2i(round(x_start + dx * (i - y_start) / dy), i);
                //savePixel((int)round(x_start + dx * (i - y_start) / dy), i, R, G, B);
            }
        }
    }
}
