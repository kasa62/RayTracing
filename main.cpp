#include <GL/freeglut.h>
#include <math.h>
#include <stdio.h>
#include <algorithm>

#include <vector>

#define PI 3.14159265358979323846

int width = 640;
int height = 480;

const int image_width = 640;
const int image_height = 480;

struct Sphere
{
	double c[3];//球の中心
	double r;//球の半径
};

struct Camera
{
	double view[3] = { 0.0,0.0,0.0 };			//視点
	double s[3] = { 10.0,0.0,0.0 };				//スクリーンの中心		
	double u[3] = { 0.0,float(image_width) / float(image_height),0.0 };	//スクリーンのx軸の方向ベクトル
	double v[3] = { 0.0,0.0,1.0 };				//スクリーンのy軸の方向ベクトル
};

struct Light
{
	double light[3];			//光源のベクトル
	double s_color[3];	//光源の色(RGB)
	double re_color[3];		//拡散反射の色(RGB)
};

struct Vec3
{
	double v[3];
};

std::vector<Light> g_Lights;
std::vector<Sphere> g_Spheres;
std::vector<Vec3> g_Colors;

GLuint texture = 0;

void initLights() {
	g_Lights.clear();

	Light l;
	l.light[0] = 15.0, l.light[1] = 5.0, l.light[2] = 0.0;//光源のベクトル
	l.s_color[0] = 100.0, l.s_color[1] = 100.0, l.s_color[2] = 100.0;	//光源の色(RGB)
	l.re_color[0] = 0.0, l.re_color[1] = 1.0, l.re_color[2] = 0.0;

	g_Lights.push_back(l);

	l.light[0] = 15.0, l.light[1] = 0.0, l.light[2] = 10.0;//光源のベクトル
	l.s_color[0] = 100.0, l.s_color[1] = 100.0, l.s_color[2] = 100.0;	//光源の色(RGB)
	l.re_color[0] = 1.0, l.re_color[1] = 0.0, l.re_color[2] = 0.0;

	g_Lights.push_back(l);

	l.light[0] = 15.0, l.light[1] = -5.0, l.light[2] = 0.0;//光源のベクトル
	l.s_color[0] = 0.0, l.s_color[1] = 100.0, l.s_color[2] = 100.0;	//光源の色(RGB)
	l.re_color[0] = 0.0, l.re_color[1] = 0.0, l.re_color[2] = 1.0;

	g_Lights.push_back(l);

	l.light[0] = 15.0, l.light[1] = 0.0, l.light[2] = 0.0;//光源のベクトル
	l.s_color[0] = 0.0, l.s_color[1] = 100.0, l.s_color[2] = 100.0;	//光源の色(RGB)
	l.re_color[0] = 1.0, l.re_color[1] = 0.0, l.re_color[2] = 0.0;

	g_Lights.push_back(l);
}

void initSpheres()
{
	g_Spheres.clear();

	Sphere s;
	s.c[0] = 20.0; s.c[1] = 0.0; s.c[2] = 0.0;
	s.r = 1.0;

	g_Spheres.push_back(s);

	s.c[0] = 21.0; s.c[1] = 2.0; s.c[2] = 0.0;
	s.r = 1.0;

	g_Spheres.push_back(s);

	s.c[0] = 25.0; s.c[1] = 1.0; s.c[2] = 1.0;
	s.r = 1.0;

	g_Spheres.push_back(s);
}

void initColors() {
	g_Colors.clear();
	g_Colors.resize(image_width * image_height);
}

void foo()
{
	for (int i = 0; i < g_Spheres.size(); i++)
	{
		double c[3];
		c[0] = g_Spheres[i].c[0];
		c[1] = g_Spheres[i].c[1];
		c[2] = g_Spheres[i].c[2];
	}
}

//配列の足し算 (size:3)
void plus(double a[3], double b[3], double c[3]) {
	for (int i = 0;i < 3;i++) {
		c[i] = a[i] + b[i];
	}
}

//配列の引き算 (size:3) 
void sub(double a[3], double b[3], double c[3]) {
	for (int i = 0;i < 3;i++) {
		c[i] = a[i] - b[i];
	}
}

//配列の各要素をn倍する (size:3) 
void mul(double a[3], double n, double b[3]) {
	for (int i = 0;i < 3;i++) {
		b[i] = n * a[i];
	}
}

//ベクトルa,bの内積を計算 (size:3) 
double dot(double a[3], double b[3]) {
	double c = 0;
	for (int i = 0;i < 3;i++) {
		c += a[i] * b[i];
	}
	return c;
}

//球と視線の交点の判別
double discri(double view[3], double center[3], double d[3], double r) {
	double a[3];
	sub(view, center, a);
	return dot(a, d) * dot(a, d) - dot(d, d) * (dot(a, a) - r * r);
}

//交点の計算(tM>tm)
void point(double view[3], double center[3], double d[3], double r, double* tM, double* tm) {
	double a[3];
	sub(view, center, a);
	*tm = (-dot(a, d) - sqrt(discri(view, center, d, r))) / dot(d, d);
	*tM = (-dot(a, d) + sqrt(discri(view, center, d, r))) / dot(d, d);
}

//ノルムの計算
double norm(double a[3]) {
	return sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
}

//配列の正規化 (size:3)
void seikika(double a[3], double b[3]) {
	double n = norm(a);
	for (int i = 0;i < 3;i++) {
		b[i] = a[i] / n;
	}
}

//i,jの時の方向ベクトルの計算
void cal_direction(double u[3], double v[3], double s[3], double view[3], double d[3], int i, int j) {
	double nu[3], nv[3], p[3];
	mul(u, ((double)i + 0.5) * 2 / (double)image_width - 1, nu);
	mul(v, ((double)j + 0.5) * 2 / (double)image_height - 1, nv);
	plus(s, nu, p);
	plus(p, nv, p);
	sub(p, view, d);
	seikika(d, d);
}

//交点tの決定（存在しない時は-1）
double decide_t(double view[3], double center[3], double d[3], double r) {
	double tM, tm, t = -1;
	double D = discri(view, center, d, r);
	if (D >= 0) {
		//printf("D: %e\n", D);

		point(view, center, d, r, &tM, &tm);

		if (tm > 0) {
			t = tm;
		}
		else if (tM > 0) {
			t = tM;
		}
		else {
			t = -1;
		}
	}
	return t;
}

//ベクトルxの計算
void cal_x_v(double d[3], double view[3], double t, double vector[3]) {
	double tmp[3];

	mul(d, t, tmp);
	plus(view, tmp, vector);
}

//法線ベクトルの計算
void cal_normal_v(double vector[3], double center[3], double n[3]) {
	double tmp[3];

	sub(vector, center, tmp);
	seikika(tmp, n);
}

//色の計算
void cal_color(double d[3], double view[3], double center[3], const Light& l, double t, double I[3]) {
	double vector[3], n[3], tmp[3];
	cal_x_v(d, view, t, vector);//ベクトルxの計算
	cal_normal_v(vector, center, n);//法線ベクトルの計算

	sub((double*)l.light, vector, tmp);
	double a = norm(tmp);

	//sub(vector, (double*)l.light, tmp);
	for (int k = 0;k < 3;k++) {
		I[k] += l.s_color[k] * l.re_color[k] * std::max<double>(0.0, dot(tmp, n)) / (a * a * a);
	}
}

//交点と何番目の球か(t,flag)を返す
void Intersection(double *t,int *flag,double d[3],int i,int j) {

	initSpheres();
	double center[3];
	double r;
	Camera camera;

	cal_direction(camera.u, camera.v, camera.s, camera.view, d, i, j);

	*t = 10000;
	*flag = -1;
	for (int k = 0; k < g_Spheres.size(); k++) {

		center[0] = g_Spheres[k].c[0];
		center[1] = g_Spheres[k].c[1];
		center[2] = g_Spheres[k].c[2];
		r = g_Spheres[k].r;
		double tmp = decide_t(camera.view, center, d, r);
		if (tmp > 0 && *t > tmp) {
			*t = tmp;
			*flag = k;
		}
	}
	if (*flag == -1) {
		*t = -1;
	}
}

//色の決定
void Color(double t,int flag,double d[3],double I[3]) {
	Camera camera;
	if (t > 0) {
		for (int i = 0; i < g_Lights.size();i++) {
			cal_color(d, camera.view, g_Spheres[flag].c, g_Lights[i], t, I);
		}
	}
	else {
		I[0] = 0.0;
		I[1] = 0.0;
		I[2] = 0.0;
	}
}


void set_color(){
	for (int j = 0;j < image_height;j++) {
		for (int i = 0;i < image_width;i++) {

			int flag = 0;//球の番号
			double t = 10000;//交点
			double d[3];//方向ベクトル
			Intersection(&t, &flag, d, i, j);

			double I[3] = { 0.0,0.0,0.0 };
			Color(t, flag, d, I);
			g_Colors[j * image_width + i].v[0] = I[0];
			g_Colors[j * image_width + i].v[1] = I[1];
			g_Colors[j * image_width + i].v[2] = I[2];
		}
	}
}

void plot() {
	const int height_reference = 480;
	const float half_width = 1.0 * float(width) / float(height_reference);
	const float half_height = 1.0 * float(height) / float(height_reference);

	for (int j = 0;j < image_height;j++) {
		for (int i = 0;i < image_width;i++) {

			glColor3f(g_Colors[j * image_width + i].v[0], g_Colors[j * image_width + i].v[1], g_Colors[j * image_width + i].v[2]);
			glBegin(GL_TRIANGLES);

			float h = 2.0 / height;
			float x = float(i - image_width / 2) * half_height / float(height);
			float y = float(j - image_height / 2) * half_height / float(height);

			//float _x = float(width / 2) * half_height / float(height);
			//float _y = float(height / 2) * half_height / float(height);

			glVertex3f(x, y, 0.0);
			glVertex3f(x + h, y + h, 0.0);
			glVertex3f(x, y + h, 0.0);

			glVertex3f(x + h, y + h, 0.0);
			glVertex3f(x, y, 0.0);
			glVertex3f(x + h, y, 0.0);

			glEnd();
		}

	}
}

void display()
{
	glViewport(0, 0, width, height);

	glClearColor(0.5, 0.5, 0.5, 0.0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();

	const int height_reference = 480;
	const float half_width = 1.0 * float(width) / float(height_reference);
	const float half_height = 1.0 * float(height) / float(height_reference);

	glOrtho(-half_width, half_width, -half_height, half_height, -1.0, 100.0);
	//printf("hw: %e, hh: %e\n", half_width, half_height);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	//glEnable(GL_TEXTURE_2D);
	//glBindTexture(GL_TEXTURE_2D, texture);

	/*
	テクスチャの練習（四角形）

	glBegin(GL_TRIANGLES);
	glColor3f(1.0, 1.0, 1.0);
	glTexCoord2f(0.0, 0.0);
	glVertex3f(0.0, 0.0, 0.0);
	glTexCoord2f(1.0, 0.0);
	glVertex3f(0.5, 0.0, 0.0);
	glTexCoord2f(0.0, 1.0);
	glVertex3f(0.0, 0.5, 0.0);

	glTexCoord2f(0.0, 1.0);
	glVertex3f(0.0, 0.5, 0.0);
	glTexCoord2f(1.0, 0.0);
	glVertex3f(0.5, 0.0, 0.0);
	glTexCoord2f(1.0, 1.0);
	glVertex3f(0.5, 0.5, 0.0);

	glEnd();
	*/

	//glDisable(GL_TEXTURE_2D);

	plot();

	glEnable(GL_LINE_WIDTH);
	glLineWidth(3.0);

	glutSwapBuffers();
}

void resize(int w, int h)
{
	width = w;
	height = h;
}

int main(int argc, char* argv[])
{
	initLights();
	initColors();
	set_color();
	glutInit(&argc, argv);
	glutInitWindowSize(width, height);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);

	glutCreateWindow("ex2");

	glutDisplayFunc(display);
	glutReshapeFunc(resize);

	int tex_width = 640;
	int tex_height = 480;
	float* buff = (float*)malloc(sizeof(float) * tex_width * tex_height * 3);

	for (int j = 0;j < tex_height;j++) {
		for (int i = 0;i < tex_width;i++) {
			float x = (i + 0.5) / tex_width;
			float y = (j + 0.5) / tex_height;
			int pid = j * tex_width + i;
			buff[3 * pid] = x;
			buff[3 * pid + 1] = y;
			buff[3 * pid + 2] = 0.0;

		}
	}


	glGenTextures(1, &texture);
	glBindTexture(GL_TEXTURE_2D, texture);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, tex_width, tex_height, 0, GL_RGB, GL_FLOAT, buff);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);

	glutMainLoop();
	return 0;
}

//デバック用のmain関数
//int main() {
//	double a[3] = { 5.0,3.0,7.0 };
//	double b[3] = { 2.0,4.0,6.0 };
//	double c[3];
//	double r = 5.0;							//球の半径
//	double center[3] = { 20.0,0.0,0.0 };	//球の中心
//	double view[3] = { 0.0,0.0,0.0 };		//視点
//	double s[3] = { 10.0,0.0,0.0 };
//	double u[3] = { 0.0,1.0,0.0 };
//	double v[3] = { 0.0,0.0,1.0 };
//
//	double p[3], d[3], tmp[3],D;
//	int i = 0, j = 0;
//	double result;
//
//	double nu[3], nv[3];
//	mul(u, ((double)i + 0.5) * 2 / (double)width - 1, nu);
//	//printf("nu:%4.2lf;%4.2lf;%4.2lf;\n", nu[0], nu[1], nu[2]);
//	mul(v, ((double)j + 0.5) * 2 / (double)height - 1, nv);
//	//printf("nv:%.2lf;%.2lf;%.2lf;\n", nv[0], nv[1], nv[2]);
//	plus(s, nu, p);
//	//printf("p:%.2lf;%.2lf;%.2lf;\n", p[0], p[1], p[2]);
//	plus(p, nv, p);
//	//printf("p:%.2lf;%.2lf;%.2lf;\n", p[0], p[1], p[2]);
//	sub(p, view, tmp);
//	//printf("tmp:%.2lf;%.2lf;%.2lf;\n", tmp[0], tmp[1], tmp[2]);
//	seikika(tmp, d);
//	//printf("d:%.2lf;%.2lf;%.2lf;\n", d[0], d[1], d[2]);
//	D = discri(view, center, d, r);
//	//printf("D:%.2lf;\n", D);
//	result = norm(a);
//	printf("norm:%lf\n", result);
//
//
//	//c[0] = dot(a, b);
//	//printf("%4.2lf;%4.2lf;%4.2lf;\n", c[0], c[1], c[2]);
//}
