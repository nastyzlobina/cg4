#include <stdio.h>
#include <math.h>
#include <conio.h>
#include <stdlib.h>

#include "model.h"
#include "tga.h"
#define _CRT_SECURE_NO_WARNINGS
void swap(int *a, int *b);
int iabs(int a);

/*
* Using Digital Differential Analyzer algorihm
* to draw interval connecting (x0, y0) with (x1, y1)
* on image using color
*/
void line (tgaImage *image, 
           int x0, int y0,
           int x1, int y1,
           tgaColor color);

void sortCoordsByY(int coords[3][2], double zcoords[3]) {
    int curMinInd = 0, curMaxInd = 0;
    for (int i = 1; i < 3; i++) {
        if (coords[i][1] < coords[curMinInd][1]) {
            curMinInd = i;
        }
        if (coords[i][1] > coords[curMaxInd][1]) {
            curMaxInd = i;
        }
    }
    if (curMinInd == 0 && curMaxInd == 0) {
        curMaxInd = 2;
    }
    if (curMinInd != 0) {
        int coordsTmp[2];
        coordsTmp[0] = coords[0][0];
        coordsTmp[1] = coords[0][1];
        coords[0][0] = coords[curMinInd][0];
        coords[0][1] = coords[curMinInd][1];
        coords[curMinInd][0] = coordsTmp[0];
        coords[curMinInd][1] = coordsTmp[1];
        double ztmp;
        ztmp = zcoords[0];
        zcoords[0] = zcoords[curMinInd];
        zcoords[curMinInd] = ztmp;
        if (curMaxInd == 0) {
            curMaxInd = curMinInd;
        }
    }
    if (curMaxInd != 0) {
        int coordsTmp[2];
        coordsTmp[0] = coords[2][0];
        coordsTmp[1] = coords[2][1];
        coords[2][0] = coords[curMaxInd][0];
        coords[2][1] = coords[curMaxInd][1];
        coords[curMaxInd][0] = coordsTmp[0];
        coords[curMaxInd][1] = coordsTmp[1];
        double ztmp;
        ztmp = zcoords[2];
        zcoords[2] = zcoords[curMaxInd];
        zcoords[curMaxInd] = ztmp;
    }
}

int getXOfIntersection(int y, int x0, int y0, int x1, int y1) {
    int yDiff = y1 - y0;
    if (yDiff == 0) {
        return x0;
    }
    double k = (y - y0) / (double)yDiff;
    return x0 + k*(x1 - x0);
}

int getZOfIntersection(int y, double z0, int y0, double z1, int y1) {
    int yDiff = y1 - y0;
    if (yDiff == 0) {
        return z0;
    }
    double k = (y - y0) / (double)yDiff;
    return z0 + k*(z1 - z0);
}

void triangle(int coords[3][2], double zcoords[], tgaColor color, tgaImage *image, double **zindex) {
    // ...
    // screen_coords[j][0], screen_coords[j][1] для j от 0 до 2 - это экранные координаты вершин тек. грани

    // Расположим screen_coords в порядке возрастания "y"
    sortCoordsByY(coords, zcoords);

    for (int y = coords[0][1]; y <= coords[2][1]; y++) {
        int xA = y > coords[1][1]
            ? getXOfIntersection(y, coords[2][0], coords[2][1], coords[1][0], coords[1][1])
            : getXOfIntersection(y, coords[0][0], coords[0][1], coords[1][0], coords[1][1]);
        int xB = y > coords[2][1]
            ? getXOfIntersection(y, coords[2][0], coords[2][1], coords[1][0], coords[1][1])
            : getXOfIntersection(y, coords[0][0], coords[0][1], coords[2][0], coords[2][1]);
        double zA = y > coords[1][1]
            ? getZOfIntersection(y, zcoords[2], coords[2][1], zcoords[1], coords[1][1])
            : getZOfIntersection(y, zcoords[0], coords[0][1], zcoords[1], coords[1][1]);
        double zB = y > coords[2][1]
            ? getZOfIntersection(y, zcoords[2], coords[2][1], zcoords[1], coords[1][1])
            : getZOfIntersection(y, zcoords[0], coords[0][1], zcoords[2], coords[2][1]);
        int xStart, xEnd;
        double zStart, zEnd;
        int z;
        if (xA < xB) {
            xStart = xA;
            xEnd = xB;
            zStart = zA;
            zEnd = zB;
        } else {
            xStart = xB;
            xEnd = xA;
            zStart = zB;
            zEnd = zA;
        }
        for (int x = xStart; x <= xEnd; x++) {
            if (zStart == zEnd || xStart == xEnd) {
                z = zEnd;
            } else {
                z = ((double)(x - xStart) / (xEnd - xStart)) * (zEnd - zStart) + zStart;
            }
            if (z > zindex[x][y]) {
                zindex[x][y] = z;
                tgaSetPixel(image, (unsigned int)x, (unsigned int)y, color);
            }
        }
    }
}

void getNormal_(const Vec3 planePoints[], double *normal) {
    double p1x = planePoints[0][0];
    double p1y = planePoints[0][1];
    double p1z = planePoints[0][2];
    double p2x = planePoints[1][0];
    double p2y = planePoints[1][1];
    double p2z = planePoints[1][2];
    double p3x = planePoints[2][0];
    double p3y = planePoints[2][1];
    double p3z = planePoints[2][2];
    double ax = p2x - p1x, ay = p2y - p1y, az = p2z - p1z;
    double bx = p3x - p1x, by = p3y - p1y, bz = p3z - p1z;
    normal[0] = ay * bz - az * by;
    normal[1] = az * bx - ax * bz;
    normal[2] = ax * by - ay * bx;
}

void rotate(Vec3 *point, double angleDegree) {
    /*
    double cx0 = 0, cy0 = 0, cz0 = 0;
    double del = 2;
    double lx = 
    */
    double x = (*point)[0];
    double y = (*point)[1];
    double z = (*point)[2];

    // 180 - pi
    // 30 - ?

    double angle = 3.141592653 * angleDegree / 180;

    (*point)[0] = x * cos(angle) - z * sin(angle);// -z + 0.25;
    (*point)[1] = y;// y - 0.25;
    (*point)[2] = z * cos(angle) + x * sin(angle);//  x;
}

void meshgrid(tgaImage *image, Model *model) {
    double **zindex = malloc(sizeof(void*) * image->width);
    for (int i = 0; i < image->width; i++) {
        zindex[i] = malloc(sizeof(double) * image->height);
    }
    for (int i = 0; i < image->width; i++) {
        for (int j = 0; j < image->height; j++) {
            zindex[i][j] = -10.e9;
        }
    }
    int i, j;
    // nface - число граней
    // faces[i] - i-ая грань
    // Грань задаётся 3мя точками
    // Каждая точка - это 3 координаты (x; y; z)
    // Поэтому faces[i] является массивом из 9 чисел
    //   faces[i][0] - номер вершины первой точки
    //   faces[i][1] - номер текстурной вершины первой точки
    //   faces[i][2] - номер нормали первой точки
    //   faces[i][3] - номер вершины второй точки
    //   faces[i][4] - номер текстурной вершины второй точки
    //   faces[i][5] - номер нормали второй точки
    //   faces[i][6] - номер вершины третьей точки
    //   faces[i][7] - номер текстурной вершины третьей точки
    //   faces[i][8] - номер нормали третьей точки
    for (i = 0; i < model->nface; ++i) {
        int screen_coords[3][2]; // Переводим в экранные координаты
        Vec3 facePoints[3];
        double zcoords[3];
        double N[3];
        for (j = 0; j < 3; ++j) {
            Vec3 v;
            memcpy(&v, &(model->vertices[model->faces[i][3 * j]]), sizeof(Vec3));
            rotate(&v, 30);
            memcpy(&(facePoints[j]), &v, sizeof(Vec3));
            screen_coords[j][0] = (v[0] + 1) * image->width / 2;
            screen_coords[j][1] = (1 - v[1]) * image->height / 2;
            zcoords[j] = v[2];
        }
        getNormal_(facePoints, N);
        // ld значит light direction
        double ld[] = { 0., 0., -1. };
        // (ax; ay; az) = координаты нормали
        // (bx; by; bz) = координаты направления распространения света
        // cos(teta) = (ax*bx + ay*by + az*bz) /
        //             (sqrt(ax^2+ay^2+az^2)*sqrt(bx^2+by^2+bz^2))
        double I = (N[0]*ld[0] + N[1]*ld[1] + N[2]*ld[2])
            / (sqrt(N[0]*N[0] + N[1]*N[1] + N[2]*N[2]) * sqrt(ld[0]*ld[0] + ld[1]*ld[1] + ld[2]*ld[2]));
        if (I > 0) {
            I = 0;
        } else {
            I = -I;
            tgaColor color = tgaRGB((int)(I * 255), (int)(I * 255), (int)(I * 255));
            triangle(screen_coords, zcoords, color, image, zindex);
        }
    }
    for (int i = 0; i < image->width; i++) {
        free(zindex[i]);
    }
    free(zindex);
}

int main()
{
    srand((unsigned)time(NULL));
	Model *model = loadFromObj("C:\\Users\\nasty\\source\\repos\\Project1\\Project1\\obj\\cat.obj");
  
	int height = 800;
    int width = 800;
    tgaImage * image = tgaNewImage(height, width, RGB);

	meshgrid(image, model);
	tgaSaveToFile(image, "out.tga");
    int i;

    tgaFreeImage(image);    
	freeModel(model);
    _getch();
    return 0;
}



void line(tgaImage *image, int x0, int y0, int x1, int y1, tgaColor color) {
	char steep = 0;
	if (abs(x0 - x1) < abs(y0 - y1)) {
		swap(&x0, &y0);
		swap(&x1, &y1);
		steep = 1;
	}
	if (x0 > x1) {
		swap(&x0, &x1);
		swap(&y0, &y1);
	}

	int errorAccumulation = 0;
	int deltaError = 2 * abs(y1 - y0);
	int y = y0;
	for (int x = x0; x <= x1; x++) {
		if (steep == 1) {
			tgaSetPixel(image, (unsigned int)y, (unsigned int)x, color);
		}
		else {
			tgaSetPixel(image, (unsigned int)x, (unsigned int)y, color);
		}
		
		errorAccumulation += deltaError;

		if (errorAccumulation > (x1 - x0)) {
			y += (y1 > y0 ? 1 : -1);
			errorAccumulation -= 2 * (x1 - x0);
		}
	}
}

void swap(int *a, int *b) {
    int t = *a;
    *a = *b;
    *b = t;
}

int iabs(int a) {
    return (a >= 0) ? a : -a;
}
