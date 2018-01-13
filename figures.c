#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <float.h>

#define min(a,b) ((a) < (b)) ? (a) : (b)
#define max(a,b) ((a) > (b)) ? (a) : (b)
#define PI 3.1415
#define RACCURACY 1000.
#define ACCURACY 0.001

enum FIGURE_TYPE { FIGURE_SEGMENT, FIGURE_ARC, FIGURE_RAY, FIGURE_CIRCLE, FIGURE_POLYGON };

struct figure{
    enum FIGURE_TYPE type;   // тип фигуры

    float Xc, Yc; // координаты центра фигуры
    float R;    // Радиус для окружности и дуги, радиус описаной окружности для многоугольника

    int n; // Число вершин многоугольника

    float a1, a2; // Углы ограничения дуги

    float phi; // Направляющий угол луча

    float k, b; // Коэффициенты прямой(для луча и отрезка)
    float A, B, C; // Коэффициенты уравнения прямой общего вида(для луча и отрезка)

    float Xs, Xe; // Рабочий интервал (начало и конец)

    struct figure * segment; // Отрезки составляющие правильный многоугольник
};

struct point{
    float X;
    float Y;
};

int find_line_circle_intersections(struct figure * fig1, struct figure * fig2, struct point * points);
int find_circle_circle_intersections(struct figure * fig1, struct figure * fig2, struct point * points);

void init_segment(struct figure * f){
    f->type = FIGURE_SEGMENT;

    float Ax, Ay, Bx, By;

    printf("\r\nThe segment is described by two points(A, B):\r\n");
    printf("Enter value Ax: ");
    scanf("%f", &Ax);
    printf("Enter value Ay: ");
    scanf("%f", &Ay);
    printf("Enter value Bx: ");
    scanf("%f", &Bx);
    printf("Enter value By: ");
    scanf("%f", &By);

    // Определение коэффициентов прямой
    f->k = (By-Ay)/(Bx-Ax);
    f->b = Ay - f->k*Ax;

    // коэффициенты для уравнения прямой общего вида
    f->A = Ay - By;
    f->B = Bx - Ax;
    f->C = Ax*By - Bx*Ay;

    // Определение рабочего интервала отрезка
    f->Xs = min(Ax, Bx);
    f->Xe = max(Ax, Bx);
}

void init_arc(struct figure * f){
    f->type = FIGURE_ARC;

    printf("\r\nThe arc of the circle is given by the coordinates of the center(Xc, Yc), radius(R) and bounding angles in radians{a2>a1 | a1,a2=[-2Pi;2Pi]}:\r\n");
    printf("Enter value Xc: ");
    scanf("%f", &(f->Xc));
    printf("Enter value Yc: ");
    scanf("%f", &(f->Yc));
    printf("Enter value R: ");
    scanf("%f", &(f->R));
    printf("Enter value a1: ");
    scanf("%f", &(f->a1));
    printf("Enter value a2: ");
    scanf("%f", &(f->a2));

    // Определение рабочего интервала дуги окружности
    f->Xs = f->Xc + min(cos(f->a1), cos(f->a2))*f->R;
    f->Xe = f->Xc + max(cos(f->a1), cos(f->a2))*f->R;
}

void init_ray(struct figure * f){
    f->type = FIGURE_RAY;

    printf("\r\nThe ray is given by coordinates of the center (Xc, Yc) and steering angle in radians{phi=[-2Pi;2Pi]}:\r\n");
    printf("Enter value Xc: ");
    scanf("%f", &f->Xc);
    printf("Enter value Yc: ");
    scanf("%f", &f->Yc);
    printf("Enter value phi: ");
    scanf("%f", &f->phi);

    // Определение рабочего интервала луча
    if((f->phi > -PI/2 && f->phi < PI/2) || (f->phi > 3/2*PI) || (f->phi < -3/2*PI)){   // Луч сонаправлен с осью OX
        f->Xs = f->Xc;
        f->Xe = FLT_MAX; // Примем за бесконечность
    } else if((f->phi > PI/2 && f->phi < 3/2*PI) || (f->phi > -3/2*PI && f->phi < -PI/2)){  // Луч противоположно направлен к оси OX
        f->Xs = -FLT_MAX; // Примем за минус бесконечность
        f->Xe = f->Xc;
    }else{  // Луч перпендикулярен оси OX
        f->Xs = f->Xe = f->Xc;
    }

    // Определение коэффициентов прямой
    f->k = tan(f->phi);
    f->b = f->Yc - f->k*f->Xc;

    // коэффициенты для уравнения прямой общего вида
    f->A = f->k;
    f->B = -1;
    f->C = f->b;
}

void init_circle(struct figure * f){
    f->type = FIGURE_CIRCLE;

    printf("\r\nThe circle is given by coordinates of the center(Xc, Yc) and radius(R):\r\n");
    printf("Enter value Xc: ");
    scanf("%f", &f->Xc);
    printf("Enter value Yc: ");
    scanf("%f", &f->Yc);
    printf("Enter value R: ");
    scanf("%f", &f->R);

    // Определение рабочего интервала окружности
    f->Xs = f->Xc - f->R;
    f->Xe = f->Xc + f->R;
}

void init_polygon(struct figure * f){
    f->type = FIGURE_POLYGON;

    float phi;

    printf("\r\nA regular polygon is given by coordinates of the center (Xc, Yc), radius of the circumscribed circle (R), number of vertices (n), and angular coordinate of the first vertex in radians{phi=[0;2Pi]}:\r\n");
    printf("Enter value Xc: ");
    scanf("%f", &(f->Xc));
    printf("Enter value Yc: ");
    scanf("%f", &(f->Yc));
    printf("Enter value R: ");
    scanf("%f", &(f->R));
    printf("Enter value n: ");
    scanf("%d", &(f->n));
    printf("Enter value phi: ");
    scanf("%f", &phi);

    // За правильный многоугольник принимаем n отрезков прямых
    f->segment = calloc(f->n, sizeof(struct figure));

    float Bx = f->Xc + f->R*cos(phi),  // Координаты первой вершины
        By = f->Yc + f->R*sin(phi),
        Ax,Ay;

    for(int i = 0; i < f->n; i++){
        Ax = Bx;
        Ay = By;
        Bx = f->Xc + f->R*cos(phi + 2*PI*(i+1)/f->n);
        By = f->Yc + f->R*sin(phi + 2*PI*(i+1)/f->n);

        // Определение коэффициентов прямой
        f->segment[i].k = (By-Ay)/(Bx-Ax);
        f->segment[i].b = Ay - f->segment[i].k*Ax;

        // коэффициенты для уравнения прямой общего вида
        f->segment[i].A = Ay - By;
        f->segment[i].B = Bx - Ax;
        f->segment[i].C = Ax*By - Bx*Ay;

        // Определение рабочего интервала отрезка
        f->segment[i].Xs = min(Ax, Bx);
        f->segment[i].Xe = max(Ax, Bx);
    }
}

// Поиск общего рабочега интервала двух фигур
int total_interval(struct figure * fig1, struct figure * fig2, float * Xs, float * Xe)
{
    *Xs = max(fig1->Xs, fig2->Xs);
    *Xe = min(fig1->Xe, fig2->Xe);

    if (*Xs > *Xe)
        return 1;

    return 0;
}

// Удаление дубликатов из списка
int removeCopies(struct point ** points, int len)
{
    for (int i = len-1; i > 0; i--){
		struct point p1; 
		struct point p2;
		p1.X = (*points)[i].X;
		p1.Y = (*points)[i].Y;
		p2.X = (*points)[i + 1].X;
		p2.Y = (*points)[i + 1].Y;
        if (((*points)[i].X - (*points)[(i+1)%len].X) < ACCURACY && ((*points)[i].Y - (*points)[(i+1)%len].Y) < ACCURACY){		
            for (int j = i; j < len-1; j++){
				p1.X = (*points)[i].X;
				p1.Y = (*points)[i].Y;
				p2.X = (*points)[i + 1].X;
				p2.Y = (*points)[i + 1].Y;
                (*points)[i].X = (*points)[i+1].X;
                (*points)[i].Y = (*points)[i+1].Y;
            }
            len--;
        }
    }

    (*points) = realloc((*points), len * sizeof(struct point));

    return len;
}

float angle_from_sin_cos(float sinVal, float cosVal)
{
    float x = acos(cosVal);
    float y = asin(sinVal);

    float a;
    if (y >= 0) // I и II четверти
        a = x;
    else if (x >= 0 && y < 0) // IV четверть
        a = y;
    else if (x < 0 && y < 0) // III четверть
        a = x + PI/4;

    return a;
}

int find_line_line_intersections(struct figure * fig1, struct figure * fig2, struct point * points)
{
    float Xs, Xe;

    int res = total_interval(fig1, fig2, &Xs, &Xe);
    if (res)    // Если рабочие интервалы фигур не пересекаются
        return 0;

    if (fabs(fig1->k - fig2->k) <= ACCURACY)
        if (fabs(fig1->b - fig2->b) <= ACCURACY){
            if (fabs(Xs - Xe) <= ACCURACY){  // Частный случай для двух противоположно направленных лучей с общей стартовой точкой
                points[0].X = fig1->Xc;
                points[0].Y = fig1->Yc;

                return 1; // Одна точка пересечения
            }
            return -1;  // Прямые совпадают и имеют общий интервал, точек пересечения бесконечно много
        }
        else
            return 0;   // Заданные прямые параллельны, точек пересечения нет

    // Поиск точки пересечения
    float x = round((fig2->b - fig1->b) / (fig1->k - fig2->k) * RACCURACY) / RACCURACY;

    if ((Xs - x) > ACCURACY || (x - Xe) > ACCURACY) // Точка пересечения не входит в рабочий интервал
        return 0;

    float y = round((fig1->k * x + fig1->b) * RACCURACY) / RACCURACY;

    points[0].X = x;
    points[0].Y = y;

    return 1; // Найдена одна точка пересечения
}

int find_line_arc_intersections(struct figure * fig1, struct figure * fig2, struct point * points)
{
    int len;

    len = find_line_circle_intersections(fig1, fig2, points);

    // Проверки на попадание точer на дугу
    float R = fig2->R;
    float Xc = fig2->Xc,
          Yc = fig2->Yc;

    float a1 = fig2->a1,
          a2 = fig2->a2;

    for (int i = len-1; i >= 0; i--){
        float x = points[i].X;
        float y = points[i].Y;

        float cosA = (x-Xc)/R;
        float sinA = (y-Yc)/R;

        float angle = angle_from_sin_cos(sinA, cosA);

        if ((angle - a1 > ACCURACY && a2 - angle > ACCURACY) || (fabs(angle-a1) < ACCURACY) || (fabs(angle-a2) < ACCURACY))
            continue;

        float a;
        if (angle < -ACCURACY){ // Угол меньше нуля
            a = angle + 2*PI;
            if ((a - a1 > ACCURACY && a2 - a > ACCURACY) || (fabs(a-a1) < ACCURACY) || (fabs(a-a2) < ACCURACY))
                continue;
        }else if(a > ACCURACY){ // Угол больше нуля
            a = angle - 2*PI;
            if ((a - a1 > ACCURACY && a2 - a > ACCURACY) || (fabs(a-a1) < ACCURACY) || (fabs(a-a2) < ACCURACY))
                continue;
        }else{ // Угол лежит в окрестности нуля
            a = angle + 2*PI;
            if ((a - a1 > ACCURACY && a2 - a > ACCURACY) || (fabs(a-a1) < ACCURACY) || (fabs(a-a2) < ACCURACY))
                continue;

            a = angle - 2*PI;
            if ((a - a1 > ACCURACY && a2 - a > ACCURACY) || (fabs(a-a1) < ACCURACY) || (fabs(a-a2) < ACCURACY))
                continue;
        }

        for (int j = i; j < len-1; j++){
            points[i].X = points[i+1].X;
            points[i].Y = points[i+1].Y;
        }
        len--;
    }

    return len;
}

int find_line_circle_intersections(struct figure * fig1, struct figure * fig2, struct point * points)
{
    int len;

    float r = fig2->R,
          a = fig1->A,
          b = fig1->B,
          c = fig1->C;
    float xc = fig2->Xc,
          yc = fig2->Yc;

    // Переходим к системе координат с центром окружности в (0;0)
    c = (c - a*xc - b*yc);

    float x0 = -a*c/(a*a+b*b),  y0 = -b*c/(a*a+b*b);
    if (c*c > r*r*(a*a+b*b) + ACCURACY)
        len = 0;
    else if (abs (c*c - r*r*(a*a+b*b)) < ACCURACY) {
        // Возвращаем точкy в исходную систему координат
        points[0].X = x0 + xc;
        points[0].Y = y0 + yc;

        len = 1;
    }
    else {
    	float d = r*r - c*c/(a*a+b*b);
    	float mult = sqrt (d / (a*a+b*b));
    	float ax,ay,bx,by;
    	ax = x0 + b * mult;
    	bx = x0 - b * mult;
    	ay = y0 - a * mult;
    	by = y0 + a * mult;

        // Возвращаем точки в исходную систему координат
        points[0].X = ax + xc;
        points[0].Y = ay + yc;
        points[1].X = bx + xc;
        points[1].Y = by + yc;

        len = 2;
    }

    // Проверка точек пересечения на вхождение в рабочий интервал
    for (int i = len-1; i >= 0; i--){
        float x = points[i].X;
        if ((fig1->Xs - x) > ACCURACY || (x - fig1->Xe) > ACCURACY){ // Точка пересечения не входит в рабочий интервал
            for (int j = i; j < len-1; j++){
                points[i].X = points[i+1].X;
                points[i].Y = points[i+1].Y;
            }
            len--;
        }
    }

    return len;
}

int find_line_poly_intersections(struct figure * fig1, struct figure * fig2, struct point ** points)
{
    int len = 0;

    int res;
    *points = malloc(2*sizeof(struct point));
    for(int i = 0; i < fig2->n; i++){
        res = find_line_line_intersections(fig1, &fig2->segment[i], &((*points)[len]));
        if (res == -1){
            free(*points);
            return -1;
        }
        if(res > 0){
            len += res;
            *points = realloc(*points, sizeof(struct point) * (len+res));
        }
    }

    len = removeCopies(points, len);

    return len;
}

int find_arc_arc_intersections(struct figure * fig1, struct figure * fig2, struct point * points)
{
    int len;

    len = find_circle_circle_intersections(fig1, fig2, points);

    // Проверка на принадлежность точки первой дуге окружности
    float R = fig1->R;
    float Xc = fig1->Xc,
          Yc = fig1->Yc;

    float a1 = fig1->a1,
          a2 = fig1->a2;

    // Проверки на попадание точer на дугу
    for (int i = len-1; i >= 0; i--){
        float x = points[i].X;
        float y = points[i].Y;

        float cosA = (x-Xc)/R;
        float sinA = (y-Yc)/R;

        float angle = angle_from_sin_cos(sinA, cosA);

        if ((angle - a1 > ACCURACY && a2 - angle > ACCURACY) || (fabs(angle-a1) < ACCURACY) || (fabs(angle-a2) < ACCURACY))
            continue;

        float a;
        if (angle < -ACCURACY){ // Угол меньше нуля
            a = angle + 2*PI;
            if ((a - a1 > ACCURACY && a2 - a > ACCURACY) || (fabs(a-a1) < ACCURACY) || (fabs(a-a2) < ACCURACY))
                continue;
        }else if(a > ACCURACY){ // Угол больше нуля
            a = angle - 2*PI;
            if ((a - a1 > ACCURACY && a2 - a > ACCURACY) || (fabs(a-a1) < ACCURACY) || (fabs(a-a2) < ACCURACY))
                continue;
        }else{ // Угол лежит в окрестности нуля
            a = angle + 2*PI;
            if ((a - a1 > ACCURACY && a2 - a > ACCURACY) || (fabs(a-a1) < ACCURACY) || (fabs(a-a2) < ACCURACY))
                continue;

            a = angle - 2*PI;
            if ((a - a1 > ACCURACY && a2 - a > ACCURACY) || (fabs(a-a1) < ACCURACY) || (fabs(a-a2) < ACCURACY))
                continue;
        }

        for (int j = i; j < len-1; j++){
            points[i].X = points[i+1].X;
            points[i].Y = points[i+1].Y;
        }
        len--;
    }


    // Проверка на принадлежность точки второй дуге окружности
    R = fig2->R;
    Xc = fig2->Xc;
    Yc = fig2->Yc;

    a1 = fig2->a1;
    a2 = fig2->a2;

    for (int i = len-1; i >= 0; i--){
        float x = points[i].X;
        float y = points[i].Y;

        float cosA = (x-Xc)/R;
        float sinA = (y-Yc)/R;

        float angle = angle_from_sin_cos(sinA, cosA);

        // Проверки на попадание точки на дугу
        if ((angle - a1 > ACCURACY && a2 - angle > ACCURACY) || (fabs(angle-a1) < ACCURACY) || (fabs(angle-a2) < ACCURACY))
            continue;

        float a;
        if (angle < -ACCURACY){ // Угол меньше нуля
            a = angle + 2*PI;
            if ((a - a1 > ACCURACY && a2 - a > ACCURACY) || (fabs(a-a1) < ACCURACY) || (fabs(a-a2) < ACCURACY))
                continue;
        }else if(a > ACCURACY){ // Угол больше нуля
            a = angle - 2*PI;
            if ((a - a1 > ACCURACY && a2 - a > ACCURACY) || (fabs(a-a1) < ACCURACY) || (fabs(a-a2) < ACCURACY))
                continue;
        }else{ // Угол лежит в окрестности нуля
            a = angle + 2*PI;
            if ((a - a1 > ACCURACY && a2 - a > ACCURACY) || (fabs(a-a1) < ACCURACY) || (fabs(a-a2) < ACCURACY))
                continue;

            a = angle - 2*PI;
            if ((a - a1 > ACCURACY && a2 - a > ACCURACY) || (fabs(a-a1) < ACCURACY) || (fabs(a-a2) < ACCURACY))
                continue;
        }

        for (int j = i; j < len-1; j++){
            points[i].X = points[i+1].X;
            points[i].Y = points[i+1].Y;
        }
        len--;
    }

    return len;
}

int find_arc_circle_intersections(struct figure * fig1, struct figure * fig2, struct point * points)
{
    int len;

    len = find_circle_circle_intersections(fig1, fig2, points);

    float R = fig1->R;
    float Xc = fig1->Xc,
          Yc = fig1->Yc;

    float a1 = fig1->a1,
          a2 = fig1->a2;

    for (int i = len-1; i >= 0; i--){
        float x = points[i].X;
        float y = points[i].Y;

        float cosA = (x-Xc)/R;
        float sinA = (y-Yc)/R;

        float angle = angle_from_sin_cos(sinA, cosA);

        // Проверки на попадание точки на дугу
        if ((angle - a1 > ACCURACY && a2 - angle > ACCURACY) || (fabs(angle-a1) < ACCURACY) || (fabs(angle-a2) < ACCURACY))
            continue;

        float a;
        if (angle < -ACCURACY){ // Угол меньше нуля
            a = angle + 2*PI;
            if ((a - a1 > ACCURACY && a2 - a > ACCURACY) || (fabs(a-a1) < ACCURACY) || (fabs(a-a2) < ACCURACY))
                continue;
        }else if(a > ACCURACY){ // Угол больше нуля
            a = angle - 2*PI;
            if ((a - a1 > ACCURACY && a2 - a > ACCURACY) || (fabs(a-a1) < ACCURACY) || (fabs(a-a2) < ACCURACY))
                continue;
        }else{ // Угол лежит в окрестности нуля
            a = angle + 2*PI;
            if ((a - a1 > ACCURACY && a2 - a > ACCURACY) || (fabs(a-a1) < ACCURACY) || (fabs(a-a2) < ACCURACY))
                continue;

            a = angle - 2*PI;
            if ((a - a1 > ACCURACY && a2 - a > ACCURACY) || (fabs(a-a1) < ACCURACY) || (fabs(a-a2) < ACCURACY))
                continue;
        }

        for (int j = i; j < len-1; j++){
            points[i].X = points[i+1].X;
            points[i].Y = points[i+1].Y;
        }
        len--;
    }

    return len;
}

int find_arc_poly_intersections(struct figure * fig1, struct figure * fig2, struct point ** points)
{
    int len = 0;

    int res;
    *points = malloc(2 * sizeof(struct point));
    for(int i = 0; i < fig2->n; i++){
        res = find_line_arc_intersections(&fig2->segment[i], fig1, &((*points)[len]));
        if (res == -1){
            free(*points);
            return -1;
        }
        if(res > 0){
            len += res;
            *points = realloc(*points, sizeof(struct point) * (len+res));
        }
    }

    len = removeCopies(points, len);

    return len;
}

int find_circle_circle_intersections(struct figure * fig1, struct figure * fig2, struct point * points)
{
    int len;

    float x1 = fig1->Xc,
          x2 = fig2->Xc,
          y1 = fig1->Yc,
          y2 = fig2->Yc;

    float r = fig1->R,
          R = fig2->R,
          d = sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));

    if (d < ACCURACY){
        if (R - r < ACCURACY)
            return -1;  // Окружности совпадают, число точек пересечения бесконечно

        return 0;   // Окружности концентрические, точек пересечения нет
    }

    float ax = x2-x1,
          ay = y2-y1;

    float p = (r + R + d)/2;

    float h = 2*sqrt(p*(p-r)*(p-R)*(p-d))/r;

    float d1 = sqrt(r*r - h*h);

    float mx = ax * d1/d,
          my = ay * d1/d;

    float Mx = x1 + mx,
          My = y1 + my;

    if (r+R - d > ACCURACY){
        float nx = mx / d1,
              ny = my / d1;

        float b1nx = -ny,
              b1ny = nx,
              b2nx = ny,
              b2ny = -nx;

        float b1x = b1nx * h,
              b1y = b1ny * h,
              b2x = b2nx * h,
              b2y = b2ny * h;

        float C1x = Mx + b1x,
              C1y = My + b1y,
              C2x = Mx + b2x,
              C2y = My + b2y;

        points[0].X = C1x;
        points[0].Y = C1y;
        points[1].X = C2x;
        points[1].Y = C2y;

        len = 2;
    }else if(fabs(r+R - d) < ACCURACY){
        points[0].X = Mx;
        points[0].Y = My;

        len = 1;
    }else{
        len = 0;
    }

    return len;
}

int find_circle_poly_intersections(struct figure * fig1, struct figure * fig2, struct point ** points)
{
    int len = 0;

    int res;
    *points = malloc(2 * sizeof(struct point));
    for(int i = 0; i < fig2->n; i++){
        res = find_line_circle_intersections(&fig2->segment[i], fig1, &((*points)[len]));
        if (res == -1){
            free(*points);
            return -1;
        }
        if(res > 0){
            len += res;
            *points = realloc(*points, sizeof(struct point) * (len+res));
		}
		
    }

    len = removeCopies(points, len);

    return len;
}

int find_poly_poly_intersections(struct figure * fig1, struct figure * fig2, struct point ** points)
{
    int len = 0;

    int res;
    *points = malloc(sizeof(struct point));
    for(int i = 0; i < fig1->n; i++){
        for (int j = 0; j < fig2->n; j++){
            res = find_line_line_intersections(&fig1->segment[i], &fig2->segment[j], &((*points)[len]));
            if (res == -1){
                free(*points);
                return -1;
            }
            if(res == 1){
                len++;
                *points = realloc(*points, sizeof(struct point) * (len+1));
            }
        }
    }

    len = removeCopies(points, len);

    return len;
}


int find_intersections(struct figure * fig1, struct figure * fig2, struct point ** points)
{
    int len = 0;
    switch(fig1->type){
        case FIGURE_SEGMENT:
        case FIGURE_RAY:
             switch(fig2->type){
                 case FIGURE_SEGMENT:
                 case FIGURE_RAY:
                     *points = malloc(sizeof(struct point));
                     len = find_line_line_intersections(fig1, fig2, *points);
                     if (len == -1)
                         free(*points);
                     break;

                 case FIGURE_ARC:
                     *points = malloc(2 * sizeof(struct point));

                     len = find_line_arc_intersections(fig1, fig2, *points);

                     if (len == 1)
                         *points = realloc(*points, len * sizeof(struct point));
                     else if (len == 0)
                         free(*points);
                     break;

                 case FIGURE_CIRCLE:
                     *points = malloc(2 * sizeof(struct point));

                     len = find_line_circle_intersections(fig1, fig2, *points);

                     if (len == 1)
                         *points = realloc(*points, len * sizeof(struct point));
                     else if (len == 0)
                         free(*points);
                     break;

                 case FIGURE_POLYGON:
                     len = find_line_poly_intersections(fig1, fig2, points);
                     break;
             }
            break;

        case FIGURE_ARC:
             switch(fig2->type){
                 case FIGURE_SEGMENT:
                 case FIGURE_RAY:
                     len = find_line_arc_intersections(fig2, fig1, *points);
                     break;

                 case FIGURE_ARC:
                     *points = malloc(2 * sizeof(struct point));

                     len = find_arc_arc_intersections(fig1, fig2, *points);

                     if (len == 1)
                         *points = realloc(*points, len * sizeof(struct point));
                     else if (len == 0)
                         free(*points);
                     break;

                 case FIGURE_CIRCLE:
                     *points = malloc(2 * sizeof(struct point));

                     len = find_arc_circle_intersections(fig1, fig2, *points);

                     if (len == 1)
                         *points = realloc(*points, len * sizeof(struct point));
                     else if (len == 0)
                         free(*points);
                     break;

                 case FIGURE_POLYGON:
                     len = find_arc_poly_intersections(fig1, fig2, points);
                     break;
             }
            break;

        case FIGURE_CIRCLE:
             switch(fig2->type){
                 case FIGURE_SEGMENT:
                 case FIGURE_RAY:
                     *points = malloc(2 * sizeof(struct point));

                     len = find_line_circle_intersections(fig2, fig1, *points);

                     if (len == 1)
                         *points = realloc(*points, len * sizeof(struct point));
                     else if (len == 0)
                         free(*points);
                     break;

                 case FIGURE_ARC:
                     *points = malloc(2 * sizeof(struct point));

                     len = find_arc_circle_intersections(fig2, fig1, *points);

                     if (len == 1)
                         *points = realloc(*points, len * sizeof(struct point));
                     else if (len == 0)
                         free(*points);
                     break;

                 case FIGURE_CIRCLE:
                     *points = malloc(2 * sizeof(struct point));

                     len = find_circle_circle_intersections(fig1, fig2, *points);

                     if (len == 1)
                         *points = realloc(*points, len * sizeof(struct point));
                     else if (len == 0)
                         free(*points);
                     break;

                 case FIGURE_POLYGON:
                     len = find_circle_poly_intersections(fig1, fig2, points);
                     break;
             }
            break;

        case FIGURE_POLYGON:
             switch(fig2->type){
                 case FIGURE_SEGMENT:
                 case FIGURE_RAY:
                     len = find_line_poly_intersections(fig2, fig1, points);
                     break;

                 case FIGURE_ARC:
                     len = find_arc_poly_intersections(fig2, fig1, points);
                     break;

                 case FIGURE_CIRCLE:
                     len = find_circle_poly_intersections(fig2, fig1, points);
                     break;

                 case FIGURE_POLYGON:
                     len = find_poly_poly_intersections(fig1, fig2, points);
                     break;
             }
            break;
    }
    return len;
}

void print_menu()
{
    printf("(1) Segment\r\n");
    printf("(2) Arc\r\n");
    printf("(3) Ray\r\n");
    printf("(4) Circle\r\n");
    printf("(5) Polygon\r\n");
    printf("(0) Exit\r\n");
}

int main()
{
    // Переменная для выбора пункта меню пользователем
    char choice = 0;
    // Две фигуры
    struct figure fig1, fig2;

    // Выбор первой фигуры
    while (choice < '0' || choice > '5'){ // Пока пользователь не выберет пункт меню
        printf("\r\nSelect first figure:\r\n");
        print_menu();

        scanf("%c", &choice);
    }

    switch(choice){
    case '1': // Отрезок
        init_segment(&fig1);
        break;
    case '2': // Дуга
        init_arc(&fig1);
        break;
    case '3': // Луч
        init_ray(&fig1);
        break;
    case '4': // Окружность
        init_circle(&fig1);
        break;
    case '5': // Правильный многоугольник
        init_polygon(&fig1);
        break;
    case '0': // Выход
        exit(0);
        break;
    }

    // Выбор второй фигуры
    scanf("%c", &choice); // Кушаем enter, оставшийся с прошлого ввода
    while (choice < '0' || choice > '5'){ // Пока пользователь не выберет пункт меню
        printf("\r\nSelect second figure:\r\n");
        print_menu();

        scanf("%c", &choice);
    }

    switch(choice){
    case '1': // Отрезок
        init_segment(&fig2);
        break;
    case '2': // Дуга
        init_arc(&fig2);
        break;
    case '3': // Луч
        init_ray(&fig2);
        break;
    case '4': // Окружность
        init_circle(&fig2);
        break;
    case '5': // Правильный многоугольник
        init_polygon(&fig2);
        break;
    case '0': // Выход
        exit(0);
        break;
    }

    // Приступаем к поиску точек пересечения
    struct point * points;
    int len; // Количество найденных точек пересечения

    len = find_intersections(&fig1, &fig2, &points);

    // Вывод точек пересечения
    if (len == -1)
        printf("\r\nThe intersection points are infinitely many\r\n");
    else if (len == 0)
        printf("\r\nNo intersection points\r\n");
    else{
        printf("\r\nFound %d points of intersection:\r\n", len);
        for (int i = 0; i < len; i++)
            printf("%d) %.3f : %.3f\r\n", i+1, points[i].X, points[i].Y);
    }

    return 0;
}
