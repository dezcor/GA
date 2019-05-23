// SGA.cpp: define el punto de entrada de la aplicación de consola.
//

#include "GA.h"
#include <opencv2/opencv.hpp>
#include <string>
#include <iostream>

bool compar(GA::Individuo& i,GA::Individuo& j)
{ 
	return i.GetObjetivo()<j.GetObjetivo();
}
struct Points
{
    float Xi;
    float Yi;
    Points()
    {
        Xi = 0;
        Yi = 0;
    }
    Points(float x, float y)
    {
        Xi = x;
        Yi = y;
    }
};
class DetecteCircle : public GA::Poblacion
{
  public:
    DetecteCircle();
    DetecteCircle(const char *path,bool numeroReal);
    ~DetecteCircle();
    //Funcion Virtual
    void FuncionObjetivo(GA::Individuo &individuo);
    //No virtuales
    void PutResul(const string &path);
    void SetImage(const char *path);
    void PutImageallCircle(const string & path);

  private:
    void drawcircle(int x0, int y0, int radius, int gray);
    float isCircle(float x0, float y0, float radius);
    double getX0(float xi, float yi, float xj, float yj, float xk, float yk);
    double getY0(float xi, float yi, float xj, float yj, float xk, float yk);
    double getR(float x, float y, float x0, float y0);
    cv::Mat img;
    vector<Points> points;
};

int main(int argc, char *argv[])
{
    string in, out,is;
    bool isReal=true;
    if (argc >= 3)
    {
        in = argv[1];
        printf("%s\n", in.c_str());
        int x = in.find(".");
        out = in.substr(0, x);
        //out += "Result.bmp";
        is = argv[2];
        if(is == "true")
            isReal = true;
        else
        {
            isReal = false;
            is = "false";
        }
        std::cout << is << std::endl;
    }
    else
        return 1;
    DetecteCircle poblacion(in.c_str(),isReal);
    poblacion.resizePoblacion(50u);
    poblacion.setProbabilidadCruza(0.5);
    poblacion.setProbabilidadMuta(0.1);
    poblacion.Run(1000, 40);
    //poblacion.MostrarPoblacion();
    printf("Mejor Individuo: %d\n", poblacion.GetidMejor());
    poblacion.MostrarIndividuo(poblacion.GetidMejor());
    poblacion.PutResul(out + is + "Result.bmp");
    poblacion.PutImageallCircle(out + is + "allCircle.png");
    cv::waitKey(0);
    return 0;
}

DetecteCircle::DetecteCircle() : Poblacion()
{
    //img = NULL;
    MaxError = 100;
}

DetecteCircle::DetecteCircle(const char *path,bool numeroReal) : Poblacion()
{
    MaxError = 100;
    cv::Mat image,gris,blur,borde;
	gris = cv::imread(path, cv::IMREAD_REDUCED_COLOR_4);
    cv::imwrite("fabi.bmp",gris);
    gris = cv::imread(path, cv::IMREAD_GRAYSCALE);
    bool no=false;
    for (int i = 0; i < gris.rows; i++)
	{
		for (int j = 0; j < gris.cols; j++)
		{
			if (gris.at<uchar>(i, j)<225 && gris.at<uchar>(i, j)>0)
			{
                //Suavizamos la imagen (desenfoque Gaussiano) utilizando la función Gaussian Blur
                cv::GaussianBlur(gris,blur,cv::Size(7,7),1.5,1.5);
                //Utilizamos la función Canny para detectar bordes
                cv::Canny(blur,borde,0,30,3);
                img =borde;
                no=true;
                break;
			}
		}
        if(no) break;
	}
    if(!no) img=gris;
	for (int i = 0; i < borde.rows; i++)
	{
		for (int j = 0; j < borde.cols; j++)
		{
			if (img.at<uchar>(i, j)>225)
			{
				img.at<uchar>(i,j)=0;
			}
			else
			img.at<uchar>(i,j)=255;
		}
	}
    for (int i = 0; i < img.rows; i++)
    {
        for (int j = 0; j < img.cols; j++)
        {
            if (img.at<uchar>(i, j) == 0)
            {
                points.push_back(Points(j, i));
            }
        }
    }
    isNumberReal = numeroReal;
    if(isNumberReal)
    {
        this->NumeroBitGet.resize(3, 16);
        Limites.resize(2);
        Limites[0] = 0;
        Limites[1] = points.size()-1;
    }
    else
    {
        this->NumeroBitGet.resize(3, (unsigned int)ceil(log2(points.size())));
        Limites.resize(2);
        Limites[0] = 0;
        Limites[1] = points.size()-1;
    }
    NumeroDeGenes = 3;
    SizeCromosoma = 3 * NumeroBitGet[0];
}

void DetecteCircle::SetImage(const char *path)
{
    img = cv::imread(path, cv::IMREAD_GRAYSCALE);
    points.clear();
    for (int i = 0; i < img.rows; i++)
    {
        for (int j = 0; j < img.cols; j++)
        {
            if (img.at<uchar>(i, j) == 0)
            {
                points.push_back(Points(j, i));
            }
        }
    }
    this->NumeroBitGet.resize(3, (unsigned int)ceil(log2(points.size())));
    NumeroDeGenes = 3;
    SizeCromosoma = 3 * NumeroBitGet[0];
}

DetecteCircle::~DetecteCircle()
{
}

void DetecteCircle::FuncionObjetivo(GA::Individuo &individuo)
{
    this->DecodificarToInt(individuo);
    int i = individuo.GetValor(0);
    int j = individuo.GetValor(1);
    int k = individuo.GetValor(2);
    if (i >= points.size())
    {
        individuo.SetObjetivo(0.1);
        return;
    }
    if (j >= points.size())
    {
        individuo.SetObjetivo(0.1);
        return;
    }
    if (k >= points.size())
    {
        individuo.SetObjetivo(0.1);
        return;
    }
    float x0 = getX0(points[i].Xi, points[i].Yi, points[j].Xi, points[j].Yi, points[k].Xi, points[k].Yi);
    float y0 = getY0(points[i].Xi, points[i].Yi, points[j].Xi, points[j].Yi, points[k].Xi, points[k].Yi);
    float r = getR(points[i].Xi, points[i].Yi, x0, y0);
    if (r < 3)
    {
        individuo.SetObjetivo(0.1);
        return;
    }
    int num;
    individuo.SetObjetivo(isCircle(x0, y0, r) * 100);
}

double DetecteCircle::getX0(float xi, float yi, float xj, float yj, float xk, float yk)
{
    float a = (xj * xj + yj * yj - (xi * xi + yi * yi));
    float b = 2 * (yj - yi);
    float c = (xk * xk + yk * yk - (xi * xi + yi * yi));
    float d = 2 * (yk - yi);
    float dt = a * d - b * c;
    float mt = (4 * ((xj - xi) * (yk - yi) - (xk - xi) * (yj - yi)));
    if (mt == 0)
        return 0;
    return dt / mt;
}

double DetecteCircle::getY0(float xi, float yi, float xj, float yj, float xk, float yk)
{
    float b = (xj * xj + yj * yj - (xi * xi + yi * yi));
    float a = 2 * (xj - xi);
    float d = (xk * xk + yk * yk - (xi * xi + yi * yi));
    float c = 2 * (xk - xi);
    float dt = a * d - b * c;
    float mt = (4 * ((xj - xi) * (yk - yi) - (xk - xi) * (yj - yi)));
    if (mt == 0)
        return 0;
    return dt / mt;
}

double DetecteCircle::getR(float x, float y, float x0, float y0)
{
    return sqrt((x - x0) * (x - x0) + (y - y0) * (y - y0));
}

float DetecteCircle::isCircle(float x0, float y0, float r)
{
    if (img.empty())
        return 0;
    int count = 0;
    float pi = 3.14159265358979323846;
    int i;
    int xi, yi;
    int Perimetro = 2 * r * pi;
    float dxy = 2.0 * pi / Perimetro;
    for (i = 0; i < Perimetro; i++)
    {
        xi = (r * cos(i * dxy) + x0);
        yi = (r * sin(i * dxy) + y0);
        if (xi < 0 || xi >= img.cols)
            continue;
        if (yi < 0 || yi >= img.rows)
            continue;
        count += img.at<uchar>(yi, xi) == 0;
    }

    return count / (float)Perimetro;
}

void DetecteCircle::drawcircle(int x0, int y0, int radius, int gray)
{
    int x = radius - 1;
    int count = 0, Xi = 0, Yi = 0;
    int y = 0;
    int dx = 1;
    int dy = 1;
    int err = dx - (radius << 1);

    while (x >= y)
    {
        Xi = x0 + x;
        Yi = y0 + y;
        if ((Xi >= 0 && Xi < img.cols) && (Yi >= 0 && Yi < img.rows))
            img.at<uchar>(Yi, Xi) = gray;
        Xi = x0 + y;
        Yi = y0 + x;
        if ((Xi >= 0 && Xi < img.cols) && (Yi >= 0 && Yi < img.rows))
            img.at<uchar>(Yi, Xi) = gray;
        Xi = x0 - y;
        Yi = y0 + x;
        if ((Xi >= 0 && Xi < img.cols) && (Yi >= 0 && Yi < img.rows))
            img.at<uchar>(Yi, Xi) = gray;
        Xi = x0 - x;
        Yi = y0 + y;
        if ((Xi >= 0 && Xi < img.cols) && (Yi >= 0 && Yi < img.rows))
            img.at<uchar>(Yi, Xi) = gray;
        Xi = x0 - x;
        Yi = y0 - y;
        if ((Xi >= 0 && Xi < img.cols) && (Yi >= 0 && Yi < img.rows))
            img.at<uchar>(Yi, Xi) = gray;
        Xi = x0 - y;
        Yi = y0 - x;
        if ((Xi >= 0 && Xi < img.cols) && (Yi >= 0 && Yi < img.rows))
            img.at<uchar>(Yi, Xi) = gray;
        Xi = x0 + y;
        Yi = y0 - x;
        if ((Xi >= 0 && Xi < img.cols) && (Yi >= 0 && Yi < img.rows))
            img.at<uchar>(Yi, Xi) = gray;
        Xi = x0 + x;
        Yi = y0 - y;
        if ((Xi >= 0 && Xi < img.cols) && (Yi >= 0 && Yi < img.rows))
            img.at<uchar>(Yi, Xi) = gray;

        if (err <= 0)
        {
            y++;
            err += dy;
            dy += 2;
        }

        if (err > 0)
        {
            x--;
            dx += 2;
            err += dx - (radius << 1);
        }
    }
}

void DetecteCircle::PutResul(const string &path)
{
    if (img.empty())
        return;
    int i = individuos[idMejor].GetValor(0);
    int j = individuos[idMejor].GetValor(1);
    int k = individuos[idMejor].GetValor(2);
    float x0 = getX0(points[i].Xi, points[i].Yi, points[j].Xi, points[j].Yi, points[k].Xi, points[k].Yi);
    float y0 = getY0(points[i].Xi, points[i].Yi, points[j].Xi, points[j].Yi, points[k].Xi, points[k].Yi);
    float r = getR(points[i].Xi, points[i].Yi, x0, y0);
    //drawcircle(x0, y0, r + 1, 255 / 2);
    printf("X = %f\nY = %f\nRadio = %f\n", x0, y0, r);
    cv::Scalar color(255,0,0);
    cv::Mat imagen;
    cv::cvtColor(img,imagen,cv::COLOR_GRAY2RGB);
    cv::circle(imagen, cv::Point(x0, y0), r, color);
    cv::namedWindow("Display window", cv::WINDOW_AUTOSIZE); // Create a window for display.
    cv::imshow("Display window", imagen);                      // Show our image inside it.
    cv::imwrite(path, imagen);
}


void DetecteCircle::PutImageallCircle(const string & path)
{
     if (img.empty())
        return;
        float x0,y0,r;
    printf("\n");
    cv::Mat imagen;
    cv::cvtColor(img,imagen,cv::COLOR_GRAY2RGB);
    for(size_t i1 = 0; i1 < SizePoblacion; i1++)
    {
        int i = individuos[i1].GetValor(0);
        int j = individuos[i1].GetValor(1);
        int k = individuos[i1].GetValor(2);
        x0 = getX0(points[i].Xi, points[i].Yi, points[j].Xi, points[j].Yi, points[k].Xi, points[k].Yi);
        y0 = getY0(points[i].Xi, points[i].Yi, points[j].Xi, points[j].Yi, points[k].Xi, points[k].Yi);
        r = getR(points[i].Xi, points[i].Yi, x0, y0);
        printf("\nX = %f Y = %f Radio = %f Objetivo = %f\n", x0, y0, r,individuos[i1].GetObjetivo());
        if(r<=0)continue;
        if( isnan(x0) || isnan(y0) || isnan(r)) continue;
        cv::circle(imagen, cv::Point(x0, y0), r, cv::Scalar(rand()/(float)RAND_MAX*255,(float)rand()/RAND_MAX*255,(float)rand()/RAND_MAX*255));
    }
    cv::namedWindow("Display all", cv::WINDOW_AUTOSIZE); // Create a window for display.
    cv::imshow("Display all", imagen);                      // Show our image inside it.
    cv::imwrite(path, imagen);
}