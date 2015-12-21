#include <stdio.h>
#include <stdlib.h>
#include <math.h>


double 	lamda 	= 	400,
        rho 	=	8930,
        cp		=	385,
        constantK	=	1,
        temperatureInfini	=	293,
        temperatureInitial	=	323,
        deltaTime,
        deltaRadius,
        PI	=	3.14159265358979323846,
        deltaAngular,
        commonCoeff,
        radius,
        calculationTime;
// The support functions:
void    getInput();
void    printResult();
double 	coeffA();
double 	coeffB();
double 	coeffC();
double 	coeffD();
double 	coeffE();
double 	temperatureCalcul();
double	iRadius();
double	jAngular();
double 	nTime();

int 	numberIntervalTimes,
        numberNodeRadius,
        numberNodeAngular;

int N;
int M;

double 	** temperature, ** temperatureFlag, A, B, C, D, E;

int main()
{
    int n,i,j;
    double temp1, temp2, temp3, temp4, temp5;

    lamda = 400;
    rho = 8930;
    cp = 385;
    constantK = 100;
    radius = 0.1;
    temperatureInfini = 293;
    temperatureInitial = 353;

    calculationTime = 100;
    numberIntervalTimes = 1000;
    numberNodeRadius = 5;
    numberNodeAngular = 7;
    //getInput();

    deltaAngular 	=	PI / (numberNodeAngular - 1);
    deltaRadius 	=	radius / (numberNodeRadius - 1);
    deltaTime 	    =	calculationTime / (numberIntervalTimes - 1);

    commonCoeff		=	(lamda*deltaTime)/(rho*cp);
    N = numberNodeRadius-1;
    M = numberNodeAngular-1;
    /*	Init temperature array */
    temperature = (double **) malloc(numberNodeRadius*sizeof(double *));
    temperatureFlag     =   (double **) malloc(numberNodeRadius*sizeof(double *));
    for (i = 0; i < numberNodeRadius; i++) {
        temperature[i] = (double *)malloc(numberNodeAngular*sizeof(double));
        temperatureFlag[i] = (double *)malloc(numberNodeAngular*sizeof(double));
    }

    FILE *datFile;
    datFile=fopen("result/data.txt","w");
    FILE *timeFile;
    timeFile=fopen("result/time.txt","w");

    FILE *radiusFile;
    radiusFile=fopen("result/radius.txt","w");

    FILE *angularFile;
    angularFile=fopen("result/angular.txt","w");

    FILE *inputFile;
    inputFile=fopen("result/input.txt","w");

    // Start calculation
    for (n = 0; n< numberIntervalTimes; n++) {
        //Calculation at any time interval
        for(i=0; i<numberNodeRadius; i++){
            for(j=0; j<numberNodeAngular; j++){
                if(n == 0){
                    temperature[i][j]	=	temperatureInitial;
                } else {
                    if(i == 0) {
                        //At center of the sphere
                        if (j == 0){
                            temperature[0][0]	=	((3*commonCoeff*2)/pow(deltaRadius, 2))*(temperatureFlag[1][0] - temperatureFlag[0][0]) + temperatureFlag[0][0];
                        } else {
                            temperature[0][j] 	=	temperatureFlag[0][0];
                        }
                    } else if(i < N){
                        //At any point between center and boundary
                        if (j == 0){
                            //At any points on the right radius
                            temp1 	=	temperatureFlag[i-1][0];
                            //apply condition
                            temp2	=	temperatureFlag[i][1];
                            temp3	=	temperatureFlag[i][0];
                            temp4	=	temperatureFlag[i][1];
                            temp5	=	temperatureFlag[i+1][0];

                            temperature[i][0]	=	temperatureCalcul(i, j, temp1, temp2, temp3, temp4, temp5);
                        } else if(j == M){
                            //At any points on the left radius
                            temp1 	=	temperatureFlag[i-1][M];
                            temp2	=	temperatureFlag[i][M-1];
                            temp3	=	temperatureFlag[i][M];
                            //apply condition
                            temp4	=	temperatureFlag[i][M-1];
                            temp5	=	temperatureFlag[i+1][M];

                            temperature[i][M]	=	temperatureCalcul(i, j, temp1, temp2, temp3, temp4, temp5);
                        } else {
                            //At any other points inside
                            temp1 	=	temperatureFlag[i-1][j];
                            temp2	=	temperatureFlag[i][j-1];
                            temp3	=	temperatureFlag[i][j];
                            temp4	=	temperatureFlag[i][j+1];
                            temp5	=	temperatureFlag[i+1][j];

                            temperature[i][j]	=	temperatureCalcul(i, j, temp1, temp2, temp3, temp4, temp5);

                        }
                    } else {
                        // At the boundary of the sphere
                        double flag = ((-2)*constantK*sin(jAngular(j))*deltaRadius)/lamda;
                        if (j == 0){
                            //At the extreme limit point on the right
                            temp1 	=	temperatureFlag[N-1][0];
                            //apply condition
                            temp2	=	temperatureFlag[N][1];
                            temp3	=	temperatureFlag[N][0];
                            temp4	=	temperatureFlag[N][1];
                            //apply condition
                            temp5	=	flag*(temperatureFlag[N][j] - temperatureInfini) + temperatureFlag[N-1][j];

                            temperature[N][0]	=	temperatureCalcul(i, j, temp1, temp2, temp3, temp4, temp5);
                        } else if(j == M){
                            //At the extreme limit point on the left
                            temp1 	=	temperatureFlag[N-1][M];
                            temp2	=	temperatureFlag[N][M-1];
                            temp3	=	temperatureFlag[N][M];
                            //Apply condition
                            temp4	=	temperatureFlag[N][M-1];
                            //Apply condition
                            //temp5	=	(-2)*constantK*sin(jAngular(j))*deltaRadius/lamda*(temperatureFlag[N][j] - temperatureInfini) + temperatureFlag[N-1][j];
                            temp5	=	flag*(temperatureFlag[N][j] - temperatureInfini) + temperatureFlag[N-1][j];

                            temperature[N][M]	=	temperatureCalcul(i, j, temp1, temp2, temp3, temp4, temp5);
                        } else {
                            //At any other points on the boundary
                            temp1 	=	temperatureFlag[N-1][j];
                            temp2	=	temperatureFlag[N][j-1];
                            temp3	=	temperatureFlag[N][j];
                            temp4	=	temperatureFlag[N][j+1];
                            //Apply condition
                            //temp5	=	(-2)*constantK*sin(jAngular(j))*deltaRadius/lamda*(temperatureFlag[N][j] - temperatureInfini) + temperatureFlag[N-1][j];
                            temp5	=	flag*(temperatureFlag[N][j] - temperatureInfini) + temperatureFlag[N-1][j];

                            temperature[N][j]	=	temperatureCalcul(i, j, temp1, temp2, temp3, temp4, temp5);
                        }
                    }
                }
                fprintf(datFile,"%lf\t",temperature[i][j]);
                if(j == M){
                    fprintf(datFile,"\n");
                }
            }
        }
        for (i=0; i<numberNodeRadius; i++) {
            for(j=0; j<numberNodeAngular; j++){
                temperatureFlag[i][j]   =   temperature[i][j];
            }
        }
    }
    for (n = 0; n< numberIntervalTimes; n++) {
        for(i=0; i<numberNodeRadius; i++){
            fprintf(timeFile,"%d\t%lf\n",n,nTime(n));
            fprintf(radiusFile,"%lf\n",iRadius(i));
        }
    }
    for(j=0; j<numberNodeAngular; j++){
        fprintf(angularFile,"%lf\t",jAngular(j));
    }

    fprintf( inputFile, "Calculation time \t %lf \t [s]\n", calculationTime);
    fprintf( inputFile, "Number of interval times  \t %d\n", numberIntervalTimes);
    fprintf( inputFile, "Time interval \t %lf \t[s]\n", deltaTime);
    fprintf( inputFile, "Number of interval radius  \t %d\n", numberNodeRadius);
    fprintf( inputFile, "Radius interval \t %lf [m]\n", deltaRadius);
    fprintf( inputFile, "Number of interval angular  \t %d\n", numberNodeAngular);
    fprintf( inputFile, "Angular interval \t %lf \t[rad]\n", deltaAngular);

    fprintf( inputFile, "Thermal conductivity \t %lf \t[J/K]\n", lamda);
    fprintf( inputFile, "Density \t %lf \t[kg/m3]\n", rho);
    fprintf( inputFile, "Heat specific pressure constant \t %lf \t[J/(kg.K)]\n", cp);
    fprintf( inputFile, "Constant K \t %lf\n", constantK);
    fprintf( inputFile, "Radius \t %lf \t[m]\n", radius);
    fprintf( inputFile, "Temperature infinitive  \t %lf \t[K]\n", temperatureInfini);
    fprintf( inputFile, "Temperature initial \t %lf \t[K]\n", temperatureInitial);

    fclose(inputFile);
    fclose(datFile);
    fclose(timeFile);
    fclose(radiusFile);
    fclose(angularFile);
    printResult();
    printf("\nEnd of program!\n");
    system("pause");
    return 0;
}
void    getInput(){
    double flag;
    printf("\n|--------------Inputs of problems------------------|\n");
    printf( " Thermal conductivity [J/K]:");
    scanf("%lf", &flag);
    lamda = flag;

    printf( "\n Density [kg/m3]:");
    scanf("%lf", &flag);
    rho = flag;

    printf( "\n Heat specific pressure constant [J/(kg.K)]:");
    scanf("%lf", &flag);
    cp = flag;

    printf( "\n Constant K :");
    scanf("%lf", &flag);
    constantK = flag;

    printf( "\n Radius [m]:");
    scanf("%lf", &flag);
    radius = flag;

    printf( "\n Temperature infinitive [K]:");
    scanf("%lf", &flag);
    temperatureInfini = flag;

    printf( "\n Temperature initial [K]:");
    scanf("%lf", &flag);
    temperatureInitial = flag;

    printf("\n|--------------------------------------------------|\n");
    printf("\n|--------------------------------------------------|\n");
    printf("\n|--------Inputs of finite method difference--------|\n");

    printf( "\n Calculation time [s] :");
    scanf("%lf", &flag);
    calculationTime = flag;

    int d;
    printf( "\n Number of interval times :");
    scanf("%d", &d);
    numberIntervalTimes = d;

    printf( "\n Number of nodes by radius :");
    scanf("%d", &d);
    numberNodeRadius = d;

    printf( "\n Number of nodes by angular :");
    scanf("%d", &d);
    numberNodeAngular = d;

    printf("\n|--------------------------------------------------|\n");
    printf("\n|--------------------------------------------------|\n");
}

void printResult(){
    printf("\n|--------------------Calculation-------------------|\n");

    printf( " Thermal conductivity: \t %lf [J/K]\n", lamda);
    printf( "\n Density: \t %lf [kg/m3]\n", rho);
    printf( "\n Heat specific pressure constant: \t %lf [J/(kg.K)]\n", cp);
    printf( "\n Constant K: \t %lf\n", constantK);
    printf( "\n Radius: \t %lf [m]\n", radius);
    printf( "\n Temperature infinitive : \t %lf [K]\n", temperatureInfini);
    printf( "\n Temperature initial: \t %lf [K]\n", temperatureInitial);
    printf( "\n Calculation time: \t %lf [s]\n", calculationTime);
    printf( "\n Number of interval times : \t %d\n", numberIntervalTimes);
    printf( "\n Time interval: \t %lf [s]\n", deltaTime);
    printf( "\n Number of interval radius : \t %d\n", numberNodeRadius);
    printf( "\n Radius interval: \t %lf [m]\n", deltaRadius);
    printf( "\n Number of interval angular : \t %d\n", numberNodeAngular);
    printf( "\n Angular interval: \t %lf [rad]\n", deltaAngular);
}
/*
Function to calculate the coeff A
*/
double temperatureCalcul(int i, int j, double temp1, double temp2, double temp3, double temp4, double temp5)
{
    double result;
    result 	=	coeffA(i, j) * temp1;
    result	+=	coeffB(i, j) * temp2;
    result	+=	coeffC(i, j) * temp3;
    result	+=	coeffD(i, j) * temp4;
    result	+=	coeffE(i, j) * temp5;

    return result;
}
/*
Function to calculate the time from beginning.
*/
double nTime(int n)
{
    return n*deltaTime;
}
/*
Function to calculate radius at node [i]
*/
double iRadius(int i)
{
    //return (i-1)*deltaRadius;
    return i*deltaRadius;
}
/*
Function to calculate angular at node [j]
*/
double jAngular(int j)
{
    //return (j-1)*deltaAngular;
    return j*deltaAngular;
}
/*
Function to calculate the coeff A
corresponding with T[n][i-1][j]
*/
double coeffA(int i, int j)
{
    double 	radius_i 	=	iRadius(i),
            x = 1 / pow(deltaRadius, 2),
            y =	1 / (radius_i * deltaRadius);
    return commonCoeff*(x - y);
}
/*
Function to calculate the coeff B
corresponding with T[n][i][j-1]
*/
double coeffB(int i, int j)
{
    double 	angular_j	=	jAngular(j),
            radius_i	=	iRadius(i),
            x = 1 / pow((radius_i * deltaAngular), 2),
            y =	cos(angular_j) / (2*pow(radius_i, 2)*sin(angular_j)*deltaAngular);
    if (j == 0 || j == M) {
        return commonCoeff*2*x;
    } else {
        return commonCoeff*(x - y);
    }

}
/*
Function to calculate the coeff C
corresponding with T[n][i][j]
*/
double coeffC(int i, int j)
{
    double 	radius_i	=	iRadius(i),
            x = 1 / pow(deltaRadius, 2),
            y = 1 / pow((radius_i * deltaAngular), 2);
    if (j == 0 || j == M) {
        return 	1 - 2*commonCoeff*(x + 2*y);
    } else {
        return 	1 - 2*commonCoeff*(x + y);
    }

}
/*
Function to calculate the coeff D
corresponding with T[n][i][j+1]
*/
double coeffD(int i, int j)
{
    double 	angular_j	=	jAngular(j),
            radius_i	=	iRadius(i),
            x = 1 / pow((radius_i * deltaAngular), 2),
            y =	cos(angular_j) / (2*pow(radius_i, 2)*sin(angular_j)*deltaAngular);
    //return commonCoeff*(x + y);
    if (j == 0 || j == M) {
        return commonCoeff*2*x;
    } else {
        return commonCoeff*(x + y);
    }
}
/*
Function to calculate the coeff E
corresponding with T[n][i+1][j]
*/
double coeffE(int i, int j)
{
    double 	radius_i 	=	iRadius(i),
            x = 1 / pow(deltaRadius, 2),
            y =	1 / (radius_i * deltaRadius);

    return commonCoeff*(x + y);
}
