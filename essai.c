#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

const double rc=2.5;
typedef struct 
{
    double x; 
    double y;
} Particule;

double Energie(Particule tab_particules[],int N)
{
    double rij,U=0.0;
    double U_decale=4*(pow(1/rc,12) - pow(1/rc,6));

    for (int i = 0; i < N; i++)
    {
        for (int j = i+1; j < N; j++)
        {
            rij=sqrt((tab_particules[i].x-tab_particules[j].x)*(tab_particules[i].x-tab_particules[j].x)\
            +(tab_particules[i].y-tab_particules[j].y)*(tab_particules[i].y-tab_particules[j].y));
            if( rij<rc)
            {
                U += 4 *((pow(1/rij,12) - pow(1/rij,6)) )-U_decale;
            }
            else
            {
                U+=0;
            }        
        }
    }
    return U;
}

int Tirage_changement(Particule tab_particules[], int N) 
{
    // Tirer une particule au hasard
    int index = rand() % N;
    double R_MAX=3.0;
    Particule particule = tab_particules[index];

    // Tirer un rayon aléatoire entre 0 et R_MAX
    double r = ((double)rand() / RAND_MAX) * R_MAX;

    // Tirer un angle aléatoire entre 0 et 2π
    double theta = ((double)rand() / RAND_MAX) * 2.0 * M_PI;

    // Calculer les nouveaux déplacements en x et y
    double dx = r * cos(theta);
    double dy = r * sin(theta);

    tab_particules[index].x+=dx;
    tab_particules[index].y+=dy;
    return index;
}

int Metropolis(Particule tab_particules_lamda[],Particule tab_particules_lamda_prime[], int N, double T) 
{    
    double E_lamda=Energie(tab_particules_lamda,N);
    double E_lamda_prime=Energie(tab_particules_lamda_prime,N);
    double facteur_boltzmann=exp(-(E_lamda_prime-E_lamda)/T);
    double A=fmin(1,facteur_boltzmann);
    printf("%g\n",A);
    double nombre = rand();

    if(nombre<A)
    {
        for (int i=0; i<N;i++)
        {
        tab_particules_lamda[i].x=tab_particules_lamda_prime[i].x;
        tab_particules_lamda[i].y=tab_particules_lamda_prime[i].y;
        }
    }
}

int main() 
{
    double L=45.0;
    double T=1;
    int N=100;
    srand((unsigned int)time(NULL));
    Particule tab_particules_lamda[N];
    Particule tab_particules_lamda_prime[N];
    FILE *pos = fopen("Positions_xy.txt", "w");
    int nx = (int)sqrt(N);
    int k = 0;
    //((float)rand() / RAND_MAX)
    for (int i = 0; i < nx && k < N; i++) {
        for (int j = 0; j < nx && k < N; j++) {
            tab_particules_lamda[k].x = 4.25 * (j+1) - L / 2.;
            tab_particules_lamda[k].y = 4.25 * (i+1) - L / 2.;
            tab_particules_lamda_prime[k].x=tab_particules_lamda[k].x;
            tab_particules_lamda_prime[k].y=tab_particules_lamda[k].y;
            //fprintf(pos, "%.2f %.2f\n", tab_particules[k].x, tab_particules[k].y);
            k++;
        }
    }
    for (int step=0; step<10;step++)
    {
        Tirage_changement(tab_particules_lamda_prime,N);
        Metropolis(tab_particules_lamda,tab_particules_lamda_prime,N,T);
    }

    for (int i=0; i<N;i++)
    {
        fprintf(pos, "%.2f %.2f\n", tab_particules_lamda[i].x, tab_particules_lamda[i].y);
    }

    fclose(pos);
    
   /* for (int j = 0; j <N; j++) 
   {
    printf("%f %f \n",tab_particules_lamda_prime[j].x,tab_particules_lamda_prime[j].y);
            
    } */
   /*  int index=Tirage_changement(tab_particules_lamda_prime,N);
    printf("%f %f \n",tab_particules_lamda[index].x,tab_particules_lamda[index].y);
    printf("%f %f \n",tab_particules_lamda_prime[index].x,tab_particules_lamda_prime[index].y);
 */
   
}