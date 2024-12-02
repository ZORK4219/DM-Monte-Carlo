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


void configuration_cristalline(Particule tab_particules[], int N, double L) {
    int nx = (int)ceil(sqrt(N));  // Nombre de particules par ligne
    int ny = (int)ceil((double)N / nx);  // Nombre de particules par colonne
    double spacing_x = L / nx;  // Espacement horizontal
    double spacing_y = L / ny;  // Espacement vertical
    int k = 0;

    for (int i = 0; i < ny && k < N; i++) {
        for (int j = 0; j < nx && k < N; j++) {
            // Position des particules, centrées dans la boîte
            tab_particules[k].x = spacing_x * j - L / 2. + spacing_x / 2.;
            tab_particules[k].y = spacing_y * i - L / 2. + spacing_y / 2.;
            k++;
        }
    }
}



void configuration_aleatoire(Particule tab_particules[],int N, double L)
{
    for (int k = 0; k < N; k++) 
    {
    // Génère des coordonnées aléatoires dans la plage [-L/2, L/2]
    tab_particules[k].x = ((double)rand() / RAND_MAX) * L - L / 2.;
    tab_particules[k].y = ((double)rand() / RAND_MAX) * L - L / 2.;
    }
}
void copy_configuration(Particule tab_particules_lamda[],Particule tab_particules_lamda_prime[],int N)
{
    for (int k = 0; k < N; k++) 
    {
    // Copie des positions dans tab_particules_lamda_prime
    tab_particules_lamda_prime[k].x = tab_particules_lamda[k].x;
    tab_particules_lamda_prime[k].y = tab_particules_lamda[k].y;
    }
}
double Energie(Particule tab_particules[],int N,double L)
{
    double rij,U=0.0;
    double U_decale=4*(pow(1/rc,12) - pow(1/rc,6)),dx,dy;

    for (int i = 0; i < N; i++)
    {
        for (int j = i+1; j < N; j++)
        {
            dx=tab_particules[i].x-tab_particules[j].x;
            dy=tab_particules[i].y-tab_particules[j].y;
            dx=dx-round(dx/L)*L;
            dy=dy-round(dy/L)*L;
            rij=sqrt(dx*dx+dy*dy);
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
void conditions_periodiques(Particule tab_particules[],int N,double L)
{
    for (int i = 0; i < N; i++) 
    {
        if(tab_particules[i].x>L/2) tab_particules[i].x-=L;
        else if(tab_particules[i].x<-L/2) tab_particules[i].x+=L;
        if(tab_particules[i].y>L/2) tab_particules[i].y-=L;
        else if(tab_particules[i].y<-L/2) tab_particules[i].y+=L;
        
    }
}
int Tirage_changement(Particule tab_particules[], int N) 
{
    // Tirer une particule au hasard
    int index = rand() % N;
    double Dr=0.2;
    Particule particule = tab_particules[index];

    // Tirer un rayon aléatoire entre 0 et R_MAX
    double r = ((double)rand() / RAND_MAX) * Dr;

    // Tirer un angle aléatoire entre 0 et 2π
    double theta = ((double)rand() / RAND_MAX) * 2.0 * M_PI;

    // Calculer les nouveaux déplacements en x et y
    double dx = r * cos(theta);
    double dy = r * sin(theta);

    tab_particules[index].x+=dx;
    tab_particules[index].y+=dy;
    
    return index;
}

void Metropolis(Particule tab_particules_lamda[],Particule tab_particules_lamda_prime[], int N, double T,double L) 
{    
    double E_lamda=Energie(tab_particules_lamda,N,L);
    double E_lamda_prime=Energie(tab_particules_lamda_prime,N,L);
    double facteur_boltzmann=exp(-(E_lamda_prime-E_lamda)/T);
    double A=fmin(1,facteur_boltzmann);
    printf("%g\n",A);
    double nombre = ((double)rand() / RAND_MAX);

    if(nombre<A)
    {
        for (int i=0; i<N;i++)
        {
        tab_particules_lamda[i].x=tab_particules_lamda_prime[i].x;
        tab_particules_lamda[i].y=tab_particules_lamda_prime[i].y;
        }
    }
  
}
void write_xy (FILE *fp,Particule tab_particules[],int N) {
  int i;
  fprintf(fp, "%d\n", N);
  fprintf(fp, "\n");
  for (i=0; i<N;i++)
  {
    fprintf(fp, "A %g %g\n",  tab_particules[i].x,tab_particules[i].y);
  }
  
}
int main() 
{
    double L=45.0;
    double T=0.1;
    int N=100;
    srand((unsigned int)time(NULL));
    Particule tab_particules_lamda[N];
    Particule tab_particules_lamda_prime[N];
    FILE *pos = fopen("Positions_xy.txt", "w");
    
    configuration_cristalline(tab_particules_lamda,N,L);
    for (int step=0; step<5000;step++)
    {
        copy_configuration(tab_particules_lamda,tab_particules_lamda_prime,N);
        Tirage_changement(tab_particules_lamda_prime,N);
        Metropolis(tab_particules_lamda,tab_particules_lamda_prime,N,T,L);
        conditions_periodiques(tab_particules_lamda,N,L);
        write_xy(pos,tab_particules_lamda,N);
    }
    fclose(pos);
   
}