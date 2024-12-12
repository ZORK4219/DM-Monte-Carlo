#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

#define PI 3.141592653589793

const double rc=2.5;
typedef struct 
{
    double x; 
    double y;
} Particule;

double Potentiel(Particule tab_particules[],int N,double L)
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

double Energie_1part(Particule tab_particules[], Particule part_concernee, int index_concernee, int taille_tableau, double taille_boite)
{
    double rij,U=0.0;
    double U_decale=4*(pow(1/rc,12) - pow(1/rc,6)),dx,dy;
    for (int j = 0; j < taille_tableau; j++)
    {
        if (j == index_concernee)
        {
            continue;
        }
        dx=part_concernee.x-tab_particules[j].x;
        dy=part_concernee.y-tab_particules[j].y;
        dx=dx-round(dx/taille_boite)*taille_boite;
        dy=dy-round(dy/taille_boite)*taille_boite;
        rij=sqrt(dx*dx+dy*dy);
        if( rij<rc)
        {
            U += 4 *((pow(1/rij,12) - pow(1/rij,6)) )-U_decale;
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

double Metropolis(Particule tab_particules_lamda_i[], int N, double T, double L) 
{    
    // Tirer une particule au hasard
    int index = rand() % (N);
    printf("index tiré = %d\n", index);
    Particule particule_i = tab_particules_lamda_i[index];
    Particule particule_ip1 = tab_particules_lamda_i[index];
    double Dr=0.01;
    // Tirer un rayon aléatoire entre 0 et R_MAX
    double r = ((double)rand() / RAND_MAX) * Dr;

    // Tirer un angle aléatoire entre 0 et 2π
    double theta = ((double)rand() / RAND_MAX) * 2.0 * PI;

    // Calculer les nouveaux déplacements en x et y
    double dx = r * cos(theta);
    double dy = r * sin(theta);

    particule_ip1.x+=dx;
    particule_ip1.y+=dy;

    double E_lamda_i=Energie_1part(tab_particules_lamda_i, particule_i, index, N, L);
    double E_lamda_ip1=Energie_1part(tab_particules_lamda_i, particule_ip1, index, N, L);
    double delta_E = E_lamda_ip1 - E_lamda_i;

    double facteur_boltzmann=exp(-delta_E/T);
    double A=fmin(1,facteur_boltzmann);
    printf("A = %f | boltzmann = %f\n", A, facteur_boltzmann);
    double nombre = ((double)rand() / RAND_MAX);

    if(nombre<A)
    {
        tab_particules_lamda_i[index].x = particule_ip1.x;
        tab_particules_lamda_i[index].y = particule_ip1.y;
        return 1.;
    }
    else
    {
        return 0.;
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
    double L = 10.0;
    double T = 0.01;
    int N = 128;
    srand((unsigned int)time(NULL));
    Particule tab_particules_lamda_i[N];
    FILE *pos = fopen("Positions_xy.txt", "w");
    FILE *pot = fopen("Potentiel.txt", "w");
    double somme = 0;
    int nstep = 100000;
    
    configuration_cristalline(tab_particules_lamda_i, N, L);
    
    for (int step = 0; step < nstep; step++)
    {
        somme += Metropolis(tab_particules_lamda_i, N, T, L);
        conditions_periodiques(tab_particules_lamda_i, N, L);
        
        // Écrit le potentiel dans le fichier tous les 10 pas pour réduire le volume de données
        if (step % 10 == 0)
        {
            fprintf(pot, "%d %g\n", step, Potentiel(tab_particules_lamda_i, N, L));
        }
    }

    printf("Moyenne : %.2g\n", somme / nstep);
    fclose(pos);
    fclose(pot);

    // Lancer gnuplot pour tracer le potentiel
    system("gnuplot -persist -e \"plot 'Potentiel.txt' with lines title 'Potentiel'\"");

    return 0;
}