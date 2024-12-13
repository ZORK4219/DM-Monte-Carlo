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

double Energie_totale(Particule tab_particules[], int taille_tableau, double taille_boite)
{
    double rij,U=0.0;
    double U_decale=4*(pow(1/rc,12) - pow(1/rc,6)),dx,dy;
    for (int i = 0; i < taille_tableau; i++)
    {
        for (int j = i + 1; j < taille_tableau; j++)
        {
            dx=tab_particules[i].x-tab_particules[j].x;
            dy=tab_particules[i].y-tab_particules[j].y;
            dx=dx-round(dx/taille_boite)*taille_boite;
            dy=dy-round(dy/taille_boite)*taille_boite;
            rij=sqrt(dx*dx+dy*dy);
            if( rij<rc)
            {
                U += 4 *((pow(1/rij,12) - pow(1/rij,6)) )-U_decale;
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

void Metropolis(Particule tab_particules_lamda_i[], int N, double T, double L) 
{    
    // Tirer une particule au hasard
    int index = rand() % (N);
    double Dr=0.2;
    Particule particule_i = tab_particules_lamda_i[index];
    Particule particule_ip1 = tab_particules_lamda_i[index];

    // Tirer un rayon aléatoire entre 0 et R_MAX
    double r = ((double)rand() / RAND_MAX) * Dr;

    // Tirer un angle aléatoire entre 0 et 2π
    double theta = ((double)rand() / RAND_MAX) * 2.0 * PI;

    // Calculer les nouveaux déplacements en x et y
    double dx = r * cos(theta);
    double dy = r * sin(theta);

    particule_ip1.x+=dx;
    particule_ip1.y+=dy;

    double E_lamda_i=0;//Energie_1part(tab_particules_lamda_i, particule_i, index, N, L);
    double E_lamda_ip1=Energie_1part(tab_particules_lamda_i, particule_ip1, index, N, L);
    double delta_E = E_lamda_ip1 - E_lamda_i;

    double facteur_boltzmann=exp(-delta_E/T);
    double A=fmin(1,facteur_boltzmann);
    double nombre = ((double)rand() / RAND_MAX);

    if(nombre<A)
    {
        tab_particules_lamda_i[index].x = particule_ip1.x;
        tab_particules_lamda_i[index].y = particule_ip1.y;
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

void write_parametre_valeur (FILE* fichier, double parametre, double valeur, int type) {
    if (type == 0) {
        fprintf(fichier, "%-10s | %-10s\n", "Paramètre", "Fonction");
    }

    if (type == 1) {
        fprintf(fichier, "%-10.6f | %-10.6f\n", parametre, valeur);
    }
}

void afficher_barre_progression(int etape_actuelle, int total_etapes) {
    const int largeur_barre = 50; // Largeur de la barre en caractères
    int progression = (etape_actuelle * largeur_barre) / total_etapes;
    printf("\r["); // Retour au début de la ligne
    for (int i = 0; i < largeur_barre; i++) {
        if (i < progression)
            printf("=");
        else
            printf(" ");
    }
    printf("] %d%%", (etape_actuelle * 100) / total_etapes);
    fflush(stdout); // Forcer l'affichage immédiat
}

int main() 
{
    double L=10.0;
    double T=0.05;
    int N=100;
    srand((unsigned int)time(NULL));
    Particule tab_particules_lamda_i[N];

    FILE* fichier_energie = fopen("Energie_rho.txt", "w");
    double densite_min = 0.001; double densite_max = 1.3;
    double energie; double densite = densite_min;
    int nbr_mesure = 1000;
    int nbr_step = 500000;
    double delta_densite = (densite_max - densite_min) / nbr_mesure;

    write_parametre_valeur(fichier_energie, densite, energie/N, 0);
    for (int i = 0; i < nbr_mesure; i++) {
        double Volume = N / densite;
        L = sqrt(Volume);

        Particule tab_particules_lamda_i[N];
        configuration_cristalline(tab_particules_lamda_i,N,L);
        for (int step=0; step<=nbr_step; step++)
        {
            Metropolis(tab_particules_lamda_i,N,T,L);
            conditions_periodiques(tab_particules_lamda_i,N,L);
            if (step = nbr_step) 
            {
                energie = Energie_totale(tab_particules_lamda_i, N, L);
            }
        }
        write_parametre_valeur(fichier_energie, densite, energie/N, 1);

        densite = densite + delta_densite;
        
    }



    fclose(fichier_energie);
   
}