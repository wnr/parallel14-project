#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <vector>
#include "common.h"


using namespace std;

//
//  benchmarking program
//
int main( int argc, char **argv )
{    
    if( find_option( argc, argv, "-h" ) >= 0 )
    {
        printf( "Options:\n" );
        printf( "-h to see this help\n" );
        printf( "-n <int> to set the number of particles\n" );
        printf( "-o <filename> to specify the output file name\n" );
        printf( "-s <int> to set the number of steps\n" );
        printf( "-f <int> to set the save frequency\n" );
        return 0;
    }
    
    const int n = read_int( argc, argv, "-n", 1000 );

    const int s = read_int( argc, argv, "-s", NSTEPS);

    const int f = read_int( argc, argv, "-f", SAVEFREQ);

    char *savename = read_string( argc, argv, "-o", NULL );
    
    FILE *fsave = savename ? fopen( savename, "w" ) : NULL;
    particle_t *particles = (particle_t*) malloc( n * sizeof(particle_t) );
    set_size( n );
    init_particles( n, particles );
    
    int num_cells_side = ((sqrt(n)) + 0.5);
    double cell_size = (get_size()/num_cells_side);
    cell_size *= 1.01; // Fix so that the cells cover 101 % of the map, to make sure that we don't miss any coordinates

    vector<vector<Cell*> > *area = new vector<vector<Cell*> >();
    vector<vector<vector<Cell*>*> > *collissionCells = new vector<vector<vector<Cell*>*> >;

    for(int i = 0; i < num_cells_side; i++) {
        area->push_back(vector<Cell*>());
        collissionCells->push_back(vector<vector<Cell*>*>());

        for(int j = 0; j < num_cells_side; j++) {
            (*area)[i].push_back(new Cell());
            (*collissionCells)[i].push_back(new vector<Cell*>());
        }
    }

    for(int i = 0; i < num_cells_side; i++) {
        for(int j = 0; j < num_cells_side; j++) {
            for(int row = i-1; row <= i+1; row++) {
                for(int col = j-1; col <= j+1; col++) {
                    if(row >= 0 && col >= 0 && row < num_cells_side && col < num_cells_side) {
                        (*collissionCells)[i][j]->push_back((*area)[row][col]);
                    }
                }
            }
        }
    }    

    for(int i = 0; i < n; i++) {
        double x = particles[i].x;
        double y = particles[i].y;

        int row = x / cell_size;
        int col = y / cell_size;
        (*area)[row][col]->add(&particles[i]);
    }

    double before = 0;

    double t = 0;
    double tc2 = 0;
    double ta = 0;
    double tm = 0;

    //
    //  simulate a number of time steps
    //
    double simulation_time = read_timer( );
    for( int step = 0; step < s; step++ )
    {
        before = read_timer();

        for(int i = 0; i < num_cells_side; i++) {
            for(int j = 0; j < num_cells_side; j++) {
                for(auto it = (*area)[i][j]->begin(); it != (*area)[i][j]->end();) {
                    particle_t *particle = *it;

                    double x = particle->x;
                    double y = particle->y;

                    int row = x / cell_size;
                    int col = y / cell_size;

                    if(row != i || col != j) {
                        it = (*area)[i][j]->particles.erase(it);
                        (*area)[row][col]->add(particle);
                    } else {
                        ++it;
                    }
                }
            }   
        }

        ta += read_timer() - before;

        //
        //  compute forces
        //
        
        before = read_timer();
        for(int i = 0; i < n; i++) {
            particle_t *p = &particles[i];

            p->ax = p->ay = 0;


            double x = p->x;
            double y = p->y;

            int r = x / cell_size;
            int c = y / cell_size;

            for(int row = r-1; row <= r+1; row++) {
                for(int col = c-1; col <= c+1; col++) {
                    if(row >= 0 && col >= 0 && row < num_cells_side && col < num_cells_side) {
                        Cell *cell = (*area)[row][col];
                        for(auto p_it = cell->begin(); p_it != cell->end(); p_it++) {
                            particle_t *other_p = *p_it;
                            apply_force(*p, *other_p);
                        }
                    }
                }
            }
        }
        t += read_timer() - before;

        //
        //  move particles
        //
        before = read_timer();
        for( int i = 0; i < n; i++ ) {
            //printf("serial: %d, x: %f, y: %f, vx: %f, vy: %f, ax: %f, ay: %f\n", i, particles[i].x, particles[i].y, particles[i].vx, particles[i].vy, particles[i].ax, particles[i].ay);
            fflush(stdout);
            move( particles[i] );
        }
        tm += read_timer() - before;
        
        //
        //  save if necessary
        //
        if( fsave && (step%f) == 0 )
            save( fsave, n, particles );
    }
    simulation_time = read_timer( ) - simulation_time;
    
    printf("\nclear: %f\nadd: %f\nalgo: %f\nmove: %f", tc2, ta, t, tm);

    printf( "\nn = %d, simulation time = %g seconds\n", n, simulation_time );
    
    free( particles );
    if( fsave )
        fclose( fsave );
    
    return 0;
}
