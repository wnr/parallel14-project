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
    vector<vector<list<Cell*>*> > *collissionCells = new vector<vector<list<Cell*>*> >;

    for(int i = 0; i < num_cells_side; i++) {
        area->push_back(vector<Cell*>());
        collissionCells->push_back(vector<list<Cell*>*>());

        for(int j = 0; j < num_cells_side; j++) {
            (*area)[i].push_back(new Cell());
            (*collissionCells)[i].push_back(new list<Cell*>());
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

    //
    //  simulate a number of time steps
    //
    double simulation_time = read_timer( );
    for( int step = 0; step < s; step++ )
    {
        for(int i = 0; i < num_cells_side; i++) {
            for(int j = 0; j < num_cells_side; j++) {
                (*area)[i][j]->clear();
            }
        }

        for(int i = 0; i < n; i++) {
            double x = particles[i].x;
            double y = particles[i].y;


            int row = x / cell_size;
            int col = y / cell_size;
            (*area)[row][col]->add(&particles[i]);

            // double min_x = x - PARTICLE_SIZE;
            // double min_y = y - PARTICLE_SIZE;
            // double max_x = x + PARTICLE_SIZE;
            // double max_y = y + PARTICLE_SIZE;


            // int min_row = min_x / cell_size;
            // int max_row = max_x / cell_size;
            // int min_col = min_y / cell_size;
            // int max_col = max_y / cell_size;
            
            // for(int row = min_row; row <= max_row; row++) {
            //     for(int col = min_col; col <= max_col; col++) {
            //         if(row > 0 && col > 0 && row < num_cells_side && col < num_cells_side) {
            //             area[row][col].add(particles[i]);
            //         }
            //     }
            // }
        }

        //
        //  compute forces
        //
        
        for(int i = 0; i < num_cells_side; i++) {
            for(int j = 0; j < num_cells_side; j++) {
                Cell *cell = (*area)[i][j];
                list<Cell*> *cells = (*collissionCells)[i][j];

                auto it = cell->begin();
                for(; it != cell->end(); it++) {
                    particle_t *p = *it;
                    p->ax = p->ay = 0;

                    auto cells_it = cells->begin();
                    for(; cells_it != cells->end(); cells_it++) {
                        auto c = *cells_it;

                        auto p_it = c->begin();
                        for(; p_it != c->end(); p_it++) {
                            particle_t *other_p = *p_it;
                            apply_force(*p, *other_p);
                        }
                    }
                }
            }
        }

        //
        //  move particles
        //
        for( int i = 0; i < n; i++ ) 
            move( particles[i] );
        
        //
        //  save if necessary
        //
        if( fsave && (step%f) == 0 )
            save( fsave, n, particles );
    }
    simulation_time = read_timer( ) - simulation_time;
    
    printf( "n = %d, simulation time = %g seconds\n", n, simulation_time );
    
    free( particles );
    if( fsave )
        fclose( fsave );
    
    return 0;
}
