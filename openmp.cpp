#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <vector>
#include <omp.h>
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
        printf( "-p <int> to set the number of worker threads\n");
        printf( "-f <int> to set the save frequency\n" );
        return 0;
    }
    
    const int n = read_int( argc, argv, "-n", 1000 );

    const int s = read_int( argc, argv, "-s", NSTEPS);

    const int f = read_int( argc, argv, "-f", SAVEFREQ);

    const int num_threads = read_int( argc, argv, "-p", omp_get_num_procs() );


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

    //
    //  simulate a number of time steps
    //
    double simulation_time = read_timer( );

    double before = 0;

    double t = 0;
    double tc2 = 0;
    double ta = 0;
    double tm = 0;

    int count = 0;

    int threads_used = 0;

    vector<vector<particle_t*> > to_move;
    for(int i = 0; i < omp_get_num_procs(); i++) {
        to_move.push_back(vector<particle_t*>());
    }

    omp_set_num_threads(num_threads);

    #pragma omp parallel
    {
        threads_used = omp_get_num_threads();

        for(int step = 0; step < s; step++ )
        {
            before = read_timer();
            
            #pragma omp for
            for(int i = 0; i < num_cells_side; i++) {
                int thread = omp_get_thread_num();
                for(int j = 0; j < num_cells_side; j++) {
                    for(auto it = (*area)[i][j]->begin(); it != (*area)[i][j]->end();) {
                        particle_t *particle = *it;

                        double x = particle->x;
                        double y = particle->y;

                        int row = x / cell_size;
                        int col = y / cell_size;

                        if(row != i || col != j) {
                            it = (*area)[i][j]->particles.erase(it);
                            to_move[thread].push_back(particle);
                        } else {
                            ++it;
                        }
                    }
                }   
            }

            #pragma omp single
            {
                for(int i = 0; i < to_move.size(); i++) {
                    for(auto it = to_move[i].begin(); it != to_move[i].end(); it++) {
                        particle_t *particle = *it;

                        double x = particle->x;
                        double y = particle->y;

                        int row = x / cell_size;
                        int col = y / cell_size;

                        (*area)[row][col]->add(particle);
                    }

                    count += to_move[i].size();

                    to_move[i].clear();
                }
                
                ta += read_timer() - before;   
            }


            //
            //  compute forces
            //
            
            double before = read_timer();
            #pragma omp for
            for(int i = 0; i < num_cells_side; i++) {
                for(int j = 0; j < num_cells_side; j++) {
                    Cell *cell = (*area)[i][j];
                    vector<Cell*> *cells = (*collissionCells)[i][j];

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

            #pragma omp barrier

            #pragma omp master
            t += read_timer() - before;


            //
            //  move particles
            //
            before = read_timer();
            #pragma omp for
            for( int i = 0; i < n; i++ ) {
                move( particles[i] );
            }

            #pragma omp barrier

            #pragma omp master
            tm += read_timer() - before;
            //
            //  save if necessary
            //
            #pragma omp master
            if( fsave && (step%f) == 0 )
                save( fsave, n, particles );
        }
    }

    simulation_time = read_timer( ) - simulation_time;
    
    printf("\nclear: %f\nadd: %f\nalgo: %f\nmove: %f", tc2, ta, t, tm);
    printf("\ncount: %d", count);

    printf( "\nn = %d, threads = %d, simulation time = %g seconds\n", n, threads_used, simulation_time );
    
    free( particles );
    if( fsave )
        fclose( fsave );
    
    return 0;
}
