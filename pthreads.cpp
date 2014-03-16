#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <pthread.h>
#include <vector>
#include "common.h"

using namespace std;

//
//  global variables
//
int n, n_threads;
particle_t *particles;
FILE *fsave;
pthread_mutex_t lock;

int f;
int num_cells_side;
double cell_size;
vector<vector<Cell*> > *area;
vector<vector<vector<Cell*>*> > *collissionCells;
double t = 0;
double tc2 = 0;
double ta = 0;
double tm = 0;

//Barrier
pthread_mutex_t barrier_mutex;
pthread_cond_t go;
int num_arrived = 0;
void barrier() {
    pthread_mutex_lock(&barrier_mutex);
    num_arrived++;
    if(num_arrived == n_threads) {
        num_arrived = 0;
        pthread_cond_broadcast(&go);
    } else {
        pthread_cond_wait(&go, &barrier_mutex);
    }
    pthread_mutex_unlock(&barrier_mutex);
}

//
//  check that pthreads routine call was successful
//
#define P( condition ) {if( (condition) != 0 ) { printf( "\n FAILURE in %s, line %d\n", __FILE__, __LINE__ );exit( 1 );}}

//
//  This is where the action happens
//
void *thread_routine( void *pthread_id )
{
    int thread_id = *(int*)pthread_id;

    int particles_per_thread = (n + n_threads - 1) / n_threads;
    int first_particle = min(  thread_id    * particles_per_thread, n );
    int last_particle  = min( (thread_id+1) * particles_per_thread, n );

    int cells_per_thread = (num_cells_side + n_threads - 1) / n_threads;
    int first_cell = min(thread_id * cells_per_thread, num_cells_side);
    int last_cell = min((thread_id+1) * cells_per_thread, num_cells_side);

    double before;

    vector<particle_t*> to_move;
    
    //
    //  simulate a number of time steps
    //
    for( int step = 0; step < NSTEPS; step++ )
    {
        if(thread_id == 0) {
            before = read_timer();
        }
        
        //for(int i = 0; i < num_cells_side; i++) {
        for( int i = first_cell; i < last_cell; i++ ) {
            for(int j = 0; j < num_cells_side; j++) {
                for(auto it = (*area)[i][j]->begin(); it != (*area)[i][j]->end();) {
                    particle_t *particle = *it;

                    double x = particle->x;
                    double y = particle->y;

                    int row = x / cell_size;
                    int col = y / cell_size;

                    if(row != i || col != j) {
                        it = (*area)[i][j]->particles.erase(it);
                        to_move.push_back(particle);
                    } else {
                        ++it;
                    }
                }
            }   
        }

        barrier();

        pthread_mutex_lock(&lock);

        for(auto it = to_move.begin(); it != to_move.end(); it++) {
            particle_t *particle = *it;

            double x = particle->x;
            double y = particle->y;

            int row = x / cell_size;
            int col = y / cell_size;

            (*area)[row][col]->add(particle);
        }

        pthread_mutex_unlock(&lock);
        
        to_move.clear();

        barrier();

        if(thread_id == 0) {
            ta += read_timer() - before;
        }

        //
        //  compute forces
        //
        
        if(thread_id == 0) {
            before = read_timer();
        }

        for(int i = first_cell; i < last_cell; i++) {
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

        barrier();

        if(thread_id == 0) {
            t += read_timer() - before;
        }

        //
        //  move particles
        //
        if(thread_id == 0) {
            before = read_timer();
        }

        for( int i = first_particle; i < last_particle; i++ ) {
            move( particles[i] );
        }

        barrier();

        if(thread_id == 0) {
            tm += read_timer() - before;

            if( fsave && (step % f) == 0 ) {
                save( fsave, n, particles );
            }
        }
    }
    
    return NULL;
}

//
//  benchmarking program
//
int main( int argc, char **argv )
{    
    //
    //  process command line
    //
    if( find_option( argc, argv, "-h" ) >= 0 )
    {
        printf( "Options:\n" );
        printf( "-h to see this help\n" );
        printf( "-n <int> to set the number of particles\n" );
        printf( "-o <filename> to specify the output file name\n" );
        printf( "-s <int> to set the number of steps\n" );
        printf( "-f <int> to set the save frequency\n" );
        printf( "-p <int> to set the number of worker threads\n");
        return 0;
    }
    
    n = read_int( argc, argv, "-n", 1000 );

    const int s = read_int( argc, argv, "-s", NSTEPS);

    f = read_int( argc, argv, "-f", SAVEFREQ);

    char *savename = read_string( argc, argv, "-o", NULL );

    n_threads = read_int( argc, argv, "-p", 2 );
    
    //
    //  allocate resources
    //
    fsave = savename ? fopen( savename, "w" ) : NULL;

    particles = (particle_t*) malloc( n * sizeof(particle_t) );
    set_size( n );
    init_particles( n, particles );

    num_cells_side = ((sqrt(n)) + 0.5);
    cell_size = (get_size()/num_cells_side);
    cell_size *= 1.01; // Fix so that the cells cover 101 % of the map, to make sure that we don't miss any coordinates

    area = new vector<vector<Cell*> >();
    collissionCells = new vector<vector<vector<Cell*>*> >;

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

    pthread_attr_t attr;
    P( pthread_attr_init( &attr ) );
    P( pthread_cond_init( &go, NULL) );
    P( pthread_mutex_init( &barrier_mutex, NULL) );


    pthread_mutex_init(&lock, NULL);

    int *thread_ids = (int *) malloc( n_threads * sizeof( int ) );
    for( int i = 0; i < n_threads; i++ ) 
        thread_ids[i] = i;

    pthread_t *threads = (pthread_t *) malloc( n_threads * sizeof( pthread_t ) );
    
    //
    //  do the parallel work
    //
    double simulation_time = read_timer( );

    t = 0;
    ta = 0;
    tm = 0;

    for( int i = 1; i < n_threads; i++ ) 
        P( pthread_create( &threads[i], &attr, thread_routine, &thread_ids[i] ) );
    
    thread_routine( &thread_ids[0] );
    
    for( int i = 1; i < n_threads; i++ ) 
        P( pthread_join( threads[i], NULL ) );
    simulation_time = read_timer( ) - simulation_time;

    printf("\nclear: %f\nadd: %f\nalgo: %f\nmove: %f", tc2, ta, t, tm);
    
    printf( "n = %d, n_threads = %d, simulation time = %g seconds\n", n, n_threads, simulation_time );
    
    //
    //  release resources
    //
    P( pthread_cond_destroy( &go ) );
    P( pthread_mutex_destroy(&barrier_mutex));
    P( pthread_mutex_destroy(&lock));
    P( pthread_attr_destroy( &attr ) );
    free( thread_ids );
    free( threads );
    free( particles );
    if( fsave )
        fclose( fsave );
    
    return 0;
}
