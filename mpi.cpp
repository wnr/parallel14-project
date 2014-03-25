#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <vector>
#include <math.h>
#include <vector>
#include "common.h"
#include <algorithm>

using namespace std;

typedef vector<indexed_particle*> cell_t;

struct to_move_element {
    int cell_index;
    int current_row;
    int current_col;

    to_move_element() : cell_index(-1), current_row(-1), current_col(-1) {}
    to_move_element(int ci, int cr, int cc) : cell_index(ci), current_row(cr), current_col(cc) {}
};

//
//  benchmarking program
//
int main( int argc, char **argv )
{    
    //
    //  process command line parameters
    //
    if( find_option( argc, argv, "-h" ) >= 0 )
    {
        printf( "Options:\n" );
        printf( "-h to see this help\n" );
        printf( "-n <int> to set the number of particles\n" );
        printf( "-f <int> to set the save frequency\n" );
        printf( "-s <int> to set the number of steps\n" );
        printf( "-o <filename> to specify the output file name\n" );
        return 0;
    }
    
    int n = read_int( argc, argv, "-n", 1000 );
    char *savename = read_string( argc, argv, "-o", NULL );

    const int s = read_int( argc, argv, "-s", NSTEPS);

    int f = read_int( argc, argv, "-f", SAVEFREQ);
    
    //
    //  set up MPI
    //
    int n_proc, rank;
    MPI_Init( &argc, &argv );
    MPI_Comm_size( MPI_COMM_WORLD, &n_proc );
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    
    //
    //  allocate generic resources
    //
    FILE *fsave = savename && rank == 0 ? fopen( savename, "w" ) : NULL;
    indexed_particle *particles = (indexed_particle*) malloc( n * sizeof(indexed_particle) );
    
    MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI::ERRORS_THROW_EXCEPTIONS);

    MPI_Datatype PARTICLE;
    {
        indexed_particle p;

        int blen[2] = {6, 1};
        MPI_Aint offsets[2] = {
            (int)((char*)&p.x - (char*)&p),
            (int)((char*)&p.index - (char*)&p)
        };
        MPI_Datatype types[2] = {MPI_DOUBLE, MPI_INT};
        if(MPI_Type_create_struct(2, blen, offsets, types, &PARTICLE) != MPI_SUCCESS) {
            exit(3);
        }
        MPI_Type_commit( &PARTICLE );
    }
    
    
    MPI_Datatype TO_MOVE_ELEMENT;
    MPI_Type_contiguous(3, MPI_INTEGER, &TO_MOVE_ELEMENT);
    MPI_Type_commit(&TO_MOVE_ELEMENT);

    //
    //  set up the data partitioning across processors
    //
    int particle_per_proc = (n + n_proc - 1) / n_proc;
    int *partition_offsets = (int*) malloc( (n_proc+1) * sizeof(int) );
    for( int i = 0; i < n_proc+1; i++ )
        partition_offsets[i] = min( i * particle_per_proc, n );
    
    int *partition_sizes = (int*) malloc( n_proc * sizeof(int) );
    for( int i = 0; i < n_proc; i++ )
        partition_sizes[i] = partition_offsets[i+1] - partition_offsets[i];
    
    double before;

    set_size( n );

    //
    //  allocate storage for local partition
    //
    int nlocal = partition_sizes[rank];
    indexed_particle *local = (indexed_particle*) malloc( nlocal * sizeof(indexed_particle) );

    int num_cells_side = ((sqrt(n)) + 0.5);
    double cell_size = (get_size()/num_cells_side);
    cell_size *= 1.01; // Fix so that the cells cover 101 % of the map, to make sure that we don't miss any coordinates

    cell_t* area[num_cells_side * num_cells_side];

    vector<to_move_element> to_move;
    vector<indexed_particle*> to_move_particles;
    vector<indexed_particle> moved_particles;

    int n_threads = n_proc;
    int cells_per_thread = (num_cells_side + n_threads - 1) / n_threads;
    int first_cell = min(rank * cells_per_thread, num_cells_side);
    int last_cell = min((rank+1) * cells_per_thread, num_cells_side);

    MPI_Status status;
    
    //
    //  initialize and distribute the particles (that's fine to leave it unoptimized)
    //
    
    if( rank == 0 ) {
        srand48( 0 ); //TODO time(NULL)
        
        int sx = (int)ceil(sqrt((double)n));
        int sy = (n+sx-1)/sx;
        
        int *shuffle = (int*)malloc( n * sizeof(int) );
        for( int i = 0; i < n; i++ )
            shuffle[i] = i;
        
        for( int i = 0; i < n; i++ ) 
        {
            //
            //  make sure particles are not spatially sorted
            //
            int j = lrand48()%(n-i);
            int k = shuffle[j];
            shuffle[j] = shuffle[n-i-1];
            
            //
            //  distribute particles evenly to ensure proper spacing
            //
            particles[i].x = get_size()*(1.+(k%sx))/(1+sx);
            particles[i].y = get_size()*(1.+(k/sx))/(1+sy);

            //
            //  assign random velocities within a bound
            //
            particles[i].vx = drand48()*2-1;
            particles[i].vy = drand48()*2-1;
        }
        free( shuffle );
    }

    MPI_Bcast( particles, n, PARTICLE, 0, MPI_COMM_WORLD );

    for(int i = 0; i < num_cells_side; i++) {
        for(int j = 0; j < num_cells_side; j++) {
            area[i * num_cells_side + j] = new cell_t;
        }
    }

    for(int i = 0; i < n; i++) {
        particles[i].index = i;
        double x = particles[i].x;
        double y = particles[i].y;

        int row = x / cell_size;
        int col = y / cell_size;

        area[row * num_cells_side + col]->push_back(&particles[i]);
    }

    double tfm = 0;
    double tsm = 0;
    double tcm = 0;
    double tf = 0;
    double tsf = 0;
    double tm = 0;

    //
    //  simulate a number of time steps
    //
    double simulation_time = read_timer( );
    for( int step = 0; step < s; step++ )
    {
        before = read_timer();

        //for(int i = 0; i < num_cells_side; i++) {
        for( int i = first_cell; i < last_cell; i++ ) {
            for(int j = 0; j < num_cells_side; j++) {
                int cnt = 0;
                for(auto it = area[i * num_cells_side + j]->begin(); it != area[i * num_cells_side + j]->end(); it++, cnt++) {
                    indexed_particle *particle = *it;

                    double x = particle->x;
                    double y = particle->y;

                    int row = x / cell_size;
                    int col = y / cell_size;

                    if(row != i || col != j) {
                        int cell_index = distance(area[i * num_cells_side + j]->begin(), it);
                        to_move.push_back(to_move_element(cell_index, i, j));
                    }
                }
            }
        }

        tfm += read_timer() - before;

        const int PARTICLE_CELL_MOVES_TAG = 2;

        MPI_Request request;
        before = read_timer();

        if(rank == 0) {
            int to_move_n;
            for(int r = 1; r < n_proc; r++) {
                MPI_Probe(r, PARTICLE_CELL_MOVES_TAG, MPI_COMM_WORLD, &status);
                MPI_Get_count(&status, TO_MOVE_ELEMENT, &to_move_n);

                if(to_move_n == MPI_UNDEFINED) {
                    printf("error, invalid count sent."); fflush(stdout);
                    exit(1);
                }

                to_move_element *buffer = (to_move_element*)malloc(to_move_n * sizeof(to_move_element));

                MPI_Recv(buffer, to_move_n, TO_MOVE_ELEMENT, r, PARTICLE_CELL_MOVES_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                //to_move.reserve(to_move.size() + to_move_n);
                copy(buffer, buffer + to_move_n, back_inserter(to_move));
            }

            for(int r = 1; r < n_proc; r++) {
                MPI_Ssend((void*)&to_move[0], to_move.size(), TO_MOVE_ELEMENT, r, PARTICLE_CELL_MOVES_TAG, MPI_COMM_WORLD);
            }

        } else {
            int to_move_n;
            MPI_Ssend((void*)&to_move[0], to_move.size(), TO_MOVE_ELEMENT, 0, PARTICLE_CELL_MOVES_TAG, MPI_COMM_WORLD);
            MPI_Probe(0, PARTICLE_CELL_MOVES_TAG, MPI_COMM_WORLD, &status);
            MPI_Get_count(&status, TO_MOVE_ELEMENT, &to_move_n);

            if(to_move_n == MPI_UNDEFINED) {
                printf("error, invalid count sent 2."); fflush(stdout);
                exit(1);
            }

            //to_move.resize(to_move.size() + to_move_n);
            to_move_element *buffer = (to_move_element*)malloc(to_move_n * sizeof(to_move_element));
            MPI_Recv(buffer, to_move_n, TO_MOVE_ELEMENT, 0, PARTICLE_CELL_MOVES_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            to_move = vector<to_move_element>(buffer, buffer + to_move_n);

        }

        tsm += read_timer() - before;

        before = read_timer();

        for(auto it = to_move.begin(); it != to_move.end(); it++) {
            to_move_element e = *it;

            indexed_particle *p = area[e.current_row * num_cells_side + e.current_col]->at(e.cell_index);
            to_move_particles.push_back(p);
        }

        for(int i = 0; i < to_move_particles.size(); i++) {
            indexed_particle *particle = to_move_particles[i];

            double x = particle->x;
            double y = particle->y;

            int row = x / cell_size;
            int col = y / cell_size;
    
            area[row * num_cells_side + col]->push_back(particle);

            cell_t *c = area[to_move[i].current_row * num_cells_side + to_move[i].current_col];
            c->erase(remove(c->begin(), c->end(), particle));
        }
        
        to_move.clear();
        to_move_particles.clear();

        moved_particles.clear();

        tcm += read_timer() - before;

        before = read_timer();
        for(int i = first_cell; i < last_cell; i++) {
            for(int j = 0; j < num_cells_side; j++) {
                cell_t *cell = area[i * num_cells_side + j];

                auto it = cell->begin();
                for(; it != cell->end(); it++) {
                    indexed_particle *p = *it;
                    p->ax = p->ay = 0;

                    for(int row = i-1; row <= i+1; row++) {
                        for(int col = j-1; col <= j+1; col++) {
                            if(row >= 0 && col >= 0 && row < num_cells_side && col < num_cells_side) {
                                cell_t *c = area[row * num_cells_side + col];
                                auto p_it = c->begin();
                                for(; p_it != c->end(); p_it++) {
                                    indexed_particle *other_p = *p_it;
                                    apply_force_mpi(*p, *other_p);
                                }
                            }
                        }
                    }

                    moved_particles.push_back(*p);
                }
            }
        }

        tf += read_timer() - before;

        const int FORCE_APPLIED_SYNC_TAG = 1;

        before = read_timer();
        if(rank == 0) {
            int moved_n;
            for(int r = 1; r < n_proc; r++) {
                MPI_Probe(r, FORCE_APPLIED_SYNC_TAG, MPI_COMM_WORLD, &status);
                MPI_Get_count(&status, PARTICLE, &moved_n);

                if(moved_n == MPI_UNDEFINED) {
                    printf("error, invalid count sent."); fflush(stdout);
                    exit(2);
                }

                indexed_particle *buffer = (indexed_particle*)malloc(moved_n * sizeof(indexed_particle));

                MPI_Recv(buffer, moved_n, PARTICLE, r, FORCE_APPLIED_SYNC_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                //moved_particles.reserve(moved_particles.size() + moved_n);
                copy(buffer, buffer + moved_n, back_inserter(moved_particles));
            }

            assert(moved_particles.size() == n);

            for(int r = 1; r < n_proc; r++) {
                MPI_Ssend(&moved_particles[0], moved_particles.size(), PARTICLE, r, FORCE_APPLIED_SYNC_TAG, MPI_COMM_WORLD);
            }
        } else {
            int moved_n;
            MPI_Ssend(&moved_particles[0], moved_particles.size(), PARTICLE, 0, FORCE_APPLIED_SYNC_TAG, MPI_COMM_WORLD);
            MPI_Probe(0, FORCE_APPLIED_SYNC_TAG, MPI_COMM_WORLD, &status);
            MPI_Get_count(&status, PARTICLE, &moved_n);

            if(moved_n == MPI_UNDEFINED) {
                printf("error, invalid count sent."); fflush(stdout);
                exit(2);
            }

            //moved_particles.reserve(moved_particles.size() + moved_n);
            indexed_particle *buffer = (indexed_particle*)malloc(moved_n * sizeof(indexed_particle));
            MPI_Recv(buffer, moved_n, PARTICLE, 0, FORCE_APPLIED_SYNC_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            moved_particles = vector<indexed_particle>(buffer, buffer + moved_n);
        }

        tsf += read_timer() - before;

        before = read_timer();
        for(int i = 0; i < moved_particles.size(); i++) {

            //printf("mpi: %d, x: %f, y: %f, vx: %f, vy: %f, ax: %f, ay: %f, index: %d\n", rank, moved_particles[i].x, moved_particles[i].y, moved_particles[i].vx, moved_particles[i].vy, moved_particles[i].ax, moved_particles[i].ay, moved_particles[i].index);
                indexed_particle *mp = &moved_particles[i];
                indexed_particle *p = &particles[mp->index];
                p->ax = mp->ax;
                p->ay = mp->ay;
                move_mpi(particles[mp->index]);
            
        }

        tm += read_timer() - before;

        if(rank == 0 && fsave && (step%f) == 0) {
            save_mpi( fsave, n, particles );
        }
    }
    simulation_time = read_timer( ) - simulation_time;

    printf("rank: %d\ntfm: %f\ntsm: %f\ntcm: %f\ntf: %f\ntsf: %f\ntm: %f\n", rank, tfm, tsm, tcm, tf, tsf, tm);
    if( rank == 0 ) {
        printf( "n = %d, n_procs = %d, simulation time = %g s\n", n, n_proc, simulation_time );
    }
    
    //
    //  release resources
    //
    free( partition_offsets );
    free( partition_sizes );
    free( local );
    free( particles );
    if( fsave )
        fclose( fsave );
    
    MPI_Finalize( );
    
    return 0;
}
