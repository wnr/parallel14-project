#ifndef __CS267_COMMON_H__
#define __CS267_COMMON_H__

#include <vector>

inline int min( int a, int b ) { return a < b ? a : b; }
inline int max( int a, int b ) { return a > b ? a : b; }

#define density 0.0005
#define mass    0.01
#define cutoff  0.01
#define min_r   (cutoff/100)
#define dt      0.0005

//
//  saving parameters
//
const int NSTEPS = 1000;
const int SAVEFREQ = 10;

//
// particle data structure
//
typedef struct 
{
  double x;
  double y;
  double vx;
  double vy;
  double ax;
  double ay;
} particle_t;


struct indexed_particle : public particle_t {
    int index;
};

//
//  timing routines
//
double read_timer( );

//
//  simulation routines
//
void set_size( int n );
double get_size();
void init_particles( int n, particle_t *p );
void apply_force( particle_t &particle, particle_t &neighbor );
void apply_force_mpi( indexed_particle &particle, indexed_particle &neighbor );
void move( particle_t &p );
void move_mpi( indexed_particle &p );

//
//  I/O routines
//
FILE *open_save( char *filename, int n );
void save( FILE *f, int n, particle_t *p );
void save_mpi( FILE *f, int n, indexed_particle *p );


//
//  argument processing routines
//
int find_option( int argc, char **argv, const char *option );
int read_int( int argc, char **argv, const char *option, const int default_value );
char *read_string( int argc, char **argv, const char *option, char *default_value );

#define PARTICLE_SIZE 0.001

class Cell {
public:
	std::vector<particle_t*> particles;

public:
	Cell() {}

	void add(particle_t *p) {
		particles.push_back(p);
	}

	void clear() {
		particles.clear();
	}

	void remove(std::vector<particle_t*>::iterator it) {
		particles.erase(it);
	}

	std::vector<particle_t*>::iterator begin() {
		return particles.begin();
	}

	std::vector<particle_t*>::iterator end() {
		return particles.end();
	}
};

#endif
