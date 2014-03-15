#ifndef __CS267_COMMON_H__
#define __CS267_COMMON_H__

#include <list>

inline int min( int a, int b ) { return a < b ? a : b; }
inline int max( int a, int b ) { return a > b ? a : b; }

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
void move( particle_t &p );

//
//  I/O routines
//
FILE *open_save( char *filename, int n );
void save( FILE *f, int n, particle_t *p );

//
//  argument processing routines
//
int find_option( int argc, char **argv, const char *option );
int read_int( int argc, char **argv, const char *option, const int default_value );
char *read_string( int argc, char **argv, const char *option, char *default_value );

#define PARTICLE_SIZE 0.001

class Cell {
public:
	std::list<particle_t*> particles;

public:
	Cell() {}

	void add(particle_t *p) {
		particles.push_back(p);
	}

	void clear() {
		particles.clear();
	}

	std::list<particle_t*>::iterator begin() {
		return particles.begin();
	}

	std::list<particle_t*>::iterator end() {
		return particles.end();
	}
};

#endif
