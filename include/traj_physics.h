#define trilinon_
#define Uon_  // whether to calculate the electric (V) potential and potential energy (U). Needs Eon to be enabled.
#define Eon_  // whether to calculate the electric (E) field
#define Bon_  // whether to calculate the magnetic (B) field
#define EFon_ // whether to apply electric force
#define BFon_ // whether to apply magnetic force
#define printDensity
#define printParticles
// #define printV //print out V
#define printB // print out B field
#define printE // print out E field
// #define FileIn //whether to load from input file (unused)
#define RamDisk // whether to use RamDisk
#define maxcells 32
#define cldevice 1

#define Temp_e 1e8 // in Kelvin
#define Temp_d 1e8 // in Kelvin
#define Hist_n 256
#define Hist_max 100000 // in eV Kelvin to eV is divide by 11600
                   

constexpr float a0 = 200e-3; // typical dimensions of a cell
// technical parameters
constexpr int n_space = 64;                              // must be 2 to power of n
constexpr int n_partd = n_space * n_space * n_space * 8; // must be 2 to power of n
constexpr int n_parte = n_partd;
constexpr unsigned int ncoeff = 8;

constexpr int n_output_part = (n_partd > 8192) ? 8192 : n_partd; // maximum number of particles to output to file
// const int nprtd=floor(n_partd/n_output_part);

constexpr int ndatapoints = 300; // total number of time steps to calculate
constexpr int nc = 10;           // number of times to calculate E and B between printouts
constexpr int md_me = 60;        // ratio of electron speed/deuteron speed at the same KE. Used to calculate electron motion more often than deuteron motion
constexpr int ncalc0[2] = {md_me, 1};
constexpr int nthreads = 8; // match with your CPU

// The maximum expected E and B fields. If fields go beyond this, the the time step, cell size etc will be wrong. Should adjust and recalculate.
//  maximum expected magnetic field
constexpr float Bmax0 = 0.1;
constexpr float Emax0 = 1e7;

constexpr float target_part = 1e12;
constexpr float r_part_spart = target_part / n_partd; // 1e12 / n_partd; // ratio of particles per tracked "super" particle
// ie. the field of N particles will be multiplied by (1e12/N), as if there were 1e12 particles

constexpr int n_space_divx = n_space;
constexpr int n_space_divy = n_space;
constexpr int n_space_divz = n_space;
constexpr int n_space_divx2 = n_space_divx * 2;
constexpr int n_space_divy2 = n_space_divy * 2;
constexpr int n_space_divz2 = n_space_divz * 2;
constexpr int n_cells = n_space_divx * n_space_divy * n_space_divz;
constexpr int n_cells8 = n_cells * 8;
// physical "constants"
constexpr float kb = 1.38064852e-23;       // m^2kss^-2K-1
constexpr float e_charge = 1.60217662e-19; // C
constexpr float ev_to_j = e_charge;
constexpr float e_mass = 9.10938356e-31;
constexpr float e_charge_mass = e_charge / e_mass;
constexpr float kc = 8.9875517923e9;         // kg m3 s-2 C-2
constexpr float epsilon0 = 8.8541878128e-12; // F m-1
constexpr float pi = 3.1415926536;
constexpr float u0 = 4e-7 * pi;
