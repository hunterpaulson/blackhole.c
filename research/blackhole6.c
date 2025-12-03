// bh_ascii.c — Simplified Schwarzschild ASCII rings (3 rotating hotspots)
// build: cc -O3 -lm bh_ascii.c -o bh_ascii
// run:   ./bh_ascii [inc_deg] [FOVx_deg] [robs]

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#define W 120
#define HEIGHT 50

// Units: G=c=1
const double Mbh = 1.0;
static inline double A(double r){ return 1.0 - 2.0*Mbh/r; }

// Scene params (fixed; only observer is configurable)
const double rin = 6.0, rout = 40.0;
double robs = 50.0;                // observer radius
double inc_deg = 10.0;             // inclination (deg)
double theta_obs = 0.0;            // set in main from inc_deg
const double phi_obs = 0.0;
double FOVx = 60.0*M_PI/180.0;     // horizontal FOV
double FOVy = 0.0;                 // set in main from FOVx
const double emiss_p = 2.0;        // emissivity ~ r^-p

// Simple Saturn-like rings: bright/dark radial bands, no azimuthal texture
static inline double ring_mul(double r){
  if(r < rin) r = rin; if(r > rout) r = rout;
  double s = (r - rin) / (rout - rin); // 0..1 across radial extent
  // Fewer bands with wider dark gaps and smooth edges
  const double Nbands    = 8.0;   // fewer concentric rings
  const double fill_frac = 0.30;  // bright fraction per band (smaller => wider gaps)
  const double edge_soft = 0.02;  // soft edge width (in band fraction)
  const double band_floor= 0.12;  // darkness of gaps
  const double peak      = 1.45;  // brightness of bands
  double pos = Nbands * s;
  double f   = pos - floor(pos);  // position within current band [0,1)
  double w   = edge_soft + 1e-6;
  // Smooth step: near 1 inside bright band (f < fill_frac), near 0 in gap
  double t   = 0.5 + 0.5 * tanh((fill_frac - f) / w);
  return band_floor + (peak - band_floor) * t;
}

// Minimal azimuthal modulation for animation (rotates with phase)
// Three large circular hotspots (disks with soft edges) rotating clockwise.
// Centers sit at radius ~ rout/2; diameter ~ rout (radius ~ rout/2).
static inline double hotspots_mul(double r, double phi, double phase){
  const int    N      = 1;           // number of hotspots
  const double amp    = 3.0;        // strength per hotspot (added to 1)
  const double rc     = 0.5 * rout;  // center radius (halfway to edge)
  const double Rh     = 0.5 * rout;  // hotspot radius (diameter = rout)
  const double edge   = 0.1 * rout; // soft edge width
  // pixel position in disk plane
  double x = r * cos(phi), y = r * sin(phi);
  double m = 1.0;
  for(int k=0;k<N;k++){
    double ang = -phase + 2.0*M_PI*(double)k/(double)N; // clockwise
    double cx = rc * cos(ang);
    double cy = rc * sin(ang);
    double dx = x - cx, dy = y - cy;
    double d  = sqrt(dx*dx + dy*dy);
    // soft disk: ~1 inside Rh, falls to 0 outside over 'edge'
    double t  = 0.5 + 0.5 * tanh((Rh - d) / (edge + 1e-9));
    m += amp * t;
  }
  return m;
}

// (azimuthal modulation removed to keep focus on tricorn bumps)

typedef struct { double r, phi, g, emiss; int hit; } Hit; // hit=1 if disk

// Metric g_{\mu\nu} at (r,th)
static void metric(double r, double th, double g[4][4]){
  double Ar = A(r), s = sin(th), s2 = s*s;
  memset(g,0,sizeof(double)*16);
  g[0][0] = -Ar;
  g[1][1] = 1.0/Ar;
  g[2][2] = r*r;
  g[3][3] = r*r*s2;
}

// Christoffels needed and acceleration a^\mu
static void accel(double x[4], double v[4], double a[4]){
  double r=x[1], th=x[2], s=sin(th), c=cos(th), Ar=A(r);
  (void)Ar; // used only for readability
  double Gttr = Mbh/(r*(r-2.0*Mbh));               // \Gamma^t_{tr}=\Gamma^t_{rt}
  double Grtt = Ar*Mbh/(r*r);                      // \Gamma^r_{tt}
  double Grrr = -Mbh/(r*(r-2.0*Mbh));              // \Gamma^r_{rr}
  double Grthth = -(r-2.0*Mbh);                    // \Gamma^r_{\theta\theta}
  double Grphph = -(r-2.0*Mbh)*s*s;                // \Gamma^r_{\phi\phi}
  double Gthrth = 1.0/r;                           // \Gamma^\theta_{r\theta}=\Gamma^\theta_{\theta r}
  double Gthphph = -s*c;                           // \Gamma^\theta_{\phi\phi}
  double Gphrph = 1.0/r;                           // \Gamma^\phi_{r\phi}=\Gamma^\phi_{\phi r}
  double Gphthph = (c/(s+1e-12));                  // \Gamma^\phi_{\theta\phi}=\Gamma^\phi_{\phi\theta}

  double vt=v[0], vr=v[1], vth=v[2], vph=v[3];
  a[0] = -2.0*Gttr*vt*vr;
  a[1] = - (Grtt*vt*vt + Grrr*vr*vr + Grthth*vth*vth + Grphph*vph*vph);
  a[2] = - (2.0*Gthrth*vr*vth + Gthphph*vph*vph);
  a[3] = - (2.0*Gphrph*vr*vph + 2.0*Gphthph*vth*vph);
}

static void rk4(double x[4], double v[4], double h){
  double k1x[4],k2x[4],k3x[4],k4x[4];
  double k1v[4],k2v[4],k3v[4],k4v[4];
  double a[4], xt[4], vt[4];
  // k1
  accel(x,v,a);
  for(int i=0;i<4;i++){ k1x[i]=h*v[i]; k1v[i]=h*a[i]; xt[i]=x[i]+0.5*k1x[i]; vt[i]=v[i]+0.5*k1v[i]; }
  // k2
  accel(xt,vt,a);
  for(int i=0;i<4;i++){ k2x[i]=h*vt[i]; k2v[i]=h*a[i]; xt[i]=x[i]+0.5*k2x[i]; vt[i]=v[i]+0.5*k2v[i]; }
  // k3
  accel(xt,vt,a);
  for(int i=0;i<4;i++){ k3x[i]=h*vt[i]; k3v[i]=h*a[i]; xt[i]=x[i]+k3x[i]; vt[i]=v[i]+k3v[i]; }
  // k4
  accel(xt,vt,a);
  for(int i=0;i<4;i++){ k4x[i]=h*vt[i]; k4v[i]=h*a[i]; }
  for(int i=0;i<4;i++){
    x[i]+= (k1x[i]+2*k2x[i]+2*k3x[i]+k4x[i])/6.0;
    v[i]+= (k1v[i]+2*k2v[i]+2*k3v[i]+k4v[i])/6.0;
  }
  if(x[2]<1e-6) x[2]=1e-6; if(x[2]>M_PI-1e-6) x[2]=M_PI-1e-6;
}

static void pix_ray(int px,int py,double x0[4], double v0[4]){
  double u = (px+0.5)/(double)W - 0.5;
  double v = (py+0.5)/(double)HEIGHT - 0.5;
  double ax = u*FOVx;
  double ay = v*FOVy;  // flip vertical so image is not upside-down
  double nr=-1.0, nth=tan(ay), nph=tan(ax);
  double norm = sqrt(nr*nr+nth*nth+nph*nph);
  nr/=norm; nth/=norm; nph/=norm;
  double Ar=A(robs), s=sin(theta_obs);
  x0[0]=0.0; x0[1]=robs; x0[2]=theta_obs; x0[3]=phi_obs;
  v0[0]=1.0/sqrt(Ar);
  v0[1]=nr*sqrt(Ar);
  v0[2]=nth/robs;
  v0[3]=nph/(robs*(s>1e-12?s:1e-12));
}

static Hit trace_pixel(int px,int py){
  Hit H; H.hit=0;
  double x[4], v[4]; pix_ray(px,py,x,v);
  double th_prev=x[2], x_prev[4], v_prev[4];
  for(int i=0;i<4;i++){ x_prev[i]=x[i]; v_prev[i]=v[i]; }
  const double h0=0.5, rh=2.0*Mbh;
  for(int step=0; step<5000; ++step){
    double h=h0;
    if(x[1]<10.0) h=0.25*h0;
    if(x[1]<6.0)  h=0.125*h0;
    rk4(x,v,h);
    if(x[1]<=1.001*rh) return H;               // captured
    if(x[1]>1.2*robs && step>10) return H;     // escaped
    // equatorial crossing
    if((th_prev-M_PI/2.0)*(x[2]-M_PI/2.0)<=0.0){
      double f=(M_PI/2.0 - th_prev) / (x[2]-th_prev + 1e-15);
      double rhit = x_prev[1] + f*(x[1]-x_prev[1]);
      double phit = x_prev[3] + f*(x[3]-x_prev[3]);
      if(rhit>=rin && rhit<=rout){
        // p_mu at hit (linear interp of v)
        double vh[4]; for(int i=0;i<4;i++) vh[i]=v_prev[i]+f*(v[i]-v_prev[i]);
        double gmn[4][4]; metric(rhit,M_PI/2.0,gmn);
        double pmu[4]={0};
        for(int a=0;a<4;a++) for(int b=0;b<4;b++) pmu[a]+=gmn[a][b]*vh[b];
        // E_obs
        double ut_obs = 1.0/sqrt(A(robs));
        double Eobs = -(pmu[0]*ut_obs);
        // E_em (Keplerian)
        double denom = sqrt(1.0-3.0*Mbh/rhit);
        double ut = 1.0/denom;
        double uphi = sqrt(Mbh/(rhit*rhit*rhit))/denom;
        double Eem = -(pmu[0]*ut + pmu[3]*uphi);
        double g = (Eobs/(Eem>1e-15?Eem:1e-15));
        H.hit=1; H.r=rhit; H.phi=fmod(phit+1000.0*M_PI*2, 2*M_PI);
        H.g = g>0?g:0; H.emiss = pow(rhit, -emiss_p);
        return H;
      }
    }
    th_prev=x[2]; for(int i=0;i<4;i++){ x_prev[i]=x[i]; v_prev[i]=v[i]; }
  }
  return H;
}

// ASCII ramp (dark -> bright)
static char RAMP[] = " `.,-:'_;~*/\\^\"<>!=()+?{}|[]#%$&@";
static double gamma_c = 0.30; // display gamma

int main(int argc, char** argv){
  // Observer-only CLI controls: inc_deg, FOVx (deg), robs
  if(argc > 1){ double val = atof(argv[1]); if(val > 0.0 && val < 89.0) inc_deg = val; }
  if(argc > 2){ double fdeg = atof(argv[2]); if(fdeg > 5.0 && fdeg < 170.0) FOVx = fdeg * M_PI/180.0; }
  if(argc > 3){ double r = atof(argv[3]); if(r > 10.0 && r < 2000.0) robs = r; }
  theta_obs = M_PI/2.0 - (inc_deg*M_PI/180.0);
  FOVy = FOVx * ((double)HEIGHT/W);
  static Hit map[HEIGHT][W];
  // precompute lens map
  for(int y=0;y<HEIGHT;y++){
    for(int x=0;x<W;x++) map[y][x]=trace_pixel(x,y);
  }
  // Precompute a fixed normalization scale from the base (no hotspot) field
  double norm_scale = 1e-12;
  for(int y=0;y<HEIGHT;y++){
    for(int x=0;x<W;x++) if(map[y][x].hit){
      double g3 = pow(map[y][x].g,3.0);
      double base = map[y][x].emiss * g3 * ring_mul(map[y][x].r);
      if(base > norm_scale) norm_scale = base;
    }
  }
  // Animation: rotate three hotspots around the disk
  printf("\x1b[2J"); // clear once
  double phase = 0.0;
  const double dphase = 2*M_PI / 180.0; // 180 frames per full rotation
  for(;;){
    printf("\x1b[H"); // move cursor home
    // accumulate intensity per frame (fixed normalization baseline)
    static double I[HEIGHT][W];
    for(int y=0;y<HEIGHT;y++) for(int x=0;x<W;x++){
      double val=0.0;
      if(map[y][x].hit){
        double g3 = pow(map[y][x].g,3.0);
        double base = map[y][x].emiss * g3 * ring_mul(map[y][x].r);
        // Multiply by rotating hotspot field (clockwise), fixed radius and size
        val = base * hotspots_mul(map[y][x].r, map[y][x].phi, phase);
      }
      I[y][x]=val;
    }
    // draw frame
    int ramp_len = (int)sizeof(RAMP) - 1; // exclude null
    for(int y=0;y<HEIGHT;y++){
      for(int x=0;x<W;x++){
        double v = I[y][x] / norm_scale;
        if(v<0.0) v=0.0; if(v>1.0) v=1.0;
        double q = pow(v, gamma_c);
        int idx_base = (int)(q * (ramp_len-1));
        if(idx_base<0) idx_base=0; if(idx_base>ramp_len-1) idx_base=ramp_len-1;
        char c = RAMP[idx_base];
        putchar(c);
      }
      putchar('\n');
    }
    fflush(stdout);
    usleep(40000); // ~25 fps
    phase += dphase;
    if(phase > 2*M_PI) phase -= 2*M_PI;
  }
  return 0;
}
