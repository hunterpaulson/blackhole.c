// bh_ascii.c — Schwarzschild ASCII disk (precompute lens map + animate)
// build: cc -O3 -lm bh_ascii.c -o bh_ascii
// run:   ./bh_ascii

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

// Params
const double rin = 6.0, rout = 40.0;
double robs = 50.0;                // observer radius (zoom out default)
double inc_deg = 5.0;              // default inclination (deg), can override via argv
double theta_obs = 0.0;            // set in main from inc_deg
const double phi_obs = 0.0;
double FOVx = 70.0*M_PI/180.0;     // wider FOV for zoomed-out default
double FOVy = 0.0;                 // set in main from FOVx
const double emiss_p = 2.0;

// Hotspot
const double hs_sigma = 0.18;
const double hs_gain  = 0.0;  // disable hotspot to keep mass uniform
// Uniformity + precomputed rotating ring textures
#define NBINS 64
static double bin_mean_g3[NBINS];
static int    bin_counts[NBINS];
static inline int r_to_bin(double r){
  if(r<rin) r=rin; if(r>rout) r=rout;
  int b = (int)((r - rin) / (rout - rin) * (NBINS-1) + 0.5);
  if(b<0) b=0; if(b>=NBINS) b=NBINS-1; return b;
}
const double uniform_alpha = 0.8;   // 0 => fully uniform per-ring, 1 => original g^3
// Ring textures: per-bin azimuthal pattern that rotates one turn per loop
#define NPHI 256
static double ring_tex[NBINS][NPHI];
static double radial_tex[NBINS];
// Saturn-like radial band controls (CLI-overridable)
static int    ring_count = 8;     // number of bright bands across [rin,rout]
static double ring_fill  = 0.35;  // fraction of each band that is bright
static double ring_peak  = 2.0;   // multiplier for bright rings
static double ring_floor = 0.10;  // multiplier for gaps
static double edge_soft  = 0.02;  // softness of band edges (in band fraction)
static inline double frand(){ return rand()/(double)RAND_MAX; }
static void gen_ring_textures(void){
  for(int b=0;b<NBINS;b++){
    // initialize with random values
    for(int i=0;i<NPHI;i++) ring_tex[b][i] = frand();
    // smooth a few passes (circular)
    for(int pass=0; pass<2; pass++){
      double tmp[NPHI];
      for(int i=0;i<NPHI;i++){
        int iL=(i-1+NPHI)%NPHI, iR=(i+1)%NPHI;
        tmp[i] = (ring_tex[b][iL] + 2.0*ring_tex[b][i] + ring_tex[b][iR]) / 4.0;
      }
      for(int i=0;i<NPHI;i++) ring_tex[b][i]=tmp[i];
    }
    // normalize to [0,1]
    double mn=1e9,mx=-1e9;
    for(int i=0;i<NPHI;i++){ if(ring_tex[b][i]<mn) mn=ring_tex[b][i]; if(ring_tex[b][i]>mx) mx=ring_tex[b][i]; }
    double inv = (mx>mn)?(1.0/(mx-mn)):1.0;
    for(int i=0;i<NPHI;i++) ring_tex[b][i] = (ring_tex[b][i]-mn)*inv;
    // map to [floor, floor+gain] and add mild jitter (grain)
    const double tex_floor=0.20, tex_gain=1.40, jitter_amp=0.25;
    for(int i=0;i<NPHI;i++){
      double v = tex_floor + tex_gain * ring_tex[b][i] + jitter_amp * (frand() - 0.5);
      if(v < tex_floor) v = tex_floor;
      if(v > tex_floor + tex_gain) v = tex_floor + tex_gain;
      ring_tex[b][i] = v;
    }
  }
}
static inline double sample_ring_tex(int b, double phi){
  if(b<0) b=0; if(b>=NBINS) b=NBINS-1;
  double t = fmod(phi, 2*M_PI);
  if(t<0) t += 2*M_PI;
  double idx = t * (double)NPHI / (2.0*M_PI);
  int i0 = (int)floor(idx);
  int i1 = (i0+1)%NPHI;
  double f = idx - i0;
  return ring_tex[b][i0]*(1.0-f) + ring_tex[b][i1]*f;
}

// Radial rings (Saturn-like): alternating light/dark bands across radius
static void gen_radial_rings(void){
  // Gentle warp and jitter so bands vary a bit
  const double small_warp = 0.03;
  const double fill_jit   = 0.05;
  double off1 = 2*M_PI*frand();
  for(int b=0;b<NBINS;b++){
    double s = (double)b/(double)(NBINS-1); // 0..1 across [rin,rout]
    double sw = s + small_warp*sin(2*M_PI*3.0*s + off1);
    double pos = sw * (double)ring_count;
    double f = pos - floor(pos); // position within band [0,1)
    double fill = ring_fill + fill_jit*(frand()-0.5);
    if(fill < 0.05) fill = 0.05; if(fill > 0.95) fill = 0.95;
    // Smooth step from bright (f<fill) to dark (f>fill)
    double w = edge_soft + 1e-6;
    double t = 0.5 + 0.5 * tanh((fill - f)/w); // near 1 inside ring, near 0 in gap
    radial_tex[b] = ring_floor + t * (ring_peak - ring_floor);
  }
}

// Additional stylistic boost for the Doppler-bright approaching crescent
double cres_gain = 1.2;   // strength (CLI overrideable)
double cres_sigma = 0.55; // angular width in rad (CLI overrideable)

typedef struct { double r, phi, g, emiss; int bin; int hit; } Hit; // hit=1 if disk

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
        H.bin = r_to_bin(rhit);
        H.g = g>0?g:0; H.emiss = pow(rhit, -emiss_p);
        return H;
      }
    }
    th_prev=x[2]; for(int i=0;i<4;i++){ x_prev[i]=x[i]; v_prev[i]=v[i]; }
  }
  return H;
}

// High-resolution ASCII ramp for more shades
static char RAMP[] = " .'`^\",:;Il!i~+_-?][}{1)(|\\/tfjrxnuvczXYUJCLQ0OZmwqpdbkhao*#MW&8%B@$";
// gamma_c declared once above
static double norm_pct = 0.98; // percentile for normalization (CLI overrideable)
static double q_floor  = 0.05; // minimum brightness floor (CLI overrideable)
static double gamma_c = 0.40; // ASCII gamma (CLI overrideable)

int main(int argc, char** argv){
  if(argc > 1){
    double val = atof(argv[1]);
    if(val > 0.0 && val < 89.0) inc_deg = val;
  }
  if(argc > 2){
    double fdeg = atof(argv[2]);
    if(fdeg > 5.0 && fdeg < 170.0) FOVx = fdeg * M_PI/180.0;
  }
  if(argc > 3){
    double r = atof(argv[3]);
    if(r > 10.0 && r < 2000.0) robs = r;
  }
  int frames_per_cycle = 120;    // frames per full phase cycle (2π)
  int cycles = 0;                 // 0 => run forever
  if(argc > 4){
    int fpc = atoi(argv[4]);
    if(fpc >= 24 && fpc <= 2000) frames_per_cycle = fpc;
  }
  if(argc > 5){
    cycles = atoi(argv[5]);
  }
  if(argc > 6){
    cres_gain = atof(argv[6]);
  }
  if(argc > 7){
    cres_sigma = atof(argv[7]);
  }
  if(argc > 8){
    gamma_c = atof(argv[8]);
    if(gamma_c < 0.2) gamma_c = 0.2;
    if(gamma_c > 1.5) gamma_c = 1.5;
  }
  if(argc > 9){ ring_count = atoi(argv[9]); if(ring_count < 2) ring_count=2; if(ring_count>64) ring_count=64; }
  if(argc > 10){ ring_fill = atof(argv[10]); if(ring_fill<0.05) ring_fill=0.05; if(ring_fill>0.95) ring_fill=0.95; }
  if(argc > 11){ ring_peak = atof(argv[11]); if(ring_peak<1.0) ring_peak=1.0; if(ring_peak>5.0) ring_peak=5.0; }
  if(argc > 12){ ring_floor = atof(argv[12]); if(ring_floor<0.0) ring_floor=0.0; if(ring_floor>1.0) ring_floor=1.0; }
  if(argc > 13){ norm_pct = atof(argv[13]); if(norm_pct<0.50) norm_pct=0.50; if(norm_pct>0.999) norm_pct=0.999; }
  if(argc > 14){ q_floor = atof(argv[14]); if(q_floor<0.0) q_floor=0.0; if(q_floor>0.5) q_floor=0.5; }
  theta_obs = M_PI/2.0 - (inc_deg*M_PI/180.0);
  FOVy = FOVx * ((double)HEIGHT/W);
  static Hit map[HEIGHT][W];
  // precompute
  for(int y=0;y<HEIGHT;y++){
    for(int x=0;x<W;x++) map[y][x]=trace_pixel(x,y);
  }
  // precompute per-ring average g^3 (for uniformization)
  memset(bin_mean_g3, 0, sizeof(bin_mean_g3));
  memset(bin_counts, 0, sizeof(bin_counts));
  for(int y=0;y<HEIGHT;y++){
    for(int x=0;x<W;x++) if(map[y][x].hit){
      int b = map[y][x].bin;
      bin_mean_g3[b] += pow(map[y][x].g, 3.0);
      bin_counts[b]  += 1;
    }
  }
  for(int b=0;b<NBINS;b++){
    if(bin_counts[b]>0) bin_mean_g3[b] /= (double)bin_counts[b];
    else bin_mean_g3[b] = 0.0;
  }
  // generate ring textures once
  gen_ring_textures();
  gen_radial_rings();
  // animate a few seconds
  // clear screen once
  printf("\x1b[2J");
  // continuous animation with explicit phase step
  double phase = 0.0;
  double dphase = 2*M_PI / (double)frames_per_cycle;
  int k_in_cycle = 0;
  int cycle_count = 0;
  for(;;){
    // move cursor home for in-place animation
    printf("\x1b[H");
    // per-frame intensities
    double Imax=1e-12;
    static double I[HEIGHT][W];
    static double vals[HEIGHT*W];
    int idxv=0;
    for(int y=0;y<HEIGHT;y++) for(int x=0;x<W;x++){
      double val=0.0;
      if(map[y][x].hit){
        int b = map[y][x].bin;
        double g3 = pow(map[y][x].g,3.0);
        double g3u = (1.0 - uniform_alpha) * (bin_mean_g3[b] > 0 ? bin_mean_g3[b] : g3) + uniform_alpha * g3;
        double base = map[y][x].emiss * g3u * radial_tex[b];
        // Rotate a precomputed ring texture: one full turn per cycle
        double phi = map[y][x].phi;
        // For CCW orbital motion, the approaching side is on the left (phi≈π)
        // Rotate texture CCW by adding phase
        double phi_rot = phi + phase;
        double tex = sample_ring_tex(b, phi_rot);
        // Physically anchored extra boost: emphasize pixels with g^3 above ring mean
        double ring_mean = bin_mean_g3[b];
        double g3_raw = g3; // pow(map[y][x].g,3.0) already computed
        double relat = (ring_mean>1e-12) ? (g3_raw / ring_mean) : 1.0;
        double dop_boost = 1.0 + cres_gain * fmax(0.0, relat - 1.0);
        val = base * tex * dop_boost;
      }
      I[y][x]=val; if(val>Imax) Imax=val; vals[idxv++]=val;
    }
    // percentile-based normalization to use more of the ramp
    int N = idxv;
    for(int i=1;i<N;i++){ // insertion sort (small arrays)
      double key = vals[i]; int j=i-1; while(j>=0 && vals[j]>key){ vals[j+1]=vals[j]; j--; } vals[j+1]=key;
    }
    double denom = Imax;
    if(N>0){
      int pi = (int)(norm_pct*N); if(pi<1) pi=1; if(pi>=N) pi=N-1;
      double thr = vals[pi]; if(thr > 1e-12) denom = thr;
    }
    // draw
    for(int y=0;y<HEIGHT;y++){
      for(int x=0;x<W;x++){
        double v = I[y][x]/denom; if(v>1.0) v=1.0; if(v<0.0) v=0.0;
        double q = pow(v, gamma_c);
        if(q_floor>0.0) q = q_floor + (1.0 - q_floor)*q;
        int idx = (int)(q * ((int)sizeof(RAMP)-2));
        if(idx<0) idx=0; if(idx>(int)sizeof(RAMP)-2) idx=sizeof(RAMP)-2;
        putchar(RAMP[idx]);
      }
      putchar('\n');
    }
    fflush(stdout);
    usleep(40000); // ~25 fps pacing
    // advance phase, rollover at 2π
    phase += dphase;
    k_in_cycle++;
    if(k_in_cycle >= frames_per_cycle){
      k_in_cycle = 0;
      phase -= 2*M_PI;
      cycle_count++;
      if(cycles > 0 && cycle_count >= cycles) break;
    }
  }
  return 0;
}
