// bh_ascii.c — Schwarzschild ASCII disk (precompute lens map + animate)
// build: cc -O3 -lm bh_ascii.c -o bh_ascii
// run:   ./bh_ascii

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#define W 120
#define HEIGHT 36

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
const double hs_gain  = 1.7;
// Inner rim brightness boost near the disk inner edge (r ~ rin)
const double rim_sigma = 2.0;   // width in M
const double rim_gain  = 2.0;   // additional brightness multiplier strength

typedef struct { double r, phi, g, emiss, rmin; int hit; } Hit; // hit=1 if disk

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
  Hit H; H.hit=0; H.rmin=0.0;
  double x[4], v[4]; pix_ray(px,py,x,v);
  double th_prev=x[2], x_prev[4], v_prev[4];
  for(int i=0;i<4;i++){ x_prev[i]=x[i]; v_prev[i]=v[i]; }
  const double h0=0.5, rh=2.0*Mbh;
  double r_min = x[1];
  for(int step=0; step<5000; ++step){
    double h=h0;
    if(x[1]<10.0) h=0.25*h0;
    if(x[1]<6.0)  h=0.125*h0;
    rk4(x,v,h);
    if(x[1] < r_min) r_min = x[1];
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
        H.rmin = r_min;
        H.g = g>0?g:0; H.emiss = pow(rhit, -emiss_p);
        return H;
      }
    }
    th_prev=x[2]; for(int i=0;i<4;i++){ x_prev[i]=x[i]; v_prev[i]=v[i]; }
  }
  return H;
}

static char RAMP[] = " .,:;-~=*%#@";

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
  theta_obs = M_PI/2.0 - (inc_deg*M_PI/180.0);
  FOVy = FOVx * ((double)HEIGHT/W);
  static Hit map[HEIGHT][W];
  // precompute
  for(int y=0;y<HEIGHT;y++){
    for(int x=0;x<W;x++) map[y][x]=trace_pixel(x,y);
  }
  // animate a few seconds
  const int frames=72;
  // clear screen once
  printf("\x1b[2J");
  for(int k=0;k<frames;k++){
    // move cursor home for in-place animation
    printf("\x1b[H");
    double phase = 2*M_PI * (k/(double)frames);
    // per-frame max
    double Imax=1e-12;
    static double I[HEIGHT][W];
    for(int y=0;y<HEIGHT;y++) for(int x=0;x<W;x++){
      double val=0.0;
      if(map[y][x].hit){
        double dphi = atan2(sin(map[y][x].phi - phase), cos(map[y][x].phi - phase));
        double hs = 1.0 + hs_gain*exp(-(dphi*dphi)/(2*hs_sigma*hs_sigma));
        double dr = map[y][x].r - rin;
        double rim = 1.0 + rim_gain * exp(-(dr*dr)/(2*rim_sigma*rim_sigma));
        val = map[y][x].emiss * pow(map[y][x].g,3.0) * hs * rim;
      }
      I[y][x]=val; if(val>Imax) Imax=val;
    }
    // draw
    for(int y=0;y<HEIGHT;y++){
      for(int x=0;x<W;x++){
        double q = sqrt(I[y][x]/Imax); // gamma=0.5
        int idx = (int)(q * (sizeof(RAMP)-2));
        if(idx<0) idx=0; if(idx>(int)sizeof(RAMP)-2) idx=sizeof(RAMP)-2;
        putchar(RAMP[idx]);
      }
      putchar('\n');
    }
    fflush(stdout);
    usleep(40000); // ~25 fps pacing
  }
  return 0;
}
