/*
CSCI 480
Assignment 3 Raytracer

Name: Tong Wang
*/

#include <stdlib.h>
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#include <GLUT/glut.h>
#include <pic.h>
#include <string.h>
#include <iostream>
#include <vector>
#include <math.h>
using namespace std;

#define MAX_TRIANGLES 2000
#define MAX_SPHERES 10
#define MAX_LIGHTS 10

char *filename=0;
bool aa_flag = false;

//different display modes
#define MODE_DISPLAY 1
#define MODE_JPEG 2
int mode=MODE_DISPLAY;

//you may want to make these smaller for debugging purposes
#define WIDTH 640
#define HEIGHT 480

//the field of view of the camera
#define fov 60.0

#define PI 3.14159265
#define ASP_RATIO (double) WIDTH/HEIGHT
#define TANGENT_VAL (double) tan(fov*PI/180.0/2)
/* used for dealing with numeric issues. 1e-7 is chosen since it's experimented
   to be maximum small value that stil holds good rendering result */
#define APPROX_ZERO 1e-7
#define NEG_INF -1e7

unsigned char buffer[HEIGHT][WIDTH][3];

struct point {
   double x;
   double y;
   double z;
};

struct Vertex
{
  double position[3];
  double color_diffuse[3];
  double color_specular[3];
  double normal[3];
  double shininess;
};

typedef struct _Triangle
{
  struct Vertex v[3];
} Triangle;

typedef struct _Sphere
{
  double position[3];
  double color_diffuse[3];
  double color_specular[3];
  double shininess;
  double radius;
} Sphere;

typedef struct _Light
{
  double position[3];
  double color[3];
} Light;

struct color {
  double r,g,b;
  // constructor
  color() {}
  color(double R, double G, double B) : r(R), g(G), b(B) {}

  color& operator += (const color &c) {
    r += c.r;
    g += c.g;
    b += c.b;
    return *this;
  }

  color& operator *= (double d) {
    r *= d;
    g *= d;
    b *= d;
    return *this;
  }
};

class Ray {
  public: 
    point cam_origin;
    vector<double> direction;

    /*constructor*/
    Ray(){};
    Ray(const point &co, const vector<double> &dir) : cam_origin(co), direction(dir) {};

    /* intersection test functions */
    bool intersectSphere(const Sphere &sphere, point &intersection);
    bool intersectTriangle(const Triangle &triangle, point &intersection, 
                           vector<double> &t_surf_nor, color &t_diffuse,
                           color &t_specular, float &t_shiness, bool shadow_ray_flag);
};

/* rays_arr is 2d array (width*height) that stores generated rays */
Ray** rays_arr;
/* color_arr is 2D array (width*height) that stores final pixel colors*/
color** color_arr;

Triangle triangles[MAX_TRIANGLES];
Sphere spheres[MAX_SPHERES];
Light lights[MAX_LIGHTS];
double ambient_light[3];

int num_triangles=0;
int num_spheres=0;
int num_lights=0;

/* clamp color to be within range [0,1] */
color clamp_color(color c) {
  if (c.r > 1.0) c.r = 1.0;
  if (c.g > 1.0) c.g = 1.0;
  if (c.b > 1.0) c.b = 1.0;

  if (c.r < APPROX_ZERO) c.r = 0.0;
  if (c.g < APPROX_ZERO) c.g = 0.0;
  if (c.b < APPROX_ZERO) c.b = 0.0;
  return c;
}

vector<double> normalizeVec (vector<double> v) {
  vector<double> retVec;
  double x = v.at(0);
  double y = v.at(1);
  double z = v.at(2);
  double vecLength = sqrt(x*x+y*y+z*z);
  retVec.push_back(x/vecLength);
  retVec.push_back(y/vecLength);
  retVec.push_back(z/vecLength);
  return retVec;
}

/* adds a point and a vector, returns point */
point ptAddvec (point a, vector<double> b) {
  point c;
  c.x = a.x + b.at(0);
  c.y = a.y + b.at(1);
  c.z = a.z + b.at(2);
  return c;
}

/* minus two pts, returns vector */
vector<double> ptMinuspt (point a, point b) {
  vector<double> retval;
  retval.push_back(a.x-b.x);
  retval.push_back(a.y-b.y);
  retval.push_back(a.z-b.z);
  return retval;
}

/* adds two vectors, returns vector */
vector<double> vecAddvec (vector<double> a, vector<double> b) {
  vector<double> c;
  c.push_back(a.at(0) + b.at(0));
  c.push_back(a.at(1) + b.at(1));
  c.push_back(a.at(2) + b.at(2));
  return c;
}

/* multiply input vector b by scaler s */
vector<double> scalerMul (float s, vector<double> b) {
  vector<double> c;
  c.push_back(s * b.at(0));
  c.push_back(s * b.at(1));
  c.push_back(s * b.at(2));
  return c;
}

/* neagate a vector */
vector<double> negVec (vector<double> a) {
  vector<double> b;
  b.push_back( - a.at(0));
  b.push_back( - a.at(1));
  b.push_back( - a.at(2));
  return b;
}

/* calculate cross product of two vectors, returned vector is normalized */
vector<double> cross(vector<double> a, vector<double> b) {
  vector<double> retVec;
  double x,y,z;
  x = a.at(1)*b.at(2) - a.at(2)*b.at(1);
  y = a.at(2)*b.at(0) - a.at(0)*b.at(2);
  z = a.at(0)*b.at(1) - a.at(1)*b.at(0);
  retVec.push_back(x);
  retVec.push_back(y);
  retVec.push_back(z);
  return retVec;
}

/* calculate vector length */
double vecLength(vector<double> v) {
  return sqrt(v.at(0)*v.at(0) + v.at(1)*v.at(1) + v.at(2)*v.at(2));
}

/* take dot product btwn two vectors */
double dotProd(vector<double> a, vector<double> b) {
  return a.at(0)*b.at(0) + a.at(1)*b.at(1) + a.at(2)*b.at(2);
}

/* check if barycentric coordinates are valid */
bool checkTrigPRange(double value) {
  if (value < APPROX_ZERO) return false;
  if (value > 1) return false;
  return true;
}

/* transform array values to a point */
point getPtFromArr (double const arr[3]) {
  point retPt;
  retPt.x = arr[0];
  retPt.y = arr[1];
  retPt.z = arr[2];
  return retPt;
}

/* transform array values to a vector */
vector<double> getVecFromArr (double const arr[3]) {
  vector<double> retVec;
  retVec.push_back(arr[0]);
  retVec.push_back(arr[1]);
  retVec.push_back(arr[2]);
  return retVec;
}

void print_vector(vector<double> v) {
    for (auto i = v.begin(); i != v.end(); ++i)
      std::cout << *i << ' ';
    cout<<"\n";
}

/*initialize camera position (default value are 0's, which is good) */
point camera_pos;

void generate_rays(int w, int h) {
  /* create dynamic 2d array of rays */
  rays_arr = new Ray*[w];
  for(int i = 0; i < w; ++i)
    rays_arr[i] = new Ray[h];

  /* calculate boundries */
  double xmin = -ASP_RATIO*TANGENT_VAL;
  double xmax = +ASP_RATIO*TANGENT_VAL;
  double ymin = -TANGENT_VAL;
  double ymax = +TANGENT_VAL;

  /* calculate, dx, dy, which are step sizes in between rays */
  double dx,dy;
  /* assuming rays pass through center of pixels, hence -1 when dividing
     this way, ray at four corners can be properly set */
  dx = (xmax-xmin)/(w-1);
  dy = (ymax-ymin)/(h-1);

  vector<double> direction;
  /* generate rays array of viewplane size */
  for (int y=0; y<h; y++) {
    for (int x=0; x<w; x++) {
      /* get ray direction */
      direction.clear();
      direction.push_back(xmin+x*dx);
      direction.push_back(ymin+y*dy);
      direction.push_back(-1.0);
      direction = normalizeVec(direction);
      /* create ray */
      rays_arr[x][y] = Ray(camera_pos, direction);
    }
  }
}

/* sphere intersection test */
bool Ray::intersectSphere(const Sphere &sphere, point &intersection) {
  /* t0, t1 are two solutions two equation. ret=min(t0,t1), the closest intersection */
  double retT, t0, t1;

  /* get center of sphere */
  point center;
  center = getPtFromArr(sphere.position);

  /* declarations for root calculation */
  vector<double> direc;
  direc = direction; // of given ray
  double r,a,b,c;
  r = sphere.radius;

  /* calculate coefficient for the quadratic equation */
  vector<double> oriToCenterVec;
  oriToCenterVec = ptMinuspt(cam_origin, center);
  a = 1; // since |direction|=1, a=dotProd(direction, direction)=1;
  b = 2*dotProd(direc, oriToCenterVec);
  c = dotProd(oriToCenterVec, oriToCenterVec) - r*r;

  /* calculate t0, t1 */
  double discriminant = b*b-4*a*c;
  if (discriminant < 0) {             // NO SOLUTION, RETURN
    retT = -1; 
    return false;
  } else if (discriminant == 0) {     // ONE SOLUTION, SO SIMPLIFIED CALCULTION
    retT = 0.5*(-b);
  } else {                            // TWO SOLUTIONS, USE THE MIN ONE
    t0 = 0.5*(-b+sqrt(discriminant));
    t1 = 0.5*(-b-sqrt(discriminant));
    retT = min(t0,t1);
  }

  /* set closest intersection point if valid solution found */
  if(retT > 0) {
    intersection = ptAddvec(cam_origin, scalerMul(retT, direc));
    return true;
  }
  return false;
}

/* triangle intersection test */
bool Ray::intersectTriangle(const Triangle &triangle, point &intersection, vector<double> &t_surf_nor, 
                       color &t_diffuse, color &t_specular, float &t_shiness, bool shadow_ray_flag) {
  /* intersection scaler t */
  double retT;
  /* get triangle points */
  point pA, pB, pC;
  pA = getPtFromArr(triangle.v[0].position);
  pB = getPtFromArr(triangle.v[1].position);
  pC = getPtFromArr(triangle.v[2].position);

  /* calculate normal vector of triangle */
  vector<double> N = normalizeVec(cross(ptMinuspt(pB,pA), ptMinuspt(pC, pA)));

  /* parallel test */
  if (abs(dotProd(N, direction)) < APPROX_ZERO) return false;

  /* calculate retT and test if triangle is behind */
  retT = -(dotProd(N, ptMinuspt(cam_origin, pB))) / (dotProd(N, direction));
  if (retT < APPROX_ZERO) return false;

  /* at this point, retT is valid, calculate intersection to triangle plane*/
  intersection = ptAddvec(cam_origin, scalerMul(retT, direction));

  /* lastly test inside or outside triangle. reduce floating pt ops by *2 for all */
  double totalArea = vecLength(cross(ptMinuspt(pB,pA), ptMinuspt(pC, pA)));
  double a = vecLength(cross(ptMinuspt(pB,intersection), ptMinuspt(pC, intersection)))/totalArea;
  double b = vecLength(cross(ptMinuspt(pC,intersection), ptMinuspt(pA, intersection)))/totalArea;
  double c = vecLength(cross(ptMinuspt(pA,intersection), ptMinuspt(pB, intersection)))/totalArea;
  /* setting c=1-a-b and only check range of a,b,c doesn't work well due to numerical
     impresitions, which produces bigger than expect shadows. So, a better way is to 
     calculate all a,b,c, and check if they indeed add up to one with some room of impresition */
  double abcTotalCheck = fabs(1-a-b-c);

  /* check if a,b,c in range [0,1] and if they add up to one. If all do, then intersection */
  if (checkTrigPRange(a) && checkTrigPRange(b) && checkTrigPRange(c)
                         && abcTotalCheck < APPROX_ZERO) {
    /* only interpolate normal and color if the ray is not shadow ray */
    if (shadow_ray_flag == false) {
      /* interpolate surface normal with barycentric coord
        interp normal = a*pAnormal + b*pBnormal + c* pCnormal */
      vector<double> pAnormal, pBnormal, pCnormal;
      pAnormal = getVecFromArr(triangle.v[0].normal);
      pBnormal = getVecFromArr(triangle.v[1].normal);
      pCnormal = getVecFromArr(triangle.v[2].normal);

      vector<double> interpPA, interpPB, interpPC;
      interpPA = scalerMul(a, pAnormal);
      interpPB = scalerMul(b, pBnormal);
      interpPC = scalerMul(c, pCnormal);

      t_surf_nor = normalizeVec(vecAddvec(vecAddvec(interpPA, interpPB), interpPC));

      /* calculate color components with barycentric coord */
      vector<double> pAdiffuse = getVecFromArr(triangle.v[0].color_diffuse);
      vector<double> pBdiffuse = getVecFromArr(triangle.v[1].color_diffuse);
      vector<double> pCdiffuse = getVecFromArr(triangle.v[2].color_diffuse);

      vector<double> pAspecular = getVecFromArr(triangle.v[0].color_specular);
      vector<double> pBspecular = getVecFromArr(triangle.v[1].color_specular);
      vector<double> pCspecular = getVecFromArr(triangle.v[2].color_specular);


      vector<double> interpPAdif, interpPBdif, interpPCdif,
                     interpPAspe, interpPBspe, interpPCspe;

      interpPAdif = scalerMul(a, pAdiffuse);
      interpPBdif = scalerMul(b, pBdiffuse);
      interpPCdif = scalerMul(c, pCdiffuse);

      interpPAspe = scalerMul(a, pAspecular);
      interpPBspe = scalerMul(b, pBspecular);
      interpPCspe = scalerMul(c, pCspecular);

      vector<double> final_dif, final_spe;
      final_dif = vecAddvec(vecAddvec(interpPAdif, interpPBdif) , interpPCdif);
      final_spe = vecAddvec(vecAddvec(interpPAspe, interpPBspe), interpPCspe);

      t_diffuse = color(final_dif.at(0), final_dif.at(1), final_dif.at(2));
      t_specular = color(final_spe.at(0), final_spe.at(1), final_spe.at(2));
      t_shiness = a * triangle.v[0].shininess +
                  b * triangle.v[1].shininess +
                  c * triangle.v[2].shininess;
    }
    return true;
  }
  return false;
}

/* calculate lighting given a light, intersection point, surface normal,
   object diffuse color, and object shiness information */
color phongModel(Light light, point intersection, vector<double> surf_normal, color diffuse, color specular, float shiness) {
  color resultColor;

  /* get light information */
  color lightColor = color(light.color[0], light.color[1], light.color[2]);
  point lightPos;
  lightPos = getPtFromArr(light.position);

  /* calculate L,N,V,R */
  vector<double> L = normalizeVec(ptMinuspt(lightPos, intersection));
  vector<double> N = surf_normal;
  vector<double> V = normalizeVec(ptMinuspt(camera_pos, intersection));
  vector<double> R = vecAddvec(scalerMul(2*dotProd(L, N), N), negVec(L));

  /* calculate dot products and clamp to 0 if result is negative*/
  double LdotN = max(0.0, dotProd(L, N));
  double RdotV = max(0.0, dotProd(R, V));


  /* calculate and set final result color */
  resultColor.r = lightColor.r * (diffuse.r * LdotN + specular.r * pow(RdotV, shiness));
  resultColor.g = lightColor.g * (diffuse.g * LdotN + specular.g * pow(RdotV, shiness));
  resultColor.b = lightColor.b * (diffuse.b * LdotN + specular.b * pow(RdotV, shiness));

  return resultColor;
}

/* given a intersection point and current sphere index, cast shadow ray
   If shadow ray occluded by some object, color of pix should be black.
   Otherwise, use surface normal, diffuse, specular infomation to call
   Phong model to determine pixel color. */
color sphere_shad_ray_cast(point intersection, int currSphereInd, 
                           vector<double> surface_normal, color diffuse, 
                           color specular, float shiness) {
  bool shadow_ray_flag = true;
  bool visible_flag;
  color resultColor;
  point lightPos;
  vector<double> vecInterToLt;
  vector<double> vecInterToLtDirec;
  vector<double> vecShadInterToInter;
  /* define shadow intersection point */
  point shadowIntersection;

  vector<double> trig_surf_normal;
  trig_surf_normal.push_back(-1); 
  trig_surf_normal.push_back(-1); 
  trig_surf_normal.push_back(-1); 

  /* resultColor is used for all ray-obj intersections */
  resultColor = color(0,0,0);
  resultColor += color(ambient_light[0], ambient_light[1], ambient_light[2]);
  for (int i = 0; i < num_lights; i++) {
    /* set visible flag to true for each new light source */
    visible_flag = true;

    /* get light position and fire a shadow ray */
    lightPos = getPtFromArr(lights[i].position);
    vecInterToLt = ptMinuspt(lightPos, intersection);
    vecInterToLtDirec = normalizeVec(vecInterToLt);
    Ray currShadowRay (intersection, vecInterToLtDirec);

    color t_diffuse = color(0,0,0);
    color t_specular = color(0,0,0);
    float t_shiness = 0;

    /* for each shphere other than current sphere, check if any of them blocked shadow ray */
    for (int s = 0; s < num_spheres; s++) {
      /* shadow ray blocked by some sphere, go to next light */
      if(currShadowRay.intersectSphere(spheres[s], shadowIntersection) && s != currSphereInd) {
        vecShadInterToInter = ptMinuspt(shadowIntersection, intersection);
        /* additional check to avoid the situation that intersection is found 
           further away than light which causes incorrect shadows */
        if (vecLength(vecShadInterToInter) < vecLength(vecInterToLt)) visible_flag = false;
        break; 
      }
    }

    for (int t = 0; t < num_triangles; t++) {
      /* shadow ray blocked by some triangle, go to next light */
      if(currShadowRay.intersectTriangle(triangles[t], shadowIntersection, trig_surf_normal,
                                         t_diffuse, t_specular, t_shiness, shadow_ray_flag)) {
        vecShadInterToInter = ptMinuspt(shadowIntersection, intersection);
        /* additional check to avoid the situation that intersection is found 
           further away than light which causes incorrect shadows */
        if (vecLength(vecShadInterToInter) < vecLength(vecInterToLt)) visible_flag = false;
        break; 
      }
    }

    if (visible_flag) {
      resultColor+=phongModel(lights[i], intersection, surface_normal, diffuse, specular, shiness);
    }
  }
  return resultColor;
}

/* given a intersection point and current triangle index, cast shadow ray
   If shadow ray occluded by some object, color of pix should be black.
   Otherwise, use surface normal, diffuse, specular infomation to call
   Phong model to determine pixel color. */
color triangle_shad_ray_cast(point intersection, int currTriInd, 
                           vector<double> surface_normal, color diffuse, 
                           color specular, float shiness) {
  bool shadow_ray_flag = true;
  bool visible_flag;
  color resultColor;
  point lightPos;
  vector<double> vecInterToLt;
  vector<double> vecInterToLtDirec;
  vector<double> vecShadInterToInter;
  /* define shadow intersection point */
  point shadowIntersection;

  /* delare dummy components for fire shadow ray */
  vector<double> dum_normal;
  dum_normal.push_back(-1); 
  dum_normal.push_back(-1); 
  dum_normal.push_back(-1);
  color dum_dif = color(0,0,0);
  color dum_spe = color(0,0,0);
  float dum_shi = 0;

  /* resultColor is used for all ray-obj intersections */
  resultColor = color(0,0,0);
  resultColor += color(ambient_light[0], ambient_light[1], ambient_light[2]);
  for (int i = 0; i < num_lights; i++) {
    /* set visible flag to true for each new light source */
    visible_flag = true;

    /* get light position and fire a shadow ray */
    lightPos = getPtFromArr(lights[i].position);

    /* fire shadow ray */
    vecInterToLt = ptMinuspt(lightPos, intersection);
    vecInterToLtDirec = normalizeVec(vecInterToLt);
    Ray currShadowRay (intersection, vecInterToLtDirec);

    /* for each shphere other than current sphere, check if any of them blocked shadow ray */
    for (int s = 0; s < num_spheres; s++) {
      /* shadow ray blocked by some sphere, go to next light */
      if(currShadowRay.intersectSphere(spheres[s], shadowIntersection)) {
        vecShadInterToInter = ptMinuspt(shadowIntersection, intersection);
        /* additional check to avoid the situation that intersection is found 
           further away than light which causes incorrect shadows */
        if (vecLength(vecShadInterToInter) < vecLength(vecInterToLt)) visible_flag = false;
        break; 
      }
    }

    /* for each triangle, check if any of them blocked shadow ray */
    for (int t = 0; t < num_triangles; t++) {
      /* shadow ray blocked by some triangle, go to next light */
      if(currShadowRay.intersectTriangle(triangles[t], shadowIntersection, dum_normal,
                                         dum_dif, dum_spe, dum_shi, shadow_ray_flag) 
                                         && t != currTriInd) {
        vecShadInterToInter = ptMinuspt(shadowIntersection, intersection);
        /* additional check to avoid the situation that intersection is found 
           further away than light which causes incorrect shadows */
        if (vecLength(vecShadInterToInter) < vecLength(vecInterToLt)) visible_flag = false;
        break; 
      }
    }

    if (visible_flag) {
      resultColor+=phongModel(lights[i], intersection, surface_normal, diffuse, specular, shiness);
    }
  }
  
  return resultColor;
}

/* trace all rays and update 2D color array with ray tracing resulting colors */
void ray_tracer (bool aa_flag) {
  bool shadow_ray_flag = false;
  /* initialize color array */
  color_arr = new color*[WIDTH];
  for(int i = 0; i < WIDTH; ++i)
    color_arr[i] = new color[HEIGHT];

  /* make background color to be white */
  for (int y = 0; y < HEIGHT; y++) {
    for (int x = 0; x< WIDTH; x++) {
      color_arr[x][y] = color(1,1,1);
    }
  }

  /* define initial intersection point */
  point intersection;
  point center; // of sphere

  /* for sphere Phong */
  vector<double> surf_normal;
  color diffuse, specular;
  double shiness;

  /* for triangle Phong. values will be changed with calls by reference */
  color t_diffuse = color(0,0,0);
  color t_specular = color(0,0,0);
  float t_shiness = 0;
  vector<double> t_surf_normal;
  surf_normal.push_back(-10);
  surf_normal.push_back(-10);
  surf_normal.push_back(-10);

  /* zmin is used to test closest object */
  double zmin;

  /* no anti aliasing requested */
  if (aa_flag == false) {
    for (int y = 0; y < HEIGHT; y++) {
      for (int x = 0; x< WIDTH; x++) {
        /* set zmin for current ray to be negative infinite since camera facing negtive z */
        zmin = NEG_INF;
        /* for each sphere, if current ray intersects it and it is closest, find color for it */
        for (int s = 0; s < num_spheres; s++) {
          /* get current sphere color info */
          center = getPtFromArr(spheres[s].position);
          diffuse = color(spheres[s].color_diffuse[0], spheres[s].color_diffuse[1], spheres[s].color_diffuse[2]);
          specular = color(spheres[s].color_specular[0], spheres[s].color_specular[1], spheres[s].color_specular[2]);
          shiness = spheres[s].shininess;
          /* check if 1) ray-obj intersection, 2) closest*/
          if (rays_arr[x][y].intersectSphere(spheres[s],intersection)&&intersection.z>zmin) {
            /* calculate surface normal for current intersection */
            surf_normal = normalizeVec(scalerMul(1/spheres[s].radius, ptMinuspt(intersection, center)));

            /* set current color pixel to black and add illumination values */
            color_arr[x][y] = color(0,0,0);
            color_arr[x][y] += sphere_shad_ray_cast(intersection, s, surf_normal, diffuse, specular, shiness);
            zmin = intersection.z;
          }
        }

        for (int t = 0; t < num_triangles; t++) {
          if (rays_arr[x][y].intersectTriangle(triangles[t], intersection, t_surf_normal, t_diffuse, t_specular, t_shiness, shadow_ray_flag)&&intersection.z>zmin) {
            color_arr[x][y] = color(0,0,0);
            color_arr[x][y] += triangle_shad_ray_cast(intersection, t, t_surf_normal, t_diffuse, t_specular, t_shiness);
            zmin = intersection.z;
          }
        }
      }
    }
  }
  /* anti aliasing requested */
  else {
    int tripleWIDTH = WIDTH * 3;
    int tripleHEIGHT = HEIGHT * 3;

    color** sumpersample_color_arr;
    /* initialize super sampled color array */
    sumpersample_color_arr = new color*[tripleWIDTH];
    for(int i = 0; i < tripleWIDTH; ++i)
      sumpersample_color_arr[i] = new color[tripleHEIGHT];
    /* make background color to be white */
    for (int y = 0; y < tripleHEIGHT; y++) {
      for (int x = 0; x< tripleWIDTH; x++) {
        sumpersample_color_arr[x][y] = color(1,1,1);
      }
    }

    /* generate supersampled ray-traced image */
    for (int y = 0; y < tripleHEIGHT; y++) {
      for (int x = 0; x< tripleWIDTH; x++) {
        /* set zmin for current ray to be negative infinite since camera facing negtive z */
        zmin = NEG_INF;
        /* for each sphere, if current ray intersects it and it is closest, find color for it */
        for (int s = 0; s < num_spheres; s++) {
          /* get current sphere color info */
          center = getPtFromArr(spheres[s].position);
          diffuse = color(spheres[s].color_diffuse[0], spheres[s].color_diffuse[1], spheres[s].color_diffuse[2]);
          specular = color(spheres[s].color_specular[0], spheres[s].color_specular[1], spheres[s].color_specular[2]);
          shiness = spheres[s].shininess;
          /* check if 1) ray-obj intersection, 2) closest*/
          if (rays_arr[x][y].intersectSphere(spheres[s],intersection)&&intersection.z>zmin) {
            /* calculate surface normal for current intersection */
            surf_normal = normalizeVec(scalerMul(1/spheres[s].radius, ptMinuspt(intersection, center)));

            /* set current color pixel to black and add illumination values */
            sumpersample_color_arr[x][y] = color(0,0,0);
            sumpersample_color_arr[x][y] += sphere_shad_ray_cast(intersection, s, surf_normal, diffuse, specular, shiness);
            zmin = intersection.z;
          }
        }

        for (int t = 0; t < num_triangles; t++) {
          if (rays_arr[x][y].intersectTriangle(triangles[t], intersection, t_surf_normal, t_diffuse, t_specular, t_shiness, shadow_ray_flag)&&intersection.z>zmin) {
            sumpersample_color_arr[x][y] = color(0,0,0);
            sumpersample_color_arr[x][y] += triangle_shad_ray_cast(intersection, t, t_surf_normal, t_diffuse, t_specular, t_shiness);
            zmin = intersection.z;
          }
        }
      }
    }

    /* final image result is by taking average of 9 adjacent pixels in the supersamoled image */
    for (int y = 0; y < HEIGHT; y++) {
      for (int x = 0; x< WIDTH; x++) {
        color_arr[x][y] = color(0,0,0);
        color_arr[x][y] += sumpersample_color_arr[3*x  ][3*y  ];
        color_arr[x][y] += sumpersample_color_arr[3*x+1][3*y  ];
        color_arr[x][y] += sumpersample_color_arr[3*x+2][3*y  ];
        color_arr[x][y] += sumpersample_color_arr[3*x  ][3*y+1];
        color_arr[x][y] += sumpersample_color_arr[3*x+1][3*y+1];
        color_arr[x][y] += sumpersample_color_arr[3*x+2][3*y+1];
        color_arr[x][y] += sumpersample_color_arr[3*x  ][3*y+2];
        color_arr[x][y] += sumpersample_color_arr[3*x+1][3*y+2];
        color_arr[x][y] += sumpersample_color_arr[3*x+2][3*y+2];
        color_arr[x][y] *= 0.111;
      }
    }
  }


  /* all colors of pixels are now finalized, clamp them to be within [0,1] */
  for (int y = 0; y < HEIGHT; y++) {
    for (int x = 0; x< WIDTH; x++) {
      color_arr[x][y] = clamp_color(color_arr[x][y]);
    }
  }
}



void plot_pixel_display(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plot_pixel_jpeg(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plot_pixel(int x,int y,unsigned char r,unsigned char g,unsigned char b);

//MODIFY THIS FUNCTION
void draw_scene() {
  unsigned int x,y;
  //simple output
  for(x=0; x<WIDTH; x++) {
    glPointSize(2.0);  
    glBegin(GL_POINTS);
    for(y=0;y < HEIGHT;y++) {
      // plot_pixel(x,y,x%256,y%256,(x+y)%256);
      plot_pixel(x,y, color_arr[x][y].r*255, color_arr[x][y].g*255, color_arr[x][y].b*255);
    }
    glEnd();
    glFlush();
  }
  printf("Done!\n"); fflush(stdout);
}

void plot_pixel_display(int x,int y,unsigned char r,unsigned char g,unsigned char b) {
  glColor3f(((double)r)/255.f,((double)g)/255.f,((double)b)/255.f);
  glVertex2i(x,y);
}

void plot_pixel_jpeg(int x,int y,unsigned char r,unsigned char g,unsigned char b) {
  buffer[HEIGHT-y-1][x][0]=r;
  buffer[HEIGHT-y-1][x][1]=g;
  buffer[HEIGHT-y-1][x][2]=b;
}

void plot_pixel(int x,int y,unsigned char r,unsigned char g, unsigned char b) {
  plot_pixel_display(x,y,r,g,b);
  if(mode == MODE_JPEG)
      plot_pixel_jpeg(x,y,r,g,b);
}

void save_jpg() {
  Pic *in = NULL;

  in = pic_alloc(640, 480, 3, NULL);
  printf("Saving JPEG file: %s\n", filename);

  memcpy(in->pix,buffer,3*WIDTH*HEIGHT);
  if (jpeg_write(filename, in))
    printf("File saved Successfully\n");
  else
    printf("Error in Saving\n");

  pic_free(in);
}

void parse_check(char *expected,char *found) {
  if(strcasecmp(expected,found)) {
    char error[100];
    printf("Expected '%s ' found '%s '\n",expected,found);
    printf("Parse error, abnormal abortion\n");
    exit(0);
  }
}

void parse_doubles(FILE*file, char *check, double p[3])
{
  char str[100];
  fscanf(file,"%s",str);
  parse_check(check,str);
  fscanf(file,"%lf %lf %lf",&p[0],&p[1],&p[2]);
  printf("%s %lf %lf %lf\n",check,p[0],p[1],p[2]);
}

void parse_rad(FILE*file,double *r)
{
  char str[100];
  fscanf(file,"%s",str);
  parse_check("rad:",str);
  fscanf(file,"%lf",r);
  printf("rad: %f\n",*r);
}

void parse_shi(FILE*file,double *shi)
{
  char s[100];
  fscanf(file,"%s",s);
  parse_check("shi:",s);
  fscanf(file,"%lf",shi);
  printf("shi: %f\n",*shi);
}

int loadScene(char *argv)
{
  FILE *file = fopen(argv,"r");
  int number_of_objects;
  char type[50];
  int i;
  Triangle t;
  Sphere s;
  Light l;
  fscanf(file,"%i",&number_of_objects);

  printf("number of objects: %i\n",number_of_objects);
  char str[200];

  parse_doubles(file,"amb:",ambient_light);

  for(i=0;i < number_of_objects;i++) {
      fscanf(file,"%s\n",type);
      printf("%s\n",type);
      if(strcasecmp(type,"triangle")==0) {
        printf("found triangle\n");
        int j;

        for(j=0;j < 3;j++) {
            parse_doubles(file,"pos:",t.v[j].position);
            parse_doubles(file,"nor:",t.v[j].normal);
            parse_doubles(file,"dif:",t.v[j].color_diffuse);
            parse_doubles(file,"spe:",t.v[j].color_specular);
            parse_shi(file,&t.v[j].shininess);
        }

        if(num_triangles == MAX_TRIANGLES) {
            printf("too many triangles, you should increase MAX_TRIANGLES!\n");
            exit(0);
        }
        triangles[num_triangles++] = t;
      }
      else if(strcasecmp(type,"sphere")==0) {
        printf("found sphere\n");

        parse_doubles(file,"pos:",s.position);
        parse_rad(file,&s.radius);
        parse_doubles(file,"dif:",s.color_diffuse);
        parse_doubles(file,"spe:",s.color_specular);
        parse_shi(file,&s.shininess);

        if(num_spheres == MAX_SPHERES) {
            printf("too many spheres, you should increase MAX_SPHERES!\n");
            exit(0);
        }
        spheres[num_spheres++] = s;
      }
      else if(strcasecmp(type,"light")==0) {
        printf("found light\n");
        parse_doubles(file,"pos:",l.position);
        parse_doubles(file,"col:",l.color);

        if(num_lights == MAX_LIGHTS) {
            printf("too many lights, you should increase MAX_LIGHTS!\n");
            exit(0);
        }
        lights[num_lights++] = l;
      }
      else {
        printf("unknown type in scene description:\n%s\n",type);
        exit(0);
      }
  }
  return 0;
}

void display(){}

void init() {
  glMatrixMode(GL_PROJECTION);
  glOrtho(0,WIDTH,0,HEIGHT,1,-1);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  glClearColor(0,0,0,0);
  glClear(GL_COLOR_BUFFER_BIT);
}

void idle() {
  //hack to make it only draw once
  static int once=0;
  if(!once) {
    draw_scene();
    if(mode == MODE_JPEG)
    save_jpg();
  }
  once=1;
}

/* reshape function is used to work around with macOS 10.15's bug that OpenGL window
   only uses left-bottom corner of the screen.By reshaping the window again, this bug
   can be fixed */
void reshape(int width, int height) {
  glutReshapeWindow(WIDTH, HEIGHT);
}

int main (int argc, char ** argv) {
  if (argc<2 || argc > 4)
  {  
    printf ("usage: %s <scenefile> [jpeg_name] [AA_ON/AA_OFF]\n", argv[0]);
    exit(0);
  }

  string aa_string;
  if (argc == 4) {
    aa_string = argv[3];
    if ((aa_string == "AA_ON")) {
      aa_flag = true;
    }
    mode = MODE_JPEG;
    filename = argv[2];
  }

  if(argc == 3) {
      mode = MODE_JPEG;
      filename = argv[2];
  }
  else if(argc == 2)
    mode = MODE_DISPLAY;

  glutInit(&argc,argv);
  loadScene(argv[1]);

  if (aa_flag == false) generate_rays(WIDTH, HEIGHT);
  else generate_rays(WIDTH*3,HEIGHT*3);
  ray_tracer(aa_flag);

  glutInitDisplayMode(GLUT_RGBA | GLUT_SINGLE);
  glutInitWindowPosition(0,0);
  glutInitWindowSize(WIDTH,HEIGHT);
  int window = glutCreateWindow("Ray Tracer");
  glutReshapeFunc(reshape);
  glutDisplayFunc(display);
  glutIdleFunc(idle);
  init();
  glutMainLoop();
}
