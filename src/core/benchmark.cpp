// benchmark.cpp — result gathering for the NURBS renderer (throwaway).
// Sweeps surface x C (refinement) x eps (flat_threshold) and reports, vs a
// full-Newton ground truth (eps=0): time, hit%, Newton-calls/ray, skip%, and
// pixel/normal error. 
// counters are extra. CSV to stdout.  Build the `benchmark` target, run > out.csv

#include "camera.h"
#include "bvh.h"
#include "../gpu/memory.h"
#include "../gpu/intersect.h"
#include "../gpu/generateRaysKernel.h"
#include "../gpu/soa.h"
#include "../math/math_types.h"
#include "../math/math_utils.h"
#include "../math/transform.h"
#include "../math/nurbs.h"
#include "../math/bezier_patch.h"
#include <cstdio>
#include <cstring>
#include <cmath>
#include <chrono>
#include <vector>
#include <cuda_runtime.h>

using MathTypes::vector3f;
using MathTypes::vector4f;

// ------------------------------------------------------------ small helpers ---
static float h_dot(vector3f a, vector3f b){ return a.v[0]*b.v[0]+a.v[1]*b.v[1]+a.v[2]*b.v[2]; }
static vector3f h_cross(vector3f a, vector3f b){
    return {a.v[1]*b.v[2]-a.v[2]*b.v[1], a.v[2]*b.v[0]-a.v[0]*b.v[2], a.v[0]*b.v[1]-a.v[1]*b.v[0]}; }
static vector3f h_norm(vector3f a){ float l=sqrtf(h_dot(a,a)); return l>0?vector3f{a.v[0]/l,a.v[1]/l,a.v[2]/l}:a; }

static float patchDelta(const bezierPatch& P){
    auto deh=[&](int i){ vector4f c=P.control_point[i]; float iw=1.0f/c.v[3]; return vector3f{c.v[0]*iw,c.v[1]*iw,c.v[2]*iw}; };
    int p=P.p,q=P.q,row=q+1;
    vector3f c00=deh(0),c10=deh(p*row),c01=deh(q),c11=deh(p*row+q);
    float d=0;
    for(int i=0;i<=p;i++)for(int j=0;j<=q;j++){
        float u=(float)i/p,v=(float)j/q;
        vector3f B={(1-u)*(1-v)*c00.v[0]+u*(1-v)*c10.v[0]+(1-u)*v*c01.v[0]+u*v*c11.v[0],
                    (1-u)*(1-v)*c00.v[1]+u*(1-v)*c10.v[1]+(1-u)*v*c01.v[1]+u*v*c11.v[1],
                    (1-u)*(1-v)*c00.v[2]+u*(1-v)*c10.v[2]+(1-u)*v*c01.v[2]+u*v*c11.v[2]};
        vector3f Pij=deh(i*row+j);
        d=fmaxf(d,sqrtf((Pij.v[0]-B.v[0])*(Pij.v[0]-B.v[0])+(Pij.v[1]-B.v[1])*(Pij.v[1]-B.v[1])+(Pij.v[2]-B.v[2])*(Pij.v[2]-B.v[2])));
    }
    return d;
}

//  fresh NURBS surfaces 
// (heap-allocated so curvatureRefinement can own/replace the arrays)
static const float S=0.70710678f;
static const float CU_CX[9]={1,1,0,-1,-1,-1,0,1,1}, CU_CY[9]={0,1,1,1,0,-1,-1,-1,0}, CU_CW[9]={1,S,1,S,1,S,1,S,1};

static ParametricSurf::nurbsSurf makeTorus(){
    const float R=1.0f, r=0.4f, cz=6.0f;
    float* U=new float[12]{0,0,0,0.25f,0.25f,0.5f,0.5f,0.75f,0.75f,1,1,1};
    float* V=new float[12]{0,0,0,0.25f,0.25f,0.5f,0.5f,0.75f,0.75f,1,1,1};
    vector4f* cp=new vector4f[81];
    for(int i=0;i<9;i++)for(int j=0;j<9;j++){
        float radial=R+r*CU_CX[j], zz=r*CU_CY[j], w=CU_CW[i]*CU_CW[j];
        cp[i*9+j]={w*radial*CU_CX[i], w*radial*CU_CY[i], w*(zz+cz), w};
    }
    return ParametricSurf::nurbsSurf{U,V,cp,8,8,2,2};
}

static ParametricSurf::nurbsSurf makeSphere(){
    const float Rs=1.0f, cz=6.0f;
    // semicircle (radial,z) from north pole -> equator -> south pole
    const float svr[5]={0,Rs,Rs,Rs,0}, svz[5]={Rs,Rs,0,-Rs,-Rs}, svw[5]={1,S,1,S,1};
    float* U=new float[12]{0,0,0,0.25f,0.25f,0.5f,0.5f,0.75f,0.75f,1,1,1};
    float* V=new float[8]{0,0,0,0.5f,0.5f,1,1,1};
    vector4f* cp=new vector4f[45];   // 9 x 5
    for(int i=0;i<9;i++)for(int j=0;j<5;j++){
        float w=CU_CW[i]*svw[j];
        cp[i*5+j]={w*svr[j]*CU_CX[i], w*svr[j]*CU_CY[i], w*(svz[j]+cz), w};
    }
    return ParametricSurf::nurbsSurf{U,V,cp,8,4,2,2};
}

// image writers 
static unsigned char toB(float c){ c=c<0?0:(c>1?1:c); return (unsigned char)(c*255.999f); }

// normal map (0.5*n + 0.5), background dark
static void writeNormalPPM(const char* fn,int W,int H,const surfaceInteraction* si,
                           const std::vector<int>& idx,const std::vector<float>& dx,
                           const std::vector<float>& dy,const std::vector<float>& dz){
    std::vector<unsigned char> rgb(3*(size_t)W*H);
    for(long i=0;i<(long)W*H;i++){ unsigned char* p=&rgb[3*i];
        if(idx[i]==-1){ p[0]=p[1]=p[2]=15; continue; }
        vector3f n=h_norm(h_cross(si[i].dpdu, si[i].dpdv)), d={dx[i],dy[i],dz[i]};
        if(h_dot(n,d)>0){ n.v[0]=-n.v[0]; n.v[1]=-n.v[1]; n.v[2]=-n.v[2]; }
        p[0]=toB((n.v[0]+1)*0.5f); p[1]=toB((n.v[1]+1)*0.5f); p[2]=toB((n.v[2]+1)*0.5f);
    }
    FILE* f=fopen(fn,"wb"); fprintf(f,"P6\n%d %d\n255\n",W,H); fwrite(rgb.data(),1,rgb.size(),f); fclose(f);
    fprintf(stderr, "wrote %s\n", fn);
}

// per-pixel position error vs ground truth: dark=exact, red=err, green=hit/miss mismatch
static void writeErrorPPM(const char* fn,int W,int H,const surfaceInteraction* ref,const surfaceInteraction* tst,
                          const std::vector<int>& ridx,const std::vector<int>& tidx,
                          const std::vector<float>& ox,const std::vector<float>& oy,const std::vector<float>& oz,
                          const std::vector<float>& dx,const std::vector<float>& dy,const std::vector<float>& dz,
                          float cam_k,float scale){
    std::vector<unsigned char> rgb(3*(size_t)W*H);
    for(long i=0;i<(long)W*H;i++){ unsigned char* p=&rgb[3*i];
        bool rh=ridx[i]!=-1, th=tidx[i]!=-1;
        if(!rh && !th){ p[0]=p[1]=p[2]=15; continue; }
        if(rh!=th){ p[0]=0; p[1]=200; p[2]=0; continue; }   // hit/miss disagreement
        float tr=ref[i].t_hit, tt=tst[i].t_hit;
        vector3f O={ox[i],oy[i],oz[i]}, d={dx[i],dy[i],dz[i]};
        vector3f Pr={O.v[0]+tr*d.v[0],O.v[1]+tr*d.v[1],O.v[2]+tr*d.v[2]};
        vector3f Pt={O.v[0]+tt*d.v[0],O.v[1]+tt*d.v[1],O.v[2]+tt*d.v[2]};
        float wd=sqrtf((Pr.v[0]-Pt.v[0])*(Pr.v[0]-Pt.v[0])+(Pr.v[1]-Pt.v[1])*(Pr.v[1]-Pt.v[1])+(Pr.v[2]-Pt.v[2])*(Pr.v[2]-Pt.v[2]));
        float e=(tr>0)? wd*cam_k/tr : 0.0f, v=e/scale; v=v<0?0:(v>1?1:v);
        p[0]=(unsigned char)(v*255.999f); p[1]=(unsigned char)((1-v)*30); p[2]=(unsigned char)((1-v)*70);
    }
    FILE* f=fopen(fn,"wb"); fprintf(f,"P6\n%d %d\n255\n",W,H); fwrite(rgb.data(),1,rgb.size(),f); fclose(f);
    fprintf(stderr, "wrote %s\n", fn);
}

int main(){
    const int res_x=1920, res_y=1080, max_pixels=res_x*res_y;
    const float d_fov=60.0f, v_fov=MathUtils::radians(d_fov);
    const float cam_k=0.5f*(float)res_y/tanf(v_fov/2.0f);   // for host-side pixel error
    const int REP=5;

    MathTypes::matrix4x4f world_from_camera = Transform::lookAt(
        MathTypes::point3f{3,3,0}, MathTypes::point3f{0,0,6}, MathTypes::vector3f{0,1,0});
    camera_config main_cam = CameraOps::init_camera(world_from_camera, res_x, res_y);
    main_cam = CameraOps::getPerspectiveCamera(main_cam, d_fov, 0.01f, 1000.0f);

    SoAQueue::ray_queue ray_queue = {
        CudaDevAlloc<float>(max_pixels), CudaDevAlloc<float>(max_pixels), CudaDevAlloc<float>(max_pixels),
        CudaDevAlloc<float>(max_pixels), CudaDevAlloc<float>(max_pixels), CudaDevAlloc<float>(max_pixels),
        CudaDevAlloc<int>(max_pixels)};
    calculateRaysWrapper(ray_queue, perspective, main_cam);

    // rays are constant -> pull origin/direction to host once (for pixel-error)
    std::vector<float> ox(max_pixels),oy(max_pixels),oz(max_pixels),dx(max_pixels),dy(max_pixels),dz(max_pixels);
    cudaMemcpy(ox.data(),ray_queue.origin_x,   max_pixels*sizeof(float),cudaMemcpyDeviceToHost);
    cudaMemcpy(oy.data(),ray_queue.origin_y,   max_pixels*sizeof(float),cudaMemcpyDeviceToHost);
    cudaMemcpy(oz.data(),ray_queue.origin_z,   max_pixels*sizeof(float),cudaMemcpyDeviceToHost);
    cudaMemcpy(dx.data(),ray_queue.direction_x,max_pixels*sizeof(float),cudaMemcpyDeviceToHost);
    cudaMemcpy(dy.data(),ray_queue.direction_y,max_pixels*sizeof(float),cudaMemcpyDeviceToHost);
    cudaMemcpy(dz.data(),ray_queue.direction_z,max_pixels*sizeof(float),cudaMemcpyDeviceToHost);

    surfaceInteraction* interact = CudaAllocManaged<surfaceInteraction>(max_pixels);
    unsigned long long* d_newton = CudaAllocManaged<unsigned long long>(1);
    unsigned long long* d_skip   = CudaAllocManaged<unsigned long long>(1);

    std::vector<int> test_idx(max_pixels), ref_idx(max_pixels);
    std::vector<surfaceInteraction> ref_si(max_pixels);

    float Cs[]   = {0,1,2,4,8,16};
    float epss[] = {0.25f,0.5f,1.0f,2.0f,4.0f,8.0f};

    printf("surface,C,patches,eps_px,time_ms,hit_pct,newton_per_ray,skip_pct,max_pix_err,mean_pix_err,max_norm_deg\n");

    for(int s=0;s<2;s++){
        const char* sname = (s==0)?"torus":"sphere";
        for(float C : Cs){
            ParametricSurf::nurbsSurf surf = (s==0)?makeTorus():makeSphere();
            if(C>0.0f){
                NurbsOps::curvatureRefinementU(surf.u_knots, surf.p, surf.n, surf.m, C, surf.control_point);
                NurbsOps::curvatureRefinementV(surf.v_knots, surf.q, surf.m, surf.n, C, surf.control_point);
            }
            int nbU,nbV;
            bezierPatch* patches = BezierPatchOps::getPatch(surf, nbU, nbV);
            int nP = nbU*nbV;

            BVHNode* nodes = CudaAllocManaged<BVHNode>(2*nP-1);
            BvhOps::buildBVH(nodes, 0, patches, nP);

            int psz=(patches[0].p+1)*(patches[0].q+1);
            vector4f* base=patches[0].heap_pointer;
            vector4f* gpu_cp=CudaAllocManaged<vector4f>((size_t)nP*psz);
            memcpy(gpu_cp, base, (size_t)nP*psz*sizeof(vector4f));
            bezierPatch* gpu_patches=CudaAllocManaged<bezierPatch>(nP);
            float* delta=CudaAllocManaged<float>(nP);
            for(int i=0;i<nP;i++){
                gpu_patches[i]=patches[i];
                gpu_patches[i].control_point=gpu_cp+(patches[i].control_point-base);
                gpu_patches[i].heap_pointer=gpu_cp;
                delta[i]=patchDelta(gpu_patches[i]);
            }
            BezierPatchOps::freePatches(patches);
            delete[] surf.u_knots; delete[] surf.v_knots; delete[] surf.control_point;

            //  ground truth: full Newton (flat_threshold = 0 -> never skips) 
            cudaMemset(ray_queue.idx_interact, 0xFF, max_pixels*sizeof(int));
            calculateIntersectWrapper(ray_queue, nodes, gpu_patches, interact, res_x, res_y, 0.0f, v_fov, delta);
            memcpy(ref_si.data(), interact, max_pixels*sizeof(surfaceInteraction));
            cudaMemcpy(ref_idx.data(), ray_queue.idx_interact, max_pixels*sizeof(int), cudaMemcpyDeviceToHost);
            if(C==16.0f){ char fn[160]; snprintf(fn,sizeof(fn),"img_%s_C16_exact_normal.ppm",sname);   // 0% skip reference
                writeNormalPPM(fn,res_x,res_y,ref_si.data(),ref_idx,dx,dy,dz); }

            for(float eps : epss){
                // run once with counters -> populates interact/idx + skip/newton counts
                cudaMemset(ray_queue.idx_interact, 0xFF, max_pixels*sizeof(int));
                *d_newton=0; *d_skip=0;
                calculateIntersectWrapper(ray_queue, nodes, gpu_patches, interact, res_x, res_y, eps, v_fov, delta, d_newton, d_skip);
                cudaMemcpy(test_idx.data(), ray_queue.idx_interact, max_pixels*sizeof(int), cudaMemcpyDeviceToHost);
                unsigned long long nc=*d_newton, sc=*d_skip;

                // errors vs ground truth (over pixels hit in both)
                double maxpx=0,sumpx=0,maxdeg=0; long npx=0,hits=0;
                for(long i=0;i<max_pixels;i++){
                    if(test_idx[i]!=-1) ++hits;
                    if(test_idx[i]==-1 || ref_idx[i]==-1) continue;
                    vector3f dir={dx[i],dy[i],dz[i]}, O={ox[i],oy[i],oz[i]};
                    float tr=ref_si[i].t_hit, tt=interact[i].t_hit;
                    vector3f Pr={O.v[0]+tr*dir.v[0],O.v[1]+tr*dir.v[1],O.v[2]+tr*dir.v[2]};
                    vector3f Pt={O.v[0]+tt*dir.v[0],O.v[1]+tt*dir.v[1],O.v[2]+tt*dir.v[2]};
                    float wdist=sqrtf((Pr.v[0]-Pt.v[0])*(Pr.v[0]-Pt.v[0])+(Pr.v[1]-Pt.v[1])*(Pr.v[1]-Pt.v[1])+(Pr.v[2]-Pt.v[2])*(Pr.v[2]-Pt.v[2]));
                    double pxerr = (tr>0)? wdist*cam_k/tr : 0.0;
                    maxpx=fmax(maxpx,pxerr); sumpx+=pxerr; ++npx;
                    vector3f nr=h_norm(h_cross(ref_si[i].dpdu, ref_si[i].dpdv));
                    vector3f nt=h_norm(h_cross(interact[i].dpdu, interact[i].dpdv));
                    double d=h_dot(nr,nt); d=d>1?1:(d<-1?-1:d);
                    maxdeg=fmax(maxdeg, acos(fabs(d))*180.0/3.14159265); // fabs: ignore normal sign
                }

                // dump figures for a few representative configs (uses the count-run's buffers)
                {
                    char fn[160];
                    if((C==0.0f && eps==0.5f) || (C==16.0f && (eps==0.5f || eps==8.0f))){
                        snprintf(fn,sizeof(fn),"img_%s_C%.0f_eps%.2f_normal.ppm",sname,C,eps);
                        writeNormalPPM(fn,res_x,res_y,interact,test_idx,dx,dy,dz);
                    }
                    if(C==16.0f && eps==8.0f){
                        snprintf(fn,sizeof(fn),"img_%s_C%.0f_eps%.2f_error.ppm",sname,C,eps);
                        writeErrorPPM(fn,res_x,res_y,ref_si.data(),interact,ref_idx,test_idx,ox,oy,oz,dx,dy,dz,cam_k,3.0f);
                    }
                }

                // timing (no counters, pure)
                cudaMemset(ray_queue.idx_interact, 0xFF, max_pixels*sizeof(int));
                calculateIntersectWrapper(ray_queue, nodes, gpu_patches, interact, res_x, res_y, eps, v_fov, delta); // warmup
                auto t0=std::chrono::high_resolution_clock::now();
                for(int r=0;r<REP;r++)
                    calculateIntersectWrapper(ray_queue, nodes, gpu_patches, interact, res_x, res_y, eps, v_fov, delta);
                auto t1=std::chrono::high_resolution_clock::now();
                double ms=std::chrono::duration<double,std::milli>(t1-t0).count()/REP;

                unsigned long long tot=nc+sc;
                printf("%s,%.0f,%d,%.2f,%.3f,%.2f,%.4f,%.2f,%.4f,%.4f,%.4f\n",
                    sname, C, nP, eps, ms, 100.0*hits/max_pixels,
                    (double)nc/max_pixels, tot?100.0*sc/tot:0.0,
                    maxpx, npx?sumpx/npx:0.0, maxdeg);
            }
            CudaFree(nodes); CudaFree(gpu_cp); CudaFree(gpu_patches); CudaFree(delta);
        }
    }

    // dedicated FIGURE renders (illustrative only; closer camera so the
    // object fills the frame; C=16 torus; 0% / mid / 99% skip + diff maps 
    {
        MathTypes::matrix4x4f fw = Transform::lookAt(
            MathTypes::point3f{1.5f,1.5f,3.3f}, MathTypes::point3f{0,0,6}, MathTypes::vector3f{0,1,0});
        camera_config fcam = CameraOps::init_camera(fw, res_x, res_y);
        fcam = CameraOps::getPerspectiveCamera(fcam, d_fov, 0.01f, 1000.0f);
        calculateRaysWrapper(ray_queue, perspective, fcam);
        cudaMemcpy(ox.data(),ray_queue.origin_x,   max_pixels*sizeof(float),cudaMemcpyDeviceToHost);
        cudaMemcpy(oy.data(),ray_queue.origin_y,   max_pixels*sizeof(float),cudaMemcpyDeviceToHost);
        cudaMemcpy(oz.data(),ray_queue.origin_z,   max_pixels*sizeof(float),cudaMemcpyDeviceToHost);
        cudaMemcpy(dx.data(),ray_queue.direction_x,max_pixels*sizeof(float),cudaMemcpyDeviceToHost);
        cudaMemcpy(dy.data(),ray_queue.direction_y,max_pixels*sizeof(float),cudaMemcpyDeviceToHost);
        cudaMemcpy(dz.data(),ray_queue.direction_z,max_pixels*sizeof(float),cudaMemcpyDeviceToHost);

        ParametricSurf::nurbsSurf surf = makeTorus();
        NurbsOps::curvatureRefinementU(surf.u_knots, surf.p, surf.n, surf.m, 16.0f, surf.control_point);
        NurbsOps::curvatureRefinementV(surf.v_knots, surf.q, surf.m, surf.n, 16.0f, surf.control_point);
        int nbU,nbV; bezierPatch* patches=BezierPatchOps::getPatch(surf,nbU,nbV); int nP=nbU*nbV;
        BVHNode* nodes=CudaAllocManaged<BVHNode>(2*nP-1); BvhOps::buildBVH(nodes,0,patches,nP);
        int psz=(patches[0].p+1)*(patches[0].q+1); vector4f* base=patches[0].heap_pointer;
        vector4f* gcp=CudaAllocManaged<vector4f>((size_t)nP*psz); memcpy(gcp,base,(size_t)nP*psz*sizeof(vector4f));
        bezierPatch* gp=CudaAllocManaged<bezierPatch>(nP); float* del=CudaAllocManaged<float>(nP);
        for(int i=0;i<nP;i++){ gp[i]=patches[i]; gp[i].control_point=gcp+(patches[i].control_point-base);
            gp[i].heap_pointer=gcp; del[i]=patchDelta(gp[i]); }
        BezierPatchOps::freePatches(patches); delete[] surf.u_knots; delete[] surf.v_knots; delete[] surf.control_point;

        // full-Newton reference (0% skip)
        cudaMemset(ray_queue.idx_interact,0xFF,max_pixels*sizeof(int));
        calculateIntersectWrapper(ray_queue,nodes,gp,interact,res_x,res_y,0.0f,v_fov,del);
        memcpy(ref_si.data(),interact,max_pixels*sizeof(surfaceInteraction));
        cudaMemcpy(ref_idx.data(),ray_queue.idx_interact,max_pixels*sizeof(int),cudaMemcpyDeviceToHost);
        writeNormalPPM("fig_torus_eps0_normal.ppm",res_x,res_y,interact,ref_idx,dx,dy,dz);

        float figeps[]={2.0f,4.0f,8.0f};
        for(float e: figeps){
            cudaMemset(ray_queue.idx_interact,0xFF,max_pixels*sizeof(int));
            calculateIntersectWrapper(ray_queue,nodes,gp,interact,res_x,res_y,e,v_fov,del);
            cudaMemcpy(test_idx.data(),ray_queue.idx_interact,max_pixels*sizeof(int),cudaMemcpyDeviceToHost);
            char fn[160];
            snprintf(fn,sizeof(fn),"fig_torus_eps%.0f_normal.ppm",e);
            writeNormalPPM(fn,res_x,res_y,interact,test_idx,dx,dy,dz);
            snprintf(fn,sizeof(fn),"fig_torus_eps%.0f_diff.ppm",e);
            writeErrorPPM(fn,res_x,res_y,ref_si.data(),interact,ref_idx,test_idx,ox,oy,oz,dx,dy,dz,cam_k,2.5f);
        }
        CudaFree(nodes); CudaFree(gcp); CudaFree(gp); CudaFree(del);
    }

    CudaFree(interact); CudaFree(d_newton); CudaFree(d_skip);
    CudaFree(ray_queue.origin_x); CudaFree(ray_queue.origin_y); CudaFree(ray_queue.origin_z);
    CudaFree(ray_queue.direction_x); CudaFree(ray_queue.direction_y); CudaFree(ray_queue.direction_z);
    CudaFree(ray_queue.idx_interact);
    return 0;
}
