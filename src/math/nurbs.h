#ifndef NURBS_H
#define NURBS_H

#include "../core/core_defines.h"
#include "../constants.h"
#include "math_types.h"
#include "vector_ops.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <utility>
namespace ParametricSurf{



	struct nurbsSurf{
		float *u_knots, *v_knots; // these hold n+p+2/m+q+2 elements
		MathTypes::vector4f* control_point; // Row major 2d->1d array (n+1)*(m+1)
		int n,m; // Indices of the last control point
		int p,q;
	};


}



// Most algorithms came from The Nurbs Book 
namespace NurbsOps{


	// Algorithm A2.1  (The Nurbs Book)
	CUDA_FINL CORE_HD
	int findSpan(int n, int p, float u, const float* U){
		if(u == U[n+1]) return(n);
		int low = p; 
		int high = n+1;

		int mid = (low+high)/2;

		while(u<U[mid] || u >= U[mid+1]){
			if(u<U[mid]) high = mid;
			else low = mid;
			mid = (low+high)/2;
	
		}
		return mid;
	}

	// Algorithm A2.2 ( The Nurbs Book)
	CUDA_FINL CORE_HD
	void basisFuncEval(const int i, const float u, const int p, const float* U, float* N){
		N[0]=1.0;
		for(int j=1; j<=p; j++){
			float saved = 0.0f;
			for(int r=0; r<j; r++){
				float temp = N[r]/(U[i+r+1] - U[i+1-j+r]);
				N[r] = saved+(U[i+r+1]-u)*temp;
				saved = (u - U[i+1-j+r])*temp;
			}
			N[j] = saved;
		}
	}
	
	// Algorithm A4.3 (The Nurbs Book)
	CUDA_FINL CORE_HD
	MathTypes::vector4f surfacePoint(	const int n, const int p, const float* U,
						const int m, const int q, const float* V, 
						MathTypes::vector4f* P_w, 
						const float u, const float v)
	{
		float N_u[NURBS_MAX_DEGREE+1]; // Assuming p<4;
		float N_v[NURBS_MAX_DEGREE+1]; // Assuming q<4;
		int uspan = findSpan(n, p, u, U);
		basisFuncEval(uspan, u, p, U, N_u);	
		int vspan = findSpan(m, q, v, V);
		basisFuncEval(vspan, v, q, V, N_v);	

		MathTypes::vector4f temp[NURBS_MAX_DEGREE+1];
		for(int l=0; l<=q; l++){
			temp[l] = {0.0f};
			for(int k=0; k<=p; k++)
				temp[l] = VectorOps::sum(temp[l], VectorOps::multiply(N_u[k],P_w[(uspan-p+k)*(m+1) + (vspan-q+l)]));
		}
		
		MathTypes::vector4f S_w = {0.0f};
		for(int l=0; l<=q; l++){
			S_w = VectorOps::sum(S_w, VectorOps::multiply(N_v[l], temp[l]));
		}
		return VectorOps::multiply(S_w, 1/S_w.v[3]); // S_w/w
	}

	
	// Algorithm A5.1 (The Nurbs Book)
	inline CORE_H
	void curveKnotInsert(	int np, int p, const float* UP, 
						 	float u, int k, int s, 
							int r, int& nq, float* UQ, 
							const MathTypes::vector4f* P_w,
							MathTypes::vector4f* Q_w )
	{
		if(r+s>p){
			fprintf(stderr, "Multiplicity needs to be lower then degree\n");
			return;
		}
		
		if(p>NURBS_MAX_DEGREE){
			fprintf(stderr, "Degree higher than maximum allowed\n");
			return;
		}

		int mp = np + p + 1;	
		nq = np + r;
		int L=0;

		if(r<=0)
			return;

		MathTypes::vector4f R_w[NURBS_MAX_DEGREE+1]; 

		for(int i=0; i <= k; i++)	UQ[i] = UP[i];
		for(int i=1; i <= r; i++)	UQ[k+i] = u;
		for(int i=k+1; i<=mp; i++)	UQ[i+r] = UP[i];
		for(int i=0; i<=k-p; i++)	Q_w[i]	= P_w[i];
		for(int i=k-s; i<=np; i++) 	Q_w[i+r] = P_w[i];
		for(int i=0; i<=p-s; i++)	R_w[i] = P_w[k-p+i];

		for(int j=1; j<=r; j++){
			L = k-p+j;
			for(int i = 0; i<=p-j-s; i++){
				float alpha = (u-UP[L+i])/(UP[i+k+1]-UP[L+i]);
				R_w[i] = VectorOps::sum(VectorOps::multiply(alpha,R_w[i+1]), VectorOps::multiply((1.0f-alpha),R_w[i]));
			}
			Q_w[L] = R_w[0];
			Q_w[k+r-j-s] = R_w[p-j-s];
		}
		for(int i=L+1; i<k-s; i++)
			Q_w[i] = R_w[i-L];
	}
	
	//Algorithm A5.7
	inline CORE_H
	void decomposeSurface(	int n, int p, const float* U,
							int m, int q, const float* V,
							int dir, int& nb, 
							const MathTypes::vector4f* P_w,
							MathTypes::vector4f* Q_w)
	{
		
		
		if(dir == U_DIRECTION){
			const int row_stride = m+1;
			const int net_stride = (p+1)*row_stride;
			#define QW_U(b,i,r) Q_w[(b)*net_stride + (i)*row_stride+(r)]
			#define PW_U(i,r) P_w[(i)*row_stride + (r)]
			
			int a=p;
			int b=p+1;
			float alphas[NURBS_MAX_DEGREE];
			nb=0;
			for(int i=0; i<=p; i++)
				for(int row=0; row <= m; row++)
					QW_U(nb,i,row) = PW_U(i,row);
			int mu = n + p + 1;
			while(b<mu){
				int i = b;
				while(b<mu && U[b+1] == U[b]) b++;
				int mult = b-i+1;
				if(mult<p){
					float numer = U[b] - U[a];
					for(int j=p; j>mult; j--)
						alphas[j-mult-1] = numer/(U[a+j]-U[a]);
					int r = p-mult;
					for(int j=1; j<=r; j++){
						int save = r-j;
						int s = mult+j;
						for(int k=p; k>=s; k--){
							float alpha = alphas[k-s];
							for(int row=0; row<=m; row++)
								QW_U(nb,k,row) = VectorOps::sum(	VectorOps::multiply(alpha,QW_U(nb,k,row)), 
																	VectorOps::multiply((1.0f-alpha),QW_U(nb,k-1,row)));
						}
						if(b<mu)
							for(int row=0; row<=m; row++)
								QW_U(nb+1,save,row) = QW_U(nb,p,row);
					}
				}	
				nb=nb+1;
				if(b<mu){
					for(int 	i=p-mult; i<=p; i++)
						for(int row=0; row<=m; row++)
							QW_U(nb,i,row) = PW_U(b-p+i,row);
					a=b;
					b=b+1;
				}
			}
		}
		else if(dir == V_DIRECTION){
			const int in_row_stride = m+1;
			const int out_row_stride = q+1;
			const int net_stride = (n+1)*out_row_stride;
			#define QW_V(b,r,j) Q_w[(b)*net_stride + (r)*out_row_stride+(j)]
			#define PW_V(r,j) P_w[(r)*in_row_stride + (j)]
			int a=q;
			int b=q+1;
			float alphas[NURBS_MAX_DEGREE];
			nb=0;
			for(int i=0; i<=q; i++)
				for(int row=0; row <= n; row++)
					QW_V(nb,row,i) = PW_V(row,i);
			int mv = m + q + 1;
			while(b<mv){
				int i = b;
				while(b<mv && V[b+1] == V[b]) b++;
				int mult = b-i+1;
				if(mult<q){
					float numer = V[b] - V[a];
					for(int j=q; j>mult; j--)
						alphas[j-mult-1] = numer/(V[a+j]-V[a]);
					int r = q-mult;
					for(int j=1; j<=r; j++){
						int save = r-j;
						int s = mult+j;
						for(int k=q; k>=s; k--){
							float alpha = alphas[k-s];
							for(int row=0; row<=n; row++)
								QW_V(nb,row,k) = VectorOps::sum(	VectorOps::multiply(alpha,QW_V(nb,row,k)), 
																	VectorOps::multiply((1.0f-alpha),QW_V(nb,row,k-1)));
						}
						if(b<mv)
							for(int row=0; row<=n; row++)
								QW_V(nb+1,row,save) = QW_V(nb,row,q);
					}
				}	
				nb=nb+1;
				if(b<mv){
					for(int 	i=q-mult; i<=q; i++)
						for(int row=0; row<=n; row++)
							QW_V(nb,row,i) = PW_V(row,b-q+i);
					a=b;
					b=b+1;
				}
			}
		}
		#undef QW_V
		#undef QW_U
		#undef PW_V
		#undef PW_U
	}
	

	// Algorithms based on the formulas presented at the Practical Ray Tracing of trimmed NURBS Surfaces
	inline 
	int knotsToAdd(MathTypes::vector3f* P, float* U, int i, float C, int p){
		using namespace VectorOps;
		float V_j_sum = 0;
		MathTypes::vector3f V_j[NURBS_MAX_DEGREE+1];
		for(int j=i-p+1; j<=i; j++){
			const int base = i-p+1;
			float d = U[j+p]-U[j];
			V_j[j-base] = (d>0.0f) 	? multiply(p/d, subtract(P[j], P[j-1])) 
									: MathTypes::vector3f{0,0,0};
			V_j_sum += length(V_j[j-base]);
		}
		float A_j_mod = 0;
		for(int j=i-p+2; j<=i; j++){
			const int base = i-p+1;
			float d = U[j+p-1]-U[j];
			MathTypes::vector3f A_j = (d>0.0f) 	? multiply((p-1)/d, subtract(V_j[j-base], V_j[j-1-base]))
												: MathTypes::vector3f{0,0,0};
			A_j_mod = fmaxf(A_j_mod, length(A_j));
		}
		
		float num = A_j_mod*powf(U[i+1]-U[i],1.5f);
		float den = sqrtf(V_j_sum/p);
		if(!(den > 1e-9f)) return 0;
		float nf = C*num/den;
		if(!(nf > 1e-9f)) return 0;
		return nf;
	}
	
	inline int multiplicity(const float* U, int k, float u) {
    	int s = 0;
    	for (int i = k; i >= 0 && U[i] == u; --i) ++s;
    	return s;
	}

	inline CORE_H
	void curvatureRefinementU(float*& U, int p, int& n, int& m, float C, MathTypes::vector4f*& P){
		int* add_knots = new int[n-p+1];
		MathTypes::vector3f* curve = new MathTypes::vector3f[n+1];
		const int row_stride = m+1;
		#define PU(i,r) P[(i)*row_stride + (r)]
			

		for(int i=p; i<=n; i++){				
			int n_i = 0;
			for(int c=0; c<=m; c++){

				for(int a=0; a<=n; a++)
					curve[a] = MathTypes::vector3f{	PU(a,c).v[0]/PU(a,c).v[3], 
													PU(a,c).v[1]/PU(a,c).v[3],
													PU(a,c).v[2]/PU(a,c).v[3]}; 
				
				n_i = std::max(n_i, knotsToAdd(curve, U, i, C, p));	
			}
			add_knots[i-p] = n_i;
		}

		delete[] curve;
		#undef PU

		int total = 0;
		for(int i =0; i<=n-p; ++i) total += add_knots[i];

		int n_final = n+total;

		MathTypes::vector4f* buf_a = new MathTypes::vector4f[(n_final+1)*(m+1)];
		MathTypes::vector4f* buf_b = new MathTypes::vector4f[(n_final+1)*(m+1)];

		float* UB = new float[n_final+p+2];
		float* UA = new float[n_final+p+2];

		std::memcpy(UA, U, (n+p+2)*sizeof(float));
		std::memcpy(buf_a, P, (n+1)*(m+1)*sizeof(MathTypes::vector4f));
		int n_cur = n;
		MathTypes::vector4f *src = buf_a, *dst = buf_b;
		float *U_src = UA, *U_dst = UB;
		int nq = 0;


		#define SRC(i,r) src[(i)*row_stride + (r)]
		#define DST(i,r) dst[(i)*row_stride + (r)]

		MathTypes::vector4f* curve_in = new MathTypes::vector4f[n_final+1];
		MathTypes::vector4f* curve_out = new MathTypes::vector4f[n_final+1];

		for(int i=p; i<=n; i++){
			int cnt = add_knots[i-p];
			if(cnt<=0) continue;

			float lo = U[i], hi = U[i+1];
			for(int j=1; j<=cnt; j++){
				float u = lo+(hi-lo)*j/(float)(cnt+1);

				int k = findSpan(n_cur, p, u, U_src);
				int s = multiplicity(U_src, k, u);

				for(int c = 0; c<=m; c++){
					for(int a=0; a<=n_cur; a++)
						curve_in[a] = SRC(a,c); 
				
					curveKnotInsert(n_cur, p, U_src, u, k, s, 1, nq, U_dst, curve_in, curve_out);
					
					for(int a=0; a<=n_cur+1; a++)
						DST(a,c) = curve_out[a];
				}

				std::swap(src, dst);
				std::swap(U_src, U_dst);
				n_cur++;

			}
		} 
		
		std::swap(U, U_src);
		std::swap(P, src);

		n = n_final;
		
		delete[] U_src;    
		delete[] src;      
		delete[] U_dst;    
		delete[] dst;      
		delete[] add_knots;
		delete[] curve_in;
		delete[] curve_out;
		#undef SRC 
		#undef DST 
	}

	inline CORE_H
	void curvatureRefinementV(float*& V, int q, int& m, int n, float C, MathTypes::vector4f*& P){
		int* add_knots = new int[m-q+1];
		MathTypes::vector3f* curve = new MathTypes::vector3f[n+1];
		const int row_stride = m+1;
		#define PU(i,r) P[(i)*row_stride + (r)]
			

		for(int i=q; i<=m; i++){				
			int n_i = 0;
			for(int c=0; c<=n; c++){

				for(int a=0; a<=m; a++)
					curve[a] = MathTypes::vector3f{	PU(c,a).v[0]/PU(c,a).v[3], 
													PU(c,a).v[1]/PU(c,a).v[3],
													PU(c,a).v[2]/PU(c,a).v[3]}; 
				
				n_i = std::max(n_i, knotsToAdd(curve, V, i, C, q));	
			}
			add_knots[i-q] = n_i;
		}

		delete[] curve;
		#undef PU

		int total = 0;
		for(int i =0; i<=m-q; ++i) total += add_knots[i];

		int m_final = m+total;

		MathTypes::vector4f* buf_a = new MathTypes::vector4f[(n+1)*(m_final+1)];
		MathTypes::vector4f* buf_b = new MathTypes::vector4f[(n+1)*(m_final+1)];

		float* VB = new float[m_final+q+2];
		float* VA = new float[m_final+q+2];

		std::memcpy(VA, V, (m+q+2)*sizeof(float));
		std::memcpy(buf_a, P, (n+1)*(m+1)*sizeof(MathTypes::vector4f));

		// We need to insert into the control point directly here because the row size changed
		// when only the column size changed we could get away with it whitout restructing the array
		
		for(int c=0; c<=n; c++){
			for(int a=0; a<=m; a++){
				buf_a[c*(m_final+1) + a] = P[c*(m+1) + a];
			}
		}


		int m_cur = m;
		MathTypes::vector4f *src = buf_a, *dst = buf_b;
		float *V_src = VA, *V_dst = VB;
		int nq = 0;


		#define SRC(i,r) src[(i)*(m_cur+1) + (r)]
		#define DST(i,r) dst[(i)*(m_cur+1) + (r)]

		MathTypes::vector4f* curve_in = new MathTypes::vector4f[m_final+1];
		MathTypes::vector4f* curve_out = new MathTypes::vector4f[m_final+1];

		for(int i=q; i<=m; i++){
			int cnt = add_knots[i-q];
			if(cnt<=0) continue;

			float lo = V[i], hi = V[i+1];
			for(int j=1; j<=cnt; j++){
				float v = lo+(hi-lo)*j/(float)(cnt+1);

				int k = findSpan(m_cur, q, v, V_src);
				int s = multiplicity(V_src, k, v);

				for(int c = 0; c<=n; c++){
					for(int a=0; a<=m_cur; a++)
						curve_in[a] = src[(c)*(m_final+1) + (a)]; 
				
					curveKnotInsert(m_cur, q, V_src, v, k, s, 1, nq, V_dst, curve_in, curve_out);
					
					for(int a=0; a<=m_cur+1; a++)
						dst[(c)*(m_final+1) + (a)] = curve_out[a];
				}

				std::swap(src, dst);
				std::swap(V_src, V_dst);
				m_cur++;

			}
		} 
		
		std::swap(V, V_src);
		std::swap(P, src);

		m = m_final;
		
		delete[] V_src;    
		delete[] src;      
		delete[] V_dst;    
		delete[] dst;      
		delete[] add_knots;
		delete[] curve_in;
		delete[] curve_out;
		#undef SRC 
		#undef DST 
	}							

}


#endif
