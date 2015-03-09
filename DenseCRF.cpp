#include "stdafx.h"
#include "densecrf.h"

#ifdef __SSE__
# define SSE_DENSE_CRF
#endif

float* allocate(size_t N);
void deallocate(float*& ptr);

PairwisePotential::~PairwisePotential() {
}
SemiMetricFunction::~SemiMetricFunction() {
}
class PottsPotential: public PairwisePotential
{
protected:
	Permutohedral lattice_;
	PottsPotential( const PottsPotential&o ){}
	int N_;
	float w_;
	float *norm_;
public:
	~PottsPotential()
	{
		deallocate( norm_ );
	}
	PottsPotential(const float* features, int D, int N, float w, bool per_pixel_normalization=true) :N_(N), w_(w) 
	{
		lattice_.init( features, D, N );
		norm_ = allocate( N );
		for ( int i=0; i<N; i++ )
			norm_[i] = 1;
		// Compute the normalization factor
		lattice_.compute( norm_, norm_, 1 );
		if ( per_pixel_normalization ) {
			// use a per pixel normalization
			#pragma omp parallel for
			for ( int i=0; i<N; i++ )
				norm_[i] = 1.0 / (norm_[i]+1e-20);
		}
		else {
			float mean_norm = 0;
			for ( int i=0; i<N; i++ )
				mean_norm += norm_[i];
			mean_norm = N / mean_norm;
			// use a per pixel normalization
			for ( int i=0; i<N; i++ )
				norm_[i] = mean_norm;
		}
	}

	void apply(float* out_values, const float* in_values, float* tmp, int value_size) const {
		lattice_.compute( tmp, in_values, value_size );
#pragma omp parallel for
		for ( int i=0; i<N_; i++ ){
			for ( int j=0; j<value_size; j++)	{
				const int k = i*value_size + j;
				out_values[k] += w_*norm_[i]*tmp[k];
			}
		}
	}
};
class SemiMetricPotential: public PottsPotential{
protected:
	const SemiMetricFunction * function_;
public:
	void apply(float* out_values, const float* in_values, float* tmp, int value_size) const {
		lattice_.compute( tmp, in_values, value_size );

		// To the metric transform
		float * tmp2 = new float[value_size];
		for ( int i=0; i<N_; i++ ) {
			float * out = out_values + i*value_size;
			float * t1  = tmp  + i*value_size;
			function_->apply( tmp2, t1, value_size );
			for ( int j=0; j<value_size; j++ )
				out[j] -= w_*norm_[i]*tmp2[j];
		}
		delete[] tmp2;
	}
	SemiMetricPotential(const float* features, int D, int N, float w, const SemiMetricFunction* function, bool per_pixel_normalization=true) :PottsPotential( features, D, N, w, per_pixel_normalization ),function_(function) {
	}
};



/////////////////////////////
/////  Alloc / Dealloc  /////
/////////////////////////////
DenseCRF::DenseCRF(int N, int M) : N_(N), M_(M) {

	// initialize higher order terms
	addHO = 0;
	addDet = 0; 
	addCooc = 0; 
	addAtt = 0;
	addAttDet=0;

	unary_ = allocate( N_*M_ );
	additional_unary_ = allocate( N_*M_ );
	current_ = allocate( N_*M_ ); 
	next_ = allocate( N_*M_ ); 
	tmp_ = allocate( 2*N_*M_ ); 
	memset( additional_unary_, 0, sizeof(float)*N_*M_ );
}
//DenseCRF::~DenseCRF() 
//{
//	deallocate( unary_ ); 
//	deallocate( additional_unary_ ); 
//	deallocate( current_ ); 
//	deallocate( next_ ); 
//	deallocate( tmp_ );
//	for( unsigned int i=0; i<pairwise_.size(); i++ )
//		delete pairwise_[i];
//}

DenseCRF::DenseCRF(int N, int M, int A) : N_(N), M_(M), A_(A) {

	// initialize higher order terms
	addHO = 0;
	addDet = 0; 
	addCooc = 0; 
	addAtt = 0;
	addAttDet=0;

	unary_ = allocate( N_*M_ );
	additional_unary_ = allocate( N_*M_ );
	current_ = allocate( N_*M_ ); 
	next_ = allocate( N_*M_ ); 
	tmp_ = allocate( 2*N_*M_ ); 
	memset( additional_unary_, 0, sizeof(float)*N_*M_ );

	unary_att_ = allocate( N_*A_*2 );
	additional_unary_att_ = allocate( N_*A_*2 );
	current_att_ = allocate( N_*A_*2 ); 
	next_att_ = allocate( N_*A_*2 ); 
	tmp_att_ = allocate( 2*N_*A_*2 );
	memset( additional_unary_att_, 0, sizeof(float)*N_*A_*2 );
	pair_att_ = allocate( A_*A_);
	pair_join_att_ = allocate(A_*M_);
}

DenseCRF::DenseCRF(int N, int M, int A, int NS, int AD) : N_(N), M_(M), A_(A),NS_(NS),AD_(AD) {

	// initialize higher order terms
	addHO = 0;
	addDet = 0; 
	addCooc = 0;
	addAttDet = 0;
	addAtt = 0;

	unary_ = allocate( N_*M_ );
	additional_unary_ = allocate( N_*M_ );
	current_ = allocate( N_*M_ ); 
	next_ = allocate( N_*M_ ); 
	tmp_ = allocate( 2*N_*M_ );
	memset( additional_unary_, 0, sizeof(float)*N_*M_ );

	unary_att_ = allocate( N_*A_*2 );
	additional_unary_att_ = allocate( N_*A_*2 );
	current_att_ = allocate( N_*A_*2 ); 
	next_att_ = allocate( N_*A_*2 ); 
	tmp_att_ = allocate( 2*N_*A_*2 );
	memset( additional_unary_att_, 0, sizeof(float)*N_*A_*2 );
	pair_att_ = allocate( A_*A_);
	pair_join_att_ = allocate(A_*M_);


	unary_objdet_ = allocate( NS_*2);
	current_det_ = allocate( NS*2 );
	next_det_ = allocate( NS*2 );
	tmp_det_ = allocate( NS*2 );

	unary_attdet_ = allocate( NS_*2*AD_);
	current_attdet_ = allocate( NS_*AD_*2 ); 
	next_attdet_ = allocate( NS_*AD_*2 ); 
	tmp_attdet_ = allocate( 2*NS_*AD_*2 );
	pair_attdet_ = allocate( AD_*AD_);
	pair_join_attdet_ = allocate(AD_*M_);

}

DenseCRF::~DenseCRF() 
{
	deallocate( unary_ ); 
	deallocate( additional_unary_ ); 
	deallocate( current_ ); 
	deallocate( next_ ); 
	deallocate( tmp_ );
	for( unsigned int i=0; i<pairwise_.size(); i++ )
		delete pairwise_[i];

	if (addAtt==1)
	{
		deallocate( unary_att_ ); 
		deallocate( additional_unary_att_ ); 
		deallocate( current_att_ ); 
		deallocate( next_att_ ); 
		deallocate( tmp_att_ );
		deallocate( pair_att_);
		deallocate( pair_join_att_);

	}
	if (addAttDet==1)
	{
		deallocate( unary_attdet_ ); 
		deallocate( pair_attdet_ ); 
		deallocate( pair_join_attdet_ ); 
		deallocate( unary_objdet_ ); 
		deallocate( current_det_ );
		deallocate( next_det_);
		deallocate( tmp_det_);
	}
}



DenseCRF2D::DenseCRF2D(int W, int H, int M) : DenseCRF(W*H,M), W_(W), H_(H)
{
}



DenseCRF2D::DenseCRF2D(int W, int H, int M, int A) : DenseCRF(W*H,M,A), W_(W), H_(H)
{
}

DenseCRF2D::DenseCRF2D(int W, int H, int M, int A, int NS, int AD) : DenseCRF(W*H,M,A,NS,AD), W_(W), H_(H)
{
}

DenseCRF2D::~DenseCRF2D()
{
}
/////////////////////////////////
/////  Pairwise Potentials  /////
/////////////////////////////////
void DenseCRF::addPairwiseEnergy (const float* features, int D, float w, const SemiMetricFunction * function) 
{
	if (function)
		addPairwiseEnergy( new SemiMetricPotential( features, D, N_, w, function ) );
	else
		addPairwiseEnergy( new PottsPotential( features, D, N_, w ) );
}

void DenseCRF::addPairwiseEnergy ( PairwisePotential* potential )
{
	pairwise_.push_back( potential );
}

void DenseCRF2D::addPairwiseGaussian ( float sx, float sy, float w, const SemiMetricFunction * function ) 
{
	float * feature = new float [N_*2];
#pragma omp parallel for
	for( int j=0; j<H_; j++ ){
		for( int i=0; i<W_; i++ ){
			feature[(j*W_+i)*2+0] = i / sx;
			feature[(j*W_+i)*2+1] = j / sy;
		}
	}
	addPairwiseEnergy( feature, 2, w, function );
	delete [] feature;
}

void DenseCRF2D::addPairwiseBilateral ( float sx, float sy, float sr, float sg, float sb, const unsigned char* im, float w, const SemiMetricFunction * function ) {
	float * feature = new float [N_*5];
#pragma omp parallel for
	for( int j=0; j<H_; j++ ){
		float* fj = feature + j*W_*5;
		const byte* imj = im + j*W_*3;
		for( int i=0; i<W_; i++, fj+= 5, imj+=3 ){
			fj[0] = i/sx;
			fj[1] = j/sy;
			fj[2] = imj[0]/sr;
			fj[3] = imj[1]/sg;
			fj[4] = imj[2]/sb;
			//feature[(j*W_+i)*5+0] = i / sx;
			//feature[(j*W_+i)*5+1] = j / sy;
			//feature[(j*W_+i)*5+2] = im[(i+j*W_)*3+0] / sr;
			//feature[(j*W_+i)*5+3] = im[(i+j*W_)*3+1] / sg;
			//feature[(j*W_+i)*5+4] = im[(i+j*W_)*3+2] / sb;
		}
	}
	addPairwiseEnergy( feature, 5, w, function );
	delete [] feature;
}
//////////////////////////////
/////  Unary Potentials  /////
//////////////////////////////
void DenseCRF::setUnaryEnergy ( const float* unary/*,  float *cooc_unary, float *cooc_pairwise*/) 
{
	memcpy( unary_, unary, N_*M_*sizeof(float) );
}
void DenseCRF::setUnaryEnergy ( const float* unary, const float* att_unary,float* att_pairwise,float* att_joint_pairwise) 
{
	memcpy( unary_, unary, N_*M_*sizeof(float) );
	//start attributes
	memcpy( unary_att_, att_unary, N_*A_*2*sizeof(float) );
	memcpy( pair_att_, att_pairwise,A_*A_*sizeof(float) );
	memcpy( pair_join_att_,att_joint_pairwise,A_*M_*sizeof(float) );
	//end attributes
}
void DenseCRF::setUnaryEnergy ( const float* unary, const float* unary_objdet, const float* att_unary, const float* attdet_unary,float* att_pairwise,float* att_joint_pairwise,  float *attdet_pairwise, float *attdet_joint_pairwise) 
{
	memcpy( unary_, unary, N_*M_*sizeof(float) );
	memcpy( unary_att_, att_unary, N_*A_*2*sizeof(float) );
	memcpy( pair_att_, att_pairwise,A_*A_*sizeof(float) );
	memcpy( pair_join_att_,att_joint_pairwise,A_*M_*sizeof(float) );

	//start attdet
	memcpy( unary_attdet_, att_unary, NS_*AD_*2*sizeof(float) );
	memcpy( pair_attdet_, att_pairwise,AD_*AD_*sizeof(float) );
	memcpy( pair_join_attdet_,att_joint_pairwise,AD_*M_*sizeof(float) );
	memcpy( unary_objdet_,unary_objdet,NS_*2*sizeof(float));
	//end attdet
}


void DenseCRF::setUnaryEnergy ( int n, const float* unary ) 
{
	memcpy( unary_+n*M_, unary, M_*sizeof(float) );
}
void DenseCRF2D::setUnaryEnergy ( int x, int y, const float* unary ) {
	memcpy( unary_+(x+y*W_)*M_, unary, M_*sizeof(float) );
}
///////////////////////
/////  Inference  /////
///////////////////////

Mat DenseCRF2D::map(int nIterations, vecM &resAtts1u, float relax)
{
	// Find the map
	const float*prob = runInference(nIterations, relax);
	Mat obj1u(H_, W_, CV_8U);
	byte* objLabel = (byte*)(obj1u.data);
#pragma omp parallel for
	for( int i=0; i<N_; i++ ){
		const float * p = prob + i*M_;
		float mx = p[0];
		int imx = 0;
		for( int j=1; j<M_; j++)
			if( mx < p[j] )	
				mx = p[j], imx = j;
		objLabel[i] = imx;
	}

	// Find attributes results if needed
	if (resAtts1u.size() != A_)
		return obj1u;
	for (int a = 0; a < A_; a++)
		resAtts1u[a] = Mat(obj1u.size(), CV_8U);
#pragma omp parallel for
	for( int i=0; i<N_; i++ )for( int a=0; a<A_; a++ ){
		const float * p = current_att_ + i*A_*2+a*2;
		resAtts1u[a].data[i] = p[0] > p[1] ? 0 : 255;
	}
	return obj1u;
}


void DenseCRF::map ( int n_iterations, short* result, short* result_att, short* result_attdet, float relax) 
{
	//start attdet
	runInference( n_iterations, relax);
	const float* prob = current_; 
	const float* prob_att = current_att_; 
	const float* prob_attdet_ = current_attdet_;
	const float* prob_det = current_det_;

	//find the map
	for( int i=0; i<N_; i++ ){
		const float * p = prob + i*M_;
		float mx = p[0];
		int imx = 0;
		for( int j=1; j<M_; j++)
			if( mx < p[j] )
				mx = p[j], imx = j;
		result[i] = imx;
	}

	if (result_att == NULL)
		return;
	for( int i=0; i<N_; i++ )for( int a=0; a<A_; a++ ){
		const float * p = prob_att + i*A_*2+a*2;
		result_att[i*A_ + a] = p[0] > p[1] ? 0 : 255;
	}

	if(result_attdet == NULL)
		return;
	for( int i=0; i<NS_; i++ )for( int a=0; a<AD_; a++ ){
		const float * p = prob_attdet_ + i*AD_*2+a*2;
		int imx = p[0] > p[1] ? 0 : 255;
		result_attdet[i*AD_+a] = prob_det[i*2+1] < prob_det[i*2] ? 0 : imx;
	}
}

float* DenseCRF::runInference( int n_iterations, float relax ) 
{
	startInference();
	for( int it=0; it<n_iterations; it++)
		stepInference(relax);
	return current_;
}

void DenseCRF::expAndNormalize_ojbdet(float *objdet_out, float *objdet_in, float scale, float relax )
{
	//float *V = new float [NS_+10 ];
#pragma omp parallel for
	for (int i = 0; i < NS_; i++ ){
		//(num_of_attributes*2)*(i*imagewidth+j)+(l+a*2)
		const float * b = objdet_in + (2*i);
		// Find the max and subtract it so that the exp doesn't explode
		float mx = scale*max(b[0], b[1]);
		float tt = 0;
		float V[2];
		for( int j=0; j<2; j++ ){
			V[j] = fast_exp( scale*b[j]-mx );
			tt += V[j];
		}
		// Make it a probability
		for( int j=0; j<2; j++ )
			V[j] /= tt;

		float * aa = objdet_out + (2*i);
		for( int j=0; j<2; j++ )
			aa[j] = relax == 1 ? V[j] : (1-relax)*aa[j] + relax*V[j];
	}
	//delete[] V;
}

void DenseCRF::expAndNormalize_attdet( float *attdet_out, float *attdet_in, float scale, float relax )
{
#pragma omp parallel for
	for (int i = 0; i < NS_; i++ ){
		for( int a = 0; a < AD_; a++){
			//(num_of_attributes*2)*(i*imagewidth+j)+(l+a*2)

			const float * b = attdet_in + (AD_*2)*(i)+(a*2);
			// Find the max and subtract it so that the exp doesn't explode
			float mx = scale* max(b[0], b[1]);
			float tt = 0;
			float V[2];
			for( int j=0; j<2; j++ ){
				V[j] = fast_exp(scale*b[j]-mx );
				tt += V[j];
			}
			// Make it a probability
			for( int j=0; j<2; j++ )
				V[j] /= tt;

			float * aa = attdet_out + (AD_*2)*(i)+(a*2);
			for( int j=0; j<2; j++ )
				aa[j] = relax == 1 ? V[j] : (1-relax)*aa[j] + relax*V[j];
		}
	}
}

void DenseCRF::expAndNormalize_att( float *att_out, float *att_in, float scale, float relax )
{
#pragma omp parallel for
	for (int i = 0; i < N_; i++ ){
		for( int a = 0; a < A_; a++){
			//(num_of_attributes*2)*(i*imagewidth+j)+(l+a*2)
			const float * b = att_in + (A_*2)*(i)+(a*2);			
			float mx = scale * max(b[0], b[1]);// Find the max and subtract it so that the exp doesn't explode
			float tt = 0;
			float V[2];
			for( int j=0; j<2; j++ ){
				V[j] = fast_exp( scale*b[j]-mx );
				tt += V[j];
			}
			// Make it a probability
			for( int j=0; j<2; j++ )
				V[j] /= tt;

			float * aa = att_out + (A_*2)*(i)+(a*2);
			for( int j=0; j<2; j++ )
				aa[j] = relax == 1 ? V[j] : (1-relax)*aa[j] + relax*V[j];
		}
	}
}

void DenseCRF::expAndNormalize ( float* out, const float* in, float scale, float relax ) 
{
#pragma omp parallel for
	for( int i=0; i<N_; i++ ){
		const float * b = in + i*M_;
		// Find the max and subtract it so that the exp doesn't explode
		float mx = scale*b[0];
		for( int j=1; j<M_; j++ )
			if( mx < scale*b[j] )
				mx = scale*b[j];
		float tt = 0;
		float *V = new float[M_];
		for( int j=0; j<M_; j++ ){
			V[j] = fast_exp( scale*b[j]-mx );
			tt += V[j];
		}
		// Make it a probability
		for( int j=0; j<M_; j++ )
			V[j] /= tt;

		float * a = out + i*M_;
		for( int j=0; j<M_; j++ )
			a[j] = relax == 1 ? V[j] : (1-relax)*a[j] + relax*V[j];
		delete[] V;
	}
}



void DenseCRF::expAndNormalize_cooc ( float *cooc_out, float *cooc_in, float scale, float relax ) 
{
#pragma omp parallel for
	for( int i=0; i<M_; i++ ){
		const float * b_cooc = cooc_in + i*2;
		// Find the max and subtract it so that the exp doesn't explode
		float mx_cooc = scale*b_cooc[0];
		for( int j=1; j < 2; j++ )
			if( mx_cooc < scale * b_cooc[j] )
				mx_cooc = scale*b_cooc[j];

		float tt = 0;
		float V_cooc[2];
		for( int j=0; j<2; j++ ){
			V_cooc[j] = fast_exp( scale*b_cooc[j]-mx_cooc );
			tt += V_cooc[j];
		}
		// Make it a probability
		for( int j=0; j<2; j++ )
			V_cooc[j] /= tt;

		float * a_cooc = cooc_out + i*2;
		for( int j=0; j<2; j++ )
			a_cooc[j] = V_cooc[j];		
	}
}



///////////////////
/////  Debug  /////
///////////////////

void DenseCRF::unaryEnergy(const short* ass, float* result) {
	for( int i=0; i<N_; i++ )
		if ( 0 <= ass[i] && ass[i] < M_ )
			result[i] = unary_[ M_*i + ass[i] ];
		else
			result[i] = 0;
}
void DenseCRF::pairwiseEnergy(const short* ass, float* result, int term) 
{
	float * current = allocate( N_*M_ );
	// Build the current belief [binary assignment]
	for( int i=0,k=0; i<N_; i++ )
		for( int j=0; j<M_; j++, k++ )
			current[k] = (ass[i] == j);

	for( int i=0; i<N_*M_; i++ )
		next_[i] = 0;
	if (term == -1)
		for( unsigned int i=0; i<pairwise_.size(); i++ )
			pairwise_[i]->apply( next_, current, tmp_, M_ );
	else
		pairwise_[ term ]->apply( next_, current, tmp_, M_ );
	for( int i=0; i<N_; i++ )
		if ( 0 <= ass[i] && ass[i] < M_ )
			result[i] =-next_[ i*M_ + ass[i] ];
		else
			result[i] = 0;
	deallocate( current );
}


void DenseCRF::startInference()
{
	if(addCooc)
	{
		int *total_num_labels = new int[M_];
		memset(total_num_labels, 0, M_);
		for(int i = 0; i < N_; i++)	{
			int class_label = 0; 
			float temp_unary_cost = unary_[i*M_]; 
			for(int j = 1; j < M_; j++)	{
				if(temp_unary_cost < unary_[i*M_+j]){
					temp_unary_cost = unary_[i*M_+j]; 
					class_label = j; 
				}
			}
			total_num_labels[class_label]++; 
		}

		float pairwise_cooc = 0.0; // float p1, p2, p, p12; 
		for(int i = 0; i < M_; i++)
		{
			if(total_num_labels[i] > 0)
			{
				next_cooc_[2*i+1] = total_num_labels[i];
				next_cooc_[2*i] = 1;
			}
			else
			{
				next_cooc_[2*i+1] = 1;
				next_cooc_[2*i] = 100;
			}
		}

		delete []total_num_labels; 
	}
	// Initialize using the unary energies
	expAndNormalize( current_, unary_, -1 );

	if(addCooc)
	{
		expAndNormalize_cooc ( current_cooc_, next_cooc_) ; 
	}


	//start attributes
	//float *unary_att_, *pair_att_, *pair_join_att_, *current_att_, *next_att_, *tmp_att_, *additional_unary_att_;
	//float att_factor_joint,att_fac
	if(addAtt)
		expAndNormalize_att( current_att_,unary_att_, -1);

	//end attributes
	if(addAttDet){
		expAndNormalize_ojbdet( current_det_,unary_objdet_, -1);
		expAndNormalize_attdet( current_attdet_,unary_attdet_, -1);
	}
}


void DenseCRF::stepInference( float relax )
{
	// Set the unary potential
#pragma omp parallel for
	for( int i=0; i<N_*M_; i++ )
		next_[i] = -unary_[i] - additional_unary_[i];

	//start PN Potts
	if(addHO){
		set_mem_init_higherorder(); 
		int num_of_layers = num_layers; 
		for(int i = 0; i < num_of_layers; i++)
			calculate_higherorder_pot(i); 

		//HO->calculate_higherorder_pot(i, current, next); 
		float ho_weight_ = 1.0;
#pragma omp parallel for
		for(int i = 0; i < N_*M_; i++)
			next_[i] = next_[i] - ho_weight_ * higher_order[i]; 
	}

	// start add co-occurrence terms
	if(addCooc){
		int *higher_labels = new int[M_];
		for(int i = 0; i < M_; i++)
			higher_labels[i] = current_cooc_[2*i] < current_cooc_[2*i+1] ? 1 : 0;
		float *temp_prob_mult = new float[2*M_];
#pragma omp parallel for
		for(int i = 0; i < M_; i++){
			float mult_prob = 0.0, mult_prob1 = 0.0; 
			for(int j = 0; j < N_; j++){
				mult_prob += 1.0-current_[j*M_+i];
				mult_prob1 += current_[j*M_+i]; 
			}
			temp_prob_mult[2*i] = mult_prob; 
			temp_prob_mult[2*i+1] = mult_prob1;
			if(temp_prob_mult[2*i] < 1e-4) 
				temp_prob_mult[2*i] = 1e-4; 
			if(temp_prob_mult[2*i+1] < 1e-4) 
				temp_prob_mult[2*i+1] = 1e-4; 
		}

		float pairwise_cooc = 0.0;
		for(int i = 0; i < M_; i++){
			float p1, p2, p, p12; 
			pairwise_cooc = 0.0; 
			p1 = unary_cooc_[i]; 
			for(int j = 0; j < M_; j++){
				p2 = unary_cooc_[j]; 
				p12 = pair_cooc_[i*M_+j]; 
				p = 1 - (1 - p12 / p2) * (1 - p12 / p1);
				if(p > 1) 
					p = 1;
				if(p < 1e-6) 
					p = 1e-6;
				if(i != j)
					pairwise_cooc = pairwise_cooc - (0.005*N_) * log(p) * current_cooc_[j*2+1];
			}
			next_cooc_[2*i+1] = -pairwise_cooc; 		
		}

		for(int i = 0; i < M_; i++)
			next_cooc_[2*i] = -1.0*cooc_factor*(temp_prob_mult[2*i+1]); 
#pragma omp parallel for
		for(int i = 0; i < N_; i++)
			for(int j = 0; j < M_; j++)
				next_[i*M_+j] -= cooc_factor*current_cooc_[2*j];

		delete []temp_prob_mult; 
		delete []higher_labels; 
	}

	// start det
	if(addDet){
		//det_set_mem_init_higherorder(); 
		det_calculate_higherorder_pot(); 
#pragma omp parallel for
		for(int i = 0; i < N_*M_; i++)
			next_[i] = next_[i] - det_higher_order[i]; 
	}

	// start attributes
	if (addAtt){
#pragma omp parallel for
		for( int i=0; i<N_*A_*2; i++ )
			next_att_[i] = -unary_att_[i];//- additional_unary_att_[i];

		// joint unary, objs -> att
#pragma omp parallel for
		for( int i=0; i<N_; i++ ) 
			for( int a=0; a<A_; a++ )
				for( int l=0; l<M_; l++ )
					next_att_[i*A_*2 + a*2] -= current_[i*M_+l] * pair_join_att_[a*M_ + l] * att_factor_joint;

		// joint unary, Att -> att
#pragma omp parallel for
		for( int i=0; i<N_; i++ ) 
			for( int a=0; a<A_; a++ ) 
				for( int b=0; b<2; b++ )
					for( int a1=0; a1<A_; a1++ )
						for( int b1=0; b1<2; b1++ )
							if(a1!=a && b1!= b)
								next_att_[i*A_*2 + a*2 + b] -= current_att_[i*A_*2 + a1*2 + b1] * pair_att_[a*A_ + a1] * att_factor_correlation;

		// pairwise, att
		for( unsigned int i=0; i<pairwise_.size(); i++)
			pairwise_[i]->apply( next_att_, current_att_, tmp_att_, A_*2);

		// joint unary, attr -> objs
#pragma omp parallel for
		for( int i=0; i<N_; i++ ) 
			for( int l=0; l<M_; l++ ) 
				for( int a=0; a<A_; a++ )
					next_[i*M_+l] -= current_att_[i*A_*2 + a*2 + 0] * pair_join_att_[a*M_ + l] * att_factor_joint;

		expAndNormalize_att(current_att_,next_att_,1.0,relax);
	}

	// start attdet
	if(addAttDet){
#pragma omp parallel for
		for( int i=0; i<NS_*AD_*2; i++ )
			next_attdet_[i] = -unary_attdet_[i];//- additional_unary_att_[i];
#pragma omp parallel for
		for( int i=0; i<NS_*2;i++)
			next_det_[i] = -unary_objdet_[i];

		// joint unary, objdets -> attdet
#pragma omp parallel for
		for( int i=0; i<NS_; i++ ){
			for( int a=0; a<AD_; a++ ){
				int l = det_stats_potential[i];
				next_attdet_[i*AD_*2 + a*2] -= current_det_[i*2+1] * pair_join_attdet_[a*M_ + l] * attdet_factor_joint;
			}
		}


		// joint unary, Attdet -> attdet
#pragma omp parallel for
		for( int i=0; i<NS_; i++ ) 
			for( int a=0; a<AD_; a++ ) 
				for( int b=0; b<2; b++ )
					for( int a1=0; a1<AD_; a1++ )
						for( int b1=0; b1<2; b1++ )
							if((!(a1==a))&&(!(b1==b)))
								next_attdet_[i*AD_*2 + a*2 + b] -= current_attdet_[i*AD_*2 + a1*2 + b1] * pair_attdet_[a*AD_ + a1] * attdet_factor_correlation;

		// joint unary, attd -> objdet
#pragma omp parallel for
		for( int i=0; i<NS_; i++ ){
			int l = det_stats_potential[i];
			for( int a=0; a<AD_; a++ )
				next_det_[i*2+1] -= current_attdet_[i*AD_*2 + a*2 + 0] * pair_join_attdet_[a*2 + l] * attdet_factor_joint;
		}

		// joint unary, objectdet - > object
#pragma omp parallel for
		for (int i = 0; i<NS_; i++){
			int l = det_stats_potential[i];
			int basesegmentcounts = det_baseSegmentCounts[i];
			for (int j=0; j<basesegmentcounts; j++){
				int curr_pix_index = det_baseSegmentIndexes[i][j];
				next_det_[i*2+1] -=  current_[curr_pix_index*M_+l]*(-obj_det_factor);
				next_[curr_pix_index*M_+l] -= current_det_[i*2+1]*(-obj_det_factor);
			}
		}

		expAndNormalize_attdet(current_attdet_, next_attdet_, 1.0, relax);
		expAndNormalize_ojbdet(current_det_,next_det_,1.0,relax);
	}

	// pairwise potentials
	for( unsigned int i=0; i<pairwise_.size(); i++)
		pairwise_[i]->apply( next_, current_, tmp_, M_);


	// Exponentiate and normalize
	expAndNormalize( current_, next_, 1.0, relax );
	if(addCooc)
		expAndNormalize_cooc ( current_cooc_, next_cooc_) ; 
}
void DenseCRF::currentMap( short * result )
{
	// Find the map
	for( int i=0; i<N_; i++ ){
		const float * p = current_ + i*M_;
		// Find the max and subtract it so that the exp doesn't explode
		float mx = p[0];
		int imx = 0;
		for( int j=1; j<M_; j++ )
			if( mx < p[j] ){
				mx = p[j];
				imx = j;
			}
			result[i] = imx;
	}
}


// start co-occurrence

void DenseCRF::set_hocooc(int add_Cooc)
{
	addCooc = 0;
	if(add_Cooc) 
	{
		addCooc = 1;	

		unary_cooc_ = new float[M_];
		memset(unary_cooc_, 0, sizeof(float)*M_);
		pair_cooc_ = new float[M_*M_];
		memset(pair_cooc_, 0, sizeof(float)*M_*M_);
		current_cooc_ = allocate(2*M_); 
		next_cooc_ = allocate(2*M_);	
	}
}

void DenseCRF::setco_occurrence(float *cooc_unary, float *cooc_pairwise,  float coocFactor) 
{
	if(addCooc)
	{
		memcpy( unary_cooc_, cooc_unary, M_*sizeof(float) );
		for(int i = 0; i < M_; i++)
			unary_cooc_[i] = unary_cooc_[i]/600.0;

		memcpy(pair_cooc_, cooc_pairwise, M_*M_*sizeof(float));

		for(int i = 0; i < M_*M_; i++)
			pair_cooc_[i] = pair_cooc_[i]/600.0;

		cooc_factor = coocFactor; 
	}
}

// end co-occurrence
// start attributes
void DenseCRF::set_attributes(float *att_unary, float* att_pairwise, float* att_joint_pairwise,float weight_joint, float weight_correlation) 
{
	addAtt = 1;
	memcpy( unary_att_, att_unary, N_*2*A_*sizeof(float) );
	memcpy(pair_att_, att_pairwise, A_*A_*sizeof(float));
	memcpy(pair_join_att_, att_joint_pairwise, A_*M_*sizeof(float));

	att_factor_joint = weight_joint;
	att_factor_correlation = weight_correlation;
}

// end atttributes
// start attdet
void DenseCRF::set_attdet(float *obj_unary,float *att_unary, float* att_pairwise, float* att_joint_pairwise,float weight_joint, float weight_correlation, float weight_objdet)
{
	addAttDet=1;
	memcpy(unary_objdet_,obj_unary,NS_*2*sizeof(float));

	memcpy(unary_attdet_, att_unary, NS_*2*AD_*sizeof(float) );

	memcpy(pair_attdet_, att_pairwise, AD_*AD_*sizeof(float));

	memcpy(pair_join_attdet_, att_joint_pairwise, AD_*M_*sizeof(float));


	attdet_factor_joint = weight_joint;
	attdet_factor_correlation = weight_correlation;
	obj_det_factor = weight_objdet;
}
// end attdet

// start adding higher order potential related stuffs

void DenseCRF::set_ho(int add_ho)
{
	addHO = 0;
	if(add_ho) 
		addHO = 1; 		
}

void DenseCRF::ho_mem_init(int imagewidth, int imageheight, char **layers_name, int num_of_layers, char *ho_stats, char *ho_seg_ext, char *ho_sts_ext, float hoParam1, float hoParam2)
{
	if(addHO)
	{
		num_layers = num_of_layers; 
		layers_dir_name = new char*[num_layers]; 
		for(int i = 0; i < num_layers; i++)
			layers_dir_name[i] = layers_name[i]; 

		stats_pot = ho_stats; 

		ho_segment_ext = ho_seg_ext; 
		ho_stats_ext = ho_sts_ext; 

		ho_width = imagewidth; ho_height = imageheight; 
		higher_order = new float[ho_width*ho_height*M_];

		for(int i = 0; i < ho_width*ho_height*M_; i++)
		{
			higher_order[i] = 0.0; 
		}

		temp_segmentimagedata = new int[N_];
		segmentimagedata = new int[N_];
		baseSegmentIndexes = new int**[num_layers]; 
		segmentCount = new int[num_layers];
		baseSegmentCounts = new int*[num_layers];
		ho_param1 = hoParam1; 
		ho_param2 = hoParam2; 
	}
}

void DenseCRF::readSegments(char *filename)
{
	if (addHO)
	{
		char *seg_file_name = new char[200];

		// read segments
		for(int i = 0; i < num_layers; i++)
		{
			//if(i > 5)
			//ho_segment_ext = "msh"; 
			sprintf(seg_file_name, "%s%s.%s", layers_dir_name[i], filename, ho_segment_ext);
			readSegmentIndex(seg_file_name, i);
		}

		// read stats potential
		sprintf(seg_file_name, "%s%s.%s", stats_pot, filename, ho_stats_ext);
		readStatsPotential(seg_file_name, num_layers); 

		delete []seg_file_name; 
	}
}

void DenseCRF::readSegmentIndex(char *file_name, int layer)
{
	if(addHO)
	{
		FILE *fin = fopen(file_name, "rb");
		int width, height, bands; 
		fread(&width, 4, 1, fin); 
		fread(&height, 4, 1, fin); 
		fread(&bands, 4, 1, fin);
		fread(temp_segmentimagedata, sizeof(int), width * height * bands, fin);
		fclose(fin);

		for(int i = 0; i < height; i++)
			for(int j = 0; j < width; j++)
				segmentimagedata[i*width+j] = temp_segmentimagedata[(height-(i+1))*width+j];

		segmentCount[layer] = 0;
		int points = width*height; 
		for(int j = 0; j < points; j++) 
			if(segmentimagedata[j]+1 > segmentCount[layer]) 
				segmentCount[layer] = segmentimagedata[j] + 1;

		baseSegmentCounts[layer] = new int[segmentCount[layer]];
		memset(baseSegmentCounts[layer], 0, segmentCount[layer] * sizeof(int));

		for(int j = 0; j < points; j++) 
			baseSegmentCounts[layer][segmentimagedata[j]]++;

		int temp_count = 0; 

		for(int i = 0; i < segmentCount[layer]; i++)
			temp_count += baseSegmentCounts[layer][i]; 

		baseSegmentIndexes[layer] = new int *[segmentCount[layer]];
		for(int j = 0; j < segmentCount[layer]; j++) 
			baseSegmentIndexes[layer][j] = new int[baseSegmentCounts[layer][j]];

		segmentationIndex = new int[segmentCount[layer]];
		memset(segmentationIndex, 0, segmentCount[layer] * sizeof(int));

		for(int j = 0; j < points; j++)
		{
			baseSegmentIndexes[layer][segmentimagedata[j]][segmentationIndex[segmentimagedata[j]]] = j;
			segmentationIndex[segmentimagedata[j]]++;
		}
		delete []segmentationIndex;
	}
}

void DenseCRF::readStatsPotential(char *file_name, int num_of_layers)
{
	if(addHO)
	{
		int totalsegmentCount = 0; 

		for(int i = 0; i < num_of_layers; i++)
		{
			totalsegmentCount += segmentCount[i]; 
		}

		stats_potential = new double[totalsegmentCount*M_];
		for(int i = 0; i < totalsegmentCount*M_; i++)
			stats_potential[i] = 0.0; 

		FILE *fin1 = fopen(file_name, "rb");
		fread(stats_potential, sizeof(double), totalsegmentCount*M_, fin1);
		fclose(fin1); 

		float sum = 0; int start_loc = 0; 
		for(int i1 = 0; i1 < num_of_layers; i1++)
		{
			for(int i = 0; i < segmentCount[i1]; i++)
			{
				sum = 0.0; 
				for(int k = 0; k < M_; k++) 
				{
					sum += exp(stats_potential[start_loc+k]);
				}
				for(int k = 0; k < M_; k++) 
					stats_potential[start_loc+k] = exp(stats_potential[start_loc+k]) / sum; 
				start_loc = start_loc + M_; 
			}
		}
	}
}


void DenseCRF::calculate_higherorder_pot(int layer)
{
	if(addHO)
	{
		h_norm = new double[segmentCount[layer]*M_];

		float norm_val = 0.0; 
		int basesegmentcounts = 0; 
		int curr_pix_label = 0, curr_pix_index; // int x, y; 
		//	int neigh_pix_index, neigh_pix_label; 

		double higher_order_prob; 

		for(int i = 0; i < segmentCount[layer]; i++)
			for(int j = 0; j < M_; j++)
				h_norm[i*M_+j] = 1.0; 

		for(int i = 0; i < segmentCount[layer]; i++)
		{
			basesegmentcounts = baseSegmentCounts[layer][i];
			higher_order_prob = 1.0;
			for(int j = 0; j < M_; j++)
			{
				higher_order_prob = 1.0; 
				for(int k = 0; k < basesegmentcounts; k++)
				{
					curr_pix_index = baseSegmentIndexes[layer][i][k];
					higher_order_prob = higher_order_prob * current_[curr_pix_index*M_+j]; 
				}
				h_norm[i*M_+j] = higher_order_prob; 
			}
		}


		double alpha = 0.5, maxcost, weight, costdata = 0.0; int start_loc = 0; 

		for(int i = 0; i < layer; i++)
			start_loc = start_loc + segmentCount[i]; 

		start_loc = start_loc * M_; 

		for(int i = 0; i < segmentCount[layer]; i++)
		{
			basesegmentcounts = baseSegmentCounts[layer][i];

			weight = 0.3 * basesegmentcounts; 
			maxcost = -weight * log(alpha); 

			for(int j = 0; j < basesegmentcounts; j++)
			{
				curr_pix_index = baseSegmentIndexes[layer][i][j]; 
				for(int k = 0; k < M_; k++)
				{
					higher_order_prob = h_norm[i*M_+k]/(current_[curr_pix_index*M_+k]+0.0001);
					costdata = - weight * log(stats_potential[start_loc+k]); 
					higher_order[curr_pix_index*M_+k] += (ho_param1*costdata - ho_param2*higher_order_prob);
				}		
			}	
			start_loc = start_loc+M_; 
		}
		delete []h_norm; 	
	}
}

void DenseCRF::set_mem_init_higherorder()
{
	if(addHO)
	{
		for(int i = 0; i < N_*M_; i++)
			higher_order[i] = 0.0; 
	}
}


void DenseCRF::del_mem_higherorder()
{
	if(addHO)
	{
		delete []higher_order;
		delete []segmentimagedata;
		delete []temp_segmentimagedata;
		for(int i = 0; i < num_layers; i++)
		{
			for(int j = 0; j < segmentCount[i]; j++)
			{
				delete []baseSegmentIndexes[i][j];
			}
			delete []baseSegmentIndexes[i]; 
		}

		delete []baseSegmentIndexes; 

		for(int i = 0; i < num_layers; i++)
		{
			delete []baseSegmentCounts[i]; 
		}
		delete []baseSegmentCounts;
		delete[]segmentCount; 
		delete []stats_potential; 		
	}
	if(addCooc)
	{
		delete []unary_cooc_; 
		delete []pair_cooc_; 	
	}
	if(addDet)
	{
		delete []det_higher_order; 
		delete []det_baseSegmentCounts;
		delete []det_baseSegmentIndexes;
		delete []det_stats_potential; 
		delete []det_h_norm; 
	}
}

void DenseCRF::add_higher_order()
{

}

// end higher order potential related stuffs


// start det 

void DenseCRF::set_hodet(int add_det)
{
	addDet = 0; 
	if(add_det) addDet = 1; 
}

void DenseCRF::det_ho_mem_init(int imagewidth, int imageheight, char *det_seg, char *det_bb, char *det_seg_ext, char *det_bb_ext, float detParam1, float detParam2)
{
	if(addDet)
	{
		det_seg_dir = det_seg; 
		det_bb_dir = det_bb; 
		det_segment_ext = det_seg_ext; 
		det_bbox_ext = det_bb_ext; 
		det_ho_width = imagewidth; 
		det_ho_height = imageheight; 
		det_higher_order = new float[det_ho_width*det_ho_height*M_];
		for(int i = 0; i < det_ho_width*det_ho_height*M_; i++)
		{
			det_higher_order[i] = 0.0; 
		}

		det_param1 = detParam1; 
		det_param2 = detParam2; 
	}	
}

void DenseCRF::det_readSegmentIndex(char *fileName)
{
	if(addDet)
	{
		char *file_name = new char[200]; 
		sprintf(file_name, "%s%s.%s",det_seg_dir, fileName, det_segment_ext);
		FILE *fin = fopen(file_name, "rb");
		fread(&det_segmentCount, sizeof(int), 1, fin);
		det_baseSegmentCounts = new int[det_segmentCount];
		det_baseSegmentIndexes = new int *[det_segmentCount];
		int **det_tempbaseSegmentIndexes = new int*[det_segmentCount];
		for(int i = 0; i < det_segmentCount; i++)
		{
			fread(&det_baseSegmentCounts[i], sizeof(int), 1, fin);
			det_baseSegmentIndexes[i] = new int[det_baseSegmentCounts[i]];
			det_tempbaseSegmentIndexes[i] = new int[det_baseSegmentCounts[i]];
			if(det_baseSegmentCounts[i] != 0) fread(det_tempbaseSegmentIndexes[i], sizeof(int), det_baseSegmentCounts[i], fin);

			// for apascal
			/*		for(int j = 0; j < det_baseSegmentCounts[i]; j++)
			{
			int temp_id1 = det_tempbaseSegmentIndexes[i][j], x1, y1; 
			x1 = temp_id1%det_ho_width; y1 = temp_id1/det_ho_width; 
			det_baseSegmentIndexes[i][j] = (det_ho_height-(y1+1))*det_ho_width+x1; 
			}*/

			// for our segments not for apascal
			for(int j = 0; j < det_baseSegmentCounts[i]; j++)
			{
				det_baseSegmentIndexes[i][j] = det_tempbaseSegmentIndexes[i][j]; 
			}
		}
		fclose(fin);

		sprintf(file_name, "%s%s.%s", det_bb_dir, fileName, det_bbox_ext);
		unsigned char det_segmented; int det_objects; 
		FILE *fin1 = fopen(file_name, "rb");

		fread(&det_segmented, sizeof(unsigned char), 1, fin1);
		fread(&det_objects, sizeof(int), 1, fin1);

		det_stats_potential = new int[det_objects];

		int x1, y1, x2, y2; unsigned char type; 

		for(int i = 0; i < det_objects; i++)
		{
			fread(&type, sizeof(unsigned char), 1, fin1); 
			fread(&x1, sizeof(int), 1, fin1); fread(&y1, sizeof(int), 1, fin1); fread(&x2, sizeof(int), 1, fin1); fread(&y2, sizeof(int), 1, fin1);
			det_stats_potential[i] = (int)type; 
		}	

		det_resp = new double[det_objects];
		if(det_objects > 0) fread(det_resp, sizeof(double), det_objects, fin1);
		fclose(fin1);	

		det_h_norm = new double[det_segmentCount*M_];
	}
}

void DenseCRF::det_calculate_higherorder_pot()
{
	if(addDet)
	{
		float norm_val = 0.0; 

		int basesegmentcounts = 0; 
		int curr_pix_label = 0, curr_pix_index; //int x, y; 
		//int neigh_pix_index, neigh_pix_label; 

		double higher_order_prob; 

		for(int i = 0; i < det_segmentCount; i++)
			for(int j = 0; j < M_; j++)
				det_h_norm[i*M_+j] = 1.0; 

		for(int i = 0; i < det_segmentCount; i++)
		{
			basesegmentcounts = det_baseSegmentCounts[i];
			higher_order_prob = 1.0;
			for(int j = 0; j < M_; j++)
			{
				higher_order_prob = 1.0; 
				for(int k = 0; k < basesegmentcounts; k++)
				{
					curr_pix_index = det_baseSegmentIndexes[i][k];
					higher_order_prob = higher_order_prob * current_[curr_pix_index*M_+j]; 
				}
				det_h_norm[i*M_+j] = higher_order_prob; 
			}
		}


		for(int i = 0; i < det_ho_width*det_ho_height*M_; i++)
		{
			det_higher_order[i] = 0.0; 
		}

		double alpha = 0.5, maxcost, weight, costdata = 0.0;

		for(int i = 0; i < det_segmentCount; i++)
		{
			basesegmentcounts = det_baseSegmentCounts[i];

			weight = 0.3 * basesegmentcounts; 
			maxcost = -weight * log(alpha); 
			costdata = 6.0 * basesegmentcounts * (det_resp[i]+1.2);

			if(costdata < 0) costdata = 0; 

			for(int j = 0; j < basesegmentcounts; j++)
			{
				curr_pix_index = det_baseSegmentIndexes[i][j]; 

				for(int k = 0; k < M_; k++)
				{
					higher_order_prob = det_h_norm[i*M_+k]/(current_[curr_pix_index*M_+k]+0.0001);
					if(k != det_stats_potential[i])
					{
						det_higher_order[curr_pix_index*M_+k] += det_param1*costdata-det_param2*higher_order_prob; 
					}				
				}
			}
		}
	}
}

// end det

float* allocate(size_t N) {
	float * r = NULL;
	if (N>0)
#ifdef SSE_DENSE_CRF
		r = (float*)_mm_malloc( N*sizeof(float)+16, 16 );
#else
		r = new float[N];
#endif
	memset( r, 0, sizeof(float)*N);
	return r;
}
void deallocate(float*& ptr) {
	if (ptr)
#ifdef SSE_DENSE_CRF
		_mm_free( ptr );
#else
		delete[] ptr;
#endif
	ptr = NULL;
}