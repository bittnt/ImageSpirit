#ifndef __densecrf_h_
#define __densecrf_h_

#pragma once


class PairwisePotential{
public:
	virtual ~PairwisePotential();
	virtual void apply( float * out_values, const float * in_values, float * tmp, int value_size ) const = 0;
};
class SemiMetricFunction{
public:
	virtual ~SemiMetricFunction();
	// For two probabilities apply the semi metric transform: v_i = sum_j mu_ij u_j
	virtual void apply( float * out_values, const float * in_values, int value_size ) const = 0;
};


class DenseCRF
{
protected:
	friend class BipartiteDenseCRF;

	// Number of variables and labels
	int N_, M_, A_, NS_,AD_;
	float *unary_, *additional_unary_, *current_, *next_, *tmp_;

	// start cooc

	float *unary_cooc_, *pair_cooc_, *current_cooc_, *next_cooc_, *tmp_cooc_; 
	float cooc_factor; 

	// end cooc

	// start attributes

	float *unary_att_, *pair_att_, *pair_join_att_, *current_att_, *next_att_, *tmp_att_, *additional_unary_att_;
	float att_factor_joint,att_factor_correlation;
	// end atttributes

	// start object attdet
	float *unary_attdet_, *pair_attdet_,*pair_join_attdet_,*unary_objdet_;
	float *current_det_,*next_det_,*tmp_det_,*current_attdet_,*next_attdet_,*tmp_attdet_;
	float attdet_factor_joint,attdet_factor_correlation,obj_det_factor;
	// end object attdet

	// Store all pairwise potentials
	std::vector<PairwisePotential*> pairwise_;

	// Run inference and return the pointer to the result
	float* runInference( int n_iterations, /*float *un_normalized_val,*/ float relax);
	void runInference(int n_iterations,float relax,float* current_, float* current_att_, 
		int hmap_objectinteraction, char* saliencyfolder,int inputobjectclass,int threshold,
		int hmap_attributeinteraction, int inputattributeclass);

	// Auxillary functions
	void expAndNormalize( float* out, const float* in, float scale = 1.0, float relax = 1.0 );
	void expAndNormalize_cooc( float *out_cooc_, float *in_cooc_, float scale = 1.0, float relax = 1.0 );
	void expAndNormalize_att ( float *att_out, float *att_in, float scale = 1.0, float relax =1.0 );
	void expAndNormalize_attdet ( float *attdet_out, float *attdet_in, float scale = 1.0, float relax =1.0 );
	void expAndNormalize_ojbdet(float *objdet_out, float *objdet_in, float scale = 1.0, float relax =1.0 );
	// Don't copy this object, bad stuff will happen
	DenseCRF( DenseCRF & o ){}
public:
	// Create a dense CRF model of size N with M labels
	DenseCRF( int N, int M );
	// Create a dense CRF model of size N with M labels and A attributes labels
	DenseCRF( int N, int M , int A);

	// Create a 2d dense CRF model of size W x H with M object labels and A attribute labels, NS number of segments, number of attdet
	DenseCRF(int N, int M, int A, int NS, int AD);


	virtual ~DenseCRF();
	// Add  a pairwise potential defined over some feature space
	// The potential will have the form:    w*exp(-0.5*|f_i - f_j|^2)
	// The kernel shape should be captured by transforming the
	// features before passing them into this function
	void addPairwiseEnergy( const float * features, int D, float w=1.0f, const SemiMetricFunction * function=NULL );

	// Add your own favorite pairwise potential (ownwership will be transfered to this class)
	void addPairwiseEnergy( PairwisePotential* potential );

	// Set the unary potential for all variables and labels (memory order is [x0l0 x0l1 x0l2 .. x1l0 x1l1 ...])
	// void setUnaryEnergy( const float * unary );
	void setUnaryEnergy( const float * unary/*, float *cooc_unary, float *cooc_pairwise*/);
	// start atttributes
	void setUnaryEnergy ( const float* unary, const float* att_unary,float* att_pairwise,float* att_joint_pairwise/*,  float *cooc_unary, float *cooc_pairwise*/);
	// end attributes
	// start attdet
	void setUnaryEnergy ( const float* unary, const float* unary_objdet, const float* att_unary, const float* attdet_unary,
		float* att_pairwise,float* att_joint_pairwise,  float *attdet_pairwise, float *attdet_joint_pairwise);
	// end attdet
	// Set the unary potential for a specific variable
	void setUnaryEnergy( int n, const float * unary );

	// Run inference and return the probabilities
	void inference( int n_iterations, float* result, float relax=1.0 );


	// Run MAP inference and return the map for each pixel
	void map( int n_iterations, short* result,short* result_att, short* result_attdet, float relax=1.0);

	// Step by step inference
	void startInference();
	void stepInference( /*float *un_normalized_val,*/ float relax = 1.0 );
	void currentMap( short * result );

	// start cooccurrence 
	void set_hocooc(int);
	void setco_occurrence(float *cooc_unary, float *cooc_pairwise, float coocFactor);
	char addCooc; 
	// end cooccurrence
	// start attributes
	char addAtt;
	//void set_hoattributes(int);
	void set_attributes(float *att_unary, float* att_pairwise, float* att_joint_pairwise,float weight_joint, float weight_correlation);
	// end attributes

	void set_attdet(float *obj_unary,float *att_unary, float* att_pairwise, float* att_joint_pairwise,float weight_joint,
					float weight_correlation, float weight_objdet);

	// start adding higher order reltated stuffs
	void set_ho(int);
	char addHO;
	int ho_width, ho_height; 
	float *higher_order;
	char *ho_segment_ext; 
	char *ho_stats_ext; 
	int *segmentIndex, *segmentCount, **baseSegmentCounts, ***baseSegmentIndexes, *segmentationIndex, *segmentimagedata, *temp_segmentimagedata; 
	double *stats_potential, *h_norm; 
	float ho_param1, ho_param2; 
	void readSegments(char *filename); 
	void calculate_higherorder_pot(int layer); 
	void readSegmentIndex(char *, int); 
	void readStatsPotential(char *file_name, int num_of_layers); 
	void add_higher_order(); 
	void mem_init_higherorder(); 
	void set_mem_init_higherorder(); 
	void del_mem_higherorder(); 
	void ho_mem_init(int, int, char **, int, char *, char *, char *, float, float); 
	char **layers_dir_name; 
	char *stats_pot; 
	int num_layers; 
	// end higher order related stuffs

	// add attributes at the region level to affect the dection term.
	// start attributes at region level
	char addAttDet;
	char num_of_attdet;

	// end attributes at region level
	// det start
	// start adding higher order reltated stuffs
	void set_hodet(int);
	char addDet; 
	char *det_seg_dir, *det_bb_dir; 
	char *det_segment_ext, *det_bbox_ext; 
	int det_ho_width, det_ho_height; 
	float *det_higher_order, det_param1, det_param2; 
	double *det_h_norm, *det_resp;
	int *det_segmentIndex, det_segmentCount, *det_baseSegmentCounts, **det_baseSegmentIndexes, *det_stats_potential; 

	void det_calculate_higherorder_pot(); 
	void det_readSegmentIndex(char *); 
	void det_set_mem_init_higherorder(); 
	void det_del_mem_higherorder(); 
	void det_ho_mem_init(int, int, char *, char *, char *, char *, float, float); 

	//void det_mem_init_higherorder(); 
	// det end

public: /* Debugging functions */
	// Compute the unary energy of an assignment
	void unaryEnergy( const short * ass, float * result );

	// Compute the pairwise energy of an assignment (half of each pairwise potential is added to each of it's endpoints)
	void pairwiseEnergy( const short * ass, float * result, int term=-1 );
};

class DenseCRF2D:public DenseCRF
{
protected:
	// Width, height of the 2d grid
	int W_, H_;
public:
	// Create a 2d dense CRF model of size W x H with M labels
	DenseCRF2D( int W, int H, int M);
	// Create a 2d dense CRF model of size W x H with M object labels and A attribute labels
	DenseCRF2D( int W, int H, int M, int A);

	// Create a 2d dense CRF model of size W x H with M object labels and A attribute labels, NS number of segments, number of attdet
	DenseCRF2D(int W, int H, int M, int A, int NS, int AD);

	virtual ~DenseCRF2D();
	// Add a Gaussian pairwise potential with standard deviation sx and sy
	void addPairwiseGaussian( float sx, float sy, float w, const SemiMetricFunction * function=NULL );

	// Add a Bilateral pairwise potential with spacial standard deviations sx, sy and color standard deviations sr,sg,sb
	void addPairwiseBilateral( float sx, float sy, float sr, float sg, float sb, const unsigned char * im, float w, const SemiMetricFunction * function=NULL );

	// Run MAP inference and return the object index map for each pixel in 1u format.
	// also return attributes if resAtts1u.size() == #Attributes
	Mat map(int nIterations, vecM &resAtts1u, float relax = 1.0);

	// Set the unary potential for a specific variable
	void setUnaryEnergy( int x, int y, const float * unary );
	using DenseCRF::setUnaryEnergy;
};


// This function defines a simplified interface to the permutohedral lattice
// We assume a filter standard deviation of 1
class Permutohedral;
class Filter{
protected:
	int n1_, o1_, n2_, o2_;
	Permutohedral * permutohedral_;
	// Don't copy
	Filter( const Filter& filter ){}
public:
	// Use different source and target features
	Filter( const float * source_features, int N_source, const float * target_features, int N_target, int feature_dim );
	// Use the same source and target features
	Filter( const float * features, int N, int feature_dim );
	//
	~Filter();
	// Filter a bunch of values
	void filter( const float * source, float * target, int value_size );
};


#endif
