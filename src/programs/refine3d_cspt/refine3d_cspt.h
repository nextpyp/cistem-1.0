


#ifdef __cplusplus
extern "C"{
#endif

//========================================
typedef struct MRCFile MRCFile; 


MRCFile* ReadStack( char * filename, bool overwrite );
MRCFile* ReadReconstruction( char * filename, bool overwrite, bool wait_for_file_to_exist );
int ReturnX( MRCFile* mrc );
int ReturnY( MRCFile* mrc );
int ReturnZ( MRCFile* mrc );
void MRCFilePrintInfo( MRCFile* mrc );
void CloseMRCFile( MRCFile* mrc );

// ========================================
typedef struct ReconstructedVolume ReconstructedVolume;

ReconstructedVolume* Init3DVolume( MRCFile* mrc, float pixel_size, char* wanted_symmetry_symbol, float molecular_mass_kDa, float outer_mask_radius, float low_resolution_limit, float high_resolution_limit );
float GetPixelSize( ReconstructedVolume* vol );

void Delete3DVolume( ReconstructedVolume* vol );
// ========================================
typedef struct ResolutionStatistics ResolutionStatistics;

ResolutionStatistics* CreateStat( int stack_x, float pixel_size, float molecular_mass_kDa );

void DeleteStat( ResolutionStatistics* stat );


// =======================================
typedef struct Image Image;

Image* CreateImage( MRCFile* stack, float index_in_stack, int stack_x, int stack_y, float padding );
Image** WeightedAverages(Image** particle_images, double ** shifts, int number_of_particles, int num_frames, int boxsize, float pixel_size, bool frame_refine, double * score_weights, double weight_width, bool write_after_average, char output_file[1000]);
void UpdateWeightedAverages(Image** particle_images, Image** weighted_average, double ** shifts, int number_of_particles, int num_frames, int boxsize, float pixel_size, double weight_width);
float GetMean( Image* image );
void NormalizeImages( Image** particle_images, int start_ind, int end_ind, int boxsize, float outer_mask_radius, float pixel_size);
void NormalizeFrames( Image** particle_images, double coordinates[][6], int start_ind, int end_ind, int boxsize, float outer_mask_radius, float pixel_size, int num_frames);
void DeleteImage( Image* image );
void DeleteImageArray( Image** images, int number_of_images );

// ========================================
typedef struct ImageProjectionComparison ImageProjectionComparison;

ImageProjectionComparison* CreateComparisonObject(ReconstructedVolume* ref_3d,
                                                      Image** image_set,
                                                      float* input_parameters,
                                                      float pixel_size, 
                                                      int boxsize,
                                                      float volt, 
                                                      float cs, 
                                                      float am, 
                                                      float mask,
                                                      float highres,
                                                      float lowres,
                                                      float mass,
                                                      float sign_cc,
                                                      bool use_statistic, 
                                                      char* statistic_file, 
                                                      bool use_priors, 
                                                      float* priors_average, 
                                                      float* priors_variance
                                                    );
void DeleteComparisonObject(ImageProjectionComparison* comparison_object);

// =======================================
typedef struct Particle Particle;
Particle* InitRefineParticle( ReconstructedVolume* vol ); 
Particle* InitProjectionImage( ReconstructedVolume* vol );
Particle* InitUnbinnedImage( MRCFile* mrc, float padding );
void SetParticleStat( Particle** set, float** par_data, int numImage, int numParameter );
void SetBasicParams( Particle* particle, Image* image, int index, int numImage, int numParameter, float outer_mask_radius, float mask_falloff, float high_resolution_limit, float low_resolution_limit, float molecular_mass_kDa, float signed_CC_limit, float ref_3d_pixel_size, float binning_factor_refine, float voltage_kV, float spherical_aberration_mm, float amplitude_contrast, ResolutionStatistics* refine_statistics );

void DeleteParticle( Particle* particle );

// =======================================
float FrealignObjectiveFunction(void *scoring_parameters, float *array_of_values);
void RunRefine3d( float** params, 
                Image** image_set, 
                ReconstructedVolume* ref_3d, 
                ImageProjectionComparison** comparison_objects,
                bool optimize_internally, 
                bool evaluateAllProjections, 
                bool only_evaluate, 
                bool frame_refine, 
                bool use_priors, 
                float* priors_average, 
                float* priors_variance, 
                int number_particles, 
                char* symm, 
                float pixelsize,
                int boxsize, 
                float volt, 
                float cs, 
                float am, 
                float mass, 
                float mask, 
                float lowres, 
                float highres, 
                float sign_cc, 
                float class_res );

// =======================================
Image** ExtractParticlesFromMRCs( char input_file[1000], bool use_frames, bool stack_avg, double coordinates[][6], int number_particles, int number_tilts, int number_frames, int binning, int boxsize, float outer_mask_radius, float pixel_size, int write_to_disk, char output_file[1000]  );
Image** ReadParticlesFromStack(char stack_file[1000], int images_to_extract[], int number_of_images, int boxsize);

void ApplyGaussianToImages(Image** particle_images, int number_of_images, int boxsize, float pixel_size, double sigma);

float FrealignObjectiveConstrainedFunction(void *scoring_parameters, float *array_of_values);
bool RunNormalRefine3d();
//void UserInput();
#ifdef __cplusplus
}
#endif


