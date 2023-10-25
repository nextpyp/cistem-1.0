#include "../../core/core_headers.h"

#include "refine3d_cspt.h"
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <string>
#include <cstdio>
#include <ctime>

class ImageProjectionComparison {
    public:
    Particle					*particle;
    ReconstructedVolume			*reference_volume;
    Image						*projection_image;

};


// wrap the functions we need in csp
extern "C" {
    
    /**
     * MRC files 
     */
    MRCFile* ReadStack( char * filename, bool overwrite ) {
        
        std::string str(filename);
        if (!DoesFileExist(str)){
            printf( "[ERROR] %s is not found.\n", str.c_str() );
            exit(1);
        }

        MRCFile* stack =  new MRCFile( str, overwrite );
        
        if ( stack->ReturnXSize() != stack->ReturnYSize() ) {
            printf( "[ERROR] %s has to be square: X = %d, Y = %d\n", filename, stack->ReturnXSize(), stack->ReturnYSize()  );
            exit(1);
        }

        return stack;
    }
    MRCFile* ReadReconstruction( char * filename, bool overwrite, bool wait_for_file_to_exist ) {
        
        std::string str(filename);
        if (!DoesFileExist(str)) {
            printf( "[ERROR] %s is not found.\n", str.c_str() );
            exit(1);
        }

        MRCFile* reconstruction =  new MRCFile( str, overwrite, wait_for_file_to_exist );
        
        if ( ( reconstruction->ReturnXSize() != reconstruction->ReturnYSize() ) || ( reconstruction->ReturnXSize() != reconstruction->ReturnZSize() ) ) {
            printf( "[ERROR] %s has to be cubic: X = %d, Y = %d, Z = %d\n", filename, reconstruction->ReturnXSize(), reconstruction->ReturnYSize(), reconstruction->ReturnZSize() );
            exit(1);
        }
        return reconstruction;
    }
    int ReturnX( MRCFile* mrc ) {
        return mrc->ReturnXSize();
    }
    int ReturnY( MRCFile* mrc ) {
        return mrc->ReturnYSize();
    }
    int ReturnZ( MRCFile* mrc ) {
        return mrc->ReturnZSize();
    }

    void MRCFilePrintInfo( MRCFile* mrc )  {
        mrc->PrintInfo();
    }
    void CloseMRCFile( MRCFile* mrc ) {
        delete mrc;
    }

    /**
     * Volume
     */
    
    ReconstructedVolume* Init3DVolume( MRCFile* mrc, float pixel_size, char* wanted_symmetry_symbol, float molecular_mass_kDa, float outer_mask_radius, float low_resolution_limit, float high_resolution_limit ) {
        
        float mask_falloff = 40.0; // in A 

        ReconstructedVolume* vol = new ReconstructedVolume;
        wxString sym(wanted_symmetry_symbol, wxConvUTF8); 
         
        vol->InitWithDimensions( mrc->ReturnXSize(), mrc->ReturnYSize(), mrc->ReturnZSize(), pixel_size, sym );
        vol->molecular_mass_in_kDa = molecular_mass_kDa;
        
        
        if (outer_mask_radius > float(mrc->ReturnXSize()) / 2.0 * pixel_size- mask_falloff) outer_mask_radius = float(mrc->ReturnXSize()) / 2.0 * pixel_size - mask_falloff;

        vol->density_map.ReadSlices(mrc,1,vol->density_map.logical_z_dimension);
        vol->density_map.CosineMask(outer_mask_radius / pixel_size, mask_falloff / pixel_size, false, true, 0.0);  
        vol->mask_radius = outer_mask_radius;

        vol->PrepareForProjections(low_resolution_limit, high_resolution_limit);
         
        return vol;
    }

    float GetPixelSize( ReconstructedVolume* vol ) {
        return vol->pixel_size;
    }

    void Delete3DVolume( ReconstructedVolume* vol ) {
        delete vol;
    }
    
    /**
     * Particle
     */
    Particle* InitRefineParticle( ReconstructedVolume* vol ) {
        Particle* particle = new Particle();

        particle->mask_center_2d_x = 0.0;
        particle->mask_center_2d_y = 0.0;
        particle->mask_center_2d_z = 0.0;
        particle->mask_radius_2d = 0.0;

        particle->parameter_map[3] = false;
        particle->parameter_map[2] = false;
        particle->parameter_map[1] = false;
        particle->parameter_map[4] = false;
        particle->parameter_map[5] = false;

        particle->apply_2D_masking = false;

        particle->constraints_used[4] = false;
        particle->constraints_used[5] = false;

        particle->Allocate(vol->density_map.logical_x_dimension, vol->density_map.logical_y_dimension);
        return particle;
    }

    
    void SetParticleStat( Particle** set, float** par_data, int numImage, int numParameter ) {
        /**
         * It seems like we don't need this method, cuz we're not imposing any constraints within frealign
         */
        
        
        float input_parameters[numParameter]; 
        float parameter_average[numParameter];
        float parameter_variance[numParameter];
        int temp_float;

        ZeroFloatArray(input_parameters, numParameter);
        ZeroFloatArray(parameter_average, numParameter);
        ZeroFloatArray(parameter_variance, numParameter);
        
        int j = 0;
        for (int current_image = 0; current_image < numImage; current_image++) {
            for ( int i = 0; i < numParameter; i++ ) {
                input_parameters[i] = par_data[current_image][i];
            }
            // swap psi and phi
            temp_float = input_parameters[1]; input_parameters[1] = input_parameters[3]; input_parameters[3] = temp_float;
            if (input_parameters[7] >= 0) {
                for (int i = 0; i < numParameter; i++) {
                        parameter_average[i] += input_parameters[i];
                        parameter_variance[i] += powf(input_parameters[i],2);
                }
                j++;
            }
        }
        for (int i = 0; i < numParameter; i++) {
            parameter_average[i] /= j;
            parameter_variance[i] /= j;
            parameter_variance[i] -= powf(parameter_average[i],2);
            
            //if (parameter_variance[i] < 0.001) refine_particle.constraints_used[i] = false;
        }

        for ( int i = 0; i < numImage; i++ ) {
            set[i]->SetParameterStatistics(parameter_average, parameter_variance);
        }
    }


    void SetBasicParams( Particle* particle, Image* image, int index, int numImage, int numParameter, float outer_mask_radius, float mask_falloff, float high_resolution_limit, float low_resolution_limit, float molecular_mass_kDa, float signed_CC_limit, float ref_3d_pixel_size, float binning_factor_refine, float voltage_kV, float spherical_aberration_mm, float amplitude_contrast, ResolutionStatistics* refine_statistics ) {
        /**
         * Set basic parameters (such as voltage, pixel size) on ONE particle object
         */
        float input_parameters[numParameter];
        float parameter_average[numParameter];
        float parameter_variance[numParameter];
        float cg_accuracy[numParameter];
        float cg_starting_point[numParameter];

        int temp_float;

        ZeroFloatArray(input_parameters, numParameter);
        ZeroFloatArray(parameter_average, numParameter);
        ZeroFloatArray(parameter_variance, numParameter);
        ZeroFloatArray(cg_accuracy, numParameter);
        ZeroFloatArray(cg_starting_point, numParameter);

        // HIGHLIGHT - actually we don't need particle parameters since we just initailize particle objects
        /**
        for ( int i = 0; i < numParameter; i++ ) {
            input_parameters[i] = par_data[index][i];
        }*/

        //swap psi and phi
        temp_float = input_parameters[1]; input_parameters[1] = input_parameters[3]; input_parameters[3] = temp_float;



        particle->ResetImageFlags();
        particle->mask_radius = outer_mask_radius;
        particle->mask_falloff = mask_falloff;
        particle->filter_radius_high = high_resolution_limit;
        particle->molecular_mass_kDa = molecular_mass_kDa;
        particle->signed_CC_limit = signed_CC_limit;
        particle->pixel_size = ref_3d_pixel_size;
        particle->is_normalized = false;
        particle->sigma_noise = input_parameters[14] / binning_factor_refine;
        particle->SetParameters(input_parameters);
        particle->MapParameterAccuracy(cg_accuracy);
        particle->SetParameterConstraints(powf(parameter_average[14],2));
       
        image->ClipInto(particle->particle_image);
        
        //particle->MapParameters(cg_starting_point);
        particle->PhaseShiftInverse();
        /**
        particle->InitCTFImage(voltage_kV, spherical_aberration_mm, amplitude_contrast, input_parameters[8], input_parameters[9], input_parameters[10], input_parameters[11]);
        
        particle->filter_radius_low = low_resolution_limit;
        particle->SetIndexForWeightedCorrelation();
        particle->WeightBySSNR(refine_statistics->part_SSNR, 0);

        particle->PhaseFlipImage();
        particle->CosineMask();
        particle->PhaseShift();
        particle->CenterInCorner();
        */
    }

    void DeleteParticle( Particle* particle ) {
        delete particle;
    }

    /**
     * Image
     */ 

    Image* CreateImage( MRCFile* stack, float index_in_stack, int stack_x, int stack_y, float padding ) {
        Image* input_image = new Image();
        input_image->Allocate( stack_x, stack_y, true);

        Image* unbinned_image = new Image();
        unbinned_image->Allocate( stack_x * padding, stack_y * padding, true);

        input_image->ReadSlice( stack, int( index_in_stack + 0.5) );
        input_image->ReplaceOutliersWithMean(5.0);
        
        input_image->ClipInto(unbinned_image);
        unbinned_image->ForwardFFT();

        delete input_image;

        return unbinned_image;
    }

    float GetMean( Image* image ) {
        return image->ReturnAverageOfRealValues();
    }

    void NormalizeImages( Image** particle_images, int start_ind, int end_ind, int boxsize, float outer_mask_radius, float pixel_size) {
            
        float mask_falloff = 40.0; // in A
        float max_samples = 2000.0;
        float mask_radius_for_noise;
        float variance;
        float average;
        int current_ind;
            
        bool use_noise_power = false;

            /**
            // output 
        MRCFile output_stack;
        std::string output(output_file);
        if (write_to_disk) { 
            output_stack.OpenFile(output,true);
        }*/
        
        Image sum_power;
        Image input_image;
        Image temp_image;
        Curve noise_power_spectrum;
        Curve number_of_terms;

        if (use_noise_power) {
            // Calculate noise power spectrum from all projection images
            input_image.Allocate(boxsize, boxsize, true);
            sum_power.Allocate(boxsize, boxsize, false);
            sum_power.SetToConstant(0.0);
            mask_radius_for_noise = outer_mask_radius / pixel_size;
            
            if (2.0 * mask_radius_for_noise + mask_falloff / pixel_size > 0.95 * particle_images[0]->logical_x_dimension)
            {
                mask_radius_for_noise = 0.95 * particle_images[0]->logical_x_dimension / 2.0 - mask_falloff / 2.0 / pixel_size;
            }
            noise_power_spectrum.SetupXAxis(0.0, 0.5 * sqrtf(2.0), int((sum_power.logical_x_dimension / 2.0 + 1.0) * sqrtf(2.0) + 1.0));
            number_of_terms.SetupXAxis(0.0, 0.5 * sqrtf(2.0), int((sum_power.logical_x_dimension / 2.0 + 1.0) * sqrtf(2.0) + 1.0));

            for (current_ind = start_ind; current_ind <= end_ind; current_ind++)
            {
                
                //if ((global_random_number_generator.GetUniformRandom() < 1.0 - 2.0 * percentage)) continue;
                //input_image.ReadSlice(&input_stack, int(input_parameters[0] + 0.5));
                input_image.CopyFrom( particle_images[current_ind] );

                variance = input_image.ReturnVarianceOfRealValues(outer_mask_radius / pixel_size, 0.0, 0.0, 0.0, true);
                if (variance == 0.0) continue;
                input_image.MultiplyByConstant(1.0 / sqrtf(variance));
                input_image.CosineMask(outer_mask_radius / pixel_size, mask_falloff / pixel_size, true);
                input_image.ForwardFFT();
                temp_image.CopyFrom(&input_image);
                temp_image.ConjugateMultiplyPixelWise(input_image);
                sum_power.AddImage(&temp_image);
            }
            sum_power.Compute1DRotationalAverage(noise_power_spectrum, number_of_terms);
            noise_power_spectrum.SquareRoot();
            noise_power_spectrum.Reciprocal();
        }
        
        float min, max;

        for (current_ind = start_ind; current_ind <= end_ind; current_ind++) {   
            particle_images[current_ind]->ReplaceOutliersWithMean(5.0);
            particle_images[current_ind]->ForwardFFT();
            if (use_noise_power) {particle_images[current_ind]->ApplyCurveFilter(&noise_power_spectrum);}
            particle_images[current_ind]->BackwardFFT();
            
            particle_images[current_ind]->GetMinMax(min, max);
            
              
            if ((min == max) || ( particle_images[current_ind]->ContainsBlankEdges(outer_mask_radius / pixel_size) ) ){
                //printf("%d - blank\n", current_ind);
                //printf("this %d\n", current_ind);
                particle_images[current_ind]->SetToConstant(0.0);
                particle_images[current_ind]->AddGaussianNoise(1);
            }
            
            variance = particle_images[current_ind]->ReturnVarianceOfRealValues(particle_images[current_ind]->physical_address_of_box_center_x - mask_falloff / pixel_size, 0.0, 0.0, 0.0, true);
            average = particle_images[current_ind]->ReturnAverageOfRealValues(particle_images[current_ind]->physical_address_of_box_center_x - mask_falloff / pixel_size, true);
            particle_images[current_ind]->AddMultiplyConstant(- average, 1.0 / sqrtf(variance));            
            
            particle_images[current_ind]->ForwardFFT();
        }
    }


    void NormalizeFrames( Image** particle_images, double coordinates[][6], int start_ind, int end_ind, int boxsize, float outer_mask_radius, float pixel_size, int num_frames) {
        
        float mask_falloff = 40.0, min, max;
        double weight, sum, scale, mean, std, average_std = 0.0;
        int counter, pixel_counter;
        bool contain_blank = false;
        std::vector< std::vector<double> > frame_weights;
        std::vector<double> stds;
        std::vector<double> curr_frame;

        // First check if the number of images is mutiple of the number of frames
        if ( end_ind - start_ind + 1 != num_frames ) {
            printf( "Number of images has to be equal to number of frames.\n" );
            printf( " - Number of particles: %d \n", end_ind - start_ind + 1 );
            printf( " - Number of frames: %d \n", num_frames );
            exit(EXIT_FAILURE);
        }
        
        Image* copies = new Image[num_frames];
        Image average;
        average.Allocate(boxsize,boxsize,1,true);
        average.SetToConstant(0.0);
        counter = 0;
        scale = 1.0 / num_frames;

        // 2. apply frame shifts to particle frames and normalize frames using same strategy from PYP
        for ( int curr_ind = start_ind; curr_ind <= end_ind; curr_ind++ ) { 
            copies[counter].Allocate(boxsize,boxsize,1,true);
            copies[counter].CopyFrom(particle_images[curr_ind]);
            copies[counter].ForwardFFT();  
            //copies[counter].PhaseShift(-coordinates[curr_ind][4] / pixel_size, -coordinates[curr_ind][5] / pixel_size);
            copies[counter].BackwardFFT();  
            
            // compute a single average from particle frames
            for (long pixel_counter = 0; pixel_counter < average.real_memory_allocated; pixel_counter++) {      
                average.real_values[pixel_counter] += copies[counter].real_values[pixel_counter] * scale;
            }

            counter ++;
        }
        
        average.ReplaceOutliersWithMean(5.0);
        average.GetMinMax(min, max);
        contain_blank = ( ( min == max ) || average.ContainsBlankEdges(outer_mask_radius/pixel_size) );
        
        if (contain_blank){    
            for (int i = 0; i < num_frames; i++) {
                copies[i].SetToConstant(0.0);
                copies[i].AddGaussianNoise(1);
                particle_images[start_ind+i]->SetToConstant(0.0);
                particle_images[start_ind+i]->AddGaussianNoise(1);
            }
        }

        //mean = average.ReturnAverageOfRealValues(average.physical_address_of_box_center_x - mask_falloff / pixel_size, true);
        mean = average.ReturnAverageOfRealValues( outer_mask_radius / pixel_size, true);
        
        if (!contain_blank) {
            // substract mean from every particle frame 
            for (int i = 0; i < num_frames; i ++) { 

                for (long pixel_counter = 0; pixel_counter < copies[i].real_memory_allocated; pixel_counter++) {      
                    copies[i].real_values[pixel_counter] -= mean;
                }

                //std = sqrtf(copies[i].ReturnVarianceOfRealValues(copies[i].physical_address_of_box_center_x - mask_falloff / pixel_size, 0.0, 0.0, 0.0, true));
                std = sqrtf(copies[i].ReturnVarianceOfRealValues( outer_mask_radius / pixel_size, 0.0, 0.0, 0.0, true));
                if (std > 0.0) {stds.push_back(std);}
            }

            for (int i = 0; i < stds.size(); i++) {
                average_std += stds[i];
            }

            if (stds.size() > 0) {
                average_std /= stds.size();
                
                // divide particle frame by average std
                for (int i = 0; i < num_frames; i ++) {
                    for (long pixel_counter = 0; pixel_counter < copies[i].real_memory_allocated; pixel_counter++) {      
                        copies[i].real_values[pixel_counter] /= average_std;
                    }
                }
            }
        }
        
        // also normalize raw particle frames
        // save particle frames without frame alignment 

        for (int curr_ind = start_ind; curr_ind <= end_ind; curr_ind++) {
            
            //particle_images[curr_ind]->BackwardFFT();
            if (!contain_blank) {
                for (long pixel_counter = 0; pixel_counter < particle_images[curr_ind]->real_memory_allocated; pixel_counter++) {      
                    particle_images[curr_ind]->real_values[pixel_counter] = (particle_images[curr_ind]->real_values[pixel_counter] - mean) / average_std;
                } 
            }
            particle_images[curr_ind]->ForwardFFT();
        }
        delete [] copies;
    }


    Image** WeightedAverages(Image** particle_images, double ** shifts, int number_of_particles, int num_frames, int boxsize, float pixel_size, bool frame_refine, double * score_weights, double weight_width, bool write_after_average, char output_file[1000] ) {
        
        double weight, sum, scale, variance, average;
        double mask_falloff = 40.0; 
        int counter;
        std::vector< std::vector<double> > frame_weights;
        std::vector<double> curr_frame;

        if (!num_frames > 0) { 
            printf("Found zero frame in your parfile. Please check");
            exit(EXIT_FAILURE); 
        }

        MRCFile output_stack;
        if (write_after_average) { 
            std::string output(output_file);
            output_stack.OpenFile(output,true);
        }

        Image** weighted_average = new Image*[number_of_particles];
        
        Image* copies = new Image[num_frames];
        for (int counter = 0; counter < num_frames; counter++) {copies[counter].Allocate(boxsize,boxsize,1,true);}

        // weight scheme
        for ( int i = 0; i < num_frames; i++ ) {
            counter = 0;
            sum = 0.0;
            curr_frame.clear();
            for ( int j = 0; j < num_frames; j++ ) {
                weight = exp( -pow( counter - i, 2 ) / weight_width );
                sum += weight;
                counter ++;
                curr_frame.push_back(weight);
            }
            for ( int j = 0; j < num_frames; j++ ) {
                curr_frame[j] = curr_frame[j] / ( sum / num_frames ) / num_frames;
                //printf("%f ", curr_frame[j]);
            }
            frame_weights.push_back(curr_frame);
            //printf("\n");
        }
        
        // generate weighted average using particle frames 
        int write_counter = 0;
        for ( int particle = 0; particle < number_of_particles; particle += num_frames) {
            
            // copy image from raw stack and apply shifts
            for ( int frame = 0; frame < num_frames; frame ++ ) {
                copies[frame].SetToConstant(0.0);
                copies[frame].CopyFrom(particle_images[particle+frame]);
                // copies[frame].ForwardFFT();  
                copies[frame].PhaseShift(-shifts[particle+frame][1] / pixel_size, -shifts[particle+frame][2] / pixel_size);
                copies[frame].BackwardFFT();  
            }
            // create weighted average from shift-applied stack
            for ( int frame = 0; frame < num_frames; frame ++ ) {
                
                weighted_average[particle+frame] = new Image();
                weighted_average[particle+frame]->Allocate(boxsize,boxsize,1,true);
                weighted_average[particle+frame]->SetToConstant(0.0);
                
                for ( int target_frame = 0; target_frame < num_frames; target_frame ++ ) {
                    
                    // having particle frames, if we intend to refine on frames -> create running weighted averages
                    // otherwise create simple average with equal weight
                    if (frame_refine || write_after_average) { scale = frame_weights[frame].at(target_frame); }
                    else { scale = 1.0 / num_frames; }//score_weights[target_frame]; }

                    // add the weighted target frame to current frame
                    for (long pixel_counter = 0; pixel_counter < particle_images[particle+frame]->real_memory_allocated; pixel_counter++) {      
                        weighted_average[particle+frame]->real_values[pixel_counter] += copies[target_frame].real_values[pixel_counter] * scale;
                    }
                }
                
                if (write_after_average) {
                    weighted_average[particle+frame]->WriteSlice(&output_stack, write_counter+1);
                    write_counter++;
                }
                weighted_average[particle+frame]->ReplaceOutliersWithMean(5.0);
                weighted_average[particle+frame]->ForwardFFT();   
                 
                if (!frame_refine && !write_after_average) { 
                    // we can actually skip other frames since one equally weighted image is used
                    // TODO: deallocate other frames to save memory
                    // break;
                }
            }
        }
        delete [] copies;
        return weighted_average;
    }


    void UpdateWeightedAverages(Image** particle_images, Image** weighted_average, double ** shifts, int number_of_particles, int num_frames, int boxsize, float pixel_size, double weight_width) {
        
        double weight, sum, scale, variance, average;
        double mask_falloff = 40.0; 
        int counter;
        std::vector< std::vector<double> > frame_weights;
        std::vector<double> curr_frame;

        MRCFile test_output;
        bool test = false;
        if (test) { 
            std::string output("test.mrc");
            test_output.OpenFile(output,true);
        }

        
        Image* copies = new Image[num_frames];
        for (int counter = 0; counter < num_frames; counter++) {copies[counter].Allocate(boxsize,boxsize,1,true);}

        // weight scheme
        for ( int i = 0; i < num_frames; i++ ) {
            counter = 0;
            sum = 0.0;
            curr_frame.clear();
            for ( int j = 0; j < num_frames; j++ ) {
                weight = exp( -pow( counter - i, 2 ) / weight_width );
                sum += weight;
                counter ++;
                curr_frame.push_back(weight);
            }
            for ( int j = 0; j < num_frames; j++ ) {
                curr_frame[j] = curr_frame[j] / ( sum / num_frames ) / num_frames;
            }
            frame_weights.push_back(curr_frame);
        }
        // generate weighted average using particle frames
        int ind_in_stack = 0;
        int write_counter = 0;
        for (int particle = 0; particle < number_of_particles; particle += num_frames) {

            // we only update one tilt at a time, so it would be followed by a break
            if (particle != shifts[ind_in_stack][0]) {continue;}

            // copy image from raw stack and apply shifts
            for ( int frame = 0; frame < num_frames; frame ++ ) {
                copies[frame].SetToConstant(0.0);
                copies[frame].CopyFrom(particle_images[particle+frame]);
                // copies[frame].ForwardFFT();  
                copies[frame].PhaseShift(-shifts[ind_in_stack+frame][1] / pixel_size, -shifts[ind_in_stack+frame][2] / pixel_size);
                copies[frame].BackwardFFT();  
            }
            // create weighted average from shift-applied stack
            for ( int frame = 0; frame < num_frames; frame ++ ) {

                weighted_average[particle+frame]->SetToConstant(0.0);
                
                for ( int target_frame = 0; target_frame < num_frames; target_frame ++ ) {
                    
                    scale = frame_weights[frame].at(target_frame);

                    // add the weighted target frame to current frame
                    for (long pixel_counter = 0; pixel_counter < particle_images[particle+frame]->real_memory_allocated; pixel_counter++) {      
                        weighted_average[particle+frame]->real_values[pixel_counter] += copies[target_frame].real_values[pixel_counter] * scale;
                    }
                }
                weighted_average[particle+frame]->ReplaceOutliersWithMean(5.0);
                weighted_average[particle+frame]->ForwardFFT();   
            }
            break;
        }
        delete [] copies;
    }



    void DeleteImage( Image* image ) {
        delete image;
    }

    void DeleteImageArray( Image** images, int number_of_images ) {
        for (int i = 0; i < number_of_images; i++) {
            delete images[i];
        }
        delete [] images;
    }

    /**
     *  Comparison object
     */
    
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
                                                    ) {

        ImageProjectionComparison* comparison_object = new ImageProjectionComparison;
        Image* projection_image = new Image();
        Particle* refine_particle = new Particle();

        wxString input_reconstruction_statistics(statistic_file, wxConvUTF8);

        float temp_float;
        temp_float = input_parameters[1]; input_parameters[1] = input_parameters[3]; input_parameters[3] = temp_float;
    
        float mask_falloff = 40.0; 
        bool normalize_particles = false;
        bool normalize_input_3d = true;
        float binning_factor_refine = ref_3d->pixel_size / pixel_size;
        
        if (mask > float(boxsize) / 2.0 * pixel_size- mask_falloff) mask = float(boxsize) / 2.0 * pixel_size - mask_falloff;

        refine_particle->mask_center_2d_x = 0.0;
        refine_particle->mask_center_2d_y = 0.0;
        refine_particle->mask_center_2d_z = 0.0;
        refine_particle->mask_radius_2d = 0.0;

        refine_particle->parameter_map[3] = false;
        refine_particle->parameter_map[2] = false;
        refine_particle->parameter_map[1] = false;
        refine_particle->parameter_map[4] = false;
        refine_particle->parameter_map[5] = false;
        refine_particle->apply_2D_masking = false;
        
        // refine_particle->constraints_used[1] = use_priors; // phi
        // refine_particle->constraints_used[2] = use_priors; // theta 
        // refine_particle->constraints_used[3] = use_priors; // psi 
        refine_particle->constraints_used[4] = use_priors; // x
        refine_particle->constraints_used[5] = use_priors; // y 

        float cg_accuracy[refine_particle->number_of_parameters];
        float cg_starting_point[refine_particle->number_of_parameters];
        ZeroFloatArray(cg_accuracy, refine_particle->number_of_parameters);
        ZeroFloatArray(cg_starting_point, refine_particle->number_of_parameters);

        // read statistic file if it exists
        ResolutionStatistics input_statistics(pixel_size, boxsize);
        ResolutionStatistics refine_statistics;

        if (! DoesFileExist(input_reconstruction_statistics) || ! use_statistic ) { 
            input_statistics.GenerateDefaultStatistics(mass);
        }
        else {
            input_statistics.ReadStatisticsFromFile(input_reconstruction_statistics);
        }
        refine_statistics = input_statistics;
        
        refine_particle->Allocate(ref_3d->density_map.logical_x_dimension, ref_3d->density_map.logical_y_dimension);
        projection_image->Allocate(ref_3d->density_map.logical_x_dimension, ref_3d->density_map.logical_y_dimension, false);

        if (use_priors) {refine_particle->SetParameterStatistics(priors_average, priors_variance);}
        
        refine_particle->ResetImageFlags();
        refine_particle->mask_radius = mask; //outer_mask_radius;
        refine_particle->mask_falloff = mask_falloff;
        refine_particle->filter_radius_high = highres; //high_resolution_limit;
        refine_particle->molecular_mass_kDa = mass; //molecular_mass_kDa;
        refine_particle->signed_CC_limit = sign_cc; //signed_CC_limit;

        refine_particle->pixel_size = ref_3d->pixel_size;
        refine_particle->is_normalized = normalize_particles;
        refine_particle->sigma_noise = input_parameters[14] / binning_factor_refine;

        refine_particle->SetParameters(input_parameters);
        refine_particle->MapParameterAccuracy(cg_accuracy);
        if (use_priors) {refine_particle->SetParameterConstraints(powf(priors_average[14],2));}
        
        int ind_local_stack = (int) (input_parameters[17]);
        image_set[ind_local_stack]->ClipInto(refine_particle->particle_image);
        
        comparison_object->reference_volume = ref_3d; //&input_3d;
        comparison_object->projection_image = projection_image;
        comparison_object->particle = refine_particle;

        refine_particle->MapParameters(cg_starting_point);
        refine_particle->PhaseShiftInverse();

        refine_particle->InitCTFImage(volt, 
                                      cs, 
                                      am, 
                                      input_parameters[8], 
                                      input_parameters[9], 
                                      input_parameters[10], 
                                      input_parameters[11]);
        //refine_particle->InitCTFImage(voltage_kV, spherical_aberration_mm, amplitude_contrast, input_parameters[8], input_parameters[9], input_parameters[10], input_parameters[11]);

        refine_particle->filter_radius_low = lowres; //low_resolution_limit;
        refine_particle->SetIndexForWeightedCorrelation();
        
        if (normalize_input_3d) refine_particle->WeightBySSNR(refine_statistics.part_SSNR, 1);
        else refine_particle->WeightBySSNR(refine_statistics.part_SSNR, 0);
        
        refine_particle->PhaseFlipImage();
        refine_particle->CosineMask();
        refine_particle->PhaseShift();
        refine_particle->CenterInCorner();
        
        return comparison_object;
    }

    void DeleteComparisonObject(ImageProjectionComparison* comparison_object) {
        delete comparison_object->projection_image;
        delete comparison_object->particle;
        delete comparison_object;
    }



    /**
     * Statistic for refinement
     */
    ResolutionStatistics* CreateStat( int stack_x, float pixel_size, float molecular_mass_kDa ) {
        
        ResolutionStatistics* stat = new ResolutionStatistics( pixel_size, stack_x );
        // use default statistics
        stat->GenerateDefaultStatistics(molecular_mass_kDa);

        return stat;
    }

    void DeleteStat( ResolutionStatistics* stat ) {
        delete stat;
    }

    /**
    *  Refine3d function
    */

    

    float FrealignObjectiveFunction(void *scoring_parameters, float *array_of_values) {

        ImageProjectionComparison *comparison_object = reinterpret_cast < ImageProjectionComparison *> (scoring_parameters);
        for (int i = 0; i < comparison_object->particle->number_of_parameters; i++) {comparison_object->particle->temp_float[i] = comparison_object->particle->current_parameters[i];}
        comparison_object->particle->UnmapParameters(array_of_values);


        if (comparison_object->particle->no_ctf_weighting) comparison_object->reference_volume->CalculateProjection(*comparison_object->projection_image,
                *comparison_object->particle->ctf_image, comparison_object->particle->alignment_parameters, 0.0, 0.0,
                comparison_object->particle->pixel_size / comparison_object->particle->filter_radius_high, false, false, false, false, false);
        // Case for normal parameter refinement with weighting applied to particle images and 3D reference
        else if (comparison_object->particle->includes_reference_ssnr_weighting) {
            comparison_object->reference_volume->CalculateProjection(*comparison_object->projection_image,
            *comparison_object->particle->ctf_image, comparison_object->particle->alignment_parameters, 0.0, 0.0,
            comparison_object->particle->pixel_size / comparison_object->particle->filter_radius_high, false, true, true, false, false);
        }
        // Case for normal parameter refinement with weighting applied only to particle images
        else {
            comparison_object->reference_volume->CalculateProjection(*comparison_object->projection_image,
            *comparison_object->particle->ctf_image, comparison_object->particle->alignment_parameters, 0.0, 0.0,
            comparison_object->particle->pixel_size / comparison_object->particle->filter_radius_high, false, true, false, true, true);
        }

    
        return 	- comparison_object->particle->particle_image->GetWeightedCorrelationWithImage(*comparison_object->projection_image, comparison_object->particle->bin_index,
                comparison_object->particle->pixel_size / comparison_object->particle->signed_CC_limit);
            //	- comparison_object->particle->ReturnParameterPenalty(comparison_object->particle->temp_float);
            // This penalty term assumes a Gaussian x,y distribution that is probably not correct in most cases. It might be better to leave it out.

    }


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
                    float class_res  ) {

        
                        
                       
                        
        Particle refine_particle;
        
        wxString wxsymm(symm, wxConvUTF8);
        
        wxString input_reconstruction_statistics 	= "/scratch/statistics_r01.txt";
        bool	 use_statistics						= false;
        int		 first_particle						= 1;
        int		 last_particle						= 0;
        float	 percent_used						= 1.0;
        float 	 pixel_size							= pixelsize;
        float    voltage_kV							= volt;
        float 	 spherical_aberration_mm			= cs;
        float    amplitude_contrast					= am;
        float	 molecular_mass_kDa					= mass;
        float    inner_mask_radius					= 0.0;
        float    outer_mask_radius					= mask;
        float    low_resolution_limit				= lowres;
        float    high_resolution_limit				= highres;
        float	 signed_CC_limit					= sign_cc;
        float	 classification_resolution_limit	= class_res;
        float    mask_radius_search					= mask;
        float	 high_resolution_limit_search		= highres;
        float	 max_search_x						= mask;
        float	 max_search_y						= mask;
        refine_particle.mask_center_2d_x			= 0.0;
        refine_particle.mask_center_2d_y			= 0.0;
        refine_particle.mask_center_2d_z			= 0.0;
        refine_particle.mask_radius_2d				= 0.0;
        bool	 global_search						= false;
        bool	 local_refinement					= true;
        
        refine_particle.parameter_map[3]			= false;
        refine_particle.parameter_map[2]			= false;
        refine_particle.parameter_map[1]			= false;
        refine_particle.parameter_map[4]			= false;
        refine_particle.parameter_map[5]			= false;
        
        if (optimize_internally && !only_evaluate) {
            refine_particle.parameter_map[4]		= true;
            refine_particle.parameter_map[5]		= true;
            if (!frame_refine) {
                refine_particle.parameter_map[1] = true;
                refine_particle.parameter_map[2] = true;
                refine_particle.parameter_map[3] = true;
            }
        }
        refine_particle.apply_2D_masking			= false;
        bool	 ctf_refinement						= false;
        bool	 normalize_particles				= false;
        bool	 normalize_input_3d					= true;
        bool	 local_global_refine				= false;
        
        
        refine_particle.constraints_used[4] = false;		// Constraint for X shifts
        refine_particle.constraints_used[5] = false;		// Constraint for Y shifts
        if (optimize_internally) {
            refine_particle.constraints_used[4] = true;
            refine_particle.constraints_used[5] = true;
        }
        
        Image projection_image;

        ReconstructedVolume			input_3d;
        ImageProjectionComparison	comparison_object;
        ConjugateGradient			conjugate_gradient_minimizer;
        
        
        int i;
        int j;
        int current_image;
        int image_counter = 0;
        float input_parameters[refine_particle.number_of_parameters];
        float output_parameters[refine_particle.number_of_parameters];
        float parameter_average[refine_particle.number_of_parameters];
        float parameter_variance[refine_particle.number_of_parameters];
        float cg_starting_point[refine_particle.number_of_parameters];
        float cg_accuracy[refine_particle.number_of_parameters];
        float binning_factor_refine;
        float mask_falloff = 40.0;	// in Angstrom
        float temp_float;

        
        ZeroFloatArray(input_parameters, refine_particle.number_of_parameters);
        ZeroFloatArray(output_parameters, refine_particle.number_of_parameters);
        ZeroFloatArray(parameter_average, refine_particle.number_of_parameters);
        ZeroFloatArray(parameter_variance, refine_particle.number_of_parameters);
        ZeroFloatArray(cg_starting_point, refine_particle.number_of_parameters);
        ZeroFloatArray(cg_accuracy, refine_particle.number_of_parameters);
        
        
        /**
         * Direct evaulation by using comparison objects initialized outside this function (faster)
         * now only supports TOMO micrograph & particle refinement (mode 5 & 6)
         */
        if (!optimize_internally) {
            for (current_image = 0; current_image < number_particles; current_image++) {
                // skip if OCC is zero
                if ( params[current_image][12] <= 0.0 ) { 
                    params[current_image][14] = 0.0;
                    params[current_image][15] = 0.0;
                    continue; 
                }
                // skip if SCANORD is not in the range for refinement
                if ( !evaluateAllProjections && params[current_image][17] < 0 ) {
                    // not in refinement range
                    continue;
                }
                for ( int i = 0; i < refine_particle.number_of_parameters; i++ ) {
                    input_parameters[i] = params[current_image][i];
                }
                temp_float = input_parameters[1]; input_parameters[1] = input_parameters[3]; input_parameters[3] = temp_float;
                for (int i = 0; i < refine_particle.number_of_parameters; i++) {output_parameters[i] = input_parameters[i];}
            
                // evaluate directly without having to repetitively initiate Image objects
                comparison_objects[current_image]->particle->SetParameters(input_parameters);
                output_parameters[15] = (use_priors) ? 
                                        - 100.0 * FrealignObjectiveConstrainedFunction(comparison_objects[current_image], cg_starting_point) : 
                                        - 100.0 * FrealignObjectiveFunction(comparison_objects[current_image], cg_starting_point);
                /**
                if (output_parameters[15] < 0.0) {
                    output_parameters[15] = 0.0;
                }
                */

                params[current_image][15] = output_parameters[15];
              
            }    
            return;
        }
        

        if (high_resolution_limit < 2.0 * pixel_size) high_resolution_limit = 2.0 * pixel_size;
        if (classification_resolution_limit < 2.0 * pixel_size) classification_resolution_limit = 2.0 * pixel_size;
        if (high_resolution_limit_search < 2.0 * pixel_size) high_resolution_limit_search = 2.0 * pixel_size;
        if (signed_CC_limit == 0.0) signed_CC_limit = pixel_size;

        if (outer_mask_radius > float(boxsize) / 2.0 * pixel_size- mask_falloff) outer_mask_radius = float(boxsize) / 2.0 * pixel_size - mask_falloff;
        if (mask_radius_search > float(boxsize) / 2.0 * pixel_size- mask_falloff) mask_radius_search = float(boxsize) / 2.0 * pixel_size - mask_falloff;
        
        input_3d = *ref_3d;

        ResolutionStatistics input_statistics(pixel_size, boxsize);
        ResolutionStatistics refine_statistics;
        if (use_statistics)
        {
            if (! DoesFileExist(input_reconstruction_statistics))
            {
                exit(-1);
            }
            input_statistics.ReadStatisticsFromFile(input_reconstruction_statistics);
        }
        else
        {
            input_statistics.GenerateDefaultStatistics(molecular_mass_kDa);
        }
        refine_statistics = input_statistics;
        if (inner_mask_radius > 0.0) input_3d.density_map.CosineMask(inner_mask_radius / pixel_size, mask_falloff / pixel_size, true);

        binning_factor_refine = input_3d.pixel_size / pixel_size;

        refine_particle.Allocate(input_3d.density_map.logical_x_dimension, input_3d.density_map.logical_y_dimension);
        projection_image.Allocate(input_3d.density_map.logical_x_dimension, input_3d.density_map.logical_y_dimension, false);

        
        /**
         * Obtain statistic (average and std) of parameters to add shift priors as constraints
         */ 
        if (optimize_internally) {
            /**
            j = 0;
            for (current_image = 0; current_image < number_particles; current_image++)
            {   
                for ( int i = 0; i < refine_particle.number_of_parameters; i++ ) {
                    input_parameters[i] = params[current_image][i];
                }
                temp_float = input_parameters[1]; input_parameters[1] = input_parameters[3]; input_parameters[3] = temp_float;
                if (input_parameters[7] >= 0)
                {
                    for (i = 1; i < refine_particle.number_of_parameters; i++)
                    {
                        parameter_average[i] += input_parameters[i];
                        parameter_variance[i] += powf(input_parameters[i],2);
                    }
                    j++;
                }
            }
            
            for (i = 0; i < refine_particle.number_of_parameters; i++)
            {
                parameter_average[i] /= j;
                parameter_variance[i] /= j;
                parameter_variance[i] -= powf(parameter_average[i],2);
                if (parameter_variance[i] < 0.001) refine_particle.constraints_used[i] = false;
            }
            refine_particle.SetParameterStatistics(parameter_average, parameter_variance);
            */
            refine_particle.SetParameterStatistics(priors_average, priors_variance);
        }
        
        /**
         * Refine each particle projection (slower)
         */
        for (current_image = 0; current_image < number_particles; current_image++)
        {  
            // skip if OCC is zero
            if ( params[current_image][12] <= 0.0 ) { 
                params[current_image][14] = 0.0;
                params[current_image][15] = 0.0;
                continue; 
            }
            
            // skip if SCANORD is not in the range for refinement
            if ( !evaluateAllProjections && params[current_image][17] < 0 ) {
                // not in refinement range
                continue;
            }

            for ( int i = 0; i < refine_particle.number_of_parameters; i++ ) {
                input_parameters[i] = params[current_image][i];
            }
            temp_float = input_parameters[1]; input_parameters[1] = input_parameters[3]; input_parameters[3] = temp_float;
            for (int i = 0; i < refine_particle.number_of_parameters; i++) {output_parameters[i] = input_parameters[i];}
            
            // 4.0
            refine_particle.ResetImageFlags();
            refine_particle.mask_radius = outer_mask_radius;
            refine_particle.mask_falloff = mask_falloff;
            refine_particle.filter_radius_high = high_resolution_limit;
            refine_particle.molecular_mass_kDa = molecular_mass_kDa;
            refine_particle.signed_CC_limit = signed_CC_limit;
            // The following line would allow using particles with different pixel sizes
            refine_particle.pixel_size = input_3d.pixel_size;
            refine_particle.is_normalized = normalize_particles;
            refine_particle.sigma_noise = input_parameters[14] / binning_factor_refine;
            
            refine_particle.SetParameters(input_parameters);
            refine_particle.MapParameterAccuracy(cg_accuracy);
            refine_particle.SetParameterConstraints(powf(parameter_average[14],2));
            
            int ind_local_stack = (int) (params[current_image][17]);
            image_set[ind_local_stack]->ClipInto(refine_particle.particle_image);
            
            // 1
            comparison_object.reference_volume = &input_3d;
            comparison_object.projection_image = &projection_image;
            comparison_object.particle = &refine_particle;
            
            // 171.0
            refine_particle.MapParameters(cg_starting_point);
            refine_particle.PhaseShiftInverse();

            // 516
            refine_particle.InitCTFImage(voltage_kV, spherical_aberration_mm, amplitude_contrast, input_parameters[8], input_parameters[9], input_parameters[10], input_parameters[11]);

            //1014
            refine_particle.filter_radius_low = low_resolution_limit;
            refine_particle.SetIndexForWeightedCorrelation();
        
            if (normalize_input_3d) refine_particle.WeightBySSNR(refine_statistics.part_SSNR, 1);
            // Apply SSNR weighting only to image since input 3D map assumed to be calculated from correctly whitened images
            else refine_particle.WeightBySSNR(refine_statistics.part_SSNR, 0);
            //particle->WeightBySSNR(refine_statistics.part_SSNR, 0);

            refine_particle.PhaseFlipImage();
            refine_particle.CosineMask();
            refine_particle.PhaseShift();
            refine_particle.CenterInCorner();

            // 790.0
            if (optimize_internally && !only_evaluate) {
                // align particle against reference with shift priors as constraints (frame refinement)
                temp_float = - 100.0 * conjugate_gradient_minimizer.Init(&FrealignObjectiveConstrainedFunction, &comparison_object, refine_particle.number_of_search_dimensions, cg_starting_point, cg_accuracy);
                output_parameters[15] = - 100.0 * conjugate_gradient_minimizer.Run();
                refine_particle.UnmapParametersToExternal(output_parameters, conjugate_gradient_minimizer.GetPointerToBestValues());
                
                // output_parameters[15] = - 100.0 * FrealignObjectiveConstrainedFunction(&comparison_object, cg_starting_point);
                
                params[current_image][1] = output_parameters[1];
                params[current_image][2] = output_parameters[2];
                params[current_image][3] = output_parameters[3];
                params[current_image][4] = output_parameters[4];
                params[current_image][5] = output_parameters[5];
                
            }
            else {
                // simply evaluate alignment (tomo)
                output_parameters[15] = - 100.0 * FrealignObjectiveFunction(&comparison_object, cg_starting_point);
            }

            if (output_parameters[15] < 0.0) {
                // output_parameters[15] = 0.0;
            }

            // update score
            params[current_image][15] = output_parameters[15];
            if (refine_particle.snr > 0.0) output_parameters[14] = sqrtf(1.0 / refine_particle.snr);
        }
    }

    /**
    Extraction 
    */

    Image** ExtractParticlesFromMRCs( char input_file[1000], bool use_frames, bool stack_avg, double coordinates[][6], int number_particles, int number_tilts, int number_frames, int binning, int boxsize, float outer_mask_radius, float pixel_size, int write_to_disk, char output_file[1000] ) {
        
        char line[1000], frame_path[1000];
        int coord_x, coord_y, idx_tilt_image, idx_frame, image_x, image_y, curr, start_ind, end_ind;
        FILE *file;

        // std::string* movie_list;
        std::vector<std::string> movie_list; 
  
        // Read filenames of all the tilted movies
        if (use_frames) {

            //movie_list = new std::string[number_tilts];
            // opne the txt 
            file = fopen(input_file, "r");
            if (!file) {
                printf("[ERROR] File %s does not exist or cannot be read.\nExiting...\n", input_file);
                exit(EXIT_FAILURE);
            }
            curr = 0;
            while (fgets(line, sizeof(line), file) != NULL) {    
                sscanf(line, "%s", frame_path);
                std::string frame_path_str(frame_path);
                if (!DoesFileExist(frame_path_str)) {
                    printf( "[ERROR] %s not found.\n", frame_path_str.c_str() );
                    exit(EXIT_FAILURE);
                }
                //movie_list[curr] = frame_path_str;
                movie_list.push_back(frame_path_str);
                curr++;        
            }      
            fclose(file);
            if (curr < number_tilts) {
                printf("[ERROR] %s has only %d lines, which does not match %d tilts in allboxes.\n", input_file, curr, number_tilts);
                exit(EXIT_FAILURE);
            }
        }
        Image** particle_images = new Image*[number_particles];
        
        // output stack 
        MRCFile output_stack;
        std::string output(output_file);
  
        
        if (use_frames) {
            /**
             * Reading all tilted movies into memory is really expensive, so we choose to extract particles from one tilted movie at a time
             * All particle frames at a given tilt are all extracted to create weighted average, but uninterested frames will be deallocated to save memory footprint
             * Image normalization and running average creation are done on a per-particle-per-tilt basis (like sliding window) to avoid saving two stacks of eqaully large size 
             */
            
            double time_writing = 0.0;

            Image* frame_images = new Image[number_frames]; 
            float image_means[number_frames];
            
            double timer_read_all_tilts = 0.0;

            for ( int tilt = 0; tilt < number_tilts; tilt++ ) {
                
                // 1. read one tilted movies at a time
                MRCFile* tilt_movie = new MRCFile(movie_list.at(tilt),false); // new MRCFile(movie_list[tilt],false);
                if (tilt_movie->ReturnZSize() < number_frames){
                    // printf("[ERROR] %s has only %d frames that does not match %d frames in allboxes\n", movie_list[tilt].c_str(), tilt_movie->ReturnZSize(), number_frames);
                    // printf("[ERROR] %s has only %d frames that does not match %d frames in allboxes\n", movie_list.at(tilt), tilt_movie->ReturnZSize(), number_frames);
                    exit(EXIT_FAILURE);
                }
                image_x = tilt_movie->ReturnXSize();
                image_y = tilt_movie->ReturnYSize();
                
                // 2. convert it to Image object
                for ( int frame = 0; frame < number_frames; frame++ ) {
                    frame_images[frame].SetToConstant(0.0);

                    // FIXME
                    std::clock_t start = std::clock();
                    frame_images[frame].ReadSlice(tilt_movie, frame+1);
                    double duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
                    timer_read_all_tilts += duration;

                    image_means[frame] = frame_images[frame].ReturnAverageOfRealValues();
                }

                delete tilt_movie;

                // 3. only extract particles at current tilt (tilts are not sorted so we need to loop over every line)
                //    but frames at every tilt has to be continuous (i.e. tilt0_frame0, tilt0_frame1, tilt0_frame2 ...)
                for ( int counter = 0; counter < number_particles; counter++ ) {
                    idx_tilt_image = coordinates[counter][2];
                    // printf("(%d, %d), tilt = %f, frame = %d\n", coordinates[counter][0], coordinates[counter][1], idx_tilt_image, (int) coordinates[counter][3]);
                    if (idx_tilt_image != tilt) {continue;}
                    else {
                        coord_x = ( coordinates[counter][0] / binning ) - (image_x/2);
                        coord_y = ( coordinates[counter][1] / binning ) - (image_y/2);
                        idx_frame = (int) coordinates[counter][3]; 
                        
                        //printf("(%f, %f) => (%d, %d)\n", coordinates[counter][0], coordinates[counter][1], coord_x, coord_y);
                        particle_images[counter] = new Image();
                        particle_images[counter]->Allocate(boxsize,boxsize,1,true);
                        frame_images[idx_frame].ClipInto(particle_images[counter],image_means[idx_frame],false,1.0,coord_x,coord_y,0); 
                        
                        // normalize one particle at a time and create running average 
                        if (idx_frame == 0) { start_ind = counter; }
                        if (idx_frame == number_frames-1) { 
                            end_ind = counter; 
                            //printf("%d - %d (%d)\n", start_ind, end_ind, idx_frame);
                            NormalizeFrames( particle_images, coordinates, start_ind, end_ind, boxsize, outer_mask_radius, pixel_size, number_frames);
                        }
                        // TODO: remove other particle frames to save memory footprint
                    }
                }
            }
            //printf("Time spent for reading all tilts: %f s\n", timer_read_all_tilts);
            //printf("Time spent for writing particle stack: %f s\n", time_writing);
            // 4. deallocate the tilt image
            delete[] frame_images;
            //delete[] movie_list; 
        }
        else { 
            // 1. read the entire tiltseries
            std::string input(input_file); 
            if (!DoesFileExist(input)) {
                printf( "[ERROR] %s not found.\n", input.c_str() );
                exit(EXIT_FAILURE);
            }
            MRCFile tiltseries(input,false);
            if (tiltseries.ReturnZSize() < number_tilts){
                printf("[ERROR] %s has only %d tilts that does not match the number of tilts %d in allboxes\n", input_file, tiltseries.ReturnZSize(), number_tilts);
                exit(EXIT_FAILURE);
            }
            image_x = tiltseries.ReturnXSize();
            image_y = tiltseries.ReturnYSize();

            // 2. convert all the tilted images into Image objects
            Image* tilted_images = new Image[number_tilts]; 
            float image_means[number_tilts];
            for ( int tilt = 0; tilt < number_tilts; tilt++ ) {
                tilted_images[tilt].ReadSlice(&tiltseries, tilt+1);
                image_means[tilt] = tilted_images[tilt].ReturnAverageOfRealValues();
            }

            // 3. extract paticles one by one            
            for ( int counter = 0; counter < number_particles; counter ++ )
            {   
                coord_x = ( coordinates[counter][0] / binning ) - (image_x/2);
                coord_y = ( coordinates[counter][1] / binning ) - (image_y/2);
                idx_tilt_image = coordinates[counter][2];
                idx_frame = coordinates[counter][3]; // it is 0 if not using frames

                particle_images[counter] = new Image();
                particle_images[counter]->Allocate(boxsize,boxsize,1,true);

                tilted_images[idx_tilt_image].ClipInto(particle_images[counter],image_means[idx_tilt_image],false,1.0,coord_x,coord_y,0);
            }  

            // 4. deallocate Images
            delete[] tilted_images;

            // 5. normalize every particle image 
            NormalizeImages(particle_images, 0, number_particles-1, boxsize, outer_mask_radius, pixel_size);
        }
        if (write_to_disk) {
            output_stack.OpenFile(output,true);
            for (int i = 0; i < number_particles; i++) {
                particle_images[i]->WriteSlice(&output_stack, i+1);
            }
            output_stack.CloseFile();
        }

        return particle_images;
    }


    Image** ReadParticlesFromStack(char stack_file[1000], int images_to_extract[], int number_of_images, int boxsize) {
        
        // FIXME
        double time_read_from_stacks = 0.0;

        std::clock_t start = std::clock(); 
        std::string input(stack_file);
        MRCFile input_stack(input, false);
        double duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
        time_read_from_stacks += duration; 


        MRCFile output_stack;
        std::string output("test.mrc");
        bool test_output = false;
        if (test_output) { 
            output_stack.OpenFile(output,true);
        }

        if (boxsize != input_stack.ReturnXSize()) {
            printf("Boxsize of input stack (%d) does not match your desired boxsize (%d)\n", input_stack.ReturnXSize(), boxsize);
            exit(EXIT_FAILURE);
        }

        Image** particle_images = new Image*[number_of_images];
        
        for (int counter = 0; counter < number_of_images; counter++) {
            particle_images[counter] = new Image();
            particle_images[counter]->Allocate(boxsize,boxsize,1,true);
            
            std::clock_t start = std::clock();    
            particle_images[counter]->ReadSlice(&input_stack, int(images_to_extract[counter]));
            double duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
            time_read_from_stacks += duration;

            if (test_output) {particle_images[counter]->WriteSlice(&output_stack, counter+1); }
            particle_images[counter]->ForwardFFT();
        }
        printf("Time spent for reading particles from stack: %f s\n", time_read_from_stacks);

        return particle_images;
    }


    void ApplyGaussianToImages(Image** particle_images, int number_of_images, int boxsize, float pixel_size, double sigma) {

        float mask_falloff = 40.0; // A

        MRCFile test_output;
        bool test = false;
        if (test) { 
            std::string output("test.mrc");
            test_output.OpenFile(output,true);
        }

        double r, s = 2.0 * sigma * sigma;
        double sum = 0.0;
        int half = boxsize / 2;

        // obtain gaussian weights 
        double gaussian_kernel[boxsize+1][boxsize+1]; 
        int start = -boxsize/2, end = boxsize/2; 
        for (int y = start; y <= end; y++) {
            for (int x = start; x <= end; x++) {        
                r = sqrt(x * x + y * y);
                gaussian_kernel[y+half][x+half] = (exp(-(r * r) / s)) / (M_PI * s);
                sum += gaussian_kernel[y+half][x+half];
            }
        }
        for (int y = start; y <= end; y++) {
            for (int x = start; x <= end; x++) {
                gaussian_kernel[y+half][x+half] /= sum;
            }
        }


        // apply shift with weights to blur image
        Image average, translated_image;
        int x_shift_in_pixel = 0, y_shift_in_pixel = 0;
        double scale = 0.0; 
        float variance, mean;

        average.Allocate(boxsize,boxsize,1,true);
        translated_image.Allocate(boxsize, boxsize,1,true);

        for (int i = 0; i < number_of_images; i++) {
            // make raw particles back to real space
            particle_images[i]->BackwardFFT();
            average.SetToConstant(0.0);
            //particle_images[i]->WriteSlice(&test_output, 1);

            for (int y = start; y <= end; y++) {
                for (int x = start; x <= end; x++) {
                    if (gaussian_kernel[y+half][x+half] > 0.000001) {
                        //printf("%d %d => %f\n", x, y, gaussian_kernel[y+half][x+half]);
                        translated_image.SetToConstant(0.0);
                        translated_image.CopyFrom(particle_images[i]);

                        scale = gaussian_kernel[y+half][x+half];
                        x_shift_in_pixel = x;
                        y_shift_in_pixel = y; 

                        translated_image.ForwardFFT();    
                        translated_image.PhaseShift(x_shift_in_pixel, y_shift_in_pixel);
                        translated_image.BackwardFFT();
                        
                        for (long pixel_counter = 0; pixel_counter < average.real_memory_allocated; pixel_counter++) {      
                            average.real_values[pixel_counter] += translated_image.real_values[pixel_counter] * scale;
                        }
                    }
                }
            }


            // replace particle images with averages
            particle_images[i]->SetToConstant(0.0);
            for (long pixel_counter = 0; pixel_counter < particle_images[i]->real_memory_allocated; pixel_counter++) {      
                particle_images[i]->real_values[pixel_counter] = average.real_values[pixel_counter];
            }

            variance = particle_images[i]->ReturnVarianceOfRealValues(particle_images[i]->physical_address_of_box_center_x - mask_falloff / pixel_size, 0.0, 0.0, 0.0, true);
            mean = particle_images[i]->ReturnAverageOfRealValues(particle_images[i]->physical_address_of_box_center_x - mask_falloff / pixel_size, true);
            particle_images[i]->AddMultiplyConstant(-mean, 1.0 / sqrtf(variance));  

            // particle_images[i]->WriteSlice(&test_output, 2);
            // test_output.CloseFile();
            // exit(1);
                   
            particle_images[i]->ForwardFFT();

        }
    }



    float FrealignObjectiveConstrainedFunction(void *scoring_parameters, float *array_of_values)
    {
        ImageProjectionComparison *comparison_object = reinterpret_cast < ImageProjectionComparison *> (scoring_parameters);
        for (int i = 0; i < comparison_object->particle->number_of_parameters; i++) {comparison_object->particle->temp_float[i] = comparison_object->particle->current_parameters[i];}
        comparison_object->particle->UnmapParameters(array_of_values);

    //	comparison_object->reference_volume->CalculateProjection(*comparison_object->projection_image, *comparison_object->particle->ctf_image,
    //			comparison_object->particle->alignment_parameters, 0.0, 0.0,
    //			comparison_object->particle->pixel_size / comparison_object->particle->filter_radius_high, false, true);

        if (comparison_object->particle->no_ctf_weighting) comparison_object->reference_volume->CalculateProjection(*comparison_object->projection_image,
                *comparison_object->particle->ctf_image, comparison_object->particle->alignment_parameters, 0.0, 0.0,
                comparison_object->particle->pixel_size / comparison_object->particle->filter_radius_high, false, false, false, false, false);
        // Case for normal parameter refinement with weighting applied to particle images and 3D reference
        else if (comparison_object->particle->includes_reference_ssnr_weighting) comparison_object->reference_volume->CalculateProjection(*comparison_object->projection_image,
                *comparison_object->particle->ctf_image, comparison_object->particle->alignment_parameters, 0.0, 0.0,
                comparison_object->particle->pixel_size / comparison_object->particle->filter_radius_high, false, true, true, false, false);
        // Case for normal parameter refinement with weighting applied only to particle images
        else comparison_object->reference_volume->CalculateProjection(*comparison_object->projection_image,
                *comparison_object->particle->ctf_image, comparison_object->particle->alignment_parameters, 0.0, 0.0,
                comparison_object->particle->pixel_size / comparison_object->particle->filter_radius_high, false, true, false, true, true);

    //	if (comparison_object->particle->origin_micrograph < 0) comparison_object->particle->origin_micrograph = 0;
    //	comparison_object->particle->origin_micrograph++;
    //	for (int i = 0; i < comparison_object->projection_image->real_memory_allocated; i++) {comparison_object->projection_image->real_values[i] *= fabs(comparison_object->projection_image->real_values[i]);}
    //	comparison_object->projection_image->ForwardFFT();
    //	comparison_object->projection_image->CalculateCrossCorrelationImageWith(comparison_object->particle->particle_image);
    //	comparison_object->projection_image->SwapRealSpaceQuadrants();
    //	comparison_object->projection_image->BackwardFFT();
    //	comparison_object->projection_image->QuickAndDirtyWriteSlice("proj.mrc", comparison_object->particle->origin_micrograph);
    //	comparison_object->projection_image->SwapRealSpaceQuadrants();
    //	comparison_object->particle->particle_image->SwapRealSpaceQuadrants();
    //	comparison_object->particle->particle_image->BackwardFFT();
    //	comparison_object->particle->particle_image->QuickAndDirtyWriteSlice("part.mrc", comparison_object->particle->origin_micrograph);
    //	comparison_object->particle->particle_image->SwapRealSpaceQuadrants();
    //	exit(0);

    //	float score =  	- comparison_object->particle->particle_image->GetWeightedCorrelationWithImage(*comparison_object->projection_image, comparison_object->particle->bin_index,
    //			  comparison_object->particle->pixel_size / comparison_object->particle->signed_CC_limit)
    //			- comparison_object->particle->ReturnParameterPenalty(comparison_object->particle->temp_float);
    //	wxPrintf("psi, theta, phi, x, y, = %g, %g, %g, %g, %g, score = %g\n",
    //			comparison_object->particle->alignment_parameters.ReturnPsiAngle(),
    //			comparison_object->particle->alignment_parameters.ReturnThetaAngle(),
    //			comparison_object->particle->alignment_parameters.ReturnPhiAngle(),
    //			comparison_object->particle->alignment_parameters.ReturnShiftX(),
    //			comparison_object->particle->alignment_parameters.ReturnShiftY(), score);
    //	return score;
    //	wxPrintf("sigma_noise, mask_volume, penalty = %g %g %g\n", comparison_object->particle->sigma_noise, comparison_object->particle->mask_volume,
    //			comparison_object->particle->ReturnParameterPenalty(comparison_object->particle->temp_float));
        return 	- comparison_object->particle->particle_image->GetWeightedCorrelationWithImage(*comparison_object->projection_image, comparison_object->particle->bin_index,
                comparison_object->particle->pixel_size / comparison_object->particle->signed_CC_limit)
                - comparison_object->particle->ReturnParameterPenalty(comparison_object->particle->temp_float);
    }



    
    bool RunNormalRefine3d(){
        Particle refine_particle;
        Particle search_particle;

        std::string input;
        bool	 use_statistics	;
        
        // wxString input_particle_images 				= "/nfs/bartesaghilab/hl325/test_stack.mrc";
        // wxString input_parameter_file 				= "/nfs/bartesaghilab/hl325/test.par";
        // wxString input_reconstruction				= "/nfs/bartesaghilab/hl325/20220305_111443_T20S_frames_00_38_r01_02.mrc";
        // wxString input_reconstruction_statistics 	= "statistics_r01.txt";
        // bool	 use_statistics						= false;
        // wxString ouput_matching_projections 		= "14sep05c_c_00003gr_00019sq_00010hl_00004es.frames_P0004_frames_r01_02_match.mrc_";
        // wxString ouput_parameter_file				= "14sep05c_c_00003gr_00019sq_00010hl_00004es.frames_P0004_frames_r01_02.par_";
        // wxString ouput_shift_file					= "/dev/null";
        // wxString my_symmetry						= "D7";
        
        
        std::cin >> input; 
        wxString input_particle_images(input.c_str());
        std::cin >> input; 
        wxString input_parameter_file(input.c_str());
        std::cin >> input; 
        wxString input_reconstruction(input.c_str());
        std::cin >> input; 
        wxString input_reconstruction_statistics(input.c_str());

        std::cin >> input;
        input = std::toupper(input.front());
        if (input.compare("N") == 0) {use_statistics = false;}

        std::cin >> input; 
        wxString ouput_matching_projections(input.c_str());
        std::cin >> input; 
        wxString ouput_parameter_file(input.c_str());
        std::cin >> input; 
        wxString ouput_shift_file(input.c_str());
        std::cin >> input; 
        wxString my_symmetry(input.c_str());

        

        int		 first_particle						= 1;
        int		 last_particle						= 38;
        float	 percent_used						= 1.0;
        float 	 pixel_size							= 0.66;
        float    voltage_kV							= 300.0;
        float 	 spherical_aberration_mm			= 2.7;
        float    amplitude_contrast					= 0.07;
        float	 molecular_mass_kDa					= 700.0;
        float    inner_mask_radius					= 0.0;
        float    outer_mask_radius					= 85.0;
        float    low_resolution_limit				= 100.0;
        float    high_resolution_limit				= 3.0;
        float	 signed_CC_limit					= 0.0;
        float	 classification_resolution_limit	= 8.0;
        float    mask_radius_search					= 127.5;
        float	 high_resolution_limit_search		= 3.0;
        float	 angular_step						= 200.0;
        int		 best_parameters_to_keep			= 20;
        float	 max_search_x						= 127.5;
        float	 max_search_y						= 127.5;
        refine_particle.mask_center_2d_x			= 0.0;
        refine_particle.mask_center_2d_y			= 0.0;
        refine_particle.mask_center_2d_z			= 0.0;
        refine_particle.mask_radius_2d				= 0.0;
        float	 defocus_search_range				= 500.0;
        float	 defocus_step						= 50.0;
        float	 padding							= 1.0;
    //	float	 filter_constant					= my_current_job.arguments[35].ReturnFloatArgument();

        std::cin >> first_particle;
        std::cin >> last_particle;
        std::cin >> percent_used;
        std::cin >> pixel_size;
        std::cin >> voltage_kV;
        std::cin >> spherical_aberration_mm;
        std::cin >> amplitude_contrast;
        std::cin >> molecular_mass_kDa;
        std::cin >> inner_mask_radius;
        std::cin >> outer_mask_radius;
        std::cin >> low_resolution_limit;
        std::cin >> high_resolution_limit;
        std::cin >> signed_CC_limit;
        std::cin >> classification_resolution_limit;
        std::cin >> mask_radius_search;
        std::cin >> high_resolution_limit_search;
        std::cin >> angular_step;
        std::cin >> best_parameters_to_keep;
        std::cin >> max_search_x;
        std::cin >> max_search_y;
        std::cin >> refine_particle.mask_center_2d_x;
        std::cin >> refine_particle.mask_center_2d_y;
        std::cin >> refine_particle.mask_center_2d_z;
        std::cin >> refine_particle.mask_radius_2d;
        std::cin >> defocus_search_range;
        std::cin >> defocus_step;
        std::cin >> padding;



        
    //     bool	 global_search						= false;
    //     bool	 local_refinement					= true;
    // // Psi, Theta, Phi, ShiftX, ShiftY
    //     refine_particle.parameter_map[3]			= false;
    //     refine_particle.parameter_map[2]			= false;
    //     refine_particle.parameter_map[1]			= false;
    //     refine_particle.parameter_map[4]			= false;
    //     refine_particle.parameter_map[5]			= false;
    //     bool 	 calculate_matching_projections		= false;
    //     refine_particle.apply_2D_masking			= false;
    //     bool	 ctf_refinement						= false;
    //     bool	 normalize_particles				= false;
    //     bool	 invert_contrast					= false;
    //     bool	 exclude_blank_edges				= true;
    //     bool	 normalize_input_3d					= true;
    //     bool	 threshold_input_3d					= false;
    //     bool	 local_global_refine				= false;
    //     int		 current_class						= false;
    //     bool	 ignore_input_angles				= false;
        

        bool	 global_search						= true;
        bool	 local_refinement					= true;
    // Psi, Theta, Phi, ShiftX, ShiftY
        refine_particle.parameter_map[3]			= true;
        refine_particle.parameter_map[2]			= true;
        refine_particle.parameter_map[1]			= true;
        refine_particle.parameter_map[4]			= true;
        refine_particle.parameter_map[5]			= true;
        bool 	 calculate_matching_projections		= true;
        refine_particle.apply_2D_masking			= true;
        bool	 ctf_refinement						= true;
        bool	 normalize_particles				= true;
        bool	 invert_contrast					= true;
        bool	 exclude_blank_edges				= true;
        bool	 normalize_input_3d					= true;
        bool	 threshold_input_3d					= true;
        bool	 local_global_refine				= true;
        int		 current_class						= true;
        bool	 ignore_input_angles				= true;


        std::cin >> input;
        input = std::toupper(input.front());
        if (input.compare("N") == 0) {global_search = false;}
        std::cin >> input;
        input = std::toupper(input.front());
        if (input.compare("N") == 0) {local_refinement = false;}
        std::cin >> input;
        input = std::toupper(input.front());
        if (input.compare("N") == 0) {refine_particle.parameter_map[3] = false;}
        std::cin >> input;
        input = std::toupper(input.front());
        if (input.compare("N") == 0) {refine_particle.parameter_map[2] = false;}
        std::cin >> input;
        input = std::toupper(input.front());
        if (input.compare("N") == 0) {refine_particle.parameter_map[1] = false;}
        std::cin >> input;
        input = std::toupper(input.front());
        if (input.compare("N") == 0) {refine_particle.parameter_map[4] = false;}
        std::cin >> input;
        input = std::toupper(input.front());
        if (input.compare("N") == 0) {refine_particle.parameter_map[5] = false;}
        std::cin >> input;
        input = std::toupper(input.front());
        if (input.compare("N") == 0) {calculate_matching_projections = false;}
        std::cin >> input;
        input = std::toupper(input.front());
        if (input.compare("N") == 0) {refine_particle.apply_2D_masking = false;}
        std::cin >> input;
        input = std::toupper(input.front());
        if (input.compare("N") == 0) {ctf_refinement = false;}
        std::cin >> input;
        input = std::toupper(input.front());
        if (input.compare("N") == 0) {normalize_particles = false;}
        std::cin >> input;
        input = std::toupper(input.front());
        if (input.compare("N") == 0) {invert_contrast = false;}
        std::cin >> input;
        input = std::toupper(input.front());
        if (input.compare("N") == 0) {exclude_blank_edges = false;}
        std::cin >> input;
        input = std::toupper(input.front());
        if (input.compare("N") == 0) {normalize_input_3d = false;}
        std::cin >> input;
        input = std::toupper(input.front());
        if (input.compare("N") == 0) {threshold_input_3d = false;}



        refine_particle.constraints_used[4] = true;		// Constraint for X shifts
        refine_particle.constraints_used[5] = true;		// Constraint for Y shifts

        Image input_image;
    //	Image ctf_input_image;
        Image projection_image;
        Image search_projection_image;
        Image unbinned_image;
        Image binned_image;
        Image final_image;
        Image temp_image;
        Image temp_image2;
        Image sum_power;
        Image *projection_cache = NULL;
    //	CTF   my_ctf;
        CTF   input_ctf;
        Image snr_image;
        ReconstructedVolume			input_3d;
        ReconstructedVolume			search_reference_3d;
        ImageProjectionComparison	comparison_object;
        ConjugateGradient			conjugate_gradient_minimizer;
        EulerSearch					global_euler_search;
    //	Kernel2D					**kernel_index = NULL;
        Curve						noise_power_spectrum;
        Curve						number_of_terms;
        RandomNumberGenerator		random_particle(true);
        ProgressBar					*my_progress;

        JobResult *intermediate_result;

        int i;
        int j;
        int fourier_size_x, fourier_size_y, fourier_size_z;
        int current_image;
        int images_to_process = 0;
        int image_counter = 0;
        int defocus_i;
        int best_defocus_i;
        int search_box_size;
        int result_parameter_counter;
        int number_of_blank_edges;
        int max_samples = 2000;
        int istart;
        float input_parameters[refine_particle.number_of_parameters];
        float output_parameters[refine_particle.number_of_parameters];
        float gui_result_parameters[refine_particle.number_of_parameters];
        float search_parameters[refine_particle.number_of_parameters];
        float parameter_average[refine_particle.number_of_parameters];
        float parameter_variance[refine_particle.number_of_parameters];
        float output_parameter_average[refine_particle.number_of_parameters];
        float output_parameter_change[refine_particle.number_of_parameters];
        float cg_starting_point[refine_particle.number_of_parameters];
        float cg_accuracy[refine_particle.number_of_parameters];
        float binning_factor_refine;
        float binning_factor_search;
        float mask_falloff = 40.0;	// in Angstrom
    //	float alpha;
    //	float sigma;
        float logp;
        float temp_float;
        float psi;
        float psi_max;
        float psi_step;
        float psi_start;
        float score;
        float best_score;
        float mask_radius_for_noise;
        float percentage;
        float variance;
        float average;
        float average_density_max;
        bool skip_local_refinement = false;

        bool take_random_best_parameter;

        if (best_parameters_to_keep < 0)
        {
            best_parameters_to_keep = -best_parameters_to_keep;
            take_random_best_parameter = true;
        }
        else
        {
            take_random_best_parameter = false;
        }

        wxDateTime my_time_in;
    //	wxDateTime my_time_out;

        ZeroFloatArray(input_parameters, refine_particle.number_of_parameters);
        ZeroFloatArray(output_parameters, refine_particle.number_of_parameters);
        ZeroFloatArray(search_parameters, refine_particle.number_of_parameters);
        ZeroFloatArray(parameter_average, refine_particle.number_of_parameters);
        ZeroFloatArray(parameter_variance, refine_particle.number_of_parameters);
        ZeroFloatArray(output_parameter_average, refine_particle.number_of_parameters);
        ZeroFloatArray(output_parameter_change, refine_particle.number_of_parameters);
        ZeroFloatArray(cg_starting_point, refine_particle.number_of_parameters);
        ZeroFloatArray(cg_accuracy, refine_particle.number_of_parameters);
        // wxPrintf("Number of parameters: %i\n", refine_particle.number_of_parameters);

        
        //	wxPrintf("\nOpening input file %s.\n", input_parameter_file);
        FrealignParameterFile input_par_file(input_parameter_file, OPEN_TO_READ);
        MRCFile input_stack(input_particle_images.ToStdString(), false);

        input_par_file.ReadFile(false, input_stack.ReturnZSize());
        random_particle.SetSeed(int(10000.0 * fabsf(input_par_file.ReturnAverage(14, true)))%10000);

        MRCFile input_file(input_reconstruction.ToStdString(), false, true);
        MRCFile *output_file;
        if (percent_used < 1.0 && calculate_matching_projections)
        {
            calculate_matching_projections = false;
            wxPrintf("\nPercent of particles used < 1, matching projections not calculated.\n");
        }
        if (calculate_matching_projections) output_file = new MRCFile(ouput_matching_projections.ToStdString(), true);
        FrealignParameterFile my_output_par_file(ouput_parameter_file, OPEN_TO_WRITE);
        FrealignParameterFile my_output_par_shifts_file(ouput_shift_file, OPEN_TO_WRITE, 16);

        if (input_stack.ReturnXSize() != input_stack.ReturnYSize())
        {
            input_stack.PrintInfo();
            //SendErrorAndCrash("Error: Particles are not square\n");
        }
        if ((input_file.ReturnXSize() != input_file.ReturnYSize()) || (input_file.ReturnXSize() != input_file.ReturnZSize()))
        {
            input_file.PrintInfo();
            //SendErrorAndCrash("Error: Input reconstruction is not cubic\n");
        }
        if (input_file.ReturnXSize() != input_stack.ReturnXSize())
        {
            input_file.PrintInfo();
            input_stack.PrintInfo();
            //SendErrorAndCrash("Error: Dimension of particles and input reconstruction differ\n");
        }
        if (last_particle < first_particle && last_particle != 0)
        {
            //SendErrorAndCrash("Error: Number of last particle to refine smaller than number of first particle to refine\n");
        }

        if (last_particle == 0) last_particle = input_stack.ReturnZSize();
        if (first_particle == 0) first_particle = 1;
        if (last_particle > input_stack.ReturnZSize()) last_particle = input_stack.ReturnZSize();

        if (max_search_x == 0.0) max_search_x = mask_radius_search;
        if (max_search_y == 0.0) max_search_y = mask_radius_search;

        my_time_in = wxDateTime::Now();
        my_output_par_file.WriteCommentLine("C Refine3D run date and time:              " + my_time_in.FormatISOCombined(' '));
        my_output_par_file.WriteCommentLine("C Input particle images:                   " + input_particle_images);
        my_output_par_file.WriteCommentLine("C Input Frealign parameter filename:       " + input_parameter_file);
        my_output_par_file.WriteCommentLine("C Input reconstruction:                    " + input_reconstruction);
        my_output_par_file.WriteCommentLine("C Input data statistics:                   " + input_reconstruction_statistics);
        my_output_par_file.WriteCommentLine("C Use statistics:                          " + BoolToYesNo(use_statistics));
        my_output_par_file.WriteCommentLine("C Output matching projections:             " + ouput_matching_projections);
        my_output_par_file.WriteCommentLine("C Output parameter file:                   " + ouput_parameter_file);
        my_output_par_file.WriteCommentLine("C Output parameter changes:                " + ouput_shift_file);
        my_output_par_file.WriteCommentLine("C Particle symmetry:                       " + my_symmetry);
        my_output_par_file.WriteCommentLine("C First particle to refine:                " + wxString::Format("%i", first_particle));
        my_output_par_file.WriteCommentLine("C Last particle to refine:                 " + wxString::Format("%i", last_particle));
        my_output_par_file.WriteCommentLine("C Percent of particles to refine:          " + wxString::Format("%f", percent_used));
        my_output_par_file.WriteCommentLine("C Pixel size of images (A):                " + wxString::Format("%f", pixel_size));
        my_output_par_file.WriteCommentLine("C Beam energy (keV):                       " + wxString::Format("%f", voltage_kV));
        my_output_par_file.WriteCommentLine("C Spherical aberration (mm):               " + wxString::Format("%f", spherical_aberration_mm));
        my_output_par_file.WriteCommentLine("C Amplitude contrast:                      " + wxString::Format("%f", amplitude_contrast));
        my_output_par_file.WriteCommentLine("C Molecular mass of particle (kDa):        " + wxString::Format("%f", molecular_mass_kDa));
        my_output_par_file.WriteCommentLine("C Inner mask radius for refinement (A):    " + wxString::Format("%f", inner_mask_radius));
        my_output_par_file.WriteCommentLine("C Outer mask radius for refinement (A):    " + wxString::Format("%f", outer_mask_radius));
        my_output_par_file.WriteCommentLine("C Low resolution limit (A):                " + wxString::Format("%f", low_resolution_limit));
        my_output_par_file.WriteCommentLine("C High resolution limit (A):               " + wxString::Format("%f", high_resolution_limit));
        my_output_par_file.WriteCommentLine("C Resolution limit for signed CC (A):      " + wxString::Format("%f", signed_CC_limit));
        my_output_par_file.WriteCommentLine("C Res limit for classification (A):        " + wxString::Format("%f", classification_resolution_limit));
        my_output_par_file.WriteCommentLine("C Mask radius for global search (A):       " + wxString::Format("%f", mask_radius_search));
        my_output_par_file.WriteCommentLine("C Approx. resolution limit for search (A): " + wxString::Format("%f", high_resolution_limit_search));
        my_output_par_file.WriteCommentLine("C Angular step:                            " + wxString::Format("%f", angular_step));
        my_output_par_file.WriteCommentLine("C Number of top hits to refine:            " + wxString::Format("%i", best_parameters_to_keep));
        my_output_par_file.WriteCommentLine("C Search range in X (A):                   " + wxString::Format("%f", max_search_x));
        my_output_par_file.WriteCommentLine("C Search range in Y (A):                   " + wxString::Format("%f", max_search_y));
        my_output_par_file.WriteCommentLine("C 2D mask X coordinate (A):                " + wxString::Format("%f", refine_particle.mask_center_2d_x));
        my_output_par_file.WriteCommentLine("C 2D mask Y coordinate (A):                " + wxString::Format("%f", refine_particle.mask_center_2d_y));
        my_output_par_file.WriteCommentLine("C 2D mask Z coordinate (A):                " + wxString::Format("%f", refine_particle.mask_center_2d_z));
        my_output_par_file.WriteCommentLine("C 2D mask radius (A):                      " + wxString::Format("%f", refine_particle.mask_radius_2d));
        my_output_par_file.WriteCommentLine("C Defocus search range (A):                " + wxString::Format("%f", defocus_search_range));
        my_output_par_file.WriteCommentLine("C Defocus step (A):                        " + wxString::Format("%f", defocus_step));
        my_output_par_file.WriteCommentLine("C Padding factor:                          " + wxString::Format("%f", padding));
    //	my_output_par_file.WriteCommentLine("C Filter constant:                         " + wxString::Format("%f", filter_constant));
        my_output_par_file.WriteCommentLine("C Global search:                           " + BoolToYesNo(global_search));
        my_output_par_file.WriteCommentLine("C Local refinement:                        " + BoolToYesNo(local_refinement));
        my_output_par_file.WriteCommentLine("C Refine Psi:                              " + BoolToYesNo(refine_particle.parameter_map[3]));
        my_output_par_file.WriteCommentLine("C Refine Theta:                            " + BoolToYesNo(refine_particle.parameter_map[2]));
        my_output_par_file.WriteCommentLine("C Refine Phi:                              " + BoolToYesNo(refine_particle.parameter_map[1]));
        my_output_par_file.WriteCommentLine("C Refine ShiftX:                           " + BoolToYesNo(refine_particle.parameter_map[4]));
        my_output_par_file.WriteCommentLine("C Refine ShiftY:                           " + BoolToYesNo(refine_particle.parameter_map[5]));
        my_output_par_file.WriteCommentLine("C Calculate matching projections:          " + BoolToYesNo(calculate_matching_projections));
        my_output_par_file.WriteCommentLine("C Apply 2D masking:                        " + BoolToYesNo(refine_particle.apply_2D_masking));
        my_output_par_file.WriteCommentLine("C Refine defocus:                          " + BoolToYesNo(ctf_refinement));
        my_output_par_file.WriteCommentLine("C Normalize particles:                     " + BoolToYesNo(normalize_particles));
        my_output_par_file.WriteCommentLine("C Invert particle contrast:                " + BoolToYesNo(invert_contrast));
        my_output_par_file.WriteCommentLine("C Exclude images with blank edges:         " + BoolToYesNo(exclude_blank_edges));
        my_output_par_file.WriteCommentLine("C Normalize input reconstruction:          " + BoolToYesNo(normalize_input_3d));
        my_output_par_file.WriteCommentLine("C Threshold input reconstruction:          " + BoolToYesNo(threshold_input_3d));
        my_output_par_file.WriteCommentLine("C");
    //	my_output_par_file.WriteCommentLine("C    Particle#            Psi          Theta            Phi         ShiftX         ShiftY            Mag     Micrograph       Defocus1       Defocus2       AstigAng      Occupancy           LogP NormSigmaNoise          Score         Change");
        my_output_par_file.WriteCommentLine("C           PSI   THETA     PHI       SHX       SHY     MAG  INCLUDE   DF1      DF2  ANGAST  PSHIFT     OCC      LogP      SIGMA   SCORE  CHANGE");
        fflush(my_output_par_file.parameter_file);
    //	my_output_par_shifts_file.WriteCommentLine("C    Particle#            Psi          Theta            Phi         ShiftX         ShiftY            Mag     Micrograph       Defocus1       Defocus2       AstigAng      Occupancy           LogP NormSigmaNoise          Score");
        my_output_par_shifts_file.WriteCommentLine("C           PSI   THETA     PHI       SHX       SHY     MAG  INCLUDE   DF1      DF2  ANGAST  PSHIFT     OCC      LogP      SIGMA   SCORE");

        if (! refine_particle.parameter_map[1] && ! refine_particle.parameter_map[2] && ! refine_particle.parameter_map[3] && ! refine_particle.parameter_map[4] && ! refine_particle.parameter_map[5])
        {
            local_refinement = false;
            global_search = false;
        }

        if (local_global_refine)
        {
            refine_particle.parameter_map[1] = true;
            refine_particle.parameter_map[2] = true;
            refine_particle.parameter_map[3] = true;
            refine_particle.parameter_map[4] = true;
            refine_particle.parameter_map[5] = true;
            local_refinement = true;
            global_search = true;
        }

        if (high_resolution_limit < 2.0 * pixel_size) high_resolution_limit = 2.0 * pixel_size;
        if (classification_resolution_limit < 2.0 * pixel_size) classification_resolution_limit = 2.0 * pixel_size;
        if (high_resolution_limit_search < 2.0 * pixel_size) high_resolution_limit_search = 2.0 * pixel_size;
        if (signed_CC_limit == 0.0) signed_CC_limit = pixel_size;

        if (outer_mask_radius > float(input_stack.ReturnXSize()) / 2.0 * pixel_size- mask_falloff) outer_mask_radius = float(input_stack.ReturnXSize()) / 2.0 * pixel_size - mask_falloff;
        if (mask_radius_search > float(input_stack.ReturnXSize()) / 2.0 * pixel_size- mask_falloff) mask_radius_search = float(input_stack.ReturnXSize()) / 2.0 * pixel_size - mask_falloff;

        input_3d.InitWithDimensions(input_file.ReturnXSize(), input_file.ReturnYSize(), input_file.ReturnZSize(), pixel_size, my_symmetry);
        input_3d.molecular_mass_in_kDa = molecular_mass_kDa;
        ResolutionStatistics input_statistics(pixel_size, input_3d.density_map.logical_y_dimension);
        ResolutionStatistics search_statistics;
        ResolutionStatistics refine_statistics;
        if (use_statistics)
        {
            if (! DoesFileExist(input_reconstruction_statistics))
            {
                //SendError(wxString::Format("Error: Input statistics %s not found\n", input_reconstruction_statistics));
                exit(-1);
            }
            input_statistics.ReadStatisticsFromFile(input_reconstruction_statistics);
        }
        else
        {
            wxPrintf("\nUsing default statistics\n");
            input_statistics.GenerateDefaultStatistics(molecular_mass_kDa);
        }
        refine_statistics = input_statistics;
        input_3d.density_map.ReadSlices(&input_file,1,input_3d.density_map.logical_z_dimension);
    //!!! This line is incompatible with ML !!!
    //	input_3d.density_map.CosineMask(outer_mask_radius / pixel_size, mask_falloff / pixel_size);
    //	input_3d.density_map.AddConstant(- input_3d.density_map.ReturnAverageOfRealValuesOnEdges());
        // Remove masking here to avoid edge artifacts later
        input_3d.density_map.CosineMask(outer_mask_radius / pixel_size, mask_falloff / pixel_size, false, true, 0.0);
        if (inner_mask_radius > 0.0) input_3d.density_map.CosineMask(inner_mask_radius / pixel_size, mask_falloff / pixel_size, true);
    //	for (i = 0; i < input_3d.density_map.real_memory_allocated; i++) if (input_3d.density_map.real_values[i] < 0.0) input_3d.density_map.real_values[i] = -log(-input_3d.density_map.real_values[i] + 1.0);
        if (threshold_input_3d)
        {
            average_density_max = input_3d.density_map.ReturnAverageOfMaxN(100, outer_mask_radius / pixel_size);
            input_3d.density_map.SetMinimumValue(-0.3 * average_density_max);
    //		input_3d.density_map.SetMinimumValue(0.0);
        }

        input_image.Allocate(input_stack.ReturnXSize(), input_stack.ReturnYSize(), true);
    //	if (outer_mask_radius > input_image.physical_address_of_box_center_x * pixel_size- mask_falloff) outer_mask_radius = input_image.physical_address_of_box_center_x * pixel_size - mask_falloff;
    //	if (mask_radius_search > input_image.physical_address_of_box_center_x * pixel_size- mask_falloff) mask_radius_search = input_image.physical_address_of_box_center_x * pixel_size - mask_falloff;
        input_3d.mask_radius = outer_mask_radius;

        if (global_search)
        {
            if (best_parameters_to_keep == 0) {best_parameters_to_keep = 1; skip_local_refinement = true;}
            // Assume square particles
            search_reference_3d = input_3d;
            search_statistics = input_statistics;
            search_box_size = ReturnClosestFactorizedUpper(myroundint(2.0 * padding * (std::max(max_search_x, max_search_y) + mask_radius_search)), 3, true);
            if (search_box_size > search_reference_3d.density_map.logical_x_dimension) search_box_size = search_reference_3d.density_map.logical_x_dimension;
            if (search_box_size != search_reference_3d.density_map.logical_x_dimension) search_reference_3d.density_map.Resize(search_box_size, search_box_size, search_box_size);
    //		search_reference_3d.PrepareForProjections(high_resolution_limit_search, true);
            search_reference_3d.PrepareForProjections(low_resolution_limit, high_resolution_limit_search, true);
    //		search_statistics.Init(search_reference_3d.pixel_size, search_reference_3d.density_map.logical_y_dimension / 2 + 1);
            search_particle.Allocate(search_reference_3d.density_map.logical_x_dimension, search_reference_3d.density_map.logical_y_dimension);
            search_projection_image.Allocate(search_reference_3d.density_map.logical_x_dimension, search_reference_3d.density_map.logical_y_dimension, false);
            temp_image2.Allocate(search_box_size, search_box_size, true);
            binning_factor_search = search_reference_3d.pixel_size / pixel_size;
            //Scale to make projections compatible with images for ML calculation
            search_reference_3d.density_map.MultiplyByConstant(powf(powf(binning_factor_search, 1.0 / 3.0), 2));
            //if (angular_step <= 0) angular_step = 360.0 * high_resolution_limit_search / PI / outer_mask_radius;
            if (angular_step <= 0) angular_step = CalculateAngularStep(high_resolution_limit_search, outer_mask_radius);
            psi_step = rad_2_deg(search_reference_3d.pixel_size / outer_mask_radius);
            psi_step = 360.0 / int(360.0 / psi_step + 0.5);
            psi_start = psi_step / 2.0 * global_random_number_generator.GetUniformRandom();
            psi_max = 0.0;
            if (refine_particle.parameter_map[3]) psi_max = 360.0;
            wxPrintf("\nBox size for search = %i, binning factor = %f, new pixel size = %f, resolution limit = %f\nAngular step size = %f, in-plane = %f\n", search_box_size, binning_factor_search, search_reference_3d.pixel_size, search_reference_3d.pixel_size * 2.0, angular_step, psi_step);
        }

        if (padding != 1.0)
        {
            input_3d.density_map.Resize(input_3d.density_map.logical_x_dimension * padding, input_3d.density_map.logical_y_dimension * padding, input_3d.density_map.logical_z_dimension * padding, input_3d.density_map.ReturnAverageOfRealValuesOnEdges());
            refine_statistics.part_SSNR.ResampleCurve(&refine_statistics.part_SSNR, refine_statistics.part_SSNR.number_of_points * padding);
        }

    //	input_3d.PrepareForProjections(high_resolution_limit);
        input_3d.PrepareForProjections(low_resolution_limit, high_resolution_limit);
        binning_factor_refine = input_3d.pixel_size / pixel_size;
        //Scale to make projections compatible with images for ML calculation
    //	input_3d.density_map.MultiplyByConstant(binning_factor_refine);
    //	input_3d.density_map.MultiplyByConstant(powf(powf(binning_factor_refine, 1.0 / 3.0), 2));
        wxPrintf("\nBinning factor for refinement = %f, new pixel size = %f\n", binning_factor_refine, input_3d.pixel_size);

        temp_image.Allocate(input_file.ReturnXSize(), input_file.ReturnYSize(), true);
        sum_power.Allocate(input_stack.ReturnXSize(), input_stack.ReturnYSize(), false);
        refine_particle.Allocate(input_3d.density_map.logical_x_dimension, input_3d.density_map.logical_y_dimension);
    //	ctf_input_image.Allocate(input_stack.ReturnXSize(), input_stack.ReturnYSize(), false);
        projection_image.Allocate(input_3d.density_map.logical_x_dimension, input_3d.density_map.logical_y_dimension, false);
        unbinned_image.Allocate(input_file.ReturnXSize() * padding, input_file.ReturnYSize() * padding, true);
        if (ctf_refinement) binned_image.Allocate(input_3d.density_map.logical_x_dimension, input_3d.density_map.logical_y_dimension, false);
        final_image.Allocate(input_file.ReturnXSize(), input_file.ReturnYSize(), true);

    // Read whole parameter file to work out average values and variances
        j = 0;
        for (current_image = 1; current_image <= input_par_file.number_of_lines; current_image++)
        {
            input_par_file.ReadLine(input_parameters); temp_float = input_parameters[1]; input_parameters[1] = input_parameters[3]; input_parameters[3] = temp_float;
            if (input_parameters[7] >= 0)
            {
                for (i = 0; i < refine_particle.number_of_parameters; i++)
                {
                    parameter_average[i] += input_parameters[i];
                    parameter_variance[i] += powf(input_parameters[i],2);
                }
                j++;
            }
            if (input_parameters[0] >= first_particle && input_parameters[0] <= last_particle) images_to_process++;
        }

        for (i = 0; i < refine_particle.number_of_parameters; i++)
        {
            parameter_average[i] /= j;
            parameter_variance[i] /= j;
            parameter_variance[i] -= powf(parameter_average[i],2);

            if (parameter_variance[i] < 0.001) refine_particle.constraints_used[i] = false;
        }
        refine_particle.SetParameterStatistics(parameter_average, parameter_variance);
        input_par_file.Rewind();

        if (normalize_particles)
        {
            wxPrintf("Calculating noise power spectrum...\n\n");
            percentage = float(max_samples) / float(images_to_process);
            sum_power.SetToConstant(0.0);
            mask_radius_for_noise = outer_mask_radius / pixel_size;
            number_of_blank_edges = 0;
            if (2.0 * mask_radius_for_noise + mask_falloff / pixel_size > 0.95 * input_image.logical_x_dimension)
            {
                mask_radius_for_noise = 0.95 * input_image.logical_x_dimension / 2.0 - mask_falloff / 2.0 / pixel_size;
            }
            noise_power_spectrum.SetupXAxis(0.0, 0.5 * sqrtf(2.0), int((sum_power.logical_x_dimension / 2.0 + 1.0) * sqrtf(2.0) + 1.0));
            number_of_terms.SetupXAxis(0.0, 0.5 * sqrtf(2.0), int((sum_power.logical_x_dimension / 2.0 + 1.0) * sqrtf(2.0) + 1.0));
            my_progress = new ProgressBar(images_to_process);
            for (current_image = 1; current_image <= input_par_file.number_of_lines; current_image++)
            {
                input_par_file.ReadLine(input_parameters);
                if (input_parameters[0] < first_particle || input_parameters[0] > last_particle) continue;
                image_counter++;
                my_progress->Update(image_counter);
                if ((global_random_number_generator.GetUniformRandom() < 1.0 - 2.0 * percentage)) continue;
                input_image.ReadSlice(&input_stack, int(input_parameters[0] + 0.5));
                if (exclude_blank_edges && input_image.ContainsBlankEdges(outer_mask_radius / pixel_size)) {number_of_blank_edges++; continue;}
                variance = input_image.ReturnVarianceOfRealValues(outer_mask_radius / pixel_size, 0.0, 0.0, 0.0, true);
                if (variance == 0.0) continue;
                input_image.MultiplyByConstant(1.0 / sqrtf(variance));
                input_image.CosineMask(outer_mask_radius / pixel_size, mask_falloff / pixel_size, true);
                input_image.ForwardFFT();
                temp_image.CopyFrom(&input_image);
                temp_image.ConjugateMultiplyPixelWise(input_image);
                sum_power.AddImage(&temp_image);
            }
            delete my_progress;
            sum_power.Compute1DRotationalAverage(noise_power_spectrum, number_of_terms);
            noise_power_spectrum.SquareRoot();
            noise_power_spectrum.Reciprocal();

            input_par_file.Rewind();
            if (exclude_blank_edges)
            {
                wxPrintf("\nImages with blank edges excluded from noise power calculation = %i\n", number_of_blank_edges);
            }
        }

        if (global_search)
        {
            for (i = 0; i < search_particle.number_of_parameters; i++) {search_particle.parameter_map[i] = refine_particle.parameter_map[i];}
    // Set parameter_map for x,y translations to true since they will always be searched and refined in a global search
    // Decided not to do this to honor user request
    //		search_particle.parameter_map[4] = true;
    //		search_particle.parameter_map[5] = true;
            for (i = 0; i < search_particle.number_of_parameters; i++) {search_particle.constraints_used[i] = refine_particle.constraints_used[i];}
            search_particle.SetParameterStatistics(parameter_average, parameter_variance);

    // Use projection_cache only if both phi and theta are searched; otherwise calculate projections on the fly
            if (search_particle.parameter_map[1] && search_particle.parameter_map[2])
            {
                global_euler_search.InitGrid(my_symmetry, angular_step, 0.0, 0.0, psi_max, psi_step, psi_start, search_reference_3d.pixel_size / high_resolution_limit_search, search_particle.parameter_map, best_parameters_to_keep);
                if (global_euler_search.best_parameters_to_keep != best_parameters_to_keep) best_parameters_to_keep = global_euler_search.best_parameters_to_keep;
                projection_cache = new Image [global_euler_search.number_of_search_positions];
                for (i = 0; i < global_euler_search.number_of_search_positions; i++)
                {
                    projection_cache[i].Allocate(search_reference_3d.density_map.logical_x_dimension, search_reference_3d.density_map.logical_y_dimension, false);
                }
                search_reference_3d.density_map.GenerateReferenceProjections(projection_cache, global_euler_search, search_reference_3d.pixel_size / high_resolution_limit_search);
                wxPrintf("\nNumber of global search views = %i (best_parameters to keep = %i)\n", global_euler_search.number_of_search_positions, global_euler_search.best_parameters_to_keep);
            }
    //		search_projection_image.RotateFourier2DGenerateIndex(kernel_index, psi_max, psi_step, psi_start);

            if (search_particle.parameter_map[4]) global_euler_search.max_search_x = max_search_x;
            else global_euler_search.max_search_x = 0.0;
            if (search_particle.parameter_map[5]) global_euler_search.max_search_y = max_search_y;
            else global_euler_search.max_search_y = 0.0;
        }

        wxPrintf("\nAverage sigma noise = %f, average LogP = %f\nAverage ShiftX = %f, average ShiftY = %f\nSigma ShiftX = %f, sigma ShiftY = %f\nNumber of particles to refine = %i\n\n",
                parameter_average[14], parameter_average[15], parameter_average[4], parameter_average[5], sqrtf(parameter_variance[4]), sqrtf(parameter_variance[5]), images_to_process);

        image_counter = 0;
        my_progress = new ProgressBar(images_to_process);
        for (current_image = 1; current_image <= input_par_file.number_of_lines; current_image++)
        {
            input_par_file.ReadLine(input_parameters);
            temp_float = random_particle.GetUniformRandom();
            if (input_parameters[0] < first_particle || input_parameters[0] > last_particle) continue;
            image_counter++;
            if (temp_float < 1.0 - 2.0 * percent_used)
            {
                input_parameters[7] = -1;//- fabsf(input_parameters[7]);
                input_parameters[16] = 0.0;
                my_output_par_file.WriteLine(input_parameters);


                // if (is_running_locally == false) // send results back to the gui..
                // {
                //     intermediate_result = new JobResult;
                //     intermediate_result->job_number = my_current_job.job_number;

                //     gui_result_parameters[0] = current_class;
                //     for (result_parameter_counter = 1; result_parameter_counter < refine_particle.number_of_parameters + 1; result_parameter_counter++)
                //     {
                //         gui_result_parameters[result_parameter_counter] = input_parameters[result_parameter_counter - 1];
                //     }

                //     intermediate_result->SetResult(refine_particle.number_of_parameters + 1, gui_result_parameters);
                //     AddJobToResultQueue(intermediate_result);
                // }

                for (i = 1; i < refine_particle.number_of_parameters; i++) output_parameters[i] = 0.0;
                output_parameters[0] = input_parameters[0];

                my_output_par_shifts_file.WriteLine(output_parameters);

                my_progress->Update(image_counter);
                continue;
            }
            else
            {
    //			input_parameters[7] = fabsf(input_parameters[7]);
                temp_float = input_parameters[1]; input_parameters[1] = input_parameters[3]; input_parameters[3] = temp_float;
            }
            for (i = 0; i < refine_particle.number_of_parameters; i++) {output_parameters[i] = input_parameters[i];}

            if (local_global_refine)
            {
                if (input_parameters[7] == 0.0) {local_refinement = false; global_search = true;}
                else {local_refinement = true; global_search = false;}
            }

    // Set up Particle object
            refine_particle.ResetImageFlags();
            refine_particle.mask_radius = outer_mask_radius;
            refine_particle.mask_falloff = mask_falloff;
    //		refine_particle.filter_radius_low = low_resolution_limit;
            refine_particle.filter_radius_high = high_resolution_limit;
            refine_particle.molecular_mass_kDa = molecular_mass_kDa;
            refine_particle.signed_CC_limit = signed_CC_limit;
            // The following line would allow using particles with different pixel sizes
            refine_particle.pixel_size = input_3d.pixel_size;
            refine_particle.is_normalized = normalize_particles;
            refine_particle.sigma_noise = input_parameters[14] / binning_factor_refine;
    //		refine_particle.logp = -std::numeric_limits<float>::max();
            refine_particle.SetParameters(input_parameters);
            refine_particle.MapParameterAccuracy(cg_accuracy);
    //		refine_particle.SetIndexForWeightedCorrelation();
            refine_particle.SetParameterConstraints(powf(parameter_average[14],2));

            input_ctf.Init(voltage_kV, spherical_aberration_mm, amplitude_contrast, input_parameters[8], input_parameters[9], input_parameters[10], 0.0, 0.0, 0.0, pixel_size, input_parameters[11]);
    //		ctf_input_image.CalculateCTFImage(input_ctf);
    //		refine_particle.is_phase_flipped = true;

            input_image.ReadSlice(&input_stack, int(input_parameters[0] + 0.5));
            input_image.ReplaceOutliersWithMean(5.0);
            if (invert_contrast) input_image.InvertRealValues();
            if (normalize_particles)
            {
                input_image.ForwardFFT();
                // Whiten noise
                input_image.ApplyCurveFilter(&noise_power_spectrum);
                // Apply cosine filter to reduce ringing
    //			input_image.CosineMask(std::max(pixel_size / high_resolution_limit, pixel_size / 7.0f + pixel_size / mask_falloff) - pixel_size / (2.0 * mask_falloff), pixel_size / mask_falloff);
                input_image.BackwardFFT();
                // Normalize background variance and average
                variance = input_image.ReturnVarianceOfRealValues(input_image.physical_address_of_box_center_x - mask_falloff / pixel_size, 0.0, 0.0, 0.0, true);
                average = input_image.ReturnAverageOfRealValues(input_image.physical_address_of_box_center_x - mask_falloff / pixel_size, true);
                input_image.AddMultiplyConstant(- average, 1.0 / sqrtf(variance));
                // At this point, input_image should have white background with a variance of 1. The variance should therefore be about 1/binning_factor^2 after binning.
            }

            // Option to add noise to images to get out of local optima
    //		input_image.AddGaussianNoise(sqrtf(2.0 * input_image.ReturnVarianceOfRealValues()));

            input_image.ClipInto(&unbinned_image);
            unbinned_image.ForwardFFT();
            unbinned_image.ClipInto(refine_particle.particle_image);
            // Multiply by binning_factor so variance after binning is close to 1.
    //		refine_particle.particle_image->MultiplyByConstant(binning_factor_refine);
            comparison_object.reference_volume = &input_3d;
            comparison_object.projection_image = &projection_image;
            comparison_object.particle = &refine_particle;
            refine_particle.MapParameters(cg_starting_point);
            refine_particle.PhaseShiftInverse();

            if (ctf_refinement && high_resolution_limit <= 20.0)
            {
    //			wxPrintf("\nRefining defocus for parameter line %i\n", current_image);
                refine_particle.filter_radius_low = 30.0;
                refine_particle.SetIndexForWeightedCorrelation();
                binned_image.CopyFrom(refine_particle.particle_image);
                refine_particle.InitCTF(voltage_kV, spherical_aberration_mm, amplitude_contrast, input_parameters[8], input_parameters[9], input_parameters[10], input_parameters[11]);
                best_score = - std::numeric_limits<float>::max();
                for (defocus_i = - myround(float(defocus_search_range)/float(defocus_step)); defocus_i <= myround(float(defocus_search_range)/float(defocus_step)); defocus_i++)
                {
                    refine_particle.SetDefocus(input_parameters[8] + defocus_i * defocus_step, input_parameters[9] + defocus_i * defocus_step, input_parameters[10], input_parameters[11]);
                    refine_particle.InitCTFImage(voltage_kV, spherical_aberration_mm, amplitude_contrast, input_parameters[8] + defocus_i * defocus_step, input_parameters[9] + defocus_i * defocus_step, input_parameters[10], input_parameters[11]);
                    if (normalize_input_3d) refine_particle.WeightBySSNR(refine_statistics.part_SSNR, 1);
    //				// Apply SSNR weighting only to image since input 3D map assumed to be calculated from correctly whitened images
                    else refine_particle.WeightBySSNR(refine_statistics.part_SSNR, 0);
                    refine_particle.PhaseFlipImage();
    //				refine_particle.CosineMask(false, true, 0.0);
                    refine_particle.CosineMask();
                    refine_particle.PhaseShift();
                    refine_particle.CenterInCorner();
    //				refine_particle.WeightBySSNR(input_3d.statistics.part_SSNR, 1);

                    score = - 100.0 * FrealignObjectiveConstrainedFunction(&comparison_object, cg_starting_point);
                    if (score > best_score)
                    {
                        best_score = score;
                        best_defocus_i = defocus_i;
    //					wxPrintf("Parameter line = %i, Defocus = %f, score = %g\n", current_image, defocus_i * defocus_step, score);
                    }
                    refine_particle.particle_image->CopyFrom(&binned_image);
                    refine_particle.is_ssnr_filtered = false;
                    refine_particle.is_masked = false;
                    refine_particle.is_centered_in_box = true;
                    refine_particle.shift_counter = 1;
                }
                output_parameters[8] = input_parameters[8] + best_defocus_i * defocus_step;
                output_parameters[9] = input_parameters[9] + best_defocus_i * defocus_step;
                refine_particle.SetDefocus(output_parameters[8], output_parameters[9], input_parameters[10], input_parameters[11]);
                refine_particle.InitCTFImage(voltage_kV, spherical_aberration_mm, amplitude_contrast, output_parameters[8], output_parameters[9], input_parameters[10], input_parameters[11]);
            }
            else
            {
                refine_particle.InitCTFImage(voltage_kV, spherical_aberration_mm, amplitude_contrast, input_parameters[8], input_parameters[9], input_parameters[10], input_parameters[11]);
            }
            refine_particle.filter_radius_low = low_resolution_limit;
            refine_particle.SetIndexForWeightedCorrelation();
            if (normalize_input_3d) refine_particle.WeightBySSNR(refine_statistics.part_SSNR, 1);
            // Apply SSNR weighting only to image since input 3D map assumed to be calculated from correctly whitened images
            else refine_particle.WeightBySSNR(refine_statistics.part_SSNR, 0);
            refine_particle.PhaseFlipImage();
    //		refine_particle.CosineMask(false, true, 0.0);
            refine_particle.CosineMask();
            refine_particle.PhaseShift();
            refine_particle.CenterInCorner();
    //		refine_particle.WeightBySSNR(input_3d.statistics.part_SSNR, 1);

    //		input_parameters[15] = 10.0;
            input_parameters[15] = - 100.0 * FrealignObjectiveConstrainedFunction(&comparison_object, cg_starting_point);

            if ((refine_particle.number_of_search_dimensions > 0) && (global_search || local_refinement))
            {
                if (global_search)
                {
    //				my_time_in = wxDateTime::UNow();
                    search_particle.ResetImageFlags();
                    search_particle.pixel_size = search_reference_3d.pixel_size;
                    if (mask_radius_search == 0.0)
                    {
                        search_particle.mask_radius = search_particle.particle_image->logical_x_dimension / 2 * search_particle.pixel_size - mask_falloff;
                    }
                    else
                    {
                        search_particle.mask_radius = mask_radius_search;
                    }
                    search_particle.mask_falloff = mask_falloff;
                    search_particle.filter_radius_low = 0.0;
                    search_particle.filter_radius_high = high_resolution_limit_search;
                    search_particle.molecular_mass_kDa = molecular_mass_kDa;
                    search_particle.signed_CC_limit = signed_CC_limit;
                    search_particle.sigma_noise = input_parameters[14] / binning_factor_search;
    //				search_particle.logp = -std::numeric_limits<float>::max();
                    search_particle.SetParameters(input_parameters);
                    search_particle.number_of_search_dimensions = refine_particle.number_of_search_dimensions;
                    search_particle.InitCTFImage(voltage_kV, spherical_aberration_mm, amplitude_contrast, input_parameters[8], input_parameters[9], input_parameters[10], input_parameters[11]);
                    temp_image.CopyFrom(&input_image);
                    // Multiply by binning_factor so variance after binning is close to 1.
    //				temp_image.MultiplyByConstant(binning_factor_search);
                    // Assume square images
                    if (search_box_size != temp_image.logical_x_dimension)
                    {
                        temp_image.ClipInto(&temp_image2);
                        temp_image2.ForwardFFT();
                        temp_image2.ClipInto(search_particle.particle_image);
                    }
                    else
                    {
                        temp_image.ForwardFFT();
                        temp_image.ClipInto(search_particle.particle_image);
                    }
                    search_particle.PhaseShiftInverse();
                    // Always apply particle SSNR weighting (i.e. whitening) reference normalization since reference
                    // projections will not have SSNR (i.e. CTF-dependent) weighting applied
                    search_particle.WeightBySSNR(search_statistics.part_SSNR, 1);
                    search_particle.PhaseFlipImage();
    //				search_particle.CosineMask(false, true, 0.0);
                    search_particle.CosineMask();
                    search_particle.PhaseShift();
    //				search_particle.CenterInCorner();
    //				search_particle.WeightBySSNR(search_reference_3d.statistics.part_SSNR);

                    if (search_particle.parameter_map[1] && ! search_particle.parameter_map[2])
                    {
                        global_euler_search.InitGrid(my_symmetry, angular_step, 0.0, input_parameters[2], psi_max, psi_step, psi_start, search_reference_3d.pixel_size / high_resolution_limit_search, search_particle.parameter_map, best_parameters_to_keep);
                        if (global_euler_search.best_parameters_to_keep != best_parameters_to_keep) best_parameters_to_keep = global_euler_search.best_parameters_to_keep;
                        if (! search_particle.parameter_map[3]) global_euler_search.psi_start = 360.0 - input_parameters[3];
                        global_euler_search.Run(search_particle, search_reference_3d.density_map, input_parameters + 1, projection_cache);
                    }
                    else
                    if (! search_particle.parameter_map[1] && search_particle.parameter_map[2])
                    {
                        global_euler_search.InitGrid(my_symmetry, angular_step, input_parameters[1], 0.0, psi_max, psi_step, psi_start, search_reference_3d.pixel_size / high_resolution_limit_search, search_particle.parameter_map, best_parameters_to_keep);
                        if (global_euler_search.best_parameters_to_keep != best_parameters_to_keep) best_parameters_to_keep = global_euler_search.best_parameters_to_keep;
                        if (! search_particle.parameter_map[3]) global_euler_search.psi_start = 360.0 - input_parameters[3];
                        global_euler_search.Run(search_particle, search_reference_3d.density_map, input_parameters + 1, projection_cache);
                    }
                    else
                    if (search_particle.parameter_map[1] && search_particle.parameter_map[2])
                    {
                        if (! search_particle.parameter_map[3]) global_euler_search.psi_start = 360.0 - input_parameters[3];
                        if (global_euler_search.best_parameters_to_keep != best_parameters_to_keep) best_parameters_to_keep = global_euler_search.best_parameters_to_keep;
                        global_euler_search.Run(search_particle, search_reference_3d.density_map, input_parameters + 1, projection_cache);
                    }
                    else
                    {
                        best_parameters_to_keep = 0;
                    }

                    // Do local refinement of the top hits to determine the best match
                    for (i = 0; i < search_particle.number_of_parameters; i++) {search_parameters[i] = input_parameters[i];}
                    for (j = 1; j < 6; j++) {global_euler_search.list_of_best_parameters[0][j - 1] = input_parameters[j];}
                    search_particle.SetParameterConstraints(powf(parameter_average[14],2));
                    comparison_object.reference_volume = &search_reference_3d;
                    comparison_object.projection_image = &search_projection_image;
                    comparison_object.particle = &search_particle;
                    search_particle.CenterInCorner();
                    search_particle.SetIndexForWeightedCorrelation();
                    search_particle.SetParameters(input_parameters);
                    search_particle.MapParameters(cg_starting_point);
                    search_particle.mask_radius = outer_mask_radius;
    //				output_parameters[15] = - 100.0 * FrealignObjectiveFunction(&comparison_object, cg_starting_point);
    //				if (! local_refinement) input_parameters[15] = output_parameters[15];
                    search_particle.UnmapParametersToExternal(output_parameters, cg_starting_point);
                    if (ignore_input_angles && best_parameters_to_keep >= 1) istart = 1;
                    else istart = 0;


                    if (take_random_best_parameter == true)
                    {
                        float best_value = global_euler_search.list_of_best_parameters[1][5];
                        float worst_value = global_euler_search.list_of_best_parameters[best_parameters_to_keep][5];;
                        float diff = best_value - worst_value;
                        float top_percent = best_value - (diff * 0.15);
                        int number_to_keep = 1;

                        for (int counter = 2; counter <= best_parameters_to_keep; counter++)
                        {
                            if (global_euler_search.list_of_best_parameters[counter][5] >= top_percent)
                            {
                                number_to_keep++;
                            }
                        }

                        int parameter_to_keep = myroundint(((global_random_number_generator.GetUniformRandom() + 1.0) / 2.0) * float(number_to_keep - 1)) + 1;
                        //wxPrintf("best_value = %f, worst_value = %f, top 10%% = %f, number_above = %i, number to take = %i\n", best_value, worst_value, top_percent, number_to_keep, parameter_to_keep);


                        for (j = 1; j < 6; j++)
                        {
                            search_parameters[j] = global_euler_search.list_of_best_parameters[parameter_to_keep][j - 1];
                        }

                        if (! search_particle.parameter_map[4]) search_parameters[4] = input_parameters[4];
                        if (! search_particle.parameter_map[5]) search_parameters[5] = input_parameters[5];
                        search_particle.SetParameters(search_parameters);
                        search_particle.MapParameters(cg_starting_point);
                        search_parameters[15] = - 100.0 * conjugate_gradient_minimizer.Init(&FrealignObjectiveConstrainedFunction, &comparison_object, search_particle.number_of_search_dimensions, cg_starting_point, cg_accuracy);
                        output_parameters[15] = search_parameters[15];
                        if (! local_refinement) input_parameters[15] = output_parameters[15];

                        temp_float = - 100.0 * conjugate_gradient_minimizer.Run();
                        search_particle.UnmapParametersToExternal(output_parameters, conjugate_gradient_minimizer.GetPointerToBestValues());
                        output_parameters[15] = temp_float;
                    }
                    else
                    {

                    for (i = istart; i <= best_parameters_to_keep; i++)
                    {
                        for (j = 1; j < 6; j++)
                        {
                            search_parameters[j] = global_euler_search.list_of_best_parameters[i][j - 1];

                        }
    //					wxPrintf("parameters in  = %i %g, %g, %g, %g, %g %g\n", i, search_parameters[3], search_parameters[2],
    //							search_parameters[1], search_parameters[4], search_parameters[5], global_euler_search.list_of_best_parameters[i][5]);
                        if (! search_particle.parameter_map[4]) search_parameters[4] = input_parameters[4];
                        if (! search_particle.parameter_map[5]) search_parameters[5] = input_parameters[5];
                        search_particle.SetParameters(search_parameters);
                        search_particle.MapParameters(cg_starting_point);
                        search_parameters[15] = - 100.0 * conjugate_gradient_minimizer.Init(&FrealignObjectiveConstrainedFunction, &comparison_object, search_particle.number_of_search_dimensions, cg_starting_point, cg_accuracy);
                        if (i == istart)
                        {
                            output_parameters[15] = search_parameters[15];
                            if (! local_refinement) input_parameters[15] = output_parameters[15];
                        }
                        if (skip_local_refinement) temp_float = search_parameters[15];
                        else temp_float = - 100.0 * conjugate_gradient_minimizer.Run();
                        // Uncomment the following line to skip local refinement.
    //					temp_float = search_parameters[15];
    //					wxPrintf("best, refine in, out, diff = %i %g %g %g %g\n", i, output_parameters[15], search_parameters[15], temp_float, temp_float - output_parameters[15]);
    //					log_diff = output_parameters[15] - temp_float;
    //					if (log_diff > log_range) log_diff = log_range;
    //					if (log_diff < - log_range) log_diff = - log_range;
                        if (temp_float > output_parameters[15])
                        // If log_diff >= 0, exp(log_diff) will always be larger than the random number and the search parameters will be kept.
                        // If log_diff < 0, there is an increasing chance that the random number is larger than exp(log_diff) and the new
                        // (worse) parameters will not be kept.
    //					if ((global_random_number_generator.GetUniformRandom() + 1.0) / 2.0 < 1.0 / (1.0 + exp(log_diff)))
                        {
    //						if (log_diff < 0.0) wxPrintf("log_diff = %g\n", log_diff);
    //						wxPrintf("image_counter = %i, i = %i, score = %g\n", image_counter, i, temp_float);
                            search_particle.UnmapParametersToExternal(output_parameters, conjugate_gradient_minimizer.GetPointerToBestValues());
                            output_parameters[15] = temp_float;
    //						wxPrintf("parameters out = %g, %g, %g, %g, %g\n", output_parameters[3], output_parameters[2],
    //								output_parameters[1], output_parameters[4], output_parameters[5]);
                        }
    //					wxPrintf("refine in, out, keep = %i %g %g %g\n", i, search_parameters[15], temp_float, output_parameters[15]);
    //					wxPrintf("parameters out = %g, %g, %g, %g, %g\n", output_parameters[3], output_parameters[2],
    //							output_parameters[1], output_parameters[4], output_parameters[5]);
                    }

                    }
                    refine_particle.SetParameters(output_parameters, true);
    //				my_time_out = wxDateTime::UNow(); wxPrintf("global search done: ms taken = %li\n", my_time_out.Subtract(my_time_in).GetMilliseconds());
                }

                if (local_refinement)
                {
    //				my_time_in = wxDateTime::UNow();
                    comparison_object.reference_volume = &input_3d;
                    comparison_object.projection_image = &projection_image;
                    comparison_object.particle = &refine_particle;
                    refine_particle.MapParameters(cg_starting_point);

                    temp_float = - 100.0 * conjugate_gradient_minimizer.Init(&FrealignObjectiveConstrainedFunction, &comparison_object, refine_particle.number_of_search_dimensions, cg_starting_point, cg_accuracy);
    //???				if (! global_search) input_parameters[15] = temp_float;
                    output_parameters[15] = - 100.0 * conjugate_gradient_minimizer.Run();

                    refine_particle.UnmapParametersToExternal(output_parameters, conjugate_gradient_minimizer.GetPointerToBestValues());

    //				my_time_out = wxDateTime::UNow(); wxPrintf("local refinement done: ms taken = %li\n", my_time_out.Subtract(my_time_in).GetMilliseconds());
                }
    //			log_diff = input_parameters[15] - output_parameters[15];
    //			wxPrintf("in = %g out = %g log_diff = %g ratio = %g\n", input_parameters[15], output_parameters[15], log_diff, 1.0 / (1.0 + exp(log_diff)));
    //			if (log_diff > log_range) log_diff = log_range;
    //			if (log_diff < - log_range) log_diff = - log_range;
                // If log_diff >= 0, exp(log_diff) will never be smaller than the random number and the new parameters will be kept.
                // If log_diff < 0 (new parameters give worse likelihood), new parameters will only be kept if random number smaller than exp(log_diff).
    //			if ((global_random_number_generator.GetUniformRandom() + 1.0) / 2.0 >= 1.0 / (1.0 + exp(log_diff))) for (i = 0; i < refine_particle.number_of_parameters; i++) {output_parameters[i] = input_parameters[i];}
    //			else output_parameters[16] = output_parameters[15] - input_parameters[15];
                output_parameters[16] = output_parameters[15] - input_parameters[15];
    //			wxPrintf("in, out, diff = %g %g %g\n", input_parameters[15], output_parameters[15], output_parameters[16]);
                if (output_parameters[16] < 0.0) for (i = 0; i < refine_particle.number_of_parameters; i++) {output_parameters[i] = input_parameters[i];}
            }
            else
            {
                input_parameters[15] = - 100.0 * FrealignObjectiveConstrainedFunction(&comparison_object, cg_starting_point);
                output_parameters[15] = input_parameters[15]; output_parameters[16] = 0.0;
            }

            refine_particle.SetParameters(output_parameters);

            refine_particle.SetAlignmentParameters(output_parameters[1], output_parameters[2], output_parameters[3], 0.0, 0.0);
    //		unbinned_image.ClipInto(refine_particle.particle_image);
    //		refine_particle.particle_image->MultiplyByConstant(binning_factor_refine);
    //		refine_particle.particle_image->QuickAndDirtyWriteSlice("part3.mrc", 1);
    //		refine_particle.PhaseFlipImage();
    //		refine_particle.CalculateProjection(projection_image, input_3d);
    //		projection_image.ClipInto(&unbinned_image);
    //		unbinned_image.BackwardFFT();
    //		unbinned_image.ClipInto(&final_image);
    //		logp = refine_particle.ReturnLogLikelihood(input_image, final_image, pixel_size, classification_resolution_limit, alpha, sigma);

    //		logp = refine_particle.ReturnLogLikelihood(input_3d, refine_statistics, classification_resolution_limit);
            output_parameters[13] = refine_particle.ReturnLogLikelihood(input_image, unbinned_image, input_ctf, input_3d, input_statistics, classification_resolution_limit);
    //		output_parameters[14] = sigma * binning_factor_refine;

    //		refine_particle.CalculateMaskedLogLikelihood(projection_image, input_3d, classification_resolution_limit);
    //		output_parameters[13] = refine_particle.logp;
            if (refine_particle.snr > 0.0) output_parameters[14] = sqrtf(1.0 / refine_particle.snr);

    //		output_parameters[14] = refine_particle.sigma_noise * binning_factor_refine;
    //		wxPrintf("logp, sigma, score = %g %g %g\n", output_parameters[13], output_parameters[14], output_parameters[15]);
    //		refine_particle.CalculateProjection(projection_image, input_3d);
    //		projection_image.BackwardFFT();
    //		wxPrintf("snr = %g mask = %g var_A = %g\n", refine_particle.snr, refine_particle.mask_volume, projection_image.ReturnVarianceOfRealValues());
    //		output_parameters[14] = sqrtf(refine_particle.snr * refine_particle.particle_image->number_of_real_space_pixels
    //				/ refine_particle.mask_volume / projection_image.ReturnVarianceOfRealValues()) * binning_factor_refine;

            if (calculate_matching_projections)
            {
                refine_particle.CalculateProjection(projection_image, input_3d);
                projection_image.ClipInto(&unbinned_image);
                unbinned_image.BackwardFFT();
                unbinned_image.ClipInto(&final_image);
                final_image.ForwardFFT();
                final_image.PhaseShift(output_parameters[4] / pixel_size, output_parameters[5] / pixel_size);
                final_image.BackwardFFT();
                final_image.WriteSlice(output_file, image_counter);
            }

            temp_float = input_parameters[1]; input_parameters[1] = input_parameters[3]; input_parameters[3] = temp_float;
            temp_float = output_parameters[1]; output_parameters[1] = output_parameters[3]; output_parameters[3] = temp_float;
            for (i = 1; i < refine_particle.number_of_parameters; i++)
            {
                if (isnanf(output_parameters[i]) != 0)
                {
    //				MyDebugAssertTrue(false, "NaN value for output parameter encountered");
                    output_parameters[i] = input_parameters[i];
                }
            }
            input_parameters[7] = 1;//fabsf(input_parameters[7]);
            output_parameters[7] = input_parameters[7];
            if (output_parameters[15] < 0.0) output_parameters[15] = 0.0;
            my_output_par_file.WriteLine(output_parameters);

            // if (is_running_locally == false) // send results back to the gui..
            // {
            //     intermediate_result = new JobResult;
            //     intermediate_result->job_number = my_current_job.job_number;

            //     gui_result_parameters[0] = current_class;
            //     for (result_parameter_counter = 1; result_parameter_counter < refine_particle.number_of_parameters + 1; result_parameter_counter++)
            //     {
            //         gui_result_parameters[result_parameter_counter] = output_parameters[result_parameter_counter - 1];
            //     }

            //     intermediate_result->SetResult(refine_particle.number_of_parameters + 1, gui_result_parameters);
            //     AddJobToResultQueue(intermediate_result);
            // }

            for (i = 1; i < refine_particle.number_of_parameters; i++) {output_parameter_average[i] += output_parameters[i];}
            for (i = 1; i < refine_particle.number_of_parameters; i++) {output_parameters[i] -= input_parameters[i];}
            for (i = 1; i < refine_particle.number_of_parameters; i++) {output_parameter_change[i] += powf(output_parameters[i],2);}
            my_output_par_shifts_file.WriteLine(output_parameters);

            fflush(my_output_par_file.parameter_file);
            fflush(my_output_par_shifts_file.parameter_file);

            my_progress->Update(image_counter);
        }

        for (i = 1; i < refine_particle.number_of_parameters; i++) {output_parameter_average[i] /= float(last_particle - first_particle + 1);}
        for (i = 1; i < refine_particle.number_of_parameters; i++) {output_parameter_change[i] /= float(last_particle - first_particle + 1);}
        for (i = 1; i < refine_particle.number_of_parameters; i++) {output_parameter_change[i] = sqrtf(output_parameter_change[i]);}
        my_output_par_file.WriteLine(output_parameter_average, true);
        my_output_par_shifts_file.WriteLine(output_parameter_change, true);
        my_output_par_file.WriteCommentLine("C  Total particles included, overall score, average occupancy "
                + wxString::Format("%11i %10.6f %10.6f", last_particle - first_particle + 1, output_parameter_average[15], output_parameter_average[12]));

        delete my_progress;
    //	delete global_euler_search;
        if (global_search)
        {
            delete [] projection_cache;
    //		search_projection_image.RotateFourier2DDeleteIndex(kernel_index, psi_max, psi_step);
        }
        if (calculate_matching_projections) delete output_file;

        wxPrintf("\nRefine3D: Normal termination\n\n");

        return true;
    }
    
}

