#ifndef __AutoRefine3DPanel__
#define __AutoRefine3DPanel__

/**
@file
Subclass of Refine3DPanel, which is generated by wxFormBuilder.
*/

#include "ProjectX_gui.h"

class AutoRefine3DPanel;

class AutoRefinementManager
{
public:
	AutoRefine3DPanel *my_parent;

	bool active_should_auto_mask;
	float active_low_resolution_limit;
	float active_mask_radius;
	float active_global_mask_radius;
	int active_number_results_to_refine;
	float active_search_range_x;
	float active_search_range_y;
	float active_inner_mask_radius;
	bool active_should_apply_blurring;
	float active_smoothing_factor;
	bool active_auto_crop;

	long current_job_starttime;
	long time_of_last_update;
	int number_of_generated_3ds;

	int running_job_type;
	long current_job_id;

	int number_of_rounds_run;

	float next_high_res_limit;

	float start_percent_used;
	float current_percent_used;

	long current_refinement_package_asset_id;
	long current_output_refinement_id;

	long number_of_received_particle_results;
	long expected_number_of_results;

	wxArrayFloat percent_used_per_round;
	wxArrayFloat resolution_per_round;
	wxArrayFloat high_res_limit_per_round;

	wxArrayFloat class_high_res_limits;
	wxArrayFloat class_next_high_res_limits;

	wxArrayInt number_of_global_alignments;
	wxArrayInt rounds_since_global_alignment;

	float classification_resolution;

	float max_percent_used;

	bool this_is_the_final_round;

	RefinementPackage *active_refinement_package;
	Refinement *input_refinement;
	Refinement *output_refinement;

	RunProfile active_refinement_run_profile;
	RunProfile active_reconstruction_run_profile;


	wxArrayString current_reference_filenames;
	wxArrayLong current_reference_asset_ids;

	void SetParent(AutoRefine3DPanel *wanted_parent);

	AutoRefinementManager();
	void BeginRefinementCycle();
	void CycleRefinement();

	void SetupRefinementJob();
	void SetupReconstructionJob();
	void SetupMerge3dJob();

	void SetupInitialReconstructionJob();
	void SetupInitialMerge3dJob();

	void RunInitialReconstructionJob();
	void RunInitialMerge3dJob();

	void RunRefinementJob();
	void RunReconstructionJob();
	void RunMerge3dJob();

	void ProcessJobResult(JobResult *result_to_process);
	void ProcessAllJobsFinished();

	void OnMaskerThreadComplete();

	void DoMasking();

//	void StartRefinement();
//	void StartReconstruction();


};


class AutoRefine3DPanel : public AutoRefine3DPanelParent
{
	friend class AutoRefinementManager;

	protected:
		// Handlers for Refine3DPanel events.
		void OnUpdateUI( wxUpdateUIEvent& event );
		void OnExpertOptionsToggle( wxCommandEvent& event );
		void OnInfoURL( wxTextUrlEvent& event );
		void TerminateButtonClick( wxCommandEvent& event );
		void FinishButtonClick( wxCommandEvent& event );
		void StartRefinementClick( wxCommandEvent& event );
		void ResetAllDefaultsClick( wxCommandEvent& event );

		void OnJobSocketEvent(wxSocketEvent& event);

		int length_of_process_number;

		AutoRefinementManager my_refinement_manager;

		int active_orth_thread_id;
		int active_mask_thread_id;
		int next_thread_id;

	public:


		long time_of_last_result_update;

		bool refinement_package_combo_is_dirty;
		bool run_profiles_are_dirty;
		bool volumes_are_dirty;

		JobResult *buffered_results;
		long my_job_id;
		long selected_refinement_package;

		//int length_of_process_number;

		JobPackage my_job_package;
		JobTracker my_job_tracker;

		bool running_job;

		AutoRefine3DPanel( wxWindow* parent );

		void Reset();
		void SetDefaults();
		void SetInfo();

		void WriteInfoText(wxString text_to_write);
		void WriteErrorText(wxString text_to_write);
		void WriteBlueText(wxString text_to_write);

		void FillRefinementPackagesComboBox();
		void FillRunProfileComboBoxes();

		void NewRefinementPackageSelected();

		void OnRefinementPackageComboBox( wxCommandEvent& event );
		void OnInputParametersComboBox( wxCommandEvent& event );

		void OnMaskerThreadComplete(wxThreadEvent& my_event);
		void OnOrthThreadComplete(ReturnProcessedImageEvent& my_event);
};

#endif
