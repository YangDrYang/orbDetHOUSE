
#ifndef ACS_CONFIG_H
#define ACS_CONFIG_H

#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>
#include <boost/filesystem.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/program_options.hpp>

#include <yaml-cpp/yaml.h>

#include "eigenIncluder.hpp"

// Need the experimental/ if we are compiling with gcc < 8, might have a namespace clash with boost/filesystem
// #include <filesystem>
#include <chrono>
#include <memory>
#include <limits>
#include <vector>
#include <array>
#include <list>
#include <map>
#include <set>

using std::vector;
using std::string;
using std::array;
using std::list;
using std::map;
using std::set;

#include "streamTrace.hpp"
#include "satSys.hpp"
#include "enums.h"


#define PI          3.141592653589793238462643383279502884197169399375105820974
#define D2R         (PI/180.0)          ///< deg to rad
#define R2D         (180.0/PI)          ///< rad to deg

/** Input source filenames and directories
*/
struct InputOptions
{
	vector<string>	atxfiles;
	vector<string>	snxfiles;
	vector<string>	blqfiles;

	vector<string>	navfiles;
	vector<string>	orbfiles;
	vector<string>	sp3files;
	vector<string>	clkfiles;

	vector<string>	erpfiles;
	vector<string>	dcbfiles;
	vector<string>	bsxfiles;
	vector<string>	ionfiles;

	string	root_stations_dir	= "./";
	string	root_input_dir		= "./";
	
	string stream_user			= "";
	string stream_pass			= "";
	string stream_root			= "";
};

/** Enabling and setting destiations of program outputs.
*/
struct OutputOptions
{
	int		trace_level					= 0;
	bool	output_trace				= false;
	string	trace_directory				= "./";
	string	trace_filename				= "<STATION><YYYY><DDD><HH>.trace";
	double	trace_rotate_period			= 60*60*24;

    bool    record_rtcm                 = false;
	string	rtcm_directory	         	= "./";
    string  obs_rtcm_filename			= "<STATION><YYYY><DDD><HH>-OBS.rtcm3";
    string  nav_rtcm_filename			= "<STATION><YYYY><DDD><HH>-NAV.rtcm3";
    double	rtcm_rotate_period			= 60*60*24;

	bool	output_log					= false;
	string	log_directory	         	= "./";
    string  log_filename				= "log<LOGTIME>.json";

    bool	output_residuals 			= false;

	bool	output_config	 			= false;

	bool	output_summary 				= false;
	string	summary_directory			= "./";
	string	summary_filename			= "pea<YYYY><DDD><HH>.summary";

	bool			output_clocks 				= false;
	E_ClkDataSrc 	clocks_receiver_source  	= E_ClkDataSrc::KALMAN;
	E_ClkDataSrc	clocks_satellite_source 	= E_ClkDataSrc::KALMAN;
	string			clocks_directory			= "./";
	string			clocks_filename				= "pea<YYYY><DDD><HH>.clk";
	bool			output_AR_clocks			= false;

	bool	output_rinex_obs			= false;
	string	rinex_obs_directory			= "./";
	string	rinex_obs_filename			= "RINEX_<STATION>_<YYYY><DDD><HH>.obs";
	bool	rinex_obs_print_C_code		= true;
	bool	rinex_obs_print_L_code		= true;
	bool	rinex_obs_print_D_code	 	= true;
	bool	rinex_obs_print_S_code 		= true;
	
	bool	output_ppp_sol 				= false;
	string	ppp_sol_directory			= "./";
	string	ppp_sol_filename			= "pea<YYYY><DDD><HH>.pppsol";

	bool	output_ionex 				= false;
	string	ionex_directory				= "./";
	string	ionex_filename				= "pea<YYYY><DDD><HH>.ionex";

	bool	output_ionstec				= false;
	string	ionstec_directory			= "./";
	string	ionstec_filename			= "pea<YYYY><DDD><HH>.STEC";

	bool	output_biasSINEX			= false;
	string	biasSINEX_directory			= "./";
	string	biasSINEX_filename			= "AUS0ACSRAP_<YYYY><DDD><HH>00_01D_30S_ABS.BIA";

	bool	output_sinex				= false;
	string 	sinex_directory				= "./";
	string	sinex_filename				= "<CONFIG><WWWW><D>.snx";

	bool	output_persistance			= false;
	bool	input_persistance			= false;
	string 	persistance_directory		= "./";
	string	persistance_filename		= "<CONFIG><WWWW><D>.persist";

	bool	output_mongo_measurements				= false;
	bool	output_mongo_states						= false;
	bool	output_intermediate_rts_mongo_states	= false;
	bool	output_mongo_metadata					= false;
	bool	output_mongo_logs						= false;
	bool	delete_mongo_history					= false;
	string	mongo_rts_suffix						= "'";
	string	mongo_uri								= "mongodb://localhost:27017";
	string	mongo_suffix							= "";
};

/** Options to be used only for debugging new features
*/
struct DebugOptions
{
	bool	debug_cs		= false;
	bool	debug_lom		= false;
	bool	debug_csacc		= false;

	int cscase = 0;         /* artificial CS case */
	int csfreq = 3;         /* cycle slip detection and repair frequency */
};

/** Options for unit testing
*/
struct TestOptions
{
	bool	output_pass		= true;
	bool	stop_on_fail	= false;
	bool	stop_on_done	= false;
	bool	output_errors	= false;
	bool	absorb_errors	= false;
	string	directory		= "";
	string	filename		= "";
};

/** Options for the general operation of the software
*/
struct GlobalOptions
{
	double	epoch_interval	= 1;
	int		max_epochs		= 0;

	boost::posix_time::ptime start_epoch	{ boost::posix_time::not_a_date_time };
	boost::posix_time::ptime end_epoch		{ boost::posix_time::not_a_date_time };

	string	config_description			= "";
	string	analysis_agency				= "GAA";
	string	analysis_center				= "Geoscience Australia";
	string	analysis_program			= "AUSACS";
	string	rinex_comment				= "AUSNETWORK1";

	bool    print_stream_statistics     = false;
	bool    caster_test                 = false;
	string  caster_stream_root          = "";

	bool	simulate_real_time			= false;

	bool	process_user				= false;
	bool	process_network 			= false;
	bool	process_minimum_constraints	= false;
	bool	process_ionosphere			= false;
	bool	process_rts					= false;
	bool	process_tests				= false;
	bool	process_ppp					= false;

	bool	process_sys	[E_Sys::NUM_SYS]	= {};
	bool	process_freq[NUM_FTYPES]		= {};

	double	elevation_mask	= 10 * D2R;

	string	pivot_station	= "<AUTO>";
	string	pivot_satellite	= "<AUTO>";

	bool	tide_solid		= false;
	bool	tide_otl		= false;
	bool	tide_pole		= false;

	bool	phase_windup	= true;
	bool	reject_eclipse	= true;
	bool	raim			= true;
	bool	antexacs		= true;
	bool 	clock_jump		= false;

	double	thres_slip   	= 0.05;
	double	max_inno     	= 30;
	double	max_gdop     	= 30;
	double	deweight_factor	= 100;
	double	ratio_limit		= 4000;

	double	wait_next_epoch		= 1;
	double	wait_all_stations	= 0;

	bool	joseph_stabilisation		= false;
	double	validity_interval_factor	= 100;

	list<string>							station_files;

	E_OffsetType ssr_input_antenna_offset = E_OffsetType::APC;

	vector<E_ObsCode>	code_priorities =
	{
		E_ObsCode::L1C,
		E_ObsCode::L1P,
		E_ObsCode::L1Y,
		E_ObsCode::L1W,
		E_ObsCode::L1M,
		E_ObsCode::L1N,
		E_ObsCode::L1S,
		E_ObsCode::L1L,
		E_ObsCode::L1X,

		E_ObsCode::L2W,
		E_ObsCode::L2P,
		E_ObsCode::L2Y,
		E_ObsCode::L2C,
		E_ObsCode::L2M,
		E_ObsCode::L2N,
		E_ObsCode::L2D,
		E_ObsCode::L2S,
		E_ObsCode::L2L,
		E_ObsCode::L2X,

		E_ObsCode::L5I,
		E_ObsCode::L5Q,
		E_ObsCode::L5X
	};

	double clock_wrap_threshold = 0.05e-3;

	//to be removed?

	double proc_noise_iono	= 0.001;		//todo aaron, these need values, or move to other type


	E_Ephemeris 	ppp_ephemeris      	= E_Ephemeris::PRECISE;
	bool 			sat_pcv				= true;
	bool 			rec_pcv				= true;
	double			predefined_fail		= 0.001;	/* pre-defined fail-rate (0.01,0.001) */
};

/** Options associated with kalman filter states
*/
struct KalmanModel
{
	vector<double>	sigma				= {0};	//{0} is very necessary
	vector<double>	apriori_val			= {0};
	vector<double>	proc_noise			= {0};
	vector<double>	tau					= {0};
	
	int				proc_noise_model;
	bool			estimate 			= false;
};

/** Options associated with the network processing mode of operation
*/
struct NetworkOptions
{
	int				phase_reject_limit	= 10;
	int				outage_reset_limit	= 10;

	E_FilterMode	filter_mode		= E_FilterMode::KALMAN;
	E_Inverter		inverter		= E_Inverter::LDLT;

	int				max_filter_iter = 2;
	int				max_prefit_remv = 2;

	int				rts_lag				= 0;
	string			rts_directory		= "./";
	string			rts_filename		= "Network-<YYYY><DDD><HH>.rts";

	KalmanModel		eop;
	KalmanModel		eop_rates;
};

/** Options associated with the ionospheric modelling processing mode of operation
*/
struct IonosphericOptions
{
	E_IonoMode 		corr_mode  			= E_IonoMode::IONO_FREE_LINEAR_COMBO;
	E_LinearCombo	iflc_freqs			= E_LinearCombo::ANY;
	
	bool			common_ionosphere	= false;
};

struct SlipOptions
{
	bool LLI	= true;
	bool GF		= true;
	bool MW		= true;
	bool EMW	= true;
	bool CJ		= true;
	bool SCDIA	= true;
};

struct IonFilterOptions
{
	E_IonoModel		model			= E_IonoModel::NONE;
	double			lat_center;
	double			lon_center;
	double			lat_width;
	double			lon_width;
	double			lat_res;
	double			lon_res;
	double			time_res;
	int				NpLayr;
	int				NBasis;
	int				func_order;
	double			model_noise;
	vector<double>	layer_heights;

	int				max_filter_iter = 2;
	int				max_prefit_remv = 2;
	E_Inverter		inverter		= E_Inverter::LDLT;

	int		rts_lag				= 0;
	string	rts_directory		= "./";
	string	rts_filename		= "Ionosphere-<YYYY><DDD><HH>.rts";

	KalmanModel	ion;
};

struct AmbROptions
{
	E_ARmode	WLmode			= E_ARmode::OFF;		///< Ambiguity resolution mode: OFF, ROUND, ITER_RND, BOOTST, LAMBDA
	E_ARmode	NLmode			= E_ARmode::OFF; 		///< Ambiguity resolution mode: OFF, ROUND, ITER_RND, BOOTST, LAMBDA
	int			lambda_set		= 2;
	int			AR_max_itr		= 1;
	double		min_el_AR		= 15;					///< minimum elevation to attempt ambigity resolution (degrees)

	double WLsuccsThres = 0.9999;	///< Thresholds for ambiguity validation: succsess rate WL
	double WLratioThres = 3;		///< Thresholds for ambiguity validation: succsess rate WL
	double WLSatPrcNois = 0.0001;	///< Process noise for WL satellite biases
	double WLStaPrcNois = 0.001;	///< Process noise for WL station biases

	double NLsuccsThres = 0.9999;	///< Thresholds for ambiguity validation: succsess rate NL
	double NLratioThres = 3;		///< Thresholds for ambiguity validation: succsess rate NL
	double NLstarttime = 3600;		///< Time before starting to calculate (and output NL zmbiguities/biases)

	bool readOSB	 = true;
	bool readDSB	 = true;
	bool readSSRbias = false;
	bool readSATbias = true;
	bool readRecBias = false;
	bool readHYBbias = false;

	bool writeOSB     = false;
	bool writeDSB     = false;
	bool writeSSRbias = false;
	bool writeSATbias = false;
	bool writeRecBias = false;

	double biasOutrate  = 0;		///< Update interval for clock update 0: no output

	bool solvGPS = false;
	bool solvGLO = false;
	bool solvGAL = false;
	bool solvBDS = false;
	bool solvQZS = false;
};

/** Options to set the tropospheric model used in calculations
*/
struct TroposphericOptions
{
	E_TropModel		model		= E_TropModel::VMF3;
	
	string			vmf3dir;
	string			orography;
	string			gpt2grid;
	
	
};	

/** Options to be applied to kalman filter states for individual satellites
*/
struct SatelliteOptions
{
	bool			_initialised	= false;
	bool			exclude			= false;

	KalmanModel		clk;
	KalmanModel		clk_rate;
	KalmanModel		keplers;
	KalmanModel		pos;
	KalmanModel		pos_rate;
	KalmanModel		orb;
	KalmanModel		code_bias;
	KalmanModel		phas_bias;
};

/** Options to be applied to kalman filter states for individual receivers
*/
struct ReceiverOptions
{
	bool			_initialised	= false;
	bool			exclude			= false;

	KalmanModel		amb;
	KalmanModel		pos;
	KalmanModel		pos_rate;
	KalmanModel		clk;
	KalmanModel		clk_rate;
	KalmanModel		keplers;
	KalmanModel		dcb;
	KalmanModel		ion;
	KalmanModel		trop;
	KalmanModel		trop_gauss_markov;
	KalmanModel		trop_grads;
	KalmanModel		trop_grads_gauss_markov;
	KalmanModel		code_bias;
	KalmanModel		phas_bias;

	E_NoiseModel	error_model	= E_NoiseModel::UNIFORM;
	vector<double>	code_sigmas	= {0.1};
	vector<double>	phas_sigmas	= {0.001};
};

/** Minimum constraint options for individual receivers
*/
struct MinimumStationOptions		//todo aaron, move to stations?
{
	bool			_initialised	= false;

	double			noise			= -1;
};

/** Options associated with the minimum constraints mode of operation
*/
struct MinimumConstraintOptions
{
	E_FilterMode	filter_mode				= E_FilterMode::LSQ;
	bool			estimate_scale			= false;
	bool			estimate_rotation		= false;
	bool			estimate_translation	= false;

	map<string,	MinimumStationOptions>		stationMap;
};

/** Options associated with the user mode (PPP) mode of operation
*/
struct PPPOptions
{
	int			phase_reject_limit	= 10;
	int			outage_reset_limit	= 10;

	int			max_filter_iter 	= 2;
	int			max_prefit_remv 	= 2;

	E_Inverter	inverter			= E_Inverter::LDLT;

	int			rts_lag				= 0;
	string		rts_directory		= "./";
	string		rts_filename		= "PPP-<Station>-<YYYY><DDD><HH>.rts";
};

/** Options associated with cycle slip detection and repair within the network filter
*/
struct CycleSlipOptions
{
	bool	enable							= false;
	bool	print_activity					= false;
	bool	debug							= false;
	bool	timer_debug						= false;
	double	freq_l1							= 0;
	double	freq_l2							= 0;
	int		new_channel_break_in_duration	= 0;
	int		d_moving_ave_win				= 0;
	string	common_mode_method				= "";
	int		slip_classify_min_pts			= 0;
	int		slip_classify_post_outlier_pts	= 0;
	double	outlier_p_alpha					= 0;
	double	t_alpha							= 0;
	double	min_size_p_alpha				= 0;
	int		outlier_detect_min_pts			= 0;
	double	int_valid_pdf_thresh			= 0;
	double	int_valid_outlier_alpha			= 0;
	double	deweight_noise					= 0;
	double	int_valid_combo_jump_alpha		= 0;
};

/** Options associated with SSR corrections and exporting RTCM messages
*/
struct SsrOptions
{
	bool			calculate_ssr			= false;
	int				update_interval			= 0;
	E_Ephemeris 	ssr_ephemeris_source	= E_Ephemeris::PRECISE;
	bool			upload_to_caster		= false;
	string			rtcm_directory			= "./";
};

/** Options associated with orbital force models
*/
struct ForceModels
{
	bool		central_force_gravity			= true;
	bool		direct_solar_radiation			= true;
	double		some_configuration_parameter	= 6;
};

/** General options object to be used throughout the software
*/
struct ACSConfig : GlobalOptions, InputOptions, OutputOptions, DebugOptions
{
	YAML::Node	yaml;
	map<string, string>	yamlDefaults;
	

	string												configFilename;
	map<string, time_t>									configModifyTimeMap;
	boost::program_options::variables_map				commandOpts;

	bool	parse(string filename, boost::program_options::variables_map& vm);
	bool	parse();
	void	info(Trace& trace);

	void	addDataFile(
		string fileName,
		string fileType,
		string dataType);

	SatelliteOptions&			getSatOpts		(SatSys&	Sat);
	ReceiverOptions&			getRecOpts		(string		id);
	MinimumStationOptions&		getMinConOpts	(string 	id);

	map<string,		SatelliteOptions>	satOptsMap;
	map<string,		ReceiverOptions>	recOptsMap;

	PPPOptions					pppOpts;
	TroposphericOptions			tropOpts;
	IonosphericOptions			ionoOpts;
	NetworkOptions				netwOpts;
	MinimumConstraintOptions	minCOpts;
	TestOptions					testOpts;
	CycleSlipOptions			csOpts;
	SsrOptions					ssrOpts;
	AmbROptions					ambrOpts;
	SlipOptions					excludeSlip;
	ForceModels					forceModels;
	IonFilterOptions			ionFilterOpts;

	KalmanModel	ion;
	list<string>	rinexFiles;
};






void replaceString(
	string&	str,
	string	subStr,
	string	replacement);

void removePath(
	string &filepath); 	// fully qualified file

void tryAddRootToPath(
	string& root,		///< Root path
	string& path);		///< Filename to prepend root path to

void replaceTimes(
	string&						str,		///< String to replace macros within
	boost::posix_time::ptime	time_time);	///< Time to use for replacements

void replaceTags(
	string&						str);		///< String to replace macros within

bool configure(
	int		argc,
	char**	argv);

bool checkValidFile(
	string&	path,				///< Filename to check
	string	description = "");	///< Description for error messages

void dumpConfig(
	Trace& trace);

extern ACSConfig acsConfig;		///< Global variable housing all options to be used throughout the software



#endif
