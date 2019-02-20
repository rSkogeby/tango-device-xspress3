static const char *RcsId = "$Header: /cvsroot/tango-ds/training/Metamotor/AcqThread.cpp,v 1.1.1.1 2007/10/04 14:20:47 jensmeyer Exp $";
//+=============================================================================
//
// file :        AcqThread.cpp
//
// description : C++ source for the AcqThread. 
//
// project :     TANGO Device Server
//
//-=============================================================================

#include <AcqThread.h>
#include <math.h>
#include <time.h>
#include "haslib.h"


#include <h5cpp/utilities/array_adapter.hpp>

namespace Xspress3_ns
{
 
  void AcqThread::append_dataset(hdf5::node::Dataset &dataset,pni::core::float32 value, int index)
  {
    dataset.extent(0,1);
    _selection.offset(0,index);
    dataset.write(value, _selection);
  }
  
  AcqThread::AcqThread (Xspress3 *xspress3_obj, omni_mutex &m) : 
	Tango::LogAdapter(xspress3_obj), omni_thread(), mutex(m)
  {
    
    xspress3 = xspress3_obj;
    
    INFO_STREAM << "AcqThread::AcqThread() \n";
    start_undetached();
  }
  
  //+------------------------------------------------------------------
  /**
   *	method:	run_undetached()
   *
   *	description:	
   *
   */
  //+------------------------------------------------------------------
  void *AcqThread::run_undetached(void*)
  {

    _selection = hdf5::dataspace::Hyperslab({0},{1});
    // Clear memory
    xsp3_histogram_clear(xspress3->m_handle, 0, xspress3->nbChannels, 0, *(xspress3->attr_NbFrames_read));

    // Initialize number of frames acquired and read out

    Tango::DevLong m_read_frame_nb = 0; // Number of frames read into buffer
    int use_dtc = 0;
    
    double *pMCA = NULL;
    u_int32_t *pMCA_RAW = NULL;
    u_int32_t *pRAW = NULL;
    u_int32_t *scalar, *mca;
    int xsp3_status = 0;
    Tango::DevLong num_frames_to_read = 0;
    Tango::DevLong last_num_frames = 0;
    
    unsigned int frame = 0;
    int running, counting, paused, finished=0;

    Tango::DevLong frames_in_file = 0;
    Tango::DevLong file_counter = 0;
    
    char filename[200];
    char new_fileprefix[200];
    int idir;
    char system_cmd[400];

    if(*(xspress3->attr_MCALength_read) > 0){
      xspress3->mca_length = *(xspress3->attr_MCALength_read);
    } else {
      xspress3->mca_length = XSP3_MAXSPECTRA;
    }

    hdf5::dataspace::Simple space_simple{{0},{hdf5::dataspace::Simple::UNLIMITED}};
    hdf5::Dimensions chunk_shape{1024};
    hdf5::dataspace::Simple space_hist{{0, xspress3->mca_length},{hdf5::dataspace::Simple::UNLIMITED,hdf5::dataspace::Simple::UNLIMITED}};
    hdf5::Dimensions chunk_hist{1024, 1024};
    
    hdf5::dataspace::Scalar space;
    auto type_float32 = hdf5::datatype::create<pni::core::float32>();
    auto type_int32 = hdf5::datatype::create<pni::core::int32>();

    // Check if directory exists
    
    sprintf(system_cmd, "test -d %s", *(xspress3->attr_FileDir_read));
    idir = system(system_cmd);
    if(idir != 0){
      sprintf(system_cmd, "mkdir -p %s", *(xspress3->attr_FileDir_read));
      system(system_cmd);
    }
    
    // Open file to store data (if only one file)
    
    if(*(xspress3->attr_FramesPerFile_read) == 0){

      
      // Check if file already exists
      sprintf(filename,"%s/%s.nxs", *(xspress3->attr_FileDir_read), *(xspress3->attr_FilePrefix_read));
      ifstream f(filename); 
      if(f.good()){
	sprintf(filename,"%s/%s_%d.nxs", *(xspress3->attr_FileDir_read), *(xspress3->attr_FilePrefix_read), int(time(NULL)));
      }
      
      hdf5::property::FileAccessList fapl;
      fapl.close_degree(hdf5::property::CloseDegree::STRONG);
      hdf5::property::FileCreationList fcpl;
      
      m_nxFile = pni::io::nexus::create_file(filename,  hdf5::file::AccessFlags::TRUNCATE, fcpl, fapl); 
    
      m_nxRootGroup = m_nxFile.root();
    
      
      m_nxGroup = pni::io::nexus::BaseClassFactory::create(m_nxRootGroup,"/entry","NXentry");
      m_nxGroup = pni::io::nexus::BaseClassFactory::create(m_nxGroup,"/entry/instrument","NXdata");
      m_nxGroup = pni::io::nexus::BaseClassFactory::create(m_nxGroup,"/entry/instrument/xspress3","NXdata");
      
      m_nxGroupCh00 = pni::io::nexus::BaseClassFactory::create(m_nxGroup,"/entry/instrument/xspress3/channel00","NXdata");
      m_nxGroupCh01 = pni::io::nexus::BaseClassFactory::create(m_nxGroup,"/entry/instrument/xspress3/channel01","NXdata");
      m_nxGroupCh02 = pni::io::nexus::BaseClassFactory::create(m_nxGroup,"/entry/instrument/xspress3/channel02","NXdata");
      m_nxGroupCh03 = pni::io::nexus::BaseClassFactory::create(m_nxGroup,"/entry/instrument/xspress3/channel03","NXdata");
      
      m_nxGroup_scaler[0] = pni::io::nexus::BaseClassFactory::create(m_nxGroupCh00, "/entry/instrument/xspress3/channel00/scaler","NXdata");
      m_nxGroup_scaler[1] = pni::io::nexus::BaseClassFactory::create(m_nxGroupCh01,"/entry/instrument/xspress3/channel01/scaler","NXdata");
      m_nxGroup_scaler[2] = pni::io::nexus::BaseClassFactory::create(m_nxGroupCh02,"/entry/instrument/xspress3/channel02/scaler","NXdata");
      m_nxGroup_scaler[3] = pni::io::nexus::BaseClassFactory::create(m_nxGroupCh03,"/entry/instrument/xspress3/channel03/scaler","NXdata");
    

      m_nxField_histogram[0] = pni::io::nexus::FieldFactory::create(m_nxGroupCh00, "histogram"  , type_int32, space_hist, chunk_hist);
      m_nxField_histogram[1] = pni::io::nexus::FieldFactory::create(m_nxGroupCh01, "histogram"  , type_int32, space_hist, chunk_hist);
      m_nxField_histogram[2] = pni::io::nexus::FieldFactory::create(m_nxGroupCh02, "histogram"  , type_int32, space_hist, chunk_hist);
      m_nxField_histogram[3] = pni::io::nexus::FieldFactory::create(m_nxGroupCh03, "histogram"  , type_int32, space_hist, chunk_hist);
      
      for(int i = 0; i < 4; i ++){
	  m_nxField_allEvent[i] = pni::io::nexus::FieldFactory::create(m_nxGroup_scaler[i], "allevent", type_float32, space_simple,chunk_shape);
	m_nxField_allGood[i] = pni::io::nexus::FieldFactory::create(m_nxGroup_scaler[i], "allgood", type_float32, space_simple,chunk_shape);
	m_nxField_inWindow0[i] = pni::io::nexus::FieldFactory::create(m_nxGroup_scaler[i], "inwindow0", type_float32, space_simple,chunk_shape);
	m_nxField_inWindow1[i] = pni::io::nexus::FieldFactory::create(m_nxGroup_scaler[i], "inwindow1", type_float32, space_simple,chunk_shape);
	m_nxField_pileup[i] = pni::io::nexus::FieldFactory::create(m_nxGroup_scaler[i], "pileup", type_float32, space_simple,chunk_shape);
	m_nxField_resetCounts[i] = pni::io::nexus::FieldFactory::create(m_nxGroup_scaler[i], "resetcounts", type_float32, space_simple,chunk_shape);
	m_nxField_resetTicks[i] = pni::io::nexus::FieldFactory::create(m_nxGroup_scaler[i], "resetticks", type_float32, space_simple,chunk_shape); 
	m_nxField_time[i] = pni::io::nexus::FieldFactory::create(m_nxGroup_scaler[i], "time", type_float32, space_simple,chunk_shape); 
	m_nxField_totalTicks[i] = pni::io::nexus::FieldFactory::create(m_nxGroup_scaler[i], "totalticks", type_float32, space_simple,chunk_shape); 
      }
      
    }
    
       
    pMCA = (double*)(calloc(XSP3_MAXFRAMES*xspress3->nbChannels*xspress3->mca_length, sizeof(double)));
    pMCA_RAW = (u_int32_t*)(calloc(XSP3_MAXFRAMES*xspress3->nbChannels*xspress3->mca_length, sizeof(u_int32_t)));    
    pRAW = (u_int32_t*)(calloc(XSP3_SW_NUM_SCALERS*xspress3->mca_length*xspress3->nbChannels, sizeof(u_int32_t)));
    
    xsp3_histogram_start(xspress3->m_handle, xspress3->m_card);

    *(xspress3->attr_LastFrame_read) = 0;

    
    while( *(xspress3->attr_LastFrame_read) < *(xspress3->attr_NbFrames_read) && (xspress3->flag_stop_acq == 0)){

      if(xspress3->debug_prints) printf("Acquisition while starts %d %d %d \n", *(xspress3->attr_LastFrame_read), *(xspress3->attr_NbFrames_read), xspress3->flag_stop_acq);
      
      usleep(10000);
      
      u_int32_t tstat;
      if ((xsp3_status = xsp3_get_glob_time_statA(xspress3->m_handle, 0, &tstat)) < 0){
	printf("ERROR: Cannot read timing status: %s\n", xsp3_get_error_message());
      }
      
      frame = XSP3_GLOB_TSTAT_A_FRAME(tstat);
      running = XSP3_GLOB_TSTAT_A_ITFG_RUNNING(tstat);
      counting = XSP3_GLOB_TSTAT_A_ITFG_COUNTING(tstat);
      paused = XSP3_GLOB_TSTAT_A_ITFG_PAUSED(tstat);
      finished = XSP3_GLOB_TSTAT_A_ITFG_FINISHED(tstat);
      
      Xsp3ErrFlag error_flags;
      int64_t cur_frame;
      
      if(xspress3->debug_prints) printf("Acquisition while 1\n");

      usleep(10000);
      
      cur_frame = xsp3_scaler_check_progress_details(xspress3->m_handle, &error_flags, 1, NULL);
      if (cur_frame < 0) {
	 if(xspress3->debug_prints) printf("ERROR calling xsp3_scaler_check_progress_details. Return Message: %s, Flags=%08X\n", xsp3_get_error_message(),  error_flags);
      }
 
      num_frames_to_read = cur_frame - last_num_frames;
      	
      *(xspress3->attr_LastFrame_read) = cur_frame;
      
      if(xspress3->debug_prints) printf("Acquisition while 2\n");
      
      usleep(10000);
      
      if(num_frames_to_read > 0){
	
	if(xspress3->debug_prints) printf("Acquisition while num_frames_to_read (%d) > 0 starts \n", num_frames_to_read);
	
	xsp3_status = xsp3_scaler_read(xspress3->m_handle, pRAW, 0, 0, last_num_frames, XSP3_SW_NUM_SCALERS, xspress3->nbChannels, num_frames_to_read);
	if (use_dtc){
	  xsp3_status = xsp3_hist_dtc_read4d(xspress3->m_handle, pMCA, NULL, 0, 0, 0, last_num_frames, xspress3->mca_length, 1, xspress3->nbChannels, num_frames_to_read);
	} else {
	  xsp3_status = xsp3_histogram_read4d(xspress3->m_handle, pMCA_RAW, 0, 0, 0, last_num_frames, xspress3->mca_length, 1, xspress3->nbChannels, num_frames_to_read);
	}
	
	if(xspress3->debug_prints) printf("Acquisition while num_frames_to_read > 0 1 \n");
	
	xsp3_histogram_circ_ack(xspress3->m_handle, 0, last_num_frames - 1,  xspress3->nbChannels, num_frames_to_read);
	
	if(xspress3->debug_prints) printf("Acquisition while num_frames_to_read > 0 2 \n");
	
	int first_to_process = 0;
	scalar = pRAW+first_to_process*xspress3->nbChannels*XSP3_SW_NUM_SCALERS;
	for(int iframe = 0; iframe < num_frames_to_read; iframe ++){

	  // Open files to store data (when the number of frames per file is reached)
      
	  if( ((*(xspress3->attr_FramesPerFile_read) != 0) &&  ( (frames_in_file >= *(xspress3->attr_FramesPerFile_read)) || file_counter ==0) ) ){
	    
	    if(file_counter == 0){
	      // Check if file already exists
	      sprintf(filename,"%s/%s_00000.nxs", *(xspress3->attr_FileDir_read), *(xspress3->attr_FilePrefix_read));
	      ifstream f(filename); 
	      if(f.good()){
		sprintf(new_fileprefix,"%s_%d", *(xspress3->attr_FilePrefix_read), int(time(NULL)));
	      } else {
		sprintf(new_fileprefix,"%s", *(xspress3->attr_FilePrefix_read));
	      }
	      
	    } else {

    
	      // Close *.nxs datafile
	      for(int i = 0; i<4; i++){
		m_nxField_histogram[i].close();
		
		m_nxField_allEvent[i].close();
		m_nxField_allGood[i].close();
		m_nxField_inWindow0[i].close();
		m_nxField_inWindow1[i].close();
		m_nxField_pileup[i].close();
		m_nxField_resetCounts[i].close();
		m_nxField_resetTicks[i].close();
		m_nxField_time[i].close();
		m_nxField_totalTicks[i].close();
	      }
	      m_nxFile.close();
	    }
	      
	    
	    sprintf(filename,"%s/%s_%05d.nxs", *(xspress3->attr_FileDir_read), new_fileprefix, int(file_counter));
	    
	    
	    frames_in_file = 0;
	    file_counter = file_counter + 1;
	    
    
	    hdf5::property::FileAccessList fapl;
	    fapl.close_degree(hdf5::property::CloseDegree::STRONG);
	    hdf5::property::FileCreationList fcpl;
	    
	    m_nxFile = pni::io::nexus::create_file(filename, hdf5::file::AccessFlags::TRUNCATE, fcpl, fapl); //  overwrite
	    
	    m_nxRootGroup = m_nxFile.root();
	    
	   
	    m_nxGroup = pni::io::nexus::BaseClassFactory::create(m_nxRootGroup,"/entry","NXentry");
	    m_nxGroup = pni::io::nexus::BaseClassFactory::create(m_nxGroup,"/entry/instrument","NXdata");
	    m_nxGroup = pni::io::nexus::BaseClassFactory::create(m_nxGroup,"/entry/instrument/xspress3","NXdata");
	    
	    m_nxGroupCh00 = pni::io::nexus::BaseClassFactory::create(m_nxGroup,"/entry/instrument/xspress3/channel00","NXdata");
	    m_nxGroupCh01 = pni::io::nexus::BaseClassFactory::create(m_nxGroup,"/entry/instrument/xspress3/channel01","NXdata");
	    m_nxGroupCh02 = pni::io::nexus::BaseClassFactory::create(m_nxGroup,"/entry/instrument/xspress3/channel02","NXdata");
	    m_nxGroupCh03 = pni::io::nexus::BaseClassFactory::create(m_nxGroup,"/entry/instrument/xspress3/channel03","NXdata");
	    
	    m_nxGroup_scaler[0] = pni::io::nexus::BaseClassFactory::create(m_nxGroupCh00, "/entry/instrument/xspress3/channel00/scaler","NXdata");
	    m_nxGroup_scaler[1] = pni::io::nexus::BaseClassFactory::create(m_nxGroupCh01,"/entry/instrument/xspress3/channel01/scaler","NXdata");
	    m_nxGroup_scaler[2] = pni::io::nexus::BaseClassFactory::create(m_nxGroupCh02,"/entry/instrument/xspress3/channel02/scaler","NXdata");
	    m_nxGroup_scaler[3] = pni::io::nexus::BaseClassFactory::create(m_nxGroupCh03,"/entry/instrument/xspress3/channel03/scaler","NXdata");
	    
	    
	    m_nxField_histogram[0] = pni::io::nexus::FieldFactory::create(m_nxGroupCh00, "histogram"  , type_int32, space);
	    m_nxField_histogram[1] = pni::io::nexus::FieldFactory::create(m_nxGroupCh01, "histogram"  , type_int32, space);
	    m_nxField_histogram[2] = pni::io::nexus::FieldFactory::create(m_nxGroupCh02, "histogram"  , type_int32, space);
	    m_nxField_histogram[3] = pni::io::nexus::FieldFactory::create(m_nxGroupCh03, "histogram"  , type_int32, space);
	    
	    for(int i = 0; i < 4; i ++){
	      m_nxField_allEvent[i] = pni::io::nexus::FieldFactory::create(m_nxGroup_scaler[i], "allevent", type_float32, space_simple,chunk_shape);
	      m_nxField_allGood[i] = pni::io::nexus::FieldFactory::create(m_nxGroup_scaler[i], "allgood", type_float32, space_simple,chunk_shape);
	      m_nxField_inWindow0[i] = pni::io::nexus::FieldFactory::create(m_nxGroup_scaler[i], "inwindow0", type_float32, space_simple,chunk_shape);
	      m_nxField_inWindow1[i] = pni::io::nexus::FieldFactory::create(m_nxGroup_scaler[i], "inwindow1", type_float32, space_simple,chunk_shape);
	      m_nxField_pileup[i] = pni::io::nexus::FieldFactory::create(m_nxGroup_scaler[i], "pileup", type_float32, space_simple,chunk_shape);
	      m_nxField_resetCounts[i] = pni::io::nexus::FieldFactory::create(m_nxGroup_scaler[i], "resetcounts", type_float32, space_simple,chunk_shape);
	      m_nxField_resetTicks[i] = pni::io::nexus::FieldFactory::create(m_nxGroup_scaler[i], "resetticks", type_float32, space_simple,chunk_shape); 
	      m_nxField_time[i] = pni::io::nexus::FieldFactory::create(m_nxGroup_scaler[i], "time", type_float32, space_simple,chunk_shape); 
	      m_nxField_totalTicks[i] = pni::io::nexus::FieldFactory::create(m_nxGroup_scaler[i], "totalticks", type_float32, space_simple,chunk_shape); 
	    }
	    

	  }
	  
	  for (int ichan=0;ichan<xspress3->nbChannels; ichan++){
	    if(use_dtc){
	      mca = (u_int32_t *)pMCA + xspress3->mca_length*xspress3->nbChannels*iframe + ichan*xspress3->mca_length; // It should not be cast if we want to use it
	    } else { 
	      mca = pMCA_RAW + xspress3->mca_length*xspress3->nbChannels*iframe + ichan*xspress3->mca_length;
	    }
	    
	    if(ichan == 0){
	      xspress3->attr_DataCh1_read = (Tango::DevLong*) mca;
	    } else if (ichan == 1){
	      xspress3->attr_DataCh2_read = (Tango::DevLong*)mca;
	    } else if (ichan == 2){
	      xspress3->attr_DataCh3_read = (Tango::DevLong*)mca;
	    } else if (ichan == 3){
	      xspress3->attr_DataCh4_read = (Tango::DevLong*)mca;
	    }
	    u_int32_t total;
	    total = 0;
	    for (int eng=0; eng<xspress3->mca_length; eng++)
	      total += mca[eng];

	    hdf5::dataspace::Hyperslab ds = hdf5::dataspace::Hyperslab({{0,0},{1,xspress3->mca_length}});

	    m_nxField_histogram[ichan].extent(0,1);
	    
	    ds.offset(0.,frames_in_file);
	    using Int32ArrayAdapter = hdf5::ArrayAdapter<int>;
	    m_nxField_histogram[ichan].write(Int32ArrayAdapter((int*)mca,xspress3->mca_length), ds);
	    
	    append_dataset(m_nxField_allEvent[ichan],(signed long)scalar[ichan*XSP3_SW_NUM_SCALERS+XSP_SW_SCALER_ALL_EVENT], frames_in_file);
	    append_dataset(m_nxField_allGood[ichan],(signed long)scalar[ichan*XSP3_SW_NUM_SCALERS+XSP_SW_SCALER_ALL_GOOD], frames_in_file);
	    append_dataset(m_nxField_inWindow0[ichan],(signed long)scalar[ichan*XSP3_SW_NUM_SCALERS+XSP_SW_SCALER_IN_WINDOW0], frames_in_file);
	    append_dataset(m_nxField_inWindow1[ichan],(signed long)scalar[ichan*XSP3_SW_NUM_SCALERS+XSP_SW_SCALER_IN_WINDOW1], frames_in_file);
	    append_dataset(m_nxField_pileup[ichan],(signed long)scalar[ichan*XSP3_SW_NUM_SCALERS+XSP_SW_SCALER_PILEUP], frames_in_file);
	    append_dataset(m_nxField_resetCounts[ichan],(signed long)scalar[ichan*XSP3_SW_NUM_SCALERS+XSP_SW_SCALER_NUM_RESETS], frames_in_file);
	    append_dataset(m_nxField_resetTicks[ichan],(signed long)scalar[ichan*XSP3_SW_NUM_SCALERS+XSP_SW_SCALER_RESET_TICKS], frames_in_file);
	    append_dataset(m_nxField_time[ichan],(signed long)scalar[ichan*XSP3_SW_NUM_SCALERS+XSP_SW_SCALER_LIVE_TICKS], frames_in_file);
	    append_dataset(m_nxField_totalTicks[ichan],(signed long)scalar[ichan*XSP3_SW_NUM_SCALERS+XSP_SW_SCALER_TOTAL_TICKS], frames_in_file);
	  }
	  frames_in_file = frames_in_file + 1;
	}
      }

      if(xspress3->debug_prints) printf("Acquisition while ends \n");
      last_num_frames = *(xspress3->attr_LastFrame_read);
      
    }

    if(xspress3->debug_prints) printf("Acquisition while close file \n");
    // Close *.nxs datafile

    if(file_counter > 0){
      for(int i = 0; i<4; i++){
	m_nxField_histogram[i].close();
	
	m_nxField_allEvent[i].close();
	m_nxField_allGood[i].close();
	m_nxField_inWindow0[i].close();
	m_nxField_inWindow1[i].close();
	m_nxField_pileup[i].close();
	m_nxField_resetCounts[i].close();
	m_nxField_resetTicks[i].close();
	m_nxField_time[i].close();
	m_nxField_totalTicks[i].close();
      }
      m_nxFile.close();
    }
    
    xspress3->set_state(Tango::ON);
    
    xspress3->flag_stop_acq = 0;
    if(xspress3->debug_prints) printf("Thread return \n");
    
    return NULL;
  }
  
} // namespace Xspress3_ns
