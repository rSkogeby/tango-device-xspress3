//=============================================================================
//
// file :         AcqThread.h
//
// description :  Include for the AcqThread class.
//                		
// project :      TANGO Device Server
//
//=============================================================================

#ifndef _ACQTHREAD_H
#define _ACQTHREAD_H

#include <tango.h>
#include <Xspress3.h>
#include <pni/core/types.hpp>
#include <pni/io/nexus.hpp>
#include <h5cpp/hdf5.hpp>

namespace Xspress3_ns

{
  class AcqThread : public omni_thread, public Tango::LogAdapter 
  {
  public :
    AcqThread (Xspress3 *xspress3_obj, omni_mutex &m); 
    
    
  private :
    Xspress3		*xspress3;
    omni_mutex 		&mutex;
    
    hdf5::file::File m_nxFile;
    hdf5::node::Group m_nxRootGroup;
    hdf5::node::Group m_nxGroup;
    
    hdf5::node::Group m_nxGroupCh00;
    hdf5::node::Group m_nxGroupCh01;
    hdf5::node::Group m_nxGroupCh02;
    hdf5::node::Group m_nxGroupCh03;
    
    hdf5::node::Group m_nxGroup_scaler[4];
    
    hdf5::node::Dataset m_nxField_histogram[4];
    
    hdf5::node::Dataset m_nxField_allEvent[4];
    hdf5::node::Dataset m_nxField_allGood[4];
    hdf5::node::Dataset m_nxField_inWindow0[4];
    hdf5::node::Dataset m_nxField_inWindow1[4];
    hdf5::node::Dataset m_nxField_pileup[4];
    hdf5::node::Dataset m_nxField_resetCounts[4];
    hdf5::node::Dataset m_nxField_resetTicks[4];
    hdf5::node::Dataset m_nxField_time[4];
    hdf5::node::Dataset m_nxField_totalTicks[4];

    void append_dataset(hdf5::node::Dataset &field,pni::core::float32 value, int index);
        
    hdf5::dataspace::Hyperslab _selection;
    void *run_undetached(void *arg);			
};



}	//	namespace Xspress3_ns

#endif // _ACQTHREAD_H
