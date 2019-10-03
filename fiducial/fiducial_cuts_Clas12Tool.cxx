// this is a Clas12Tool wrapper for Stefan Diehl's fiducial volume cuts code
// found at https://clas12-docdb.jlab.org/DocDB/0004/000444/001/fiducial_cuts.cxx
//
// `F` is an instance of `FiducialCuts`, a class encapsulating the cuts
//
// this wrapper loops over electrons in each event and prints out fiducial volume cut
// results; it was not compiled or tested, so there may be typos....


// cut level setting
F->tight=false;
F->medium=true;
F->loose=false;


// event loop
clas12::clas12reader reader("skimFile.hipo");
while(reader.next()==true) {
  printf("EVENT %d\n\n",reader.runconfig()->getEvent());

  // determine in/out bending
  auto torus = reader.runconfig()->getTorus();
  switch(torus) {
    case -1: F->inbending=true; F->outbending=false; break;
    case 1: F->inbending=false; F->outbending=true; break;
    default: fprintf(stderr,"ERROR: torus unknown\n");
  };

  // particle loop
  for(auto & part : reader.getByID(11)) {

    // set all necessary variables of F for this electron

    // -- PCAL from REC::Calorimeter
    F->pcalLayer = part->cal(clas12::PCAL)->getLayer(); // must == 1
    F->part_Cal_PCAL_sector[0] = part->cal(clas12::PCAL)->getSector();
    F->part_Cal_PCAL_lu[0] = part->cal(clas12::PCAL)->getLu();
    F->part_Cal_PCAL_lv[0] = part->cal(clas12::PCAL)->getLv();
    F->part_Cal_PCAL_lw[0] = part->cal(clas12::PCAL)->getLw();
    // NOTE: needed to add line which sets `_lw_order` to clas12::calorimeter
    //       .cpp file; otherwise calorimeter::getLw() will always return 0
    
    // -- DC from REC::Track
    F->dcTrackDetector = part->trk(clas12::DC)->getDetector(); // must == 6
    F->part_DC_sector[0] = part->trk(clas12::DC)->getSector();

    // -- DC from REC::Traj
    //    `r` loops over three regions (layers 6,18,30)
    //    and regLayer({0,1,2}) returns {6,18,30}
    for(int r=0; r<FiducialCuts::nReg; r++) {

      F->dcTrajDetector[r] = part->traj(clas12::DC,FiducialCuts::regLayer(r))->getDetector(); // must == 6
      F->dcTrajLayer[r] = part->traj(clas12::DC,FiducialCuts::regLayer(r))->getLayer(); // must == {6,18,30}

      // `dcTraj[{0,1,2}][{x,y,z}]` should map to Stefan's `part_DC_c{1,2,3}{x,y,z}[0]`
      F->dcTraj[r][FiducialCuts::x] = part->traj(clas12::DC,FiducialCuts::regLayer(r))->getX();
      F->dcTraj[r][FiducialCuts::y] = part->traj(clas12::DC,FiducialCuts::regLayer(r))->getY();
      F->dcTraj[r][FiducialCuts::z] = part->traj(clas12::DC,FiducialCuts::regLayer(r))->getZ();

    };


    // evaluate cuts
    printf("-- electron p=(%.2f,%.2f,%.2f)\n",
      part->par()->getPx(), part->par()->getPy(), part->par()->getPz() );
    printf("PCAL: %s\n",F->EC_hit_position_fiducial_cut(0) ? "yes":"no");
    printf("DC region 1: %s\n",F->DC_hit_position_region1_fiducial_cut(0) ? "yes":"no");
    printf("DC region 2: %s\n",F->DC_hit_position_region2_fiducial_cut(0) ? "yes":"no");
    printf("DC region 3: %s\n",F->DC_hit_position_region3_fiducial_cut(0) ? "yes":"no");
    printf("\n");

  };
};
