#ifndef FiducialCuts_
#define FiducialCuts_

// dihbsa implementation of Stefan Diehl's fiducial volume cut functions


class FiducialCuts
{
  public:
    FiducialCuts();
    ~FiducialCuts();


    bool EC_hit_position_fiducial_cut(int level);
    bool DC_hit_position_region1_fiducial_cut_triangle(int level);
    bool DC_hit_position_region2_fiducial_cut_triangle(int level);
    bool DC_hit_position_region3_fiducial_cut_triangle(int level);
    bool DC_hit_position_region1_fiducial_cut(int level);
    bool DC_hit_position_region2_fiducial_cut(int level);
    bool DC_hit_position_region3_fiducial_cut(int level);
    bool SetSwitches(int lev);

    
    // PCAL variables
    int pcalSec;
    int pcalLayer; // +++ should always be 1
    enum pcalLenum {u,v,w,nL};
    double pcalL[nL];

    // DC variables
    // - tracks
    int dcSec;
    int dcTrackDetector; // +++ should always be 6
    // - trajectories
    enum dcRegEnum {r1,r2,r3,nReg};
    const int regLayer[nReg] = {6,18,30};
    int dcTrajDetector[nReg]; // +++ should always be 6
    int dcTrajLayer[nReg]; // +++ should be 6,18,30
    enum dcDirEnum {x,y,z,nDir};
    float dcTraj[nRej][nDir];

    
    int torus; // torus B-field (-1=inbending, +1=outbending)

    // cut levels
    enum levelEnum { cutTight, cutMedium, cutLoose };


    bool ErrPrint(const char * str) {
      fprintf(stderr,"ERROR: FiducialCuts::%s\n",str);
      return false;
    };

  protected:
    bool inbending,outbending;
    bool tight,medium,loose;
    char msg[256];

  ClassDef(FiducialCuts,1);
};

#endif
